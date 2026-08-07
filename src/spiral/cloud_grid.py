"""Build cloud-scale SPHEREx feature maps from public L2 cutouts."""

from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from .features import (
    FEATURE_CONFIGS,
    FeatureConfig,
    TemplateFitConfig,
    add_anchored_template_feature_samples,
    add_template_fit_samples,
    add_weighted_window_samples,
    make_template_fit_accumulator,
    make_weighted_window_accumulator,
    make_anchored_template_feature_accumulator,
    measure_feature_from_anchored_template,
    measure_feature_from_template_fit,
    measure_feature_from_weighted_windows,
    measure_feature_from_windows,
)
from .flags import bad_flag_mask_from_header
from .l2 import read_l2_image, subtract_static_zodi
from .skygrid import Tile, nearest_tile_assignment, wcs_from_tile_bounds
from .wavelengths import cwave_maps_from_header


WINDOW_NAMES = ("BLUE", "FEATURE", "RED")
SAMPLE_METHODS = (
    "nearest_equal_window",
    "nearest_equal_template_window",
    "bilinear_ivar_window",
    "bilinear_ivar_template",
)
NEAREST_SAMPLE_METHODS = {"nearest_equal_window", "nearest_equal_template_window"}
TEMPLATE_SAMPLE_METHODS = {"nearest_equal_template_window", "bilinear_ivar_template"}


def load_tiles(tile_catalog: str | Path, cloud: str) -> list[Tile]:
    tiles: list[Tile] = []
    with Path(tile_catalog).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["cloud"].upper() != cloud.upper():
                continue
            tiles.append(
                Tile(
                    tile_id=row["spherex_target_id"],
                    cloud=row["cloud"].upper(),
                    ra_deg=float(row["ra_deg"]),
                    dec_deg=float(row["dec_deg"]),
                    tile_size_arcmin=float(row["tile_size_arcmin"]),
                )
            )
    if not tiles:
        raise FileNotFoundError(f"No {cloud} rows found in tile catalog {tile_catalog}")
    return tiles


def raw_cutouts(raw_root: str | Path, cloud: str, cfg: FeatureConfig) -> list[tuple[Path, str]]:
    root = Path(raw_root).expanduser().resolve()
    paths = sorted(root.glob(f"{cloud.lower()}_tile_*/data/{cfg.band}/*.fits"))
    return [(path, path.parents[2].name) for path in paths]


def make_accumulators(shape: tuple[int, int]) -> dict[str, dict[str, np.ndarray]]:
    return {
        name: {
            "sum": np.zeros(shape, dtype=float),
            "sum2": np.zeros(shape, dtype=float),
            "sumlambda": np.zeros(shape, dtype=float),
            "count": np.zeros(shape, dtype=np.int32),
        }
        for name in WINDOW_NAMES
    }


def add_values(acc: dict[str, np.ndarray], y: np.ndarray, x: np.ndarray, values: np.ndarray, lambdas: np.ndarray) -> None:
    np.add.at(acc["sum"], (y, x), values)
    np.add.at(acc["sum2"], (y, x), values * values)
    np.add.at(acc["sumlambda"], (y, x), lambdas)
    np.add.at(acc["count"], (y, x), 1)


def make_weighted_window_accumulators(shape: tuple[int, int]) -> dict[str, dict[str, np.ndarray]]:
    return {name: make_weighted_window_accumulator(shape) for name in WINDOW_NAMES}


def inverse_variance_weights(variance: np.ndarray, good: np.ndarray, *, cap_percentile: float = 99.5) -> np.ndarray:
    values = variance[good & np.isfinite(variance) & (variance > 0)]
    weights = np.divide(1.0, variance, out=np.zeros(variance.shape, dtype=float), where=np.isfinite(variance) & (variance > 0))
    if values.size:
        cap = np.nanpercentile(1.0 / values, cap_percentile)
        if np.isfinite(cap) and cap > 0:
            weights = np.minimum(weights, cap)
    return weights


def add_bilinear_weighted_window(
    acc: dict[str, np.ndarray],
    gx: np.ndarray,
    gy: np.ndarray,
    values: np.ndarray,
    lambdas: np.ndarray,
    weights: np.ndarray,
    assignment: np.ndarray,
    this_tile_index: int,
) -> int:
    x0 = np.floor(gx).astype(np.int64)
    y0 = np.floor(gy).astype(np.int64)
    dx = gx - x0
    dy = gy - y0
    shape = assignment.shape
    used = 0
    for ox, wx in ((0, 1.0 - dx), (1, dx)):
        for oy, wy in ((0, 1.0 - dy), (1, dy)):
            x = x0 + ox
            y = y0 + oy
            geom = wx * wy
            inside = (geom > 0) & (x >= 0) & (x < shape[1]) & (y >= 0) & (y < shape[0])
            if np.any(inside):
                ok = np.zeros_like(inside, dtype=bool)
                ok[inside] = assignment[y[inside], x[inside]] == this_tile_index
                inside &= ok
            if np.any(inside):
                add_weighted_window_samples(
                    acc,
                    y[inside],
                    x[inside],
                    values[inside],
                    lambdas[inside],
                    weights[inside] * geom[inside],
                    geom[inside],
                )
                used += int(np.count_nonzero(inside))
    return used


def add_bilinear_template_fit(
    acc: dict[str, np.ndarray],
    gx: np.ndarray,
    gy: np.ndarray,
    values: np.ndarray,
    lambdas: np.ndarray,
    weights: np.ndarray,
    assignment: np.ndarray,
    this_tile_index: int,
    template: TemplateFitConfig,
) -> int:
    x0 = np.floor(gx).astype(np.int64)
    y0 = np.floor(gy).astype(np.int64)
    dx = gx - x0
    dy = gy - y0
    shape = assignment.shape
    used = 0
    for ox, wx in ((0, 1.0 - dx), (1, dx)):
        for oy, wy in ((0, 1.0 - dy), (1, dy)):
            x = x0 + ox
            y = y0 + oy
            geom = wx * wy
            inside = (geom > 0) & (x >= 0) & (x < shape[1]) & (y >= 0) & (y < shape[0])
            if np.any(inside):
                ok = np.zeros_like(inside, dtype=bool)
                ok[inside] = assignment[y[inside], x[inside]] == this_tile_index
                inside &= ok
            if np.any(inside):
                add_template_fit_samples(
                    acc,
                    y[inside],
                    x[inside],
                    values[inside],
                    lambdas[inside],
                    weights[inside] * geom[inside],
                    template,
                    geom[inside],
                )
                used += int(np.count_nonzero(inside))
    return used


def build_cloud_grid(
    *,
    feature: str,
    cloud: str,
    raw_root: str | Path,
    tile_catalog: str | Path,
    cal_root: str | Path,
    pixel_scale_arcsec: float,
    flag_policy: str = "science_strict",
    sample_method: str = "nearest_equal_window",
    template_fwhm_um: float | None = None,
    ivar_cap_percentile: float = 99.5,
    max_files: int | None = None,
    progress_every: int = 1000,
) -> tuple[dict[str, np.ndarray], WCS, dict[str, object]]:
    if sample_method not in SAMPLE_METHODS:
        raise ValueError(f"Unknown sample_method={sample_method!r}; choose one of {SAMPLE_METHODS}")
    cfg = FEATURE_CONFIGS[feature]
    cloud_upper = cloud.upper()
    tiles = load_tiles(tile_catalog, cloud_upper)
    (out_wcs, shape) = wcs_from_tile_bounds(tiles, pixel_scale_arcsec)
    assignment, tile_index = nearest_tile_assignment(tiles, out_wcs, shape)
    cutouts = raw_cutouts(raw_root, cloud_upper, cfg)
    if max_files is not None:
        cutouts = cutouts[:max_files]
    if not cutouts:
        raise FileNotFoundError(f"No {cfg.band} FITS cutouts found under {raw_root} for {cloud_upper}")

    windows = {
        "BLUE": (cfg.blue_min, cfg.blue_max),
        "FEATURE": (cfg.feature_min, cfg.feature_max),
        "RED": (cfg.red_min, cfg.red_max),
    }
    if sample_method == "nearest_equal_window":
        accum = make_accumulators(shape)
        template_accum = None
        anchored_feature_accum = None
    elif sample_method == "nearest_equal_template_window":
        accum = make_accumulators(shape)
        accum.pop("FEATURE")
        anchored_feature_accum = make_anchored_template_feature_accumulator(shape)
        template_accum = None
    elif sample_method == "bilinear_ivar_window":
        accum = make_weighted_window_accumulators(shape)
        template_accum = None
        anchored_feature_accum = None
    else:
        accum = None
        anchored_feature_accum = None
        template_accum = make_template_fit_accumulator(shape)
    if template_fwhm_um is None:
        template_fwhm_um = 0.36 * cfg.width_um
    template = TemplateFitConfig(center_um=cfg.center_um, fwhm_um=float(template_fwhm_um))
    version_counts: Counter[str] = Counter()
    files_used_by_window: Counter[str] = Counter()
    pixels_used_by_window: Counter[str] = Counter()
    skipped = 0
    grid_cache: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}

    for index, (path, tile_id) in enumerate(cutouts, start=1):
        if progress_every > 0 and (index == 1 or index % progress_every == 0 or index == len(cutouts)):
            print(f"{index}/{len(cutouts)} {tile_id} {path.name}", flush=True)
        this_tile_index = tile_index.get(tile_id)
        if this_tile_index is None:
            skipped += 1
            continue
        try:
            l2 = read_l2_image(path)
            bad_flag_mask = bad_flag_mask_from_header(l2.flags_header, policy=flag_policy)
            lambda_map, _dlambda_map = cwave_maps_from_header(l2.image_header, l2.image.shape, str(cal_root))
            data, _zodi_scale = subtract_static_zodi(l2.image, l2.zodi)
        except Exception as exc:  # noqa: BLE001 - keep production runs moving over bad files.
            print(f"skip {path}: {type(exc).__name__}: {exc}", flush=True)
            skipped += 1
            continue

        base_good = (
            np.isfinite(data)
            & np.isfinite(lambda_map)
            & np.isfinite(l2.variance)
            & (l2.variance > 0)
            & ((l2.flags & bad_flag_mask) == 0)
        )
        if not np.any(base_good):
            continue

        if l2.image.shape not in grid_cache:
            yy, xx = np.indices(l2.image.shape, dtype=float)
            grid_cache[l2.image.shape] = (yy, xx)
        yy, xx = grid_cache[l2.image.shape]

        if sample_method in {"nearest_equal_window", "nearest_equal_template_window"}:
            used_any = False
            for name, (lo, hi) in windows.items():
                good = base_good & (lambda_map >= lo) & (lambda_map <= hi)
                if not np.any(good):
                    continue
                idx = np.where(good)
                sky = l2.spatial_wcs.pixel_to_world(xx[idx], yy[idx])
                gx, gy = out_wcs.world_to_pixel_values(sky.ra.deg, sky.dec.deg)
                x_index = np.floor(gx + 0.5).astype(np.int64)
                y_index = np.floor(gy + 0.5).astype(np.int64)
                inside = (x_index >= 0) & (x_index < shape[1]) & (y_index >= 0) & (y_index < shape[0])
                if np.any(inside):
                    assigned_to_this_tile = np.zeros_like(inside, dtype=bool)
                    assigned_to_this_tile[inside] = assignment[y_index[inside], x_index[inside]] == this_tile_index
                    inside &= assigned_to_this_tile
                if not np.any(inside):
                    continue
                assert accum is not None
                if sample_method == "nearest_equal_template_window" and name == "FEATURE":
                    assert anchored_feature_accum is not None
                    add_anchored_template_feature_samples(
                        anchored_feature_accum,
                        y_index[inside],
                        x_index[inside],
                        data[idx][inside],
                        lambda_map[idx][inside],
                        template,
                    )
                else:
                    add_values(accum[name], y_index[inside], x_index[inside], data[idx][inside], lambda_map[idx][inside])
                files_used_by_window[name] += 1
                pixels_used_by_window[name] += int(np.count_nonzero(inside))
                used_any = True
            if used_any:
                version_counts[str(l2.primary_header.get("VERSION", "UNKNOWN"))] += 1
            continue

        all_lo = min(lo for lo, _hi in windows.values())
        all_hi = max(hi for _lo, hi in windows.values())
        good = base_good & (lambda_map >= all_lo) & (lambda_map <= all_hi)
        if not np.any(good):
            continue
        idx = np.where(good)
        sky = l2.spatial_wcs.pixel_to_world(xx[idx], yy[idx])
        gx, gy = out_wcs.world_to_pixel_values(sky.ra.deg, sky.dec.deg)
        values = data[idx]
        lambdas = lambda_map[idx]
        ivar = inverse_variance_weights(l2.variance, base_good, cap_percentile=ivar_cap_percentile)[idx]
        finite = np.isfinite(gx) & np.isfinite(gy) & np.isfinite(values) & np.isfinite(lambdas) & np.isfinite(ivar) & (ivar > 0)
        if not np.any(finite):
            continue
        gx = gx[finite]
        gy = gy[finite]
        values = values[finite]
        lambdas = lambdas[finite]
        ivar = ivar[finite]

        used_any = False
        if sample_method == "bilinear_ivar_window":
            assert accum is not None
            for name, (lo, hi) in windows.items():
                in_window = (lambdas >= lo) & (lambdas <= hi)
                if not np.any(in_window):
                    continue
                used = add_bilinear_weighted_window(
                    accum[name],
                    gx[in_window],
                    gy[in_window],
                    values[in_window],
                    lambdas[in_window],
                    ivar[in_window],
                    assignment,
                    this_tile_index,
                )
                if used:
                    files_used_by_window[name] += 1
                    pixels_used_by_window[name] += used
                    used_any = True
        else:
            assert template_accum is not None
            used = add_bilinear_template_fit(template_accum, gx, gy, values, lambdas, ivar, assignment, this_tile_index, template)
            if used:
                files_used_by_window["TEMPLATE"] += 1
                pixels_used_by_window["TEMPLATE"] += used
                used_any = True
        if used_any:
            version_counts[str(l2.primary_header.get("VERSION", "UNKNOWN"))] += 1

    if sample_method == "nearest_equal_window":
        assert accum is not None
        maps = measure_feature_from_windows(cfg, accum["BLUE"], accum["FEATURE"], accum["RED"])
    elif sample_method == "nearest_equal_template_window":
        assert accum is not None
        assert anchored_feature_accum is not None
        maps = measure_feature_from_anchored_template(cfg, accum["BLUE"], anchored_feature_accum, accum["RED"], template)
    elif sample_method == "bilinear_ivar_window":
        assert accum is not None
        maps = measure_feature_from_weighted_windows(cfg, accum["BLUE"], accum["FEATURE"], accum["RED"])
    else:
        assert template_accum is not None
        maps = measure_feature_from_template_fit(cfg, template_accum, template)
    meta = {
        "pipeline": "SPIRAL",
        "pipeline_version": "0.1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": sample_method,
        "spatial_resampling": "nearest output pixel"
        if sample_method in NEAREST_SAMPLE_METHODS
        else "bilinear deposition to four output pixels",
        "sample_weighting": "equal weight"
        if sample_method in NEAREST_SAMPLE_METHODS
        else f"inverse L2 variance with {ivar_cap_percentile:.3g} percentile cap per cutout",
        "template_fwhm_um": float(template_fwhm_um) if sample_method in TEMPLATE_SAMPLE_METHODS else None,
        "zodi": "static ZODI subtraction with scale fixed to 1",
        "wavelength_source": "calibrated CWAVE/CBAND spectral-WCS products",
        "flag_policy": flag_policy,
        "cloud": cloud_upper,
        "feature": cfg.label,
        "feature_name": cfg.name,
        "pixel_scale_arcsec": float(pixel_scale_arcsec),
        "shape": list(shape),
        "n_tiles": len(tiles),
        "n_raw_cutouts": len(cutouts),
        "n_skipped_files": skipped,
        "version_counts": dict(version_counts),
        "files_used_by_window": dict(files_used_by_window),
        "pixels_used_by_window": dict(pixels_used_by_window),
    }
    return maps, out_wcs, meta


def write_cloud_grid(
    maps: dict[str, np.ndarray],
    out_wcs: WCS,
    meta: dict[str, object],
    output_root: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    scale = str(meta["pixel_scale_arcsec"]).replace(".", "p")
    method = str(meta.get("method", "nearest_equal_window"))
    suffix = "" if method == "nearest_equal_window" else f"_{method}"
    output = root / f"spiral_{meta['feature']}_{str(meta['cloud']).lower()}_{scale}arcsec{suffix}.fits"
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output exists; pass --overwrite to replace it: {output}")

    header = out_wcs.to_header()
    header["PIPELINE"] = "SPIRAL"
    header["GALAXY"] = str(meta["cloud"])
    header["FEATURE"] = str(meta["feature"])
    header["PIXSCALE"] = float(meta["pixel_scale_arcsec"])
    header["MAPMETH"] = method[:68]
    header["ZODI"] = "STATIC"
    hdus = [fits.PrimaryHDU(header=header)]
    units = {
        "CONTINUUM": "MJy/sr",
        "EXCESS": "MJy/sr um",
        "EXCUNC": "MJy/sr um",
        "SNR": "",
        "EW": "um",
        "BLUE": "MJy/sr",
        "FEATURE": "MJy/sr",
        "RED": "MJy/sr",
        "LAMBLUE": "um",
        "LAMFEAT": "um",
        "LAMRED": "um",
        "NTOTAL": "count",
        "NBLUE": "count",
        "NFEAT": "count",
        "NRED": "count",
        "SLOPE": "MJy/sr/um",
    }
    for name in units:
        if name not in maps:
            continue
        hdu = fits.ImageHDU(data=maps[name], header=header, name=name)
        if units[name]:
            hdu.header["BUNIT"] = units[name]
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(output, overwrite=True)
    with output.with_suffix(".json").open("w", encoding="utf-8") as handle:
        json.dump(meta | {"fits": output.name}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return output
