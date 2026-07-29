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

from .features import FEATURE_CONFIGS, FeatureConfig, measure_feature_from_windows
from .flags import bad_flag_mask_from_header
from .l2 import read_l2_image, subtract_static_zodi
from .skygrid import Tile, nearest_tile_assignment, wcs_from_tile_bounds
from .wavelengths import cwave_maps_from_header


WINDOW_NAMES = ("BLUE", "FEATURE", "RED")


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


def build_cloud_grid(
    *,
    feature: str,
    cloud: str,
    raw_root: str | Path,
    tile_catalog: str | Path,
    cal_root: str | Path,
    pixel_scale_arcsec: float,
    flag_policy: str = "science_strict",
    max_files: int | None = None,
    progress_every: int = 1000,
) -> tuple[dict[str, np.ndarray], WCS, dict[str, object]]:
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
    accum = make_accumulators(shape)
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
            add_values(accum[name], y_index[inside], x_index[inside], data[idx][inside], lambda_map[idx][inside])
            files_used_by_window[name] += 1
            pixels_used_by_window[name] += int(np.count_nonzero(inside))
            used_any = True
        if used_any:
            version_counts[str(l2.primary_header.get("VERSION", "UNKNOWN"))] += 1

    maps = measure_feature_from_windows(cfg, accum["BLUE"], accum["FEATURE"], accum["RED"])
    meta = {
        "pipeline": "SPIRAL",
        "pipeline_version": "0.1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "direct equal-weight L2 pixel binning onto cloud science grid",
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
    output = root / f"spiral_{meta['feature']}_{str(meta['cloud']).lower()}_{scale}arcsec.fits"
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output exists; pass --overwrite to replace it: {output}")

    header = out_wcs.to_header()
    header["PIPELINE"] = "SPIRAL"
    header["GALAXY"] = str(meta["cloud"])
    header["FEATURE"] = str(meta["feature"])
    header["PIXSCALE"] = float(meta["pixel_scale_arcsec"])
    header["MAPMETH"] = "EQW_L2_BIN"
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
    }
    for name in units:
        hdu = fits.ImageHDU(data=maps[name], header=header, name=name)
        if units[name]:
            hdu.header["BUNIT"] = units[name]
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(output, overwrite=True)
    with output.with_suffix(".json").open("w", encoding="utf-8") as handle:
        json.dump(meta | {"fits": output.name}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return output

