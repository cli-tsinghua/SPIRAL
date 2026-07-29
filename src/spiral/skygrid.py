"""Sky-grid helpers for SPIRAL map products."""

from __future__ import annotations

from dataclasses import dataclass

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import numpy as np


@dataclass(frozen=True)
class Tile:
    tile_id: str
    cloud: str
    ra_deg: float
    dec_deg: float
    tile_size_arcmin: float


def make_tan_wcs(center_ra: float, center_dec: float, pixel_scale_arcsec: float) -> WCS:
    pixel_scale_deg = pixel_scale_arcsec / 3600.0
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.crval = [center_ra, center_dec]
    wcs.wcs.cdelt = [-pixel_scale_deg, pixel_scale_deg]
    wcs.wcs.crpix = [1.0, 1.0]
    return wcs


def wcs_from_tile_bounds(tiles: list[Tile], pixel_scale_arcsec: float, margin_pix: int = 3) -> tuple[WCS, tuple[int, int]]:
    if not tiles:
        raise ValueError("At least one tile is required to define the output grid.")
    ra: list[float] = []
    dec: list[float] = []
    for tile in tiles:
        half = tile.tile_size_arcmin / 120.0
        ra.extend([tile.ra_deg - half, tile.ra_deg + half, tile.ra_deg - half, tile.ra_deg + half])
        dec.extend([tile.dec_deg - half, tile.dec_deg - half, tile.dec_deg + half, tile.dec_deg + half])
    ra_arr = np.asarray(ra)
    dec_arr = np.asarray(dec)
    center_ra = float(np.nanmedian(ra_arr))
    center_dec = float(np.nanmedian(dec_arr))
    trial = make_tan_wcs(center_ra, center_dec, pixel_scale_arcsec)
    x, y = trial.world_to_pixel_values(ra_arr, dec_arr)
    xmin = float(np.floor(np.nanmin(x))) - margin_pix
    xmax = float(np.ceil(np.nanmax(x))) + margin_pix
    ymin = float(np.floor(np.nanmin(y))) - margin_pix
    ymax = float(np.ceil(np.nanmax(y))) + margin_pix
    nx = int(xmax - xmin + 1)
    ny = int(ymax - ymin + 1)
    out_wcs = make_tan_wcs(center_ra, center_dec, pixel_scale_arcsec)
    out_wcs.wcs.crpix = [1.0 - xmin, 1.0 - ymin]
    out_wcs.array_shape = (ny, nx)
    return out_wcs, (ny, nx)


def nearest_tile_assignment(tiles: list[Tile], out_wcs: WCS, shape: tuple[int, int]) -> tuple[np.ndarray, dict[str, int]]:
    tile_ids = [tile.tile_id for tile in tiles]
    tile_coords = SkyCoord(
        np.asarray([tile.ra_deg for tile in tiles]) * u.deg,
        np.asarray([tile.dec_deg for tile in tiles]) * u.deg,
    )
    yy, xx = np.indices(shape, dtype=float)
    ra, dec = out_wcs.pixel_to_world_values(xx, yy)
    cell_coords = SkyCoord(ra.ravel() * u.deg, dec.ravel() * u.deg)
    nearest, _sep, _dist3d = cell_coords.match_to_catalog_sky(tile_coords)
    assignment = nearest.reshape(shape).astype(np.int16)
    return assignment, {tile_id: index for index, tile_id in enumerate(tile_ids)}

