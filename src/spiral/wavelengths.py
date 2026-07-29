"""SPHEREx wavelength calibration helpers."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from astropy.io import fits
import numpy as np


def detector_from_header(image_header: fits.Header) -> int:
    if "DETECTOR" in image_header:
        return int(image_header["DETECTOR"])
    if "DET_ID" in image_header:
        return int(image_header["DET_ID"])
    raise KeyError("Cannot determine SPHEREx detector from IMAGE header")


@lru_cache(maxsize=16)
def load_spectral_wcs_calibration(detector: int, cal_root: str) -> tuple[np.ndarray, np.ndarray]:
    """Load SPHEREx CWAVE/CBAND calibration arrays for one detector."""
    root = Path(cal_root).expanduser().resolve()
    path = root / str(detector) / f"spectral_wcs_D{detector}_spx_cal-wcs-v4-2025-254.fits"
    if not path.exists():
        raise FileNotFoundError(f"Missing SPHEREx spectral-WCS calibration product: {path}")
    with fits.open(path, memmap=False) as hdul:
        cwave = hdul["CWAVE"].data.astype(float)
        cband = hdul["CBAND"].data.astype(float)
    return cwave, cband


def cwave_maps_from_header(image_header: fits.Header, shape: tuple[int, int], cal_root: str) -> tuple[np.ndarray, np.ndarray]:
    """Map each L2 image pixel to its calibrated central wavelength and bandwidth."""
    detector = detector_from_header(image_header)
    cwave, cband = load_spectral_wcs_calibration(detector, cal_root)
    yy, xx = np.indices(shape, dtype=float)
    ix = np.rint(xx + 1.0 - float(image_header["CRPIX1W"])).astype(np.int64)
    iy = np.rint(yy + 1.0 - float(image_header["CRPIX2W"])).astype(np.int64)
    inside = (ix >= 0) & (ix < cwave.shape[1]) & (iy >= 0) & (iy < cwave.shape[0])
    lambda_map = np.full(shape, np.nan, dtype=float)
    dlambda_map = np.full(shape, np.nan, dtype=float)
    lambda_map[inside] = cwave[iy[inside], ix[inside]]
    dlambda_map[inside] = cband[iy[inside], ix[inside]]
    return lambda_map, dlambda_map

