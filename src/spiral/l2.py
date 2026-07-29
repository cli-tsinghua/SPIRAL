"""Input/output helpers for public SPHEREx L2 spectral-image cutouts."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import warnings

from astropy import log as astropy_log
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


@contextmanager
def suppress_astropy_wcs_noise():
    original_level = astropy_log.level
    astropy_log.setLevel("WARNING")
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Removed redundant SCAMP distortion parameters.*")
            warnings.filterwarnings("ignore", message=".*because SIP parameters are also present.*")
            yield
    finally:
        astropy_log.setLevel(original_level)


@dataclass(frozen=True)
class L2Image:
    image: np.ndarray
    variance: np.ndarray
    flags: np.ndarray
    zodi: np.ndarray
    spatial_wcs: WCS
    image_header: fits.Header
    flags_header: fits.Header
    primary_header: fits.Header


def read_l2_image(path: str | Path) -> L2Image:
    """Read the IMAGE, VARIANCE, FLAGS, and ZODI extensions from one L2 cutout."""
    with fits.open(path, memmap=False) as hdul:
        image_hdu = hdul["IMAGE"]
        image = image_hdu.data.astype(float)
        variance = hdul["VARIANCE"].data.astype(float)
        flags = hdul["FLAGS"].data.astype(np.uint32)
        zodi = hdul["ZODI"].data.astype(float)
        image_header = image_hdu.header.copy()
        flags_header = hdul["FLAGS"].header.copy()
        primary_header = hdul[0].header.copy()
        with suppress_astropy_wcs_noise():
            spatial_wcs = WCS(image_header)
    return L2Image(
        image=image,
        variance=variance,
        flags=flags,
        zodi=zodi,
        spatial_wcs=spatial_wcs,
        image_header=image_header,
        flags_header=flags_header,
        primary_header=primary_header,
    )


def subtract_static_zodi(image: np.ndarray, zodi: np.ndarray) -> tuple[np.ndarray, float]:
    """Subtract the supplied SPHEREx zodiacal foreground model with scale fixed to 1."""
    return image - zodi, 1.0

