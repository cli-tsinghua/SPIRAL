"""Feature-window definitions and continuum/excess measurements."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FeatureConfig:
    label: str
    name: str
    band: str
    detector: int
    center_um: float
    blue_min: float
    blue_max: float
    feature_min: float
    feature_max: float
    red_min: float
    red_max: float

    @property
    def width_um(self) -> float:
        return self.feature_max - self.feature_min


FEATURE_CONFIGS = {
    "pah33": FeatureConfig(
        label="pah33",
        name="PAH 3.3 micron",
        band="D4",
        detector=4,
        center_um=3.29,
        blue_min=2.95,
        blue_max=3.16,
        feature_min=3.18,
        feature_max=3.40,
        red_min=3.60,
        red_max=3.82,
    ),
    "bralpha": FeatureConfig(
        label="bralpha",
        name="Brackett-alpha",
        band="D5",
        detector=5,
        center_um=4.052,
        blue_min=3.84,
        blue_max=3.98,
        feature_min=4.00,
        feature_max=4.10,
        red_min=4.13,
        red_max=4.30,
    ),
}


def finalize_mean(acc: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = acc["count"].astype(float)
    mean = np.divide(acc["sum"], count, out=np.full(count.shape, np.nan), where=count > 0)
    mean2 = np.divide(acc["sum2"], count, out=np.full(count.shape, np.nan), where=count > 0)
    variance = np.maximum(mean2 - mean * mean, 0.0)
    stderr = np.sqrt(np.divide(variance, count, out=np.full(count.shape, np.nan), where=count > 1))
    wavemean = np.divide(acc["sumlambda"], count, out=np.full(count.shape, np.nan), where=count > 0)
    return mean, stderr, wavemean


def measure_feature_from_windows(
    cfg: FeatureConfig,
    blue_acc: dict[str, np.ndarray],
    feature_acc: dict[str, np.ndarray],
    red_acc: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    blue, blue_unc, blue_wave = finalize_mean(blue_acc)
    feature, feature_unc, feature_wave = finalize_mean(feature_acc)
    red, red_unc, red_wave = finalize_mean(red_acc)

    slope = np.divide(red - blue, red_wave - blue_wave, out=np.full(blue.shape, np.nan), where=np.isfinite(red_wave - blue_wave) & (red_wave != blue_wave))
    intercept = blue - slope * blue_wave
    continuum = intercept + slope * cfg.center_um
    continuum_feature = intercept + slope * feature_wave
    excess = (feature - continuum_feature) * cfg.width_um

    continuum_unc = np.sqrt(np.nan_to_num(blue_unc, nan=0.0) ** 2 + np.nan_to_num(red_unc, nan=0.0) ** 2) / np.sqrt(2.0)
    excunc = np.sqrt(np.nan_to_num(feature_unc, nan=0.0) ** 2 + np.nan_to_num(continuum_unc, nan=0.0) ** 2) * cfg.width_um
    excunc = np.where(excunc > 0, excunc, np.nan)
    snr = np.divide(excess, excunc, out=np.full(blue.shape, np.nan), where=excunc > 0)
    ew = np.divide(excess, continuum, out=np.full(blue.shape, np.nan), where=continuum != 0)
    ntotal = blue_acc["count"] + feature_acc["count"] + red_acc["count"]

    valid = (
        (blue_acc["count"] > 0)
        & (feature_acc["count"] > 0)
        & (red_acc["count"] > 0)
        & np.isfinite(continuum)
        & np.isfinite(excess)
    )
    maps = {
        "CONTINUUM": continuum,
        "EXCESS": excess,
        "EXCUNC": excunc,
        "SNR": snr,
        "EW": ew,
        "BLUE": blue,
        "FEATURE": feature,
        "RED": red,
        "LAMBLUE": blue_wave,
        "LAMFEAT": feature_wave,
        "LAMRED": red_wave,
        "NTOTAL": ntotal.astype(float),
        "NBLUE": blue_acc["count"].astype(float),
        "NFEAT": feature_acc["count"].astype(float),
        "NRED": red_acc["count"].astype(float),
    }
    for key, value in maps.items():
        if key in {"CONTINUUM", "EXCESS", "EXCUNC", "SNR", "EW"}:
            maps[key] = np.where(valid, value, np.nan).astype(np.float32)
        else:
            maps[key] = np.where(np.isfinite(value), value, np.nan).astype(np.float32)
    return maps

