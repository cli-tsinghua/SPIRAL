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


@dataclass(frozen=True)
class TemplateFitConfig:
    """Configuration for wavelength-aware feature-amplitude fitting."""

    center_um: float
    fwhm_um: float



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


def make_weighted_window_accumulator(shape: tuple[int, int]) -> dict[str, np.ndarray]:
    return {
        "sumw": np.zeros(shape, dtype=float),
        "sumwv": np.zeros(shape, dtype=float),
        "sumwlambda": np.zeros(shape, dtype=float),
        "count": np.zeros(shape, dtype=float),
    }


def add_weighted_window_samples(
    acc: dict[str, np.ndarray],
    y: np.ndarray,
    x: np.ndarray,
    values: np.ndarray,
    lambdas: np.ndarray,
    weights: np.ndarray,
    sample_weights: np.ndarray | None = None,
) -> None:
    if sample_weights is None:
        sample_weights = np.ones_like(weights, dtype=float)
    np.add.at(acc["sumw"], (y, x), weights)
    np.add.at(acc["sumwv"], (y, x), weights * values)
    np.add.at(acc["sumwlambda"], (y, x), weights * lambdas)
    np.add.at(acc["count"], (y, x), sample_weights)


def finalize_weighted_window(acc: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sumw = acc["sumw"]
    mean = np.divide(acc["sumwv"], sumw, out=np.full(sumw.shape, np.nan), where=sumw > 0)
    wavemean = np.divide(acc["sumwlambda"], sumw, out=np.full(sumw.shape, np.nan), where=sumw > 0)
    uncertainty = np.sqrt(np.divide(1.0, sumw, out=np.full(sumw.shape, np.nan), where=sumw > 0))
    return mean, uncertainty, wavemean


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


def measure_feature_from_weighted_windows(
    cfg: FeatureConfig,
    blue_acc: dict[str, np.ndarray],
    feature_acc: dict[str, np.ndarray],
    red_acc: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    blue, blue_unc, blue_wave = finalize_weighted_window(blue_acc)
    feature, feature_unc, feature_wave = finalize_weighted_window(feature_acc)
    red, red_unc, red_wave = finalize_weighted_window(red_acc)

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


def make_anchored_template_feature_accumulator(shape: tuple[int, int]) -> dict[str, np.ndarray]:
    """Accumulator for feature-window samples used by the anchored-template method."""
    return {
        "sum": np.zeros(shape, dtype=float),
        "sum2": np.zeros(shape, dtype=float),
        "sumlambda": np.zeros(shape, dtype=float),
        "count": np.zeros(shape, dtype=np.int32),
        "sum_t_value": np.zeros(shape, dtype=float),
        "sum_t": np.zeros(shape, dtype=float),
        "sum_t_lambda": np.zeros(shape, dtype=float),
        "sum_t2": np.zeros(shape, dtype=float),
    }


def add_anchored_template_feature_samples(
    acc: dict[str, np.ndarray],
    y: np.ndarray,
    x: np.ndarray,
    values: np.ndarray,
    lambdas: np.ndarray,
    template: TemplateFitConfig,
) -> None:
    """Add feature-window samples for continuum-anchored template-amplitude fitting."""
    tvalue = normalized_gaussian_template(lambdas, template.center_um, template.fwhm_um)
    np.add.at(acc["sum"], (y, x), values)
    np.add.at(acc["sum2"], (y, x), values * values)
    np.add.at(acc["sumlambda"], (y, x), lambdas)
    np.add.at(acc["count"], (y, x), 1)
    np.add.at(acc["sum_t_value"], (y, x), tvalue * values)
    np.add.at(acc["sum_t"], (y, x), tvalue)
    np.add.at(acc["sum_t_lambda"], (y, x), tvalue * lambdas)
    np.add.at(acc["sum_t2"], (y, x), tvalue * tvalue)


def measure_feature_from_anchored_template(
    cfg: FeatureConfig,
    blue_acc: dict[str, np.ndarray],
    feature_acc: dict[str, np.ndarray],
    red_acc: dict[str, np.ndarray],
    template: TemplateFitConfig,
) -> dict[str, np.ndarray]:
    """Measure feature amplitude using blue/red continuum and actual feature wavelengths.

    The blue and red windows define a local linear continuum.  The feature-window
    samples are then converted to an integrated feature amplitude by fitting only
    the template normalization.  This avoids treating all samples inside the
    feature window as equivalent when the wavelength sampling is non-uniform.
    """
    blue, blue_unc, blue_wave = finalize_mean(blue_acc)
    feature, feature_unc, feature_wave = finalize_mean(feature_acc)
    red, red_unc, red_wave = finalize_mean(red_acc)

    slope = np.divide(red - blue, red_wave - blue_wave, out=np.full(blue.shape, np.nan), where=np.isfinite(red_wave - blue_wave) & (red_wave != blue_wave))
    intercept = blue - slope * blue_wave
    continuum = intercept + slope * cfg.center_um
    numerator = feature_acc["sum_t_value"] - intercept * feature_acc["sum_t"] - slope * feature_acc["sum_t_lambda"]
    excess = np.divide(numerator, feature_acc["sum_t2"], out=np.full(blue.shape, np.nan), where=feature_acc["sum_t2"] > 0)

    continuum_unc = np.sqrt(np.nan_to_num(blue_unc, nan=0.0) ** 2 + np.nan_to_num(red_unc, nan=0.0) ** 2) / np.sqrt(2.0)
    mean_template = np.divide(feature_acc["sum_t"], feature_acc["count"], out=np.full(blue.shape, np.nan), where=feature_acc["count"] > 0)
    excunc = np.sqrt(np.nan_to_num(feature_unc, nan=0.0) ** 2 + np.nan_to_num(continuum_unc, nan=0.0) ** 2)
    excunc = np.divide(excunc, mean_template, out=np.full(blue.shape, np.nan), where=mean_template > 0)
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


def make_template_fit_accumulator(shape: tuple[int, int]) -> dict[str, np.ndarray]:
    return {
        "m00": np.zeros(shape, dtype=float),
        "m01": np.zeros(shape, dtype=float),
        "m02": np.zeros(shape, dtype=float),
        "m11": np.zeros(shape, dtype=float),
        "m12": np.zeros(shape, dtype=float),
        "m22": np.zeros(shape, dtype=float),
        "b0": np.zeros(shape, dtype=float),
        "b1": np.zeros(shape, dtype=float),
        "b2": np.zeros(shape, dtype=float),
        "count": np.zeros(shape, dtype=float),
    }


def normalized_gaussian_template(lambdas: np.ndarray, center_um: float, fwhm_um: float) -> np.ndarray:
    sigma = float(fwhm_um) / 2.354820045
    if sigma <= 0:
        raise ValueError("fwhm_um must be positive")
    return np.exp(-0.5 * ((lambdas - center_um) / sigma) ** 2) / (np.sqrt(2.0 * np.pi) * sigma)


def add_template_fit_samples(
    acc: dict[str, np.ndarray],
    y: np.ndarray,
    x: np.ndarray,
    values: np.ndarray,
    lambdas: np.ndarray,
    weights: np.ndarray,
    template: TemplateFitConfig,
    sample_weights: np.ndarray | None = None,
) -> None:
    if sample_weights is None:
        sample_weights = np.ones_like(weights, dtype=float)
    dlambda = lambdas - template.center_um
    tvalue = normalized_gaussian_template(lambdas, template.center_um, template.fwhm_um)
    terms = {
        "m00": weights,
        "m01": weights * dlambda,
        "m02": weights * tvalue,
        "m11": weights * dlambda * dlambda,
        "m12": weights * dlambda * tvalue,
        "m22": weights * tvalue * tvalue,
        "b0": weights * values,
        "b1": weights * dlambda * values,
        "b2": weights * tvalue * values,
        "count": sample_weights,
    }
    for key, term in terms.items():
        np.add.at(acc[key], (y, x), term)


def measure_feature_from_template_fit(
    cfg: FeatureConfig,
    acc: dict[str, np.ndarray],
    template: TemplateFitConfig,
    *,
    min_count: float = 6.0,
) -> dict[str, np.ndarray]:
    shape = acc["count"].shape
    continuum = np.full(shape, np.nan, dtype=float)
    slope = np.full(shape, np.nan, dtype=float)
    excess = np.full(shape, np.nan, dtype=float)
    excunc = np.full(shape, np.nan, dtype=float)
    valid = acc["count"] >= min_count
    yy, xx = np.where(valid)
    if yy.size:
        matrices = np.empty((yy.size, 3, 3), dtype=float)
        rhs = np.empty((yy.size, 3), dtype=float)
        matrices[:, 0, 0] = acc["m00"][yy, xx]
        matrices[:, 0, 1] = acc["m01"][yy, xx]
        matrices[:, 0, 2] = acc["m02"][yy, xx]
        matrices[:, 1, 0] = acc["m01"][yy, xx]
        matrices[:, 1, 1] = acc["m11"][yy, xx]
        matrices[:, 1, 2] = acc["m12"][yy, xx]
        matrices[:, 2, 0] = acc["m02"][yy, xx]
        matrices[:, 2, 1] = acc["m12"][yy, xx]
        matrices[:, 2, 2] = acc["m22"][yy, xx]
        rhs[:, 0] = acc["b0"][yy, xx]
        rhs[:, 1] = acc["b1"][yy, xx]
        rhs[:, 2] = acc["b2"][yy, xx]
        finite = np.isfinite(matrices).all(axis=(1, 2)) & np.isfinite(rhs).all(axis=1)
        det = np.full(yy.size, np.nan, dtype=float)
        det[finite] = np.linalg.det(matrices[finite])
        solve_ok = finite & np.isfinite(det) & (np.abs(det) > 0)
        if np.any(solve_ok):
            solution = np.linalg.solve(matrices[solve_ok], rhs[solve_ok, :, None])[:, :, 0]
            cov = np.linalg.inv(matrices[solve_ok])
            y_good = yy[solve_ok]
            x_good = xx[solve_ok]
            continuum[y_good, x_good] = solution[:, 0]
            slope[y_good, x_good] = solution[:, 1]
            excess[y_good, x_good] = solution[:, 2]
            excunc[y_good, x_good] = np.sqrt(np.maximum(cov[:, 2, 2], 0.0))

    snr = np.divide(excess, excunc, out=np.full(shape, np.nan), where=excunc > 0)
    ew = np.divide(excess, continuum, out=np.full(shape, np.nan), where=continuum != 0)
    maps = {
        "CONTINUUM": continuum,
        "EXCESS": excess,
        "EXCUNC": excunc,
        "SNR": snr,
        "EW": ew,
        "SLOPE": slope,
        "LAMFEAT": np.full(shape, template.center_um, dtype=float),
        "NTOTAL": acc["count"].astype(float),
    }
    for key, value in maps.items():
        maps[key] = np.where(np.isfinite(value), value, np.nan).astype(np.float32)
    return maps
