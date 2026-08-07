import numpy as np

from spiral.features import (
    FEATURE_CONFIGS,
    TemplateFitConfig,
    add_anchored_template_feature_samples,
    add_template_fit_samples,
    make_anchored_template_feature_accumulator,
    make_template_fit_accumulator,
    measure_feature_from_anchored_template,
    measure_feature_from_template_fit,
    measure_feature_from_windows,
    normalized_gaussian_template,
)


def _mean_accumulator(lambdas: np.ndarray, values: np.ndarray) -> dict[str, np.ndarray]:
    shape = (1, 1)
    return {
        "sum": np.array([[np.sum(values)]], dtype=float),
        "sum2": np.array([[np.sum(values * values)]], dtype=float),
        "sumlambda": np.array([[np.sum(lambdas)]], dtype=float),
        "count": np.array([[values.size]], dtype=np.int32),
    }


def _linear_continuum(
    lambdas: np.ndarray,
    center_um: float,
    continuum: float,
    slope: float,
) -> np.ndarray:
    return continuum + slope * (lambdas - center_um)


def test_anchored_template_recovers_integrated_excess_from_actual_wavelengths():
    cfg = FEATURE_CONFIGS["pah33"]
    template = TemplateFitConfig(center_um=cfg.center_um, fwhm_um=0.110)
    continuum = 5.0
    slope = 0.7
    true_excess = 2.5

    blue_lambda = np.array([3.00, 3.08, 3.14])
    red_lambda = np.array([3.63, 3.72, 3.80])
    feature_lambda = np.array([3.255, 3.270, 3.282, 3.294, 3.306, 3.318])
    blue_values = _linear_continuum(blue_lambda, cfg.center_um, continuum, slope)
    red_values = _linear_continuum(red_lambda, cfg.center_um, continuum, slope)
    template_values = normalized_gaussian_template(
        feature_lambda,
        template.center_um,
        template.fwhm_um,
    )
    feature_values = (
        _linear_continuum(feature_lambda, cfg.center_um, continuum, slope)
        + true_excess * template_values
    )

    feature_acc = make_anchored_template_feature_accumulator((1, 1))
    zeros = np.zeros(feature_lambda.size, dtype=np.int64)
    add_anchored_template_feature_samples(
        feature_acc,
        zeros,
        zeros,
        feature_values,
        feature_lambda,
        template,
    )

    maps = measure_feature_from_anchored_template(
        cfg,
        _mean_accumulator(blue_lambda, blue_values),
        feature_acc,
        _mean_accumulator(red_lambda, red_values),
        template,
    )
    window_maps = measure_feature_from_windows(
        cfg,
        _mean_accumulator(blue_lambda, blue_values),
        _mean_accumulator(feature_lambda, feature_values),
        _mean_accumulator(red_lambda, red_values),
    )

    np.testing.assert_allclose(maps["CONTINUUM"][0, 0], continuum, rtol=0, atol=1e-6)
    np.testing.assert_allclose(maps["EXCESS"][0, 0], true_excess, rtol=0, atol=1e-6)
    assert abs(float(window_maps["EXCESS"][0, 0]) - true_excess) > 0.5
    assert maps["NFEAT"][0, 0] == feature_lambda.size


def test_template_fit_recovers_weighted_linear_template_model():
    cfg = FEATURE_CONFIGS["pah33"]
    template = TemplateFitConfig(center_um=cfg.center_um, fwhm_um=0.110)
    continuum = 4.0
    slope = -0.4
    true_excess = 1.6
    lambdas = np.array([
        3.00,
        3.08,
        3.14,
        3.22,
        3.27,
        3.29,
        3.31,
        3.36,
        3.63,
        3.72,
        3.80,
    ])
    values = (
        _linear_continuum(lambdas, cfg.center_um, continuum, slope)
        + true_excess
        * normalized_gaussian_template(lambdas, template.center_um, template.fwhm_um)
    )
    weights = np.linspace(0.5, 1.5, lambdas.size)
    acc = make_template_fit_accumulator((1, 1))
    zeros = np.zeros(lambdas.size, dtype=np.int64)

    add_template_fit_samples(
        acc,
        zeros,
        zeros,
        values,
        lambdas,
        weights,
        template,
    )
    maps = measure_feature_from_template_fit(cfg, acc, template, min_count=6)

    np.testing.assert_allclose(maps["CONTINUUM"][0, 0], continuum, rtol=0, atol=1e-5)
    np.testing.assert_allclose(maps["SLOPE"][0, 0], slope, rtol=0, atol=1e-5)
    np.testing.assert_allclose(maps["EXCESS"][0, 0], true_excess, rtol=0, atol=1e-5)
    assert maps["NTOTAL"][0, 0] == lambdas.size
