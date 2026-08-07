from spiral import __version__
from spiral.cloud_grid import SAMPLE_METHODS
from spiral.features import FEATURE_CONFIGS


def test_version_exists():
    assert __version__


def test_feature_configs():
    assert {"pah33", "bralpha"} <= set(FEATURE_CONFIGS)


def test_sample_methods():
    assert {
        "nearest_equal_window",
        "nearest_equal_template_window",
        "bilinear_ivar_window",
        "bilinear_ivar_template",
    } <= set(SAMPLE_METHODS)
