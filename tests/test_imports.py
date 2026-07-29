from spiral import __version__
from spiral.features import FEATURE_CONFIGS


def test_version_exists():
    assert __version__


def test_feature_configs():
    assert {"pah33", "bralpha"} <= set(FEATURE_CONFIGS)

