"""Device selection is the default below the payload break-even."""
from torchgwas.reduce import SignificantPairs
from torchgwas.significance_backend import DEVICE_SELECTION_MAX_FRACTION, default_significance_backend


def test_gwas_thresholds_select_on_the_device():
    assert default_significance_backend(SignificantPairs(1e-5), 8192) == 'device'
    assert default_significance_backend(SignificantPairs(), 8192) == 'device'  # alpha / K


def test_a_threshold_passing_most_cells_selects_on_the_host():
    # 20 bytes per passing pair against 8 per cell: host copies less above 0.4.
    assert DEVICE_SELECTION_MAX_FRACTION == 0.4
    assert default_significance_backend(SignificantPairs(0.5), 64) == 'host'
    assert default_significance_backend(SignificantPairs(1.0), 64) == 'host'
    assert default_significance_backend(SignificantPairs(0.39), 64) == 'device'
