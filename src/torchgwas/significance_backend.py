"""Where significant pairs are selected by default: on the GPU or on the host.

Device selection copies 20 bytes per passing pair (row, trait, beta, t and df
packed as 32-bit values). Host selection copies the dense beta and t (8 bytes
per cell) and filters them on one core: 12.6 ms per 1024 x 8192 chunk on the
H100 host, against a 10.7 ms GPU chunk. Device selection transfers less while
the passing fraction -- about the p threshold under the null -- stays below
8/20. An explicit TORCHGWAS_SIGNIFICANCE_BACKEND overrides the choice, and a
panel with missing phenotypes always selects on the host.

Kept apart from reduce.device_significant_pairs, whose source the recorded
kernel census fingerprints.
"""
DEVICE_SELECTION_MAX_FRACTION = 8 / 20


def default_significance_backend(significance, n_traits):
    """'device' unless the per-test threshold passes a large share of cells."""
    return 'device' if significance.resolved_threshold(n_traits) < DEVICE_SELECTION_MAX_FRACTION else 'host'
