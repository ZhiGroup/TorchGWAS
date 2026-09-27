"""Version-bound work branches for flatnonzero on a contiguous boolean mask.

NumPy 2.2.6 item_selection.c, PyArray_Nonzero: an initial count is followed
by no extraction for empty output, memchr extraction at <=10% occupancy,
or branchless extraction above that threshold. Branch prices are separate;
the dense cost is not the sparse cost. Retained counts are explicit inputs.
"""
import numpy as np


def nonzero_protocol():
    if np.__version__ != '2.2.6':
        raise ValueError('Unqualified NumPy boolean nonzero implementation: '+np.__version__)
    return dict(algorithm='contiguous_bool_1d_density_v1', numpy_version=np.__version__)


def nonzero_regime(cells, retained):
    if type(cells) is not int or type(retained) is not int or not 0 <= retained <= cells:
        raise ValueError('Integer mask cells and bounded retained count required')
    if retained == 0: return 'empty'
    return 'sparse' if 10*retained <= cells else 'dense'


def validate_host_price_protocol(prices):
    from .host_significance import host_selector, PREDICATE_MAX_CELLS
    if (not isinstance(prices, dict) or prices.get('host_selector') != host_selector()
            or prices.get('predicate_max_cells') != PREDICATE_MAX_CELLS):
        raise ValueError('Independent prices must match the host selector and predicate limit')
    if prices.get('nonzero_protocol') != nonzero_protocol():
        raise ValueError('Independent prices must match the NumPy boolean nonzero protocol')
