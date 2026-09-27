"""Bounded host predicates and owned significant-pair output arrays.

The threshold remains inclusive. For FP32 statistics only, an FP64 cutoff is
rounded upward to the next representable FP32 value, preserving the original
FP64 comparison exactly. FP64 statistics retain their FP64 cutoffs.
"""
import os
import numpy as np

HOST_SELECTOR = 'bounded_flat_v2'
PREDICATE_MAX_CELLS = 1 << 20
NATIVE_HOST_SELECTOR = 'native_row_flat_v2'


def host_selector():
    """Resolve the explicit implementation choice; unknown values fail closed."""
    choice = os.getenv('TORCHGWAS_HOST_PREDICATE', 'numpy')
    if choice == 'numpy': return HOST_SELECTOR
    if choice == 'native': return NATIVE_HOST_SELECTOR
    raise ValueError('TORCHGWAS_HOST_PREDICATE must be numpy or native')


def fill_predicate_mask(values, limits, keep):
    """Shared production/probe predicate; limits are already broadcast/rounded."""
    if host_selector() == NATIVE_HOST_SELECTOR:
        from .native_host_predicate import fill_mask
        if fill_mask(values, limits, keep): return
    height, width, _ = predicate_block_shape(*values.shape)
    for row in range(0, values.shape[0], height):
        for column in range(0, values.shape[1], max(1, width)):
            sl = np.s_[row:row+height, column:column+width]
            absolute = np.abs(values[sl])
            np.greater_equal(absolute, limits[sl], out=keep[sl])
            keep[sl] &= np.isfinite(absolute)


def predicate_block_shape(markers, traits, max_cells=PREDICATE_MAX_CELLS):
    """Rectangle geometry shared by the executor and its analytical ledger."""
    if isinstance(max_cells, bool) or not isinstance(max_cells, int) or max_cells < 1:
        raise ValueError('Positive integer predicate cell limit required')
    width = min(traits, max_cells)
    height = max(1, max_cells // max(1, width))
    calls = ((markers + height - 1)//height) * ((traits + max(1, width) - 1)//max(1, width))
    return height, width, calls


def ceil_float32(critical):
    """Small cutoff-array conversion; never round a cutoff downward."""
    critical = np.asarray(critical)
    with np.errstate(over='ignore', invalid='ignore'):
        rounded = critical.astype(np.float32)
        np.nextafter(rounded, np.float32(np.inf), out=rounded,
                     where=rounded.astype(np.float64) < critical)
    return rounded


def select_host_pairs(beta, values, row_df, critical, start=0):
    """One owned result per input chunk, in variant-major coordinate order.

    The mask retains chunk shape, but abs/finite temporaries are bounded to
    PREDICATE_MAX_CELLS. Contiguous payloads use flat indices before those
    indices become row coordinates. Row-broadcast df needs only row indices.
    General broadcast df and noncontiguous inputs keep indexed gathers without
    expanding or copying the original dense statistics.
    """
    if values.dtype == np.float32:
        critical = ceil_float32(critical)
    keep = np.empty(values.shape, dtype=bool, order='C')
    limits = np.broadcast_to(critical, values.shape)
    fill_predicate_mask(values, limits, keep)
    rows = np.flatnonzero(keep).astype(np.int64, copy=False)
    flat_t = values.flags.c_contiguous
    flat_beta = beta is not None and beta.shape == values.shape and beta.flags.c_contiguous
    # These reshapes are views: never flatten a noncontiguous dense input by
    # copying it. Gather before divmod overwrites the flat index buffer.
    selected_beta = beta.reshape(-1)[rows] if flat_beta else None
    selected_t = values.reshape(-1)[rows] if flat_t else None
    columns = np.empty_like(rows)
    if values.shape[1]:
        np.divmod(rows, values.shape[1], out=(rows, columns))
    df = np.broadcast_to(row_df, values.shape)
    if not values.shape[1]:
        selected_df = np.empty(rows.shape, dtype=df.dtype)
    elif df.strides[1] == 0:
        selected_df = df[:, 0][rows]
    else:
        selected_df = df[rows, columns]
    if beta is not None and not flat_beta:
        selected_beta = beta[rows, columns]
    if not flat_t:
        selected_t = values[rows, columns]
    rows += start
    return rows, columns, selected_beta, selected_t, selected_df
