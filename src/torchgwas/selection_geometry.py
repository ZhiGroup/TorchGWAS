"""Shared bounded block geometry for device-selected significant pairs."""
from functools import lru_cache

# One device selection block is one CUDA nonzero call plus one packed copy.
# PyTorch's CUDA nonzero requires fewer than INT_MAX mask cells, so this is the
# largest block; by default a chunk is one block. Predicate temporaries are
# bounded separately (host_significance.PREDICATE_MAX_CELLS). The earlier 1M
# block cap paid six blocking round trips and ~45 Python dispatches per block:
# 47 syncs per 1024 x 8192 chunk, which held four variant shards on the GIL at
# 46% GPU busy (significant_host_profile_20260926).
DEVICE_SELECTION_MAX_CELLS = (1 << 31) - 2


@lru_cache(maxsize=256)
def device_selection_shape(rows, traits, max_cells):
    """Minimize blocking nonzero calls under the selection-cell cap.

    For a fixed variant-strip height, the widest allowed trait strip cannot
    increase the number of blocks. The two ceiling factors change only at
    quotient breakpoints, so this exact search skips constant intervals. Ties
    favor wider trait strips, then taller variant strips. It is independent
    of retained-pair counts and never expands the block list.
    """
    if any(type(value) is not int or value < 1
           for value in (rows, traits, max_cells)):
        raise ValueError('Positive device selection dimensions required')
    limit = min(rows, max_cells)
    height = 1
    best = None
    while height <= limit:
        row_quotient = (rows - 1) // height
        cell_quotient = max_cells // height
        row_end = ((rows - 1) // row_quotient if row_quotient else limit)
        cell_end = max_cells // cell_quotient
        candidate_height = min(limit, row_end, cell_end)
        width = min(traits, cell_quotient)
        blocks = ((rows + candidate_height - 1) // candidate_height) * (
            (traits + width - 1) // width)
        choice = (blocks, -width, -candidate_height)
        if best is None or choice < best[0]:
            best = (choice, width, candidate_height)
        height = candidate_height + 1
    return best[1], best[2], best[0][0]
