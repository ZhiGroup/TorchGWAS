"""JAGWAS arithmetic shapes and settings, torch-free for the planners.

The projection (jagwas_projection.JagwasReduction) cuts the lower-triangular
L^-1 into row blocks; block b multiplies columns [0, e_b).
"""
from __future__ import annotations

import math
import os


def triangular_blocks(n_traits, *, minimum=512, maximum=16):
    """Row-block bounds: about n_traits/512 blocks, at most `maximum`."""
    count = max(1, min(maximum, n_traits // minimum))
    width = math.ceil(n_traits / count)
    return [(start, min(start + width, n_traits)) for start in range(0, n_traits, width)]


def gram_rows(samples, traits):
    """Sample rows per FP64 Gram block: the FP64 copy of a block (8 * rows * K bytes)
    stays within two K x K FP64 factors (rows <= 2K) once K >= 2048."""
    return max(1, min(samples, max(2 * traits, 4096)))


def rounding_target_setting():
    """Target null rms rounding error of T, from TORCHGWAS_JAGWAS_T_ROUNDING (default 0.01).

    0.01 in T is about 0.002 in -log10 p. The kept trait set is the longest
    greedy prefix whose estimated error 2 eps_z sqrt(tr R_S^-1) meets it.
    """
    value = os.environ.get('TORCHGWAS_JAGWAS_T_ROUNDING')
    target = 0.01 if value is None else float(value)
    if not (target > 0 and math.isfinite(target)):
        raise ValueError('TORCHGWAS_JAGWAS_T_ROUNDING must be a positive number')
    return target


def projection_flops_per_variant(n_traits):
    """FP64 FLOPs one variant's projection issues (2 per multiply-add)."""
    return sum(2 * (end - start) * end for start, end in triangular_blocks(n_traits))


def projection_gemm_dimensions(n_traits, markers):
    """(inner, markers, width) per projection GEMM in execution order, as tensor_service.gemm_work takes them.

    Block b writes (w_b x markers) rows of the transposed product with inner e_b.
    """
    return [(end, end - start, markers) for start, end in triangular_blocks(n_traits)]
