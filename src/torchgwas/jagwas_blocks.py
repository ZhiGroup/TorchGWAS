"""JAGWAS arithmetic shapes and settings, torch-free for the planners.

The projection (jagwas_projection.JagwasReduction) cuts its factor into row
blocks: the default eigen factor R is upper trapezoidal (block b multiplies
columns [s_b, K)), the rounding cutoff's L^-1 lower triangular ([0, e_b)).
Either costs the same at k = K kept rows, which is what planners price: the
kept count is known only once the factor is prepared.
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


DEFAULT_RCOND = 1e-3


def rcond_setting():
    """Eigen truncation from TORCHGWAS_JAGWAS_RCOND (default 1e-3); 0 selects the rounding cutoff over traits.

    T keeps R's eigen-directions with eigenvalue above rcond x the largest
    (numpy pinv's rule), always within the rounding target as well.
    """
    value = os.environ.get('TORCHGWAS_JAGWAS_RCOND')
    if value in (None, ''):
        return DEFAULT_RCOND
    return None if float(value) == 0 else checked_rcond(float(value))


def checked_rcond(rcond):
    rcond = float(rcond)
    if not 0.0 < rcond < 1.0:
        raise ValueError('the jagwas rcond must be in (0, 1)')
    return rcond


def min_residual_setting():
    """Trait-dropping threshold from TORCHGWAS_JAGWAS_MIN_RESIDUAL; unset, none.

    Set, a trait is kept only while at least this fraction of its variance is
    not explained by the traits kept before it in the greedy pivoted order
    (its VIF given them at most 1 / min_residual).
    """
    value = os.environ.get('TORCHGWAS_JAGWAS_MIN_RESIDUAL')
    return None if value in (None, '') else checked_min_residual(float(value))


def cutoff_settings():
    """(rcond, min_residual) from the environment: one of them, or neither (the rounding cutoff alone).

    TORCHGWAS_JAGWAS_MIN_RESIDUAL alone selects trait dropping; otherwise
    TORCHGWAS_JAGWAS_RCOND decides (default 1e-3; 0 for the rounding cutoff,
    which then takes TORCHGWAS_JAGWAS_MIN_RESIDUAL if that is set too).
    """
    min_residual = min_residual_setting()
    if min_residual is not None and os.environ.get('TORCHGWAS_JAGWAS_RCOND') in (None, ''):
        return None, min_residual
    rcond = rcond_setting()
    if rcond is not None and min_residual is not None:
        raise ValueError('TORCHGWAS_JAGWAS_RCOND and TORCHGWAS_JAGWAS_MIN_RESIDUAL are exclusive')
    return rcond, min_residual


def checked_min_residual(fraction):
    fraction = float(fraction)
    if not 0.0 < fraction < 1.0:
        raise ValueError('the jagwas min_residual must be in (0, 1)')
    return fraction


def projection_flops_per_variant(n_traits):
    """FP64 FLOPs one variant's projection issues at k = K kept rows (2 per multiply-add)."""
    return sum(2 * (end - start) * (n_traits - start) for start, end in triangular_blocks(n_traits))


def projection_gemm_dimensions(n_traits, markers, *, upper=True):
    """(inner, markers, width) per projection GEMM in execution order, as tensor_service.gemm_work takes them.

    Block b writes (w_b x markers) rows of the transposed product: with the
    default eigen factor at k = K (upper) its inner size is K - s_b, with the
    rounding cutoff's L^-1 (upper=False) it is e_b.
    """
    return [((n_traits - start) if upper else end, end - start, markers)
            for start, end in triangular_blocks(n_traits)]
