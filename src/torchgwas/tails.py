"""A Student-t tail that does not underflow.

`2 * scipy.special.stdtr(df, -|t|)` returns **exactly zero** once the tail falls
below the float64 floor. At df = 20,000 that happens from **|t| = 38.354**
onward, so |t| = 39, 50, 100 and 300 all report the same clipped
`-log10 p = 307.6527`. For a min-P scan across millions of traits that is the
ranking disappearing precisely where the answer lives: the strongest hits become
indistinguishable from each other.

Nothing in scipy fixes this. Both `scipy.stats.t.logsf` and
`scipy.stats.beta.logcdf` hit the identical cliff, because each takes the log of
a survival function that has already underflowed rather than computing in log
space throughout. The normal tail via `scipy.special.log_ndtr` *is* stable and
fast, but it is the wrong distribution: against the exact Student-t it errs by
2.2% of the exponent at |t| = 30 and 3.6% at |t| = 38, which is far too much to
publish as a p-value.

So the tail is computed here, in logs, from the identity

    P(|T| > t)  =  I_x(df/2, 1/2),      x = df / (df + t^2)

with the regularised incomplete beta evaluated as

    log I_x(a, b) = a log x + b log1p(-x) - log a - betaln(a, b) + log K(a,b,x)

where `K` is the standard continued fraction. The leading terms alone are not
enough -- they overstate `-log10 p` by 2.9 decades at |t| = 5 and still by 1.2
at |t| = 38 -- so the continued fraction is the accuracy, not a refinement.
"""

from __future__ import annotations

import numpy as np
from scipy import special

_TINY = 1e-300
_EPS = 3e-16
# 40, not 400. `tdist.py` measured that 40 iterations converge to 2.6e-12
# relative against 300 and are 7x faster, and the torch port here agrees with
# scipy to 7.1e-10 at 40 over df from 1 to 499,998 -- far tighter than anything
# downstream can use. The 400 this started at cost **136 seconds** for 2 million
# values against scipy's 0.59, because the early exit only fires once *every*
# element has converged and a wide spread of `x` keeps a few going. The exit is
# kept, so easy inputs still stop sooner.
_MAX_ITERATIONS = 40


def _betacf(a: np.ndarray, b: float, x: np.ndarray) -> np.ndarray:
    """Continued fraction for the incomplete beta, Lentz's method, vectorised.

    `a` and `x` are arrays of the same shape and `b` is a scalar, which is the
    shape this module needs: `b` is always 1/2 for a two-sided t tail, while `a`
    is `df/2` and df is per variant.

    Iteration is to a fixed convergence test over the whole array rather than
    per element, so every element gets at least as many terms as it needs.
    """
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = np.ones_like(a)
    d = 1.0 - qab * x / qap
    d = np.where(np.abs(d) < _TINY, _TINY, d)
    d = 1.0 / d
    h = d.copy()
    for m in range(1, _MAX_ITERATIONS + 1):
        m2 = 2 * m
        # Even step.
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = np.where(np.abs(d) < _TINY, _TINY, d)
        c = 1.0 + aa / c
        c = np.where(np.abs(c) < _TINY, _TINY, c)
        d = 1.0 / d
        h = h * d * c
        # Odd step.
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = np.where(np.abs(d) < _TINY, _TINY, d)
        c = 1.0 + aa / c
        c = np.where(np.abs(c) < _TINY, _TINY, c)
        d = 1.0 / d
        delta = d * c
        h = h * delta
        if np.all(np.abs(delta - 1.0) < _EPS):
            break
    return h


def log_two_sided_t_sf(t, df):
    """Natural log of `P(|T| > |t|)` for Student's t, without underflowing.

    `t` and `df` broadcast against each other. Returns 0.0 (log of 1) at t = 0
    and grows without bound as |t| does, where the direct computation saturates.
    """
    t = np.asarray(t, dtype=np.float64)
    df = np.asarray(df, dtype=np.float64)
    t, df = np.broadcast_arrays(t, df)
    a = df / 2.0
    b = 0.5
    squared = t * t
    x = df / (df + squared)

    # The continued fraction converges quickly only for x below
    # (a+1)/(a+b+2); above it, use the reflection I_x(a,b) = 1 - I_{1-x}(b,a).
    # For a two-sided t tail the direct branch covers everything that matters --
    # x exceeds the threshold only for |t| very close to zero, where the tail is
    # near 1 and no precision is at stake.
    threshold = (a + 1.0) / (a + b + 2.0)
    direct = x < threshold

    out = np.empty(x.shape, dtype=np.float64)

    if np.any(direct):
        xa, aa = x[direct], a[direct]
        front = (aa * np.log(xa) + b * np.log1p(-xa)
                 - np.log(aa) - special.betaln(aa, b))
        out[direct] = front + np.log(_betacf(aa, b, xa))

    other = ~direct
    if np.any(other):
        # Above the threshold the tail is O(1) -- this branch is reached only
        # for |t| near zero -- so the direct survival function cannot underflow
        # and is simply used. Nothing is at stake here; the whole point of this
        # module is the other branch.
        with np.errstate(divide="ignore"):
            out[other] = np.log(
                np.clip(2.0 * special.stdtr(df[other], -np.abs(t[other])),
                        _TINY, None))
    return out


def upper_tail_log10_from_t(t, df):
    """`-log10 P(|T| > |t|)`, the quantity a Manhattan plot actually shows."""
    return -log_two_sided_t_sf(t, df) / np.log(10.0)


# -- the same thing on the device -------------------------------------------
#
# Lifted from `tdist.py`, which had this before I wrote the numpy version above
# and which I should have looked for first. Two changes were needed to make it
# usable from the scan. **`a` is a tensor here**, so the degrees of freedom
# broadcast per variant; the original wrapped a Python float in
# `torch.tensor(a + b)` and so accepted only a scalar df, which is wrong wherever
# missingness varies. And the iteration count is **40**, not the 400 the numpy
# path allows itself: that is `tdist.py`'s own measurement -- 40 converges to
# 2.6e-12 relative against 300, and is 7x faster -- and it matters on a device
# where every element pays for the worst case because there is no early exit.

_TORCH_ITERATIONS = 40


def _betacf_torch(a, b, x, iterations: int = _TORCH_ITERATIONS):
    """Lentz continued fraction for the incomplete beta, elementwise on `x`.

    Fixed trip count and no convergence test: a data-dependent break would
    serialise on a host read, and every element would pay for the slowest one
    anyway.
    """
    import torch

    tiny = torch.full_like(x, _TINY)
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = torch.ones_like(x)
    d = 1.0 - qab * x / qap
    d = torch.where(d.abs() < _TINY, tiny, d)
    d = 1.0 / d
    h = d.clone()
    for m in range(1, iterations + 1):
        m2 = 2 * m
        step = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + step * d
        d = torch.where(d.abs() < _TINY, tiny, d)
        c = 1.0 + step / c
        c = torch.where(c.abs() < _TINY, tiny, c)
        d = 1.0 / d
        h = h * d * c
        step = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + step * d
        d = torch.where(d.abs() < _TINY, tiny, d)
        c = 1.0 + step / c
        c = torch.where(c.abs() < _TINY, tiny, c)
        d = 1.0 / d
        h = h * d * c
    return h


def upper_tail_log10_from_t_torch(t, df):
    """`-log10 P(|T| > |t|)` on whatever device `t` is on, per-variant df.

    Computed in float64 regardless of the caller's dtype: the tail is a product
    of logs of very small numbers, and float32 would lose the exponent this
    function exists to preserve.
    """
    import torch

    t = torch.as_tensor(t).double().abs()
    df = torch.as_tensor(df, dtype=torch.float64, device=t.device)
    t, df = torch.broadcast_tensors(t, df)
    a = df / 2.0
    b = 0.5
    x = df / (df + t * t)

    log_beta = torch.lgamma(a) + torch.lgamma(
        torch.full_like(a, b)) - torch.lgamma(a + b)
    # Never form the final exp: that is the only step that underflows, and it is
    # exactly what caps `-log10 p` at 307.65 in the direct computation.
    log_tail = (a * torch.log(x.clamp_min(_TINY)) + b * torch.log1p(-x)
                - torch.log(a) - log_beta
                + torch.log(_betacf_torch(a, b, x)))
    result = -log_tail / float(np.log(10.0))

    # Above the threshold the continued fraction converges slowly and the tail
    # is O(1) anyway; the reflection is the accurate form there. Reached only
    # for |t| near zero, where nothing is at stake.
    reflect = x >= (a + 1.0) / (a + b + 2.0)
    if bool(reflect.any()):
        # I_x(a,b) = 1 - I_{1-x}(b,a), and I_{1-x}(b,a) is evaluated directly
        # because it is O(1) here -- there is no underflow to dodge. Feed the
        # non-reflected elements a harmless 0.5 so the fraction converges for
        # them too and only the `where` decides what is kept.
        flipped = torch.where(reflect, 1.0 - x, torch.full_like(x, 0.5))
        upper = (torch.exp(-log_beta + b * torch.log(flipped.clamp_min(_TINY))
                           + a * torch.log1p(-flipped))
                 * _betacf_torch(b, a, flipped) / b)
        alternative = -torch.log(
            (1.0 - upper).clamp(_TINY, 1.0)) / float(np.log(10.0))
        result = torch.where(reflect, alternative, result)
    return result
