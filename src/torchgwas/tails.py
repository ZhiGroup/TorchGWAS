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


# -- the scan's device tail ---------------------------------------------------
#
# upper_tail_log10_from_t_torch is the reference, and it is slow in a scan: 40
# iterations of eager FP64 elementwise work, each materialising a chunk-sized
# temporary, and both branches for every cell. H100, one 8.4M-cell chunk
# (1024 x 8192, per-variant df 22,238): 181 ms, against ~10.7 ms for the
# chunk's GEMM (benchmarks/logp_device_cost_20260927.py). The same identity and
# fraction are rearranged here so each cell evaluates ONE fraction -- the
# direct K(a, 1/2, x) or the reflected K(1/2, a, 1 - x), parameters swapped per
# cell -- and the iterations run as a compiled block of four reused ten times:
# 4.1 ms per chunk, 4.7 s to compile in a fresh process with a warm inductor
# cache (benchmarks/logp_blocked_probe_20260927.py). Against scipy's stdtr the
# result is within 1.3e-11 absolute at df 22,238 (2e-13 at 500, 1e-14 at 30),
# and within 4.4e-7 of the reference over a scan chunk, below float32 storage
# resolution. On the CPU, or with TORCHGWAS_COMPILE_TAILS=0 or a failed
# compile, the same stages run eagerly.

DEVICE_TAIL_BLOCK = 4
# Constants, so the exported graphs call neither scipy nor numpy.
_LGAMMA_HALF = float(special.gammaln(0.5))
_LN10 = float(np.log(10.0))
# Row strips bound the FP64 transient to ~380 MB: 747 MB was measured for one
# unstripped 8.4M-cell chunk (89 bytes per cell; the planner charges 96). Each
# strip is 12 compiled calls, so a 1024 x 8192 chunk is two strips.
DEVICE_TAIL_MAX_CELLS = 1 << 22
DEVICE_TAIL_TRANSIENT_BYTES_PER_CELL = 96

_STAGES = {}
_KINDS = {}  # device -> 'aoti' (ahead-of-time build), 'jit' (torch.compile) or 'eager'
_STAGE_LOCK = None
# The AOT build's second form (_tail_whole), and each form's measured
# (host seconds per call, GPU seconds per cell) on that device (_form_costs).
_WHOLE = {}
_FORM_COSTS = {}


def _guard(value):
    import torch

    return torch.where(value.abs() < _TINY, _TINY, value)


def _tail_prologue(t, df):
    import torch

    t = t.double().abs()
    df = df.double()
    a = df * 0.5
    squared = t * t
    total = df + squared
    x = df / total
    y = squared / total  # 1 - x, without the cancellation
    reflect = x >= (a + 1.0) / (a + 2.5)
    first = torch.where(reflect, 0.5, a)
    second = torch.where(reflect, a, 0.5)
    argument = torch.where(reflect, y, x)
    d = 1.0 / _guard(1.0 - (first + second) * argument / (first + 1.0))
    return first, second, argument, torch.ones_like(argument), d, d.clone(), x, y, reflect


def _tail_block(first, second, argument, c, d, h, start):
    qab, qap, qam = first + second, first + 1.0, first - 1.0
    for offset in range(DEVICE_TAIL_BLOCK):
        m = start + offset
        m2 = 2.0 * m
        step = m * (second - m) * argument / ((qam + m2) * (first + m2))
        d = 1.0 / _guard(1.0 + step * d)
        c = _guard(1.0 + step / c)
        h = h * d * c
        step = -(first + m) * (qab + m) * argument / ((first + m2) * (qap + m2))
        d = 1.0 / _guard(1.0 + step * d)
        c = _guard(1.0 + step / c)
        h = h * d * c
    return c, d, h


def _tail_epilogue(df, x, y, reflect, h):
    import torch

    a = df.double() * 0.5
    log_beta = torch.lgamma(a) + _LGAMMA_HALF - torch.lgamma(a + 0.5)
    log_direct = (a * torch.log(x.clamp_min(_TINY)) + 0.5 * torch.log(y.clamp_min(_TINY))
                  - torch.log(a) - log_beta + torch.log(h))
    upper = 2.0 * torch.exp(-log_beta + 0.5 * torch.log(y.clamp_min(_TINY)) + a * torch.log1p(-y)) * h
    log_reflect = torch.log((1.0 - upper).clamp(_TINY, 1.0))
    return -torch.where(reflect, log_reflect, log_direct) / _LN10


def _tail_whole(t, df):
    """Prologue, every block and the epilogue in one graph (the AOT build's second form).

    One host call per strip instead of 12, but slower on the GPU: Inductor
    fuses the 40 iterations with recomputation, 0.82-0.85 ns per cell
    against 0.36-0.42 for the blocks, which store (c, d, h) every four
    iterations (H100; the best fusion setting tried, realize_opcount_threshold
    8, reached 0.56; benchmarks/tail_graph_forms_20260927.py and
    tail_whole_build_options_20260927.py). Host time matters because a call
    holds the GIL while it launches: with the blocks, a (4096, 1) tail took
    0.59 ms on one GPU and 2.89 ms per call on each of four
    (tail_thread_scaling_20260927.py). _choose_form weighs the two.
    """
    first, second, argument, c, d, h, x, y, reflect = _tail_prologue(t, df)
    for start in range(1, _TORCH_ITERATIONS + 1, DEVICE_TAIL_BLOCK):
        c, d, h = _tail_block(first, second, argument, c, d, h, float(start))
    return _tail_epilogue(df, x, y, reflect, h)


def _eager_stages():
    return (_tail_prologue, _tail_block, _tail_epilogue)


def _starts(device):
    import torch

    return [torch.tensor(float(m), dtype=torch.float64, device=device)
            for m in range(1, _TORCH_ITERATIONS + 1, DEVICE_TAIL_BLOCK)]


def _evaluate(stages, t, df, starts):
    if callable(stages):  # the whole tail as one graph (_tail_whole)
        return stages(t, df)
    prologue, block, epilogue = stages
    first, second, argument, c, d, h, x, y, reflect = prologue(t, df)
    for start in starts:
        c, d, h = block(first, second, argument, c, d, h, start)
    return epilogue(df, x, y, reflect, h)


_FORM_PROBES = ((2048, 2), (512, 8192))


def _form_costs(device, blocks, whole, starts):
    """Each form's measured costs on `device`: host seconds at the two probe sizes, GPU seconds per cell.

    Host time is the time to issue one call (no synchronization), at a
    (2048, 2) and a (512, 8192) strip: the block form's grows with the strip
    (H100: 0.47, 0.82 and 1.41 ms at 4K, 2M and 4M cells), the whole graph's
    stays ~0.1 ms. GPU time is CUDA-event time on the larger strip. About
    40 ms per device, once.
    """
    import time

    import torch

    costs = {}
    for name, stages in (('blocks', blocks), ('whole', whole)):
        host, gpu = [], None
        for shape in _FORM_PROBES:
            t = torch.full(shape, 3.0, dtype=torch.float64, device=device)
            # Any df: the tail runs a fixed number of iterations, so its cost
            # does not depend on it.
            df = torch.full(shape, 1000.0, dtype=torch.float64, device=device)
            _evaluate(stages, t, df, starts)
            torch.cuda.synchronize(device)
            issued, elapsed = [], []
            for _ in range(3):
                start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                started = time.perf_counter()
                start.record()
                _evaluate(stages, t, df, starts)
                stop.record()
                issued.append(time.perf_counter() - started)
                stop.synchronize()
                elapsed.append(start.elapsed_time(stop) / 1e3)
            host.append(min(issued))
            gpu = min(elapsed) / (shape[0] * shape[1])
        costs[name] = (tuple(host), gpu)
    return costs


def _choose_form(device, cells):
    """'whole' or 'blocks' for a strip of `cells` on `device`, by measured cost.

    A strip holds the GIL for its host time, and every CUDA device with a
    tail in this process (a shard thread each) needs the same GIL, while its
    GPU time is its own device's. The chosen form minimizes host seconds
    (interpolated in cells between the two probe sizes) times those devices,
    plus GPU seconds per cell times `cells`: min-p's winner column is
    host-bound (whole); a dense strip takes the blocks on one GPU and the
    whole graph once several shards share the GIL.
    """
    costs = _FORM_COSTS.get(device)
    if costs is None:
        return 'blocks'
    sharing = max(1, sum(1 for key in _STAGES if getattr(key, 'type', None) == 'cuda'))
    small, large = (rows * columns for rows, columns in _FORM_PROBES)

    def host_seconds(points):
        low, high = points
        return max(low, low + (high - low) * (cells - small) / (large - small))

    price = {name: host_seconds(host) * sharing + per_cell * cells for name, (host, per_cell) in costs.items()}
    return min(price, key=price.get)


def prepare_device_tail(device):
    """Compile (or load) the device tail for `device` once; later calls are free.

    Thread-safe. Returns True when the compiled stages serve this device and
    False when the eager stages will.
    """
    import os
    import threading

    import torch

    global _STAGE_LOCK
    device = torch.device(device)
    if device.type == 'cuda' and device.index is None:
        device = torch.device('cuda', torch.cuda.current_device())
    if _STAGE_LOCK is None:
        _STAGE_LOCK = threading.Lock()
    with _STAGE_LOCK:
        if device in _STAGES:
            return _STAGES[device][0] is not None
        compiled = None
        kind = 'eager'
        # TORCHGWAS_COMPILE_TAILS: 0 eager only, jit skips the ahead-of-time build.
        mode = os.environ.get('TORCHGWAS_COMPILE_TAILS', '1')
        if device.type == 'cuda' and mode != '0':
            whole = None
            try:
                loaded = _load_aoti(device) if mode != 'jit' else None
                kind = 'aoti' if loaded is not None else 'jit'
                compiled, whole = loaded or (tuple(torch.compile(stage, dynamic=True) for stage in _eager_stages()),
                                             None)
                starts = _starts(device)
                with torch.cuda.device(device):
                    probe = torch.linspace(0.0, 40.0, 64 * 3, device=device, dtype=torch.float64).reshape(64, 3)
                    probe_df = torch.full((64, 3), 30.0, device=device, dtype=torch.float64)
                    slow = _evaluate(_eager_stages(), probe, probe_df, starts)
                    fast = _evaluate(compiled, probe, probe_df, starts)
                    if not bool(torch.allclose(fast, slow, rtol=1e-10, atol=1e-10)):
                        compiled, whole, kind = None, None, 'eager'
                    if whole is not None and not bool(torch.allclose(_evaluate(whole, probe, probe_df, starts), slow,
                                                                     rtol=1e-10, atol=1e-10)):
                        whole = None
                    if whole is not None:
                        _FORM_COSTS[device] = _form_costs(device, compiled, whole, starts)
            except Exception:  # noqa: BLE001 - no compiler, no Triton: the eager stages are exact too
                compiled, whole, kind = None, None, 'eager'
            if whole is not None:
                _WHOLE[device] = whole
        _STAGES[device] = (compiled, _starts(device))
        _KINDS[device] = kind
        return compiled is not None


def neg_log10_p_device(t, df, *, out=None, max_cells=DEVICE_TAIL_MAX_CELLS):
    """`-log10 P(|T| > |t|)` for a scan chunk on t's device, into float32 `out`.

    `t` is (rows, traits); `df` broadcasts against it, one per variant (rows,
    1), per pair (rows, traits) or per trait (1, traits). NaN t gives NaN.
    Values are computed in FP64 and stored at `out`'s dtype (float32 unless
    an `out` says otherwise).
    """
    import torch

    rows, traits = t.shape
    device = t.device
    if device.type == 'cuda' and device.index is None:
        device = torch.device('cuda', torch.cuda.current_device())
    if device not in _STAGES:
        prepare_device_tail(device)
    compiled, starts = _STAGES[device]
    stages = compiled or _eager_stages()
    df = torch.as_tensor(df, device=t.device)
    if df.dim() < 2:
        df = df.reshape(1, -1) if df.dim() == 1 and df.shape[0] == traits and traits != rows else df.reshape(-1, 1)
    result = torch.empty((rows, traits), dtype=torch.float32, device=t.device) if out is None else out
    step = max(1, int(max_cells) // max(1, traits))

    aot = _KINDS.get(device) == 'aoti'
    whole = _WHOLE.get(device)

    def form(cells):
        return whole if whole is not None and _choose_form(device, cells) == 'whole' else stages

    def run():
        for first in range(0, rows, step):
            last = min(rows, first + step)
            strip_t = t[first:last]
            strip_df = df[first:last] if df.shape[0] == rows else df
            if compiled is not None and min(strip_t.shape) < 2:
                # The tail is elementwise, so a column (a reduction's winners)
                # runs as a (cells / 2, 2) block, padded by one cell when odd:
                # the build takes at least 2 x 2, and a new shape would
                # recompile the jit stages.
                cells = strip_t.numel()
                if cells < 4:
                    result[first:last] = _evaluate(_eager_stages(), strip_t, strip_df, starts)
                    continue
                # Contiguous copies: a broadcast df reshapes to a zero-stride
                # view, which the build would read past (NaN).
                flat_t = strip_t.double().reshape(-1).contiguous()
                flat_df = strip_df.double().expand(strip_t.shape).reshape(-1).contiguous()
                if cells % 2:
                    flat_t, flat_df = torch.cat((flat_t, flat_t[-1:])), torch.cat((flat_df, flat_df[-1:]))
                values = _evaluate(form(cells), flat_t.view(-1, 2), flat_df.view(-1, 2), starts)
                result[first:last] = values.reshape(-1)[:cells].view(strip_t.shape)
                continue
            if aot:
                # The build takes FP64 (rows, traits) of at least 2 x 2.
                strip_t = strip_t.double().contiguous()
                strip_df = strip_df.double().expand(strip_t.shape).contiguous()
            result[first:last] = _evaluate(form(strip_t.numel()), strip_t, strip_df, starts)

    if compiled is None:
        run()
    elif aot:
        # Not dynamo functions, so no lock; but they launch on the current
        # device's current stream, which must be t's device.
        with torch.cuda.device(t.device):
            run()
    else:
        # PyTorch's FX-tracing flag is process-wide: a compiled call made
        # while another thread compiles (a second GPU's prepare) is refused.
        # Calls and compiles therefore share the lock; a call only queues
        # kernels, so shards hold it briefly.
        with _STAGE_LOCK:
            run()
    return result


def prepare_device_tail_async(devices):
    """Prepare the device tail for each CUDA device on a daemon thread.

    The first GPU costs ~5-7 s (compiler import, trace, cache load), each
    further one ~2 s. Started at API entry it overlaps input loading and
    preparation; a scan that needs a device first waits on the same lock.
    """
    import threading

    import torch

    targets = [torch.device(d) for d in devices if torch.device(d).type == 'cuda']
    if not targets or not torch.cuda.is_available():
        return None

    def run():
        for device in targets:
            try:
                prepare_device_tail(device)
            except Exception:  # noqa: BLE001 - the scan prepares (or falls back) itself
                pass

    worker = threading.Thread(target=run, name='torchgwas-tail-compile', daemon=True)
    worker.start()
    return worker


# -- ahead-of-time build (AOTInductor) ------------------------------------------
#
# torch.compile costs every process ~2 s to import the compiler, ~5 s to trace
# the first GPU's stages and ~1 s for each further GPU; a 10 s dense scan paid
# that before its first chunk. The same stages exported with AOTInductor load
# in milliseconds, one library serves every GPU of an architecture, and a
# four-iteration block runs 8.4M cells in 0.26 ms on the H100 (2.5e-14 from
# the eager stages). They are built once per architecture, torch and stage
# source, like the native decoders (build_device_tail.sh), into
# .build-libs/device_tail_<key>/; inputs are normalised to FP64 (rows, traits).

def _stage_source_key(capability):
    import hashlib
    import inspect

    import torch

    source = ''.join(inspect.getsource(stage) for stage in
                     (_guard, _tail_prologue, _tail_block, _tail_epilogue, _tail_whole))
    text = '|'.join((torch.__version__, str(torch.version.cuda), f'sm{capability[0]}{capability[1]}',
                     str(DEVICE_TAIL_BLOCK), source))
    return hashlib.sha256(text.encode()).hexdigest()[:24]


def device_tail_directory(capability):
    from pathlib import Path

    return Path(__file__).resolve().parents[2] / '.build-libs' / f'device_tail_{_stage_source_key(capability)}'


def build_device_tail(device='cuda'):
    """Export both forms of the tail with AOTInductor for `device`'s architecture.

    The three block stages (prologue, block, epilogue) and the one-graph
    whole tail; neg_log10_p_device picks per strip (_choose_form).
    """
    import json

    import torch
    from torch.export import Dim

    device = torch.device(device)
    capability = torch.cuda.get_device_capability(device)
    directory = device_tail_directory(capability)
    directory.mkdir(parents=True, exist_ok=True)
    rows, traits = Dim('rows', min=2, max=1 << 24), Dim('traits', min=2, max=1 << 24)
    matrix = {0: rows, 1: traits}

    class Prologue(torch.nn.Module):
        def forward(self, t, df):
            return _tail_prologue(t, df)

    class Block(torch.nn.Module):
        def forward(self, first, second, argument, c, d, h, start):
            return _tail_block(first, second, argument, c, d, h, start)

    class Epilogue(torch.nn.Module):
        def forward(self, df, x, y, reflect, h):
            return _tail_epilogue(df, x, y, reflect, h)

    class Whole(torch.nn.Module):
        def forward(self, t, df):
            return _tail_whole(t, df)

    t = torch.linspace(0.0, 40.0, 64 * 33, dtype=torch.float64, device=device).reshape(64, 33)
    df = torch.full_like(t, 30.0)
    first, second, argument, c, d, h, x, y, reflect = _tail_prologue(t, df)
    start = torch.tensor(1.0, dtype=torch.float64, device=device)
    exported = {
        'prologue': (Prologue(), (t, df), (matrix, matrix)),
        'block': (Block(), (first, second, argument, c, d, h, start), (matrix,) * 6 + (None,)),
        'epilogue': (Epilogue(), (df, x, y, reflect, h), (matrix,) * 5),
        'whole': (Whole(), (t, df), (matrix, matrix)),
    }
    for name, (module, args, shapes) in exported.items():
        torch._export.aot_compile(module, args, dynamic_shapes=shapes,
                                  options={'aot_inductor.output_path': str(directory / f'{name}.so')})
    manifest = dict(torch=torch.__version__, cuda=torch.version.cuda, capability=list(capability),
                    block_iterations=DEVICE_TAIL_BLOCK, key=_stage_source_key(capability),
                    files=[f'{name}.so' for name in exported])
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    return directory


def _keeps_subnormals():
    value = np.array([1e-40], dtype=np.float32)
    return bool((value * np.float32(1.0))[0] != 0)


def _load_aoti(device):
    """AOTInductor stages for `device`, or None when this architecture has no build."""
    import json

    import torch

    directory = device_tail_directory(torch.cuda.get_device_capability(device))
    manifest_path = directory / 'manifest.json'
    if not manifest_path.exists():
        return None
    manifest = json.loads(manifest_path.read_text())
    if any(not (directory / name).exists() for name in manifest.get('files', ())):
        return None
    # Inductor links its wrapper with -ffast-math, and loading such a library
    # (crtfastmath) switches the loading thread to flush subnormals to zero,
    # which would change NumPy results on that thread. Restore the mode.
    names = ('prologue', 'block', 'epilogue', 'whole')
    if manifest.get('files') != [f'{name}.so' for name in names]:
        return None
    kept = _keeps_subnormals()
    pro, blk, epi, whole = [torch._export.aot_load(str(directory / f'{name}.so'), device=str(device))
                            for name in names]
    if kept and not _keeps_subnormals():
        torch.set_flush_denormal(False)

    def single(graph):
        # A single-output graph returns its tensor, not a list (indexing it
        # would take the first row).
        def call(*args):
            value = graph(*args)
            return value[0] if isinstance(value, (list, tuple)) else value
        return call

    return ((lambda t, df: tuple(pro(t, df)), lambda *args: tuple(blk(*args)), single(epi)),
            single(whole))
