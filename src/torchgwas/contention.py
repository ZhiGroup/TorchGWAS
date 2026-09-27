"""Idealized equal-weight CPU-sharing scenarios, not exact Linux scheduling.

Inputs must describe the same CPU affinity/cgroup domain. A runnable count
is not load average, which includes blocked tasks and lags current activity.
Affinity overlap, CPU quotas, task weights, NUMA, SMT and wake-up delays can
invalidate a homogeneous equal-share approximation. Use independently
measured available-core capacity when it is available; retain its load/time
context. Never infer that capacity from the association runtime being tested.

These legacy helpers provide conditional sharing arithmetic. They do not
turn a host load average into a portable timing guarantee.
"""
from __future__ import annotations


def effective_cores(runnable_threads: float, foreign_load: float,
                    cores_total: int) -> float:
    """Cores under an idealized equal-weight sharing assumption.

    `runnable_threads` is how many threads the process keeps *runnable*, which
    is not always what it asked for -- see the module docstring. `foreign_load`
    is everyone else's runnable count, i.e. the load average minus our own
    contribution.
    """
    if runnable_threads <= 0:
        return 0.0
    if cores_total <= 0:
        raise ValueError("cores_total must be positive")
    foreign_load = max(float(foreign_load), 0.0)
    total = foreign_load + runnable_threads
    if total <= cores_total:
        # Nothing is waiting, so the process runs on every thread it has.
        return float(runnable_threads)
    return cores_total * runnable_threads / total


def contention_factor(runnable_threads: float, foreign_load: float,
                      cores_total: int) -> float:
    """How much slower host-bound work runs, against the same work on an idle box.

    1.0 means no contention. This is the number to divide a quiet-host rate by,
    and it is why a model carrying idle-host rates reads 3-4x fast on a loaded
    one rather than being wrong about the physics.
    """
    quiet = effective_cores(runnable_threads, 0.0, cores_total)
    busy = effective_cores(runnable_threads, foreign_load, cores_total)
    if busy <= 0:
        raise ValueError("no cores available")
    return quiet / busy


def runnable_threads_from_measurement(observed_cores: float,
                                      foreign_load: float,
                                      cores_total: int) -> float:
    """Invert the rule: how many threads was the tool actually keeping busy?

    Takes `Percent of CPU this job got` / 100 and the load at the time, and
    returns the runnable-thread count that explains it. This is how a TOOL's
    parallel behaviour gets measured separately from the scheduler's
    arithmetic, rather than the two being conflated into one fudge factor.
    """
    if observed_cores <= 0:
        return 0.0
    if observed_cores >= cores_total:
        # Saturating the box; the rule cannot resolve a thread count above it.
        return float(cores_total)
    # observed = C * T / (L + T)  =>  T = observed * L / (C - observed)
    denominator = cores_total - observed_cores
    if denominator <= 0:
        return float(cores_total)
    return observed_cores * max(float(foreign_load), 0.0) / denominator


def disk_contention_factor(our_streams: float, foreign_blocked: float,
                           achieved_bytes_per_second: float | None = None,
                           quiet_bytes_per_second: float | None = None) -> float:
    """How much slower I/O runs than on an idle box. A SEPARATE mechanism.

    CPU contention and disk contention are not the same parameter and must not
    be folded into one, because they have different physics and they move
    independently:

      CPU   CFS shares runnable threads. The slowdown is set by the RUNNABLE
            count and is roughly `(L + T) / T`, unbounded as L grows.
      DISK  a block device has a fixed queue and a bandwidth ceiling. Readers
            do NOT get an equal share of a scheduler's time slices; they get a
            share of a saturating resource, so the slowdown flattens once the
            device is at its ceiling no matter how many more readers arrive.

    They are also visible in different places. Load average counts D-state
    tasks -- blocked on I/O -- which contribute nothing to CPU contention but
    are exactly the signal for disk contention. A host with load 120 that is
    all D-state has an idle CPU and a saturated disk; one with load 120 all
    runnable has the reverse. Using load average for both, as a single
    "contention factor" would, gets one of the two wrong every time.

    The honest form is a MEASUREMENT: `achieved / quiet` from a probe run
    alongside the workload. The structural fallback below is used only when no
    probe exists, and it is deliberately crude -- it says the device is shared
    among the streams queued on it -- because a disk's behaviour under
    contention depends on the device, the scheduler and the access pattern in
    ways a one-line model has no business claiming to capture.
    """
    # `is not None`, not truthiness: a measured 0.0 means the probe failed or
    # the device stalled, and silently falling through to the structural
    # fallback would report "no contention" for a disk that delivered nothing.
    if achieved_bytes_per_second is not None and quiet_bytes_per_second is not None:
        if achieved_bytes_per_second <= 0:
            raise ValueError(
                "achieved disk bandwidth must be positive; got "
                f"{achieved_bytes_per_second!r}")
        if quiet_bytes_per_second <= 0:
            raise ValueError(
                "quiet disk bandwidth must be positive; got "
                f"{quiet_bytes_per_second!r}")
        return quiet_bytes_per_second / achieved_bytes_per_second
    if our_streams <= 0:
        raise ValueError("our_streams must be positive")
    foreign_blocked = max(float(foreign_blocked), 0.0)
    return (foreign_blocked + our_streams) / our_streams


# Which hardware each rate actually runs on. Host CPU contention slows the
# reader pool's zstd decode; it does NOT slow a GEMM executing on the GPU, and
# disk contention is a separate mechanism with its own factor.
#
# This has mattered less than it looks so far, and only by luck: the scan is
# decode-bound (78.97 GB at 20.01 GB/s is 3.95 s of a 4.68 s wall), so
# contending everything and contending only decode give nearly the same answer.
# The one validated case -- 9.01 s quiet against 35.14 s at load 160, 3.90x
# measured versus 3.667x predicted -- cannot tell the two apart. At high K the
# GEMM binds instead, and there the difference IS the prediction.
RATE_RESOURCE = {
    'decode_bytes_per_second': 'host_cpu',
    'text_parse_bytes_per_second': 'host_cpu',
    'disk_bytes_per_second': 'disk',
    'write_bytes_per_second': 'disk',
    # DMA engines move these; the host issues the copy but does not perform it.
    'h2d_bytes_per_second': 'uncontended',
    'd2h_bytes_per_second': 'uncontended',
    # On the GPU. A busy host does not slow it down.
    'gemm_flops_per_second': 'uncontended',
}


def apply_by_resource(rates: dict, host_cpu_factor: float = 1.0,
                      disk_factor: float = 1.0) -> dict:
    """Contend each rate by the factor for the hardware it actually uses.

    A rate this module does not recognise is left ALONE rather than contended
    by default: silently slowing a term whose hardware was never identified is
    how a GPU rate came to be divided by a host-CPU factor in the first place.
    """
    for factor in (host_cpu_factor, disk_factor):
        if factor <= 0:
            raise ValueError("contention factors must be positive")
    by_resource = {'host_cpu': host_cpu_factor, 'disk': disk_factor,
                   'uncontended': 1.0}
    adjusted = {}
    for key, value in rates.items():
        factor = by_resource[RATE_RESOURCE.get(key, 'uncontended')]
        if isinstance(value, dict):
            adjusted[key] = {k: v / factor for k, v in value.items()}
        elif isinstance(value, (int, float)):
            adjusted[key] = value / factor
        else:
            adjusted[key] = value
    return adjusted


def apply_to_rates(rates: dict, factor: float) -> dict:
    """Divide EVERY rate by one factor, including nested curves.

    SUPERSEDED by `apply_by_resource`; kept because callers exist. This
    docstring used to justify contending `gemm_flops_per_second` on the grounds
    that "silently skipping it would leave the one term that matters most at
    high K uncontended". That reasoning is wrong: the GEMM runs on the GPU and
    a busy host does not slow it. The one validated case cannot distinguish the
    two -- the scan is decode-bound there -- but at high K, where the GEMM
    binds, this function contends the binding term for a reason that does not
    exist.
    """
    if factor <= 0:
        raise ValueError("contention factor must be positive")
    adjusted = {}
    for key, value in rates.items():
        if isinstance(value, dict):
            adjusted[key] = {k: v / factor for k, v in value.items()}
        elif isinstance(value, (int, float)):
            adjusted[key] = value / factor
        else:
            adjusted[key] = value
    return adjusted
