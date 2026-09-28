"""Empirical autotuning: layout at startup, chunk size from the job's own chunks.

This path deliberately does not use the analytical calculator or priced
profiles. It measures the running job instead:

* Layout (GPUs, phenotype tiles or variant shards, reader workers) is chosen
  once at startup from live device availability and the existing memory-fit
  model. JAGWAS keeps the full phenotype panel on every GPU and shards
  variants. Every phenotype tile rereads the whole genotype, so tiles are only
  added when memory requires them or when the trait panel is wide enough for
  per-tile work to dominate the reread.
* Chunk size is tuned during the scan. After a short warmup, every candidate
  size runs for one segment of real (fully productive) chunks, in a forward
  then reverse order so a steady drift in server load cancels to first order.
  All partitions switch together, so the score is whole-job throughput,
  including contention for shared CPU, memory, PCIe and storage. The best
  size is then used for the rest of the job.

Changing the chunk size never changes results: sizes are multiples of the
smallest candidate, applied only to reads not yet issued (AlignedChunkSizeControl).
"""
from __future__ import annotations

import math
import os
import statistics
import subprocess
import threading
import time

import numpy as np

from .adaptive_chunks import AlignedChunkSizeControl, aligned_chunk_sizes

# 4096 is opt-in: its ring costs ~6% even when unused (lab-2080ti A/B) and it
# was the slowest size in most frozen-H100 panel repeats.
# Candidates up to 4096: sizes whose rings do not fit are dropped (plan_layout),
# and a short job drops the largest from probing (ModelChunkTuner). Full-scale
# full output (35,365 x 8.9M, K=512) ran at chunk 4096 in 20.4 s of scan where
# the tuner, limited to 2048, took 33.2 s. 8192 was never faster (1M-variant
# PGEN, fixed: warm 7.0 vs 6.0 s, cold 6.4 vs 5.1 s), and full-scale autotune
# runs that committed to it lost 60-500%: full output there is disk-bound
# (cold input read and output write on one array, 1.4 GB/s each at 4096,
# 0.7-0.9 at 8192), and sequential probes early in the run, while the page
# cache still absorbs writes, favored it.
DEFAULT_CHUNK_SIZES = (512, 1024, 2048, 4096)


def frame_aligned_sizes(sizes, frame):
    """Chunk candidates as multiples of a framed source's decode frame F.

    A framed source (zstd store, hard-call store) decodes whole frames, so a
    chunk C that is not a multiple of F decodes 1 + (F - gcd(C, F)) / C
    frames' worth per variant (boundary frames twice). Each size moves to the
    nearest multiple of F (at least F); multiples of the smallest stay so.
    Sources without frames (chunk_alignment_variants None) keep the sizes.
    """
    frame = int(frame or 0)
    if frame <= 1:
        return sorted({int(size) for size in sizes})
    return sorted({frame * max(1, round(int(size) / frame)) for size in sizes})


# Start and fallback size: best or tied-best among fixed chunks in every
# lab-2080ti layout comparison, and the frozen H100 plan's choice.
PREFERRED_CHUNK_SIZE = 1024


def _process_cpu():
    times = os.times()
    return times.user, times.system


class EmpiricalChunkTuner:
    """Chunk-size observer for one job; `control` is the matching selector.

    Pass `tuner.control` as every scan's `_chunk_size_selector` and the tuner
    itself as its `_chunk_observer` (it also answers `for_scan`/`observer`).
    The ring must be allocated for `capacity` (the largest candidate).
    It needs only start, end, device and completion time from an observation,
    so it also works on paths that report MinimalChunkObservation.
    """
    accepts_minimal_observations = True

    def __init__(self, sizes, *, total_rows, initial=None, warmup_fraction=0.03,
                 trial_fraction=0.25, repeats=2, min_gain=0.02, depth=4, concurrent=1,
                 max_trial_share=0.5, min_job_seconds=20.0, clock=time.perf_counter, cpu_clock=_process_cpu):
        """depth x concurrent: chunks already issued when the size changes (the
        prefetch ring of every concurrently scanning GPU). They drain
        unmeasured after each switch, so they count against the trial budget."""
        self.sizes = aligned_chunk_sizes(sizes)
        if isinstance(total_rows, bool) or not isinstance(total_rows, int) or total_rows < 1:
            raise ValueError('Positive total source rows required')
        if not 0 <= warmup_fraction < 1 or not 0 < trial_fraction < 1 or repeats < 1:
            raise ValueError('Invalid tuning fractions or repeats')
        self.total_rows = total_rows
        self.min_gain = float(min_gain)
        self.depth, self.concurrent = int(depth), int(concurrent)
        self._clock, self._cpu = clock, cpu_clock
        self.warmup_rows = int(warmup_fraction*total_rows)
        self.state, self.reason = 'warmup', None
        self.trial_sizes, self.schedule, self.targets, self.planned_rows = (), [], [], None
        if initial is None:
            # A fixed default, not the middle of whatever trial set fits: a
            # skipped or shortened trial must not silently land on 512.
            initial = PREFERRED_CHUNK_SIZE if PREFERRED_CHUNK_SIZE in self.sizes else self.sizes[(len(self.sizes)-1)//2]
        if len(self.sizes) == 1:
            self.state, self.reason = 'fixed', 'single_candidate'
        else:
            # Try every candidate; if the job cannot afford that, drop the
            # largest sizes (the most expensive per switch) until it can.
            for count in range(len(self.sizes), 1, -1):
                trial = self.sizes[:count]
                if initial not in trial:
                    continue  # every trial set includes the starting size
                start = initial
                schedule, targets, planned = self._plan(trial, start, trial_fraction, repeats)
                if planned <= max_trial_share*total_rows:
                    self.trial_sizes, self.schedule, self.targets, self.planned_rows = trial, schedule, targets, planned
                    initial = start
                    if count < len(self.sizes):
                        self.reason = f'largest sizes dropped from trials to fit the job: {list(self.sizes[count:])}'
                    break
            else:
                self.state, self.reason = 'skipped', 'job_too_short_for_trials'
        if initial is None or initial not in self.sizes:
            initial = self.sizes[(len(self.sizes)-1)//2]
        self.initial = initial
        self.control = AlignedChunkSizeControl(self.sizes, initial=initial)
        self.capacity = self.control.capacity
        self._lock = threading.Lock()
        self.min_job_seconds = float(min_job_seconds)
        self.max_trial_share = float(max_trial_share)
        self.check_rows = max(1, total_rows//100)  # deferred re-check every 1% of the job
        self._next_check = 0
        self.first_completed = None
        self.last_completed = None
        self.estimated_job_seconds = None
        self.completed_rows = 0
        self.devices = set()
        self.segments = []
        self._segment = None
        self._next = 0
        self.choice = initial if self.state != 'warmup' else None
        self.decided_at_rows = None
        self.started = clock()

    def _plan(self, trial, start, trial_fraction, repeats):
        order = []
        for index in range(repeats):
            order += list(trial if index % 2 == 0 else reversed(trial))
        # At the forward/reverse junction the same size runs twice in a row;
        # the second switch is a no-op and both segments count. Each segment
        # needs at least three of its own chunks per concurrent GPU.
        base = trial_fraction*self.total_rows/len(order)
        targets = [int(max(base, 3*size*self.concurrent)) for size in order]
        drain = sum(self.depth*self.concurrent*previous for previous in [start]+order[:-1])
        return order, targets, self.warmup_rows+sum(targets)+drain

    # -- observer protocol -------------------------------------------------
    def for_scan(self, device, variant_range, n_traits, **_):
        return self

    def observer(self, device, variant_range=None, trait_range=None, **_):
        return self

    def __call__(self, observation):
        rows = int(observation.end-observation.start)
        now = observation.completed
        with self._lock:
            self.completed_rows += rows
            self.last_completed = now
            device = observation.device
            self.devices.add(device)
            if self.first_completed is None:
                self.first_completed = now
            if self.state in ('warmup', 'deferred'):
                if self.completed_rows < max(self.warmup_rows, self._next_check):
                    return
                elapsed = now-self.first_completed
                remaining = self.total_rows-self.completed_rows
                rate = self.completed_rows/elapsed if elapsed > 0 else None
                estimate = None if rate is None else remaining/rate
                if self.state == 'warmup':
                    self.estimated_job_seconds = None if rate is None else self.total_rows/rate
                trials = self.planned_rows-self.warmup_rows
                if estimate is not None and estimate >= self.min_job_seconds and trials <= self.max_trial_share*remaining:
                    if self.state == 'deferred':
                        self.reason = 'trials started late: the job ran slower than its warmup suggested'
                    self._start_segment(now)
                elif trials > self.max_trial_share*remaining:
                    # Too little work left for the trial schedule.
                    self.state, self.reason = 'skipped', 'job_shorter_than_min_seconds'
                    self.choice = self.initial
                else:
                    # Trials cost switches and drains and cannot repay
                    # themselves on a short job. Early chunks can run far
                    # faster than later ones, so keep re-checking.
                    self.state, self.reason = 'deferred', 'remaining job shorter than min_job_seconds'
                    self._next_check = self.completed_rows+self.check_rows
                return
            if self.state != 'trial':
                return
            segment = self._segment
            if segment['t0'] is None:
                # Wait until every active device has delivered a chunk issued
                # at the new size; after that nothing older is in flight.
                if rows == segment['size']:
                    segment['seen'].add(device)
                if segment['seen'] >= self.devices:
                    segment['t0'] = now
                    segment['cpu0'] = self._cpu()
                elif now-segment['switched'] > 60.:
                    # A device that has finished all its work never delivers
                    # again; do not wait on it forever.
                    segment['t0'] = now
                    segment['cpu0'] = self._cpu()
                    segment['note'] = 'device_wait_timeout'
                return
            segment['rows'] += rows
            segment['other_size_rows'] += rows if rows != segment['size'] else 0
            if segment['rows'] >= segment['target']:
                cpu = self._cpu()
                seconds = now-segment['t0']
                segment.update(seconds=seconds, rows_per_second=segment['rows']/seconds if seconds > 0 else None,
                               user_seconds_per_row=(cpu[0]-segment['cpu0'][0])/segment['rows'],
                               sys_seconds_per_row=(cpu[1]-segment['cpu0'][1])/segment['rows'],
                               ended_at_rows=self.completed_rows)
                del segment['seen'], segment['cpu0']
                self.segments.append(segment)
                if self._next < len(self.schedule):
                    self._start_segment(now)
                else:
                    self._decide()

    # -- internals ---------------------------------------------------------
    def _start_segment(self, now):
        size = self.schedule[self._next]
        self._next += 1
        self.control.set_size(size)
        self.state = 'trial'
        self._segment = dict(size=size, switched=now, switched_at_rows=self.completed_rows,
                             t0=None, rows=0, other_size_rows=0, target=self.targets[self._next-1], seen=set())

    def _decide(self):
        by_size = {}
        for segment in self.segments:
            if segment['rows_per_second']:
                by_size.setdefault(segment['size'], []).append(segment['rows_per_second'])
        scores = {size: statistics.mean(v) for size, v in by_size.items()}
        best = max(scores, key=scores.get)
        # Prefer the starting size unless the winner is clearly better:
        # a switch that gains less than min_gain is within noise.
        if self.initial in scores and scores[best] < scores[self.initial]*(1+self.min_gain):
            best = self.initial
            self.reason = 'within_min_gain_of_initial'
        else:
            self.reason = 'highest_measured_throughput'
        self.control.set_size(best)
        self.choice = best
        self.state = 'committed'
        self.decided_at_rows = self.completed_rows
        self.scores = scores

    def audit(self):
        with self._lock:
            return dict(method='empirical_segments', sizes=list(self.sizes), trial_sizes=list(self.trial_sizes),
                        initial=self.initial, capacity=self.capacity, schedule=list(self.schedule),
                        state=self.state, reason=self.reason, choice=self.choice, total_rows=self.total_rows,
                        warmup_rows=self.warmup_rows, segment_targets=list(self.targets),
                        planned_trial_rows=self.planned_rows, depth=self.depth, concurrent=self.concurrent,
                        estimated_job_seconds=self.estimated_job_seconds, min_job_seconds=self.min_job_seconds,
                        # Where the executor's time went relative to the chunk flow the tuner can influence.
                        first_chunk_seconds=None if self.first_completed is None else self.first_completed-self.started,
                        last_chunk_seconds=None if self.last_completed is None else self.last_completed-self.started,
                        completed_rows=self.completed_rows, decided_at_rows=self.decided_at_rows,
                        decided_fraction=(None if self.decided_at_rows is None
                                          else self.decided_at_rows/self.total_rows),
                        scores_rows_per_second=getattr(self, 'scores', None),
                        segments=[{k: v for k, v in s.items() if k != 'switched'} for s in self.segments],
                        devices=sorted(self.devices),
                        scope='Throughput of real production chunks, all partitions switched together; '
                              'forward/reverse order cancels linear load drift, not arbitrary load changes.')


# ------------------------------------------------------------------ layout

def eligible_gpus(max_utilization=20, min_free_bytes=4 << 30):
    """CUDA devices with no foreign compute process, low utilization and free memory.

    Returns cuda:N names in this process's CUDA numbering (CUDA_VISIBLE_DEVICES
    respected) and a per-device record of what was seen (including free
    bytes). Uses nvidia-smi, so no CUDA context is created on any device;
    a context costs ~0.3-0.5 s per GPU, which added ~4 s over 7 GPUs.
    """
    import torch
    if not torch.cuda.is_available():
        return [], {}
    try:
        rows = subprocess.check_output(['nvidia-smi', '--query-gpu=uuid,utilization.gpu,memory.free',
                                        '--format=csv,noheader,nounits'], text=True, timeout=20)
        apps = subprocess.check_output(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid',
                                        '--format=csv,noheader'], text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        rows = apps = ''
    physical = {}
    for line in rows.splitlines():
        uuid, util, free = (x.strip() for x in line.split(','))
        physical[uuid.lower().removeprefix('gpu-')] = dict(utilization=int(util), free_bytes=int(free) << 20,
                                                            foreign_pids=[])
    for line in apps.splitlines():
        if line.strip():
            uuid, pid = (x.strip() for x in line.split(','))
            row = physical.get(uuid.lower().removeprefix('gpu-'))
            if row is not None and int(pid) != os.getpid():
                row['foreign_pids'].append(int(pid))
    chosen, seen = [], {}
    for index in range(torch.cuda.device_count()):
        uuid = str(torch.cuda.get_device_properties(index).uuid).lower().removeprefix('gpu-')
        row = physical.get(uuid)
        if row is None:  # nvidia-smi unavailable: fall back to CUDA's own free-memory view
            free, _ = torch.cuda.mem_get_info(index)
            row = dict(utilization=None, free_bytes=int(free), foreign_pids=[])
        seen[f'cuda:{index}'] = row
        if (not row['foreign_pids'] and (row['utilization'] is None or row['utilization'] <= max_utilization)
                and row['free_bytes'] >= min_free_bytes):
            chosen.append(f'cuda:{index}')
    return chosen, seen


def gpu_free_bytes(devices, seen=None):
    """Free memory per device from an eligible_gpus record or nvidia-smi.

    Falls back to torch.cuda.mem_get_info (which creates a context) only
    when nvidia-smi gave no answer for a device.
    """
    import torch
    if seen is None:
        _, seen = eligible_gpus(max_utilization=100, min_free_bytes=0)
    free = {}
    for device in devices:
        row = seen.get(str(device)) if seen else None
        free[str(device)] = (row['free_bytes'] if row and row.get('free_bytes') is not None
                             else int(torch.cuda.mem_get_info(torch.device(device))[0]))
    return free


def usable_cpus(sample_seconds=0.25):
    """Idle CPU capacity within this process's affinity, measured briefly.

    Shared servers often have most cores busy; counting the affinity mask
    alone would plan readers for cores other users are running on. Returns
    (usable, detail) with usable = floor of summed idle fractions, at least 1.
    """
    try:
        allowed = sorted(os.sched_getaffinity(0))
    except AttributeError:
        count = os.cpu_count() or 1
        return count, dict(affinity=count, idle_cores=None)

    def sample():
        rows = {}
        try:
            with open('/proc/stat', encoding='ascii') as stream:
                for line in stream:
                    if line.startswith('cpu') and line[3].isdigit():
                        name, *fields = line.split()
                        values = [int(x) for x in fields]
                        rows[int(name[3:])] = (sum(values), values[3]+(values[4] if len(values) > 4 else 0))
        except OSError:
            pass
        return rows
    before = sample()
    time.sleep(sample_seconds)
    after = sample()
    idle = 0.
    for cpu in allowed:
        if cpu in before and cpu in after:
            total = after[cpu][0]-before[cpu][0]
            idle += (after[cpu][1]-before[cpu][1])/total if total > 0 else 1.
        else:
            idle += 1.
    try:
        load = os.getloadavg()[0]
    except (AttributeError, OSError):
        load = None
    return max(1, int(idle)), dict(affinity=len(allowed), idle_cores=round(idle, 2),
                                   load_1min=None if load is None else round(load, 2))


def measured_gemm_rate(device, dtype, size=2048, repeats=3):
    """FLOP/s of a size^3 GEMM on `device`, median of repeats (lab-2080ti: 50 ms in FP64, 5 ms in FP32)."""
    import torch
    device = torch.device(device)
    a = torch.randn(size, size, device=device, dtype=dtype)
    b = torch.randn(size, size, device=device, dtype=dtype)
    a @ b
    torch.cuda.synchronize(device)
    seconds = []
    for _ in range(repeats):
        start = time.perf_counter()
        a @ b
        torch.cuda.synchronize(device)
        seconds.append(time.perf_counter() - start)
    return 2.0 * size ** 3 / sorted(seconds)[len(seconds) // 2]


def jagwas_projection_flops(n_traits, group_sizes=None):
    """FP64 projection FLOPs per variant: one K-wide factor, or one per JAGWAS group.

    Each group projects only its own columns (jagwas_projection.JagwasGroups),
    so a grouped job costs sum(k_g^2)-order work, not K^2: 22 groups of 100
    traits are ~22x cheaper than one 2,200-trait panel.
    """
    from .jagwas_blocks import projection_flops_per_variant
    if not group_sizes:
        return projection_flops_per_variant(int(n_traits))
    return sum(projection_flops_per_variant(int(size)) for size in group_sizes)


def jagwas_factor_bytes(n_traits, group_sizes=None):
    """FP64 factor bytes a GPU holds for JAGWAS while preparing.

    One panel: correlation, factor and inverse, 3 x K x K. Groups factor one
    at a time (reduction_tensor_work.require_jagwas_factor_capacity): the
    largest group's three k x k matrices plus every other group's retained
    k x k factor.
    """
    if not group_sizes:
        return 3 * int(n_traits) ** 2 * 8
    sizes = [int(size) for size in group_sizes]
    largest = max(sizes)
    return 3 * largest * largest * 8 + 8 * (sum(size * size for size in sizes) - largest * largest)


def gpu_seconds_per_variant(device, *, mode, n_samples, n_traits, group_sizes=None, path='torch'):
    """GPU time one variant costs: the scan's device work, plus the FP64 projection for JAGWAS.

    The larger of the FP32 GEMM alone and the scan's own per-chunk device
    work, timed (scan_device_seconds_per_variant). The GEMM alone was the
    price until 2026-09-27 and is 3.4x short at K = 512 (0.47 against 1.6 us
    per variant), which made CPU demand per GPU 3.4x too high: a K = 512
    min-p job was held to 2 GPUs (14.6 cores each) and ran 18-21 s against
    9.6-10.3 s on four. group_sizes (JAGWAS groups) prices each group's own
    projection.
    """
    import torch
    seconds = 2.0 * n_samples * n_traits / measured_gemm_rate(device, torch.float32)
    measured = scan_device_seconds_per_variant(device, mode=mode, n_samples=n_samples, n_traits=n_traits,
                                               path=path)
    if measured is not None:
        seconds = max(seconds, measured)
    if mode == 'jagwas':
        seconds += jagwas_projection_flops(n_traits, group_sizes) / measured_gemm_rate(device, torch.float64)
    return seconds


def scan_statistics_path(genotype):
    """The statistics backend a CUDA scan of `genotype` runs, as linear_scan_streaming_chunks dispatches it.

    'bed' (packed PLINK rows, iter_packed_chunks), 'native-pgen2' / 'native-dosage'
    (the fused kernels, TORCHGWAS_NATIVE_STATS=1, packed or byte transport),
    'torch' (native dosage transport, torch statistics) or 'generic' (host-decoded
    float chunks, kernels.linear_chunk_kernel).
    """
    from .scan_gpu import resolve_statistics_backend
    if hasattr(genotype, 'iter_packed_chunks'):
        return 'bed'
    if not getattr(genotype, 'supports_fused_qc', False):
        return 'generic'
    if resolve_statistics_backend() != 'native_fused':
        return 'torch'
    return 'native-pgen2' if getattr(genotype, 'native_encoding', None) == 'pgen_2bit' else 'native-dosage'


def scan_device_seconds_per_variant(device, *, mode, n_samples, n_traits, repeats=3, path='torch'):
    """GPU seconds per variant of the scan's own per-chunk device work, timed; None if not measurable.

    Times, with CUDA events on a synthetic chunk, what the scan's statistics
    backend (`path`, scan_statistics_path) does to every chunk -- decoding or
    converting the calls, the statistics (masking, centring, range, the GEMM
    against an n x (K + 1) design, residual sums) -- and the mode's own
    step: the -log10 P tail over every cell (full output), the min-p
    reduction and its winners' tail, or the significance selection. JAGWAS's
    reduction is priced by its projection FLOPs instead (gpu_seconds_per_variant).

    Panels wider than 2,048 traits are timed at 256 and 2,048 and extended
    linearly in K: the GEMM and every per-cell pass are linear in K, the
    per-variant passes constant. The chunk is 1,024 variants (fewer for large
    n, at most 2^27 genotype cells). About 50 ms.

    H100, 22,250 samples, K = 512, timed in the scan itself (TORCHGWAS_SCAN_PROFILE,
    benchmarks/gpu_pipeline_profile_20260927.py): 1.58-1.65 us per variant of
    GPU compute for min-p, 1.85 us for dense output; the GEMM alone 0.47.
    None off CUDA, or when the path's kernels are unavailable.
    """
    import torch
    device = torch.device(device)
    if device.type != 'cuda':
        return None
    if path.startswith('native'):
        from . import scan_gpu
        if not scan_gpu.available(device):
            return None
    chunk = int(max(64, min(1024, (1 << 27) // max(1, int(n_samples)))))
    widths = sorted({min(int(n_traits), 256), min(int(n_traits), 2048)})
    timed = {width: _device_chunk_seconds(device, mode, int(n_samples), width, chunk, repeats, path)
             for width in widths}
    if len(widths) == 1 or int(n_traits) <= widths[-1]:
        seconds = timed[widths[-1]]
    else:
        low, high = widths
        seconds = timed[high] + (timed[high] - timed[low]) / (high - low) * (int(n_traits) - high)
    return seconds / chunk


def _statistics_for_path(path, device, n_samples, width, chunk, generator):
    """statistics() -> (beta, t, status, variant_df) for one synthetic chunk on `path`'s kernels."""
    import torch
    design = torch.randn(n_samples, width + 1, generator=generator, device=device) / n_samples ** 0.5
    phenotype_ss = torch.ones(width, device=device)
    df = n_samples - 3
    packed_width = (n_samples + 3) // 4
    if path in ('bed', 'native-pgen2'):
        # Two-bit calls; any byte is valid (a missing code is one of the four).
        packed = torch.randint(0, 256, (chunk, packed_width), generator=generator, device=device,
                               dtype=torch.int32).to(torch.uint8)
    else:
        codes = torch.randint(0, 3, (chunk, n_samples), generator=generator, device=device, dtype=torch.uint8)
    if path == 'bed':
        from .linear import _select_packed_bed_statistics
        kernel = _select_packed_bed_statistics(None, None)
        return lambda: kernel(packed, design, phenotype_ss, n_samples, width, df, None, None)
    if path.startswith('native'):
        from . import scan_gpu

        def native():
            if path == 'native-pgen2':
                centered, ss, low, high, present = scan_gpu.prepare(packed, encoding='pgen_2bit', n_samples=n_samples)
            else:
                centered, ss, low, high, present = scan_gpu.prepare(codes, 1.0, None)
            beta, t, status = scan_gpu.finish(centered @ design, ss, low, high, phenotype_ss, present, -3.0)
            return beta, t, status, present.to(torch.float32) - 3.0
        return native
    if path == 'generic':
        from .kernels import linear_chunk_kernel
        # Host-decoded float chunks, samples x variants, as the generic scan uploads them.
        values = codes.t().float()
        phenotype = design[:, :width].contiguous()
        basis = None

        def generic():
            beta, t, variant_df = linear_chunk_kernel(values, phenotype, basis, df, covariate_rank=0)
            return beta, t, torch.zeros(t.shape[0], dtype=torch.uint8, device=device), variant_df
        return generic
    from .linear import _dosage_statistics
    return lambda: _dosage_statistics(codes.float(), design, phenotype_ss, width, df, False, covariate_rank=1)


def _device_chunk_seconds(device, mode, n_samples, width, chunk, repeats, path='torch'):
    """CUDA-event seconds of one synthetic chunk's device work (scan_device_seconds_per_variant)."""
    import torch
    with torch.cuda.device(device):
        generator = torch.Generator(device=device).manual_seed(20260927)
        statistics = _statistics_for_path(path, device, n_samples, width, chunk, generator)
        step = None
        if mode == 'full':
            from .tails import neg_log10_p_device

            def step(beta, t, status, variant_df):
                neg_log10_p_device(t, variant_df.reshape(-1, 1))
        elif mode == 'min-p':
            from .min_p import MinPReduction
            reduction = MinPReduction()

            def step(beta, t, status, variant_df):
                reduction.reduce(beta, t, status, variant_df, 1, log10_p=(None, None, torch.float32))
        elif mode == 'significant':
            from .reduce import SignificantPairs, device_significance_critical, device_significant_pairs
            critical = device_significance_critical(SignificantPairs(), n_samples, width, device)

            def step(beta, t, status, variant_df):
                for _ in device_significant_pairs(beta, t, status, variant_df.float(), critical):
                    pass

        def once():
            beta, t, status, variant_df = statistics()
            if step is not None:
                step(beta, t, status, variant_df)

        once()
        torch.cuda.synchronize(device)
        seconds = []
        for _ in range(repeats):
            start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
            start.record()
            once()
            stop.record()
            stop.synchronize()
            seconds.append(start.elapsed_time(stop) / 1e3)
        return min(seconds)


def shard_setup_seconds(device):
    """Serial cost of bringing one more GPU into a job: CUDA context plus cuBLAS handle, timed.

    Fitting T(G) = a + s*G + W/G to lab-2080ti JAGWAS shards (1, 2, 4 GPUs)
    gives s = 1.0 s; the context alone measured 0.8 s there. Call it on a
    GPU the job will use anyway (the second one), so nothing is wasted.
    """
    import torch
    device = torch.device(device)
    start = time.perf_counter()
    a = torch.ones(64, 64, device=device)
    a @ a
    torch.cuda.synchronize(device)
    return time.perf_counter() - start


def output_write_rates(directory, *, n_traits, store_beta=True, writers=4, probe_bytes=256 << 20,
                       chunk_rows=1024):
    """Dense-output bytes per second for one writer and for `writers` at once.

    Each probe writer is the scan's own BinarySumstatsWriter (staging copy,
    beta/t/df block streams, fsync at close) over probe_bytes of a zero
    panel, in a scratch directory inside `directory`, removed afterwards.
    Full-scale H100 dense output, K = 2,048: one GPU (one writer) wrote
    1.1 GB/s and took 120 s, four variant shards (four writers) 2.9 GB/s and
    45 s. The writer, not the GPU, set the one-GPU time.
    """
    import shutil
    import tempfile
    from pathlib import Path
    from .sumstats import BinarySumstatsWriter
    # beta (optional), t and -log10 P per cell, and one df per variant.
    per_variant = n_traits * (12 if store_beta else 8) + 4
    rows = max(chunk_rows, int(probe_bytes // per_variant) // chunk_rows * chunk_rows)
    beta = np.zeros((chunk_rows, n_traits), dtype=np.float32)
    df = np.full((chunk_rows, 1), 100.0, dtype=np.float32)
    names = [f't{i}' for i in range(n_traits)]
    scratch = Path(tempfile.mkdtemp(prefix='.torchgwas-write-probe-', dir=directory))

    def write(target):
        writer = BinarySumstatsWriter(target, rows, names, 102, 100, fsync=True, store_beta=store_beta,
                                      store_variant_df=True)
        try:
            for start in range(0, rows, chunk_rows):
                end = min(rows, start + chunk_rows)
                writer.write_chunk(start, end, beta[:end - start] if store_beta else None,
                                   beta[:end - start], beta[:end - start], variant_df=df[:end - start])
            writer.close()
        except BaseException:
            writer.abort()
            raise

    try:
        started = time.perf_counter()
        write(scratch / 'single')
        single = rows * per_variant / (time.perf_counter() - started)
        shutil.rmtree(scratch / 'single')
        threads = [threading.Thread(target=write, args=(scratch / f'w{i}',)) for i in range(max(1, int(writers)))]
        started = time.perf_counter()
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        aggregate = len(threads) * rows * per_variant / (time.perf_counter() - started)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return dict(bytes_per_variant=per_variant, writer_bytes_per_second=single,
                aggregate_bytes_per_second=aggregate, aggregate_writers=len(threads), probe_rows=rows)


def dense_shard_seconds(shards, *, n_variants, gpu_seconds_per_variant, setup_seconds,
                        bytes_per_variant, writer_bytes_per_second, aggregate_bytes_per_second,
                        aggregate_writers):
    """Modeled dense-output time on `shards` variant shards, each with its own writer.

    Serial setup per shard, plus the larger of the split GPU work and the
    output write at min(shards x one writer's rate, the measured aggregate);
    shards never exceed the probed writer count.
    """
    rate = min(shards * writer_bytes_per_second, aggregate_bytes_per_second)
    return (setup_seconds * shards
            + max(n_variants * gpu_seconds_per_variant / shards, n_variants * bytes_per_variant / rate))


def _attainable_cores(gpus, cpu_cores, cpu_load, readers=4):
    """Cores our threads get under fair time-slicing: T threads among L others on C cores get min(T, T*C/max(C, L+T)).

    T = readers + 1 per GPU (its readers and its shard thread) plus the
    writer; L is the 1-minute load less this process.
    """
    threads = gpus * (readers + 1) + 1
    others = max(0.0, float(cpu_load) - 1.0)
    return min(threads, threads * cpu_cores / max(cpu_cores, others + threads))


MAX_READERS_PER_GPU = 16


def readers_per_gpu(gpus, *, cpus, cpu_cores=None, cpu_load=None, decode_cpu_per_variant=None,
                    gpu_seconds_per_variant=None):
    """Decode threads per GPU, from measured demand.

    A GPU consuming a variant every g seconds offers decode/g reader-cores of
    work, more when each thread gets only a share of a core; readers are the
    smallest M/M/c server count with P(wait) <= READER_WAIT_PROBABILITY, 2 to
    MAX_READERS_PER_GPU. There is no supply cap on the count: time-slicing
    gives every runnable thread a share, so on a busy host more readers get
    the job more CPU, not less. Full-scale full output on H100 at load
    116-182: a supply cap held autotune to 2 readers, and it lost 1.7-3.2x to
    16 fixed readers (PGEN dosage, zstd). Without a demand measurement: 4.
    """
    if not (decode_cpu_per_variant and gpu_seconds_per_variant):
        return 4
    # Offered load in reader-cores: decode work per unit of GPU time, over
    # the share of a core each reader thread gets under time-slicing.
    share = 1.0
    if cpu_cores and cpu_load is not None:
        threads = gpus * 5 + 1  # the baseline of 4 readers and a shard thread per GPU, plus the writer
        share = min(1.0, cpu_cores / max(cpu_cores, max(0.0, float(cpu_load) - 1.0) + threads))
    offered = decode_cpu_per_variant / gpu_seconds_per_variant / share
    readers = MAX_READERS_PER_GPU
    for count in range(2, MAX_READERS_PER_GPU + 1):
        if _erlang_c(count, offered) <= READER_WAIT_PROBABILITY:
            readers = count
            break
    return readers


# Chunk-decode times vary, so readers matched to the mean demand leave the GPU
# waiting: size them as M/M/c servers with P(a chunk waits) <= 0.2.
READER_WAIT_PROBABILITY = 0.2


def _erlang_c(servers, offered):
    """Erlang-C probability that a request waits, `servers` servers at `offered` erlangs (1 if saturated)."""
    if offered <= 0:
        return 0.0
    if offered >= servers:
        return 1.0
    term, total = 1.0, 1.0
    for k in range(1, servers):
        term *= offered / k
        total += term
    top = term * offered / servers * servers / (servers - offered)
    return top / (total + top)


def best_shard_count(work_seconds, setup_seconds, limit):
    """argmin over G in 1..limit of setup*G + work/G (serial setup per shard, work split evenly)."""
    return min(range(1, max(1, int(limit)) + 1), key=lambda g: setup_seconds * g + work_seconds / g)


def decode_cpu_seconds_per_variant(genotype, first, last, variants=512):
    """CPU seconds the host spends per variant on the scan's own transport, or None if not measurable.

    Reads `variants` rows at the start of the range the way the scan does:
    a native fill into the transfer form (PGEN, zstd store), or packed 2-bit
    rows for PLINK sources, which the GPU unpacks (read_packed_into; with a
    hard-call store this is its zstd decompression). The first half warms
    the reader (and the page cache), the second half is timed in process CPU
    time, so disk waits do not count. Other sources (BGEN, whose decode may
    run on the GPU) return None.

    A framed source (chunk_alignment_variants F: zstd store, hard-call store)
    decodes whole frames, so it is timed on one whole aligned frame: 256 rows
    inside a 2,048-2,500-variant frame charged the entire frame to them
    (8-10x). The buffer is faulted in before timing, as the scan's pinned
    ring is.
    """
    frame = int(getattr(genotype, 'chunk_alignment_variants', None) or 0)
    if frame > 1:
        first = -(-int(first) // frame) * frame
        half = frame if first + 2 * frame <= int(last) else 0
    else:
        half = min(int(variants), int(last) - int(first)) // 2
    if half < 8:
        return None
    if (not getattr(genotype, 'allows_direct_native_fill', False) and hasattr(genotype, 'read_packed_into')
            and getattr(genotype, '_bytes_per_variant', None)):
        width = int(genotype._bytes_per_variant)
        packed = np.empty(2 * half * width, dtype=np.uint8)
        packed.fill(0)  # fault the pages in (np.zeros may map lazily)
        view = memoryview(packed).cast('B')
        try:
            genotype.read_packed_into(view[:half * width], first, first + half)
            start = time.process_time()
            genotype.read_packed_into(view[half * width:], first + half, first + 2 * half)
            seconds = time.process_time() - start
        except Exception:
            return None
        return seconds / half
    if not (getattr(genotype, 'allows_direct_native_fill', False) and hasattr(genotype, 'native_reader_session')):
        return None
    dtype = np.dtype(getattr(genotype, 'native_transfer_dtype', getattr(genotype, 'native_dtype', np.int8)))
    width = int(getattr(genotype, 'native_row_width', 0) or genotype.shape[0])
    raw = np.empty(2 * half * width * dtype.itemsize + 64, dtype=np.uint8)
    raw.fill(0)
    offset = (-raw.ctypes.data) % 64  # packed PGEN rows must be 64-byte aligned
    out = raw[offset:offset + 2 * half * width * dtype.itemsize].view(dtype).reshape(2 * half, width)
    try:
        with genotype.native_reader_session() as fill:
            fill(first, first + half, out[:half])
            start = time.process_time()
            fill(first + half, first + 2 * half, out[half:])
            seconds = time.process_time() - start
    except Exception:
        return None
    return seconds / half


def plan_layout(*, mode, n_samples, n_traits, covariate_rank, n_variants, devices, cpus,
                capacity, depth, transfer_bytes_per_variant, device_free_bytes, host_free_bytes,
                min_tile_traits=2048, max_tile_traits=32768, cpus_per_device=None, reduction_width=None,
                allow_partitions=True, chunk_sizes=None, decode_cpu_per_variant=None,
                gpu_seconds_per_variant=None, host_cores_per_device=0.25, shard_setup_seconds=None,
                cpu_cores=None, cpu_load=None, output_rates=None, group_sizes=None):
    """Choose GPUs and the phenotype/variant partition for one job.

    mode is 'full', 'significant', 'min-p' or 'jagwas'. Returns a dict with
    trait_block, trait_devices, variant_devices, device, reader_workers and a
    'why' list. The memory model is pipeline_model.device_ring_bytes (via
    auto_trait_block) evaluated at the largest chunk candidate.
    allow_partitions=False (filtered full output, which only the single-device
    writer supports) keeps full and JAGWAS output on one GPU.

    chunk_sizes: the tuner's candidates (default: capacity alone). Sizes whose
    rings cannot hold even a minimal tile (for JAGWAS: the whole panel) on
    the device or in pinned host memory are dropped, largest first, and the
    rest are planned at the largest survivor; if none fits, ValueError.

    CPU cap on GPUs: with decode_cpu_per_variant and gpu_seconds_per_variant
    (measured), one GPU needs decode/gpu cores of decode plus
    host_cores_per_device for its shard thread, and GPUs are capped at
    (cpus - 1) / that (one core for the writer and main thread). With the
    host's cpu_cores and 1-minute cpu_load, the supply is what our threads can
    attain under fair time-slicing rather than the idle cores alone: T
    runnable threads among L others on C cores get min(T, T*C/max(C, L+T)).
    Without the measurements, or with an explicit cpus_per_device, the cap is
    cpus // cpus_per_device (default 3).

    Variant-shard count: with gpu_seconds_per_variant and shard_setup_seconds
    (measured), G minimizes setup*G + n_variants*gpu_seconds_per_variant/G;
    otherwise at most n_variants // (8*capacity) shards.

    group_sizes (JAGWAS groups): the factor memory is priced per group
    (jagwas_factor_bytes); the panel itself still sits whole on every GPU.

    output_rates (dense output, from output_write_rates): each shard writes its
    own store, so a shard serves a variant in max(GPU time, output bytes / one
    writer's rate), and G minimizes dense_shard_seconds. That service time,
    not the GPU time alone, also sizes the readers and the CPU cap.
    """
    if not devices:
        raise ValueError('At least one usable GPU is required')
    sizes = sorted({int(s) for s in (chunk_sizes or [capacity])})

    def feasible(depth, shards, need):
        """Drop sizes whose rings cannot hold `need` traits; returns (kept, dropped)."""
        from .pipeline_model import device_ring_bytes, host_pinned_bytes
        kept, dropped = list(sizes), []
        def fits(chunk):
            if device_ring_bytes(chunk_variants=chunk, depth=depth, n_samples=n_samples, n_traits=need,
                                 covariate_rank=covariate_rank,
                                 transfer_bytes_per_variant=transfer_bytes_per_variant) > 0.85*device_free_bytes:
                return False
            return not host_free_bytes or host_pinned_bytes(
                chunk_variants=chunk, depth=depth, n_traits=need,
                transfer_bytes_per_variant=transfer_bytes_per_variant,
                reduction_width=reduction_width)*shards <= 0.85*host_free_bytes
        while kept and not fits(kept[-1]):
            dropped.append(kept.pop())
        if not kept:
            what = (f'the full phenotype panel ({n_traits} traits); JAGWAS cannot tile phenotypes'
                    if mode == 'jagwas' and need == n_traits else f'{need} traits')
            raise ValueError(f'Not even chunk {sizes[0]} fits {what}: {n_samples} samples, '
                             f'{device_free_bytes/2**30:.1f} GiB free on the GPU'
                             + (f', {host_free_bytes/2**30:.1f} GiB host' if host_free_bytes else ''))
        return kept, dropped

    if mode in ('full', 'jagwas') and not allow_partitions:
        readers = max(2, min(4, cpus))
        depth = max(int(depth), readers)
        kept, dropped = feasible(depth, 1, n_traits if mode == 'jagwas' else min(n_traits, 64))
        why = ['row-filtered output uses the single-device writer']
        if dropped:
            why.append(f'chunk sizes {sorted(dropped)} dropped: rings do not fit memory')
        return dict(trait_block=None, trait_devices=None, variant_devices=None, device=devices[0],
                    devices_considered=list(devices), devices_used=1, chunk_sizes=kept,
                    chunk_sizes_dropped=sorted(dropped),
                    reader_workers=readers, prefetch_chunks=depth, why=why)
    from .pipeline_model import auto_trait_block, MIN_AUTO_TRAIT_BLOCK
    why = []
    # A brief idle-CPU sample on a saturated host can read ~1 core; two GPUs
    # still beat one in every tiled comparison on both hosts (our threads get
    # a fair share), so the CPU check never goes below two.
    cpu_demand = None
    # What one GPU's pipeline spends per variant: its GEMM, or for dense output
    # its writer when that is slower (K = 512: 0.47 us GEMM, 3.2 us of writing).
    writer_seconds = (output_rates['bytes_per_variant'] / output_rates['writer_bytes_per_second']
                      if output_rates else None)
    service_seconds = (max(gpu_seconds_per_variant, writer_seconds)
                       if gpu_seconds_per_variant and writer_seconds else gpu_seconds_per_variant)
    if cpus_per_device is None and decode_cpu_per_variant is not None and service_seconds:
        # Measured demand: GEMM-heavy scans (wide panels, JAGWAS FP64) need a
        # fraction of a core per GPU; 4 JAGWAS shards on lab-2080ti used 1.5
        # cores in total, 4 significant-pair shards at K=8,192 about 2.
        # The writer thread is busy writer_seconds of every service interval.
        per_gpu = (decode_cpu_per_variant / service_seconds + host_cores_per_device
                   + (writer_seconds / service_seconds if writer_seconds else 0.0))
        if cpu_cores and cpu_load is not None:
            # Busy is not full: other users' threads share cores with ours.
            # Full-scale A100 (96 cores, load ~82, 6-9 idle): 4 significant-
            # pair shards 230.8 s where the idle-core cap took 2 (348.1 s).
            fits = [g for g in range(1, len(devices) + 1)
                    if g * per_gpu + 1 <= _attainable_cores(g, cpu_cores, cpu_load)]
            by_cpu = max(min(2, len(devices)), max(fits, default=1))
        else:
            by_cpu = max(min(2, len(devices)), int(max(cpus - 1, 1) / per_gpu))
        cpu_demand = dict(decode_cpu_per_variant=decode_cpu_per_variant, gpu_seconds_per_variant=gpu_seconds_per_variant,
                          writer_seconds_per_variant=writer_seconds,
                          cores_per_gpu=round(per_gpu, 3), host_cores_per_device=host_cores_per_device, gpus_by_cpu=by_cpu,
                          cpu_cores=cpu_cores, cpu_load=cpu_load,
                          supply='fair_share' if cpu_cores and cpu_load is not None else 'idle_cores')
        rule = (f'{per_gpu:.2f} cores per GPU: decode {decode_cpu_per_variant*1e6:.1f} us CPU per '
                f'{service_seconds*1e6:.1f} us per variant '
                + (f'(writer {writer_seconds*1e6:.1f} us, GPU {gpu_seconds_per_variant*1e6:.1f} us)'
                   if writer_seconds else 'GPU')
                + f', plus {host_cores_per_device} host')
    else:
        cpus_per_device = cpus_per_device or 3
        by_cpu = max(min(2, len(devices)), cpus//cpus_per_device)
        rule = f'{cpus_per_device} CPUs per GPU'
    if len(devices) > by_cpu:
        why.append(f'{len(devices)} idle GPUs but {cpus} CPUs; using {by_cpu} ({rule})')
        devices = devices[:by_cpu]
    # Decode concurrency per GPU is capped by its prefetch depth, so the ring
    # is as deep as the readers it feeds; memory is checked at that depth.
    # Measured on lab-2080ti (same 2-GPU layout): 4 readers per GPU at depth 4
    # beat 2 per GPU (26.8 s vs 31.9 s API) and 8 per GPU at depth 8 (scan
    # 17.9 s vs 19.5 s). Fewer only when the host has almost no idle CPU.
    # Readers from measured demand within the fair-share supply: full-scale
    # A100 shards with 2 readers per GPU (idle-core sizing) ran 276 s, with 4
    # ran 218-231 s.
    per_device_readers = readers_per_gpu(len(devices), cpus=cpus, cpu_cores=cpu_cores, cpu_load=cpu_load,
                                         decode_cpu_per_variant=decode_cpu_per_variant,
                                         gpu_seconds_per_variant=service_seconds)
    depth = max(int(depth), per_device_readers)
    # Every candidate the tuner may switch to must fit: plan at the largest
    # size that can hold a minimal tile (JAGWAS: the whole panel).
    kept, dropped = feasible(depth, len(devices),
                             n_traits if mode == 'jagwas' else min(n_traits, MIN_AUTO_TRAIT_BLOCK))
    if mode == 'jagwas':
        # JAGWAS runs only when the full panel and its K x K float64 factor
        # (correlation, factor and inverse while preparing) fit every GPU.
        from .pipeline_model import device_ring_bytes
        ring = device_ring_bytes(chunk_variants=kept[0], depth=depth, n_samples=n_samples, n_traits=n_traits,
                                 covariate_rank=covariate_rank, transfer_bytes_per_variant=transfer_bytes_per_variant)
        if ring + jagwas_factor_bytes(n_traits, group_sizes) > 0.85 * device_free_bytes:
            raise ValueError(f'the JAGWAS factor for {n_traits} traits does not fit one GPU '
                             f'({device_free_bytes/2**30:.1f} GiB free); JAGWAS needs the full panel on every GPU')
    if dropped:
        why.append(f'chunk sizes {sorted(dropped)} dropped: their rings do not fit memory')
    capacity = kept[-1]
    fit = auto_trait_block(n_samples=n_samples, n_traits=n_traits, covariate_rank=covariate_rank,
                           chunk_variants=capacity, depth=depth,
                           transfer_bytes_per_variant=transfer_bytes_per_variant,
                           device_memory_bytes=int(device_free_bytes), host_memory_bytes=host_free_bytes or None,
                           reduction_width=reduction_width, trait_devices=len(devices))
    result = dict(trait_block=None, trait_devices=None, variant_devices=None, device=devices[0],
                  fit_traits=int(fit), devices_considered=list(devices), chunk_sizes=kept,
                  chunk_sizes_dropped=sorted(dropped))
    if cpu_demand is not None:
        result['cpu_demand'] = cpu_demand
    min_shard_variants = 8*capacity
    # min-p writes one row per variant: like JAGWAS its shards need no
    # writer of their own, and like full output it tiles phenotypes (merged
    # per variant) only when the panel does not fit.
    if mode == 'jagwas' or (mode in ('full', 'min-p') and fit >= n_traits):
        if mode == 'jagwas' and fit < n_traits:
            raise ValueError(f'The full phenotype panel ({n_traits} traits) does not fit one GPU at chunk {capacity}; '
                             'JAGWAS cannot tile phenotypes')
        # Shards split the variants, so nothing is reread, but each GPU adds
        # setup and a narrow panel gives it little work. lab-2080ti, 3 rounds
        # of 3 repeats: JAGWAS at K=256 was 4-13% faster on 2 shards than on 1
        # and slower on 4 (its per-variant K x K projection is real GPU work);
        # dense output at K=64 was 3-6% slower on 2 (writer-bound). Narrow
        # JAGWAS uses at most 2 GPUs, narrow dense output 1; wide panels all.
        usable = max(1, min(len(devices), n_variants//min_shard_variants))
        if mode == 'full' and output_rates and gpu_seconds_per_variant and shard_setup_seconds:
            # Dense output: each shard adds a writer as well as a GPU.
            modeled = {g: dense_shard_seconds(g, n_variants=n_variants, gpu_seconds_per_variant=gpu_seconds_per_variant,
                                              setup_seconds=shard_setup_seconds,
                                              **{k: output_rates[k] for k in ('bytes_per_variant', 'writer_bytes_per_second',
                                                                              'aggregate_bytes_per_second', 'aggregate_writers')})
                       for g in range(1, usable + 1)}
            best = min(modeled, key=modeled.get)
            result['shard_model'] = dict(dense_seconds={str(g): round(s, 3) for g, s in modeled.items()},
                                         setup_seconds=round(shard_setup_seconds, 3), best=best, limit=usable,
                                         output_rates=dict(output_rates))
            why.append(f'dense output: {best} shard(s) model {modeled[best]:.1f} s '
                       f'(one writer {output_rates["writer_bytes_per_second"]/1e9:.2f} GB/s, '
                       f'{output_rates["aggregate_writers"]} writers {output_rates["aggregate_bytes_per_second"]/1e9:.2f} GB/s, '
                       f'GPU {gpu_seconds_per_variant*1e6:.2f} us per variant)')
            usable = best
        elif gpu_seconds_per_variant and shard_setup_seconds:
            # Each shard adds serial setup; the work splits. H100 JAGWAS at
            # K=8,192, 200k variants: 7 shards 9-15 s, 2 shards 2.4 s.
            work = n_variants * gpu_seconds_per_variant
            best = best_shard_count(work, shard_setup_seconds, usable)
            result['shard_model'] = dict(work_seconds=round(work, 3), setup_seconds=round(shard_setup_seconds, 3),
                                         best=best, limit=usable)
            if best < usable:
                why.append(f'{best} shards minimize {shard_setup_seconds:.2f} s setup per GPU + {work:.1f} s GPU work / GPUs')
            usable = best
        measured = (bool(output_rates) if mode == 'full'
                    else bool(gpu_seconds_per_variant and shard_setup_seconds))
        if n_traits < 2*min_tile_traits and not measured:
            # Without measurements (for dense output, of the writer): the
            # lab-2080ti rule above; min-p output is as small as JAGWAS's.
            # With them the setup/work model decides: at K = 512 on the H100,
            # min-p ran 9.6-10.3 s on four shards against 13.5 on two.
            usable = min(usable, 2) if mode in ('jagwas', 'min-p') else 1
        if usable > 1:
            result['variant_devices'] = devices[:usable]
            why.append(f'full panel fits one GPU: {usable} variant shards, genotype read once')
        elif len(devices) > 1:
            why.append('single GPU: ' + (f'{n_traits} traits is too narrow for sharding to pay'
                                         if n_traits < 2*min_tile_traits else 'too few variants to shard'))
        else:
            why.append('single GPU available')
        used = usable
    else:
        # Phenotype tiles. Each tile rereads and decodes the whole genotype,
        # so tiles cost CPU; they pay only when the per-tile trait work is
        # large. Measured: K=8,192 on lab-2080ti, 2 tiles 26.7 s, 1 tile 37.1 s,
        # 4 tiles 44.4 s; K=16,385 on H100, 2 tiles beat 4 and 1. At voxel
        # scale (K ~ 10^5-10^6) the GEMM dominates and more tiles win, hence
        # one more tile per max_tile_traits. Memory can force more tiles.
        min_tiles = 1 if n_traits < 2*min_tile_traits else 2
        by_width = max(min_tiles, -(-n_traits//max_tile_traits))
        used = max(1, min(len(devices), by_width))
        tile = min(int(fit), -(-n_traits//used))
        if tile >= n_traits:
            why.append('whole panel on one GPU' + (f' (only {n_traits} traits; tiles below '
                       f'{min_tile_traits} would reread the genotype for little work)' if len(devices) > 1 else ''))
        else:
            limited = tile == int(fit) and tile < -(-n_traits//used)
            tiles = -(-n_traits//tile)
            used = min(len(devices), tiles)
            if tiles > used and tiles % used:
                # Memory-limited: each GPU scans several tiles in turn. Round
                # the count up to a multiple of the GPUs so no GPU runs an
                # extra pass while the others idle (17 tiles on 4 GPUs is 5
                # rounds; 20 narrower tiles are 5 full rounds).
                tiles = -(-tiles//used)*used
                tile = -(-n_traits//tiles)
                tiles = -(-n_traits//tile)
            result['trait_block'] = int(tile)
            result['trait_devices'] = devices[:used]
            rounds = -(-tiles//used)
            why.append(f'{tiles} phenotype tiles of {tile} traits over {used} GPUs'
                       + (f', {rounds} tiles per GPU in turn' if rounds > 1 else '')
                       + (' (memory-limited)' if limited else ''))
    per_used = readers_per_gpu(used, cpus=cpus, cpu_cores=cpu_cores, cpu_load=cpu_load,
                               decode_cpu_per_variant=decode_cpu_per_variant,
                               gpu_seconds_per_variant=service_seconds)
    from .pipeline_model import device_ring_bytes, host_pinned_bytes
    need = int(result['trait_block']) if result.get('trait_block') else n_traits

    def ring_fits(slots):
        chunk = max(kept)
        device = device_ring_bytes(chunk_variants=chunk, depth=slots, n_samples=n_samples, n_traits=need,
                                   covariate_rank=covariate_rank, transfer_bytes_per_variant=transfer_bytes_per_variant)
        if mode == 'jagwas':
            device += jagwas_factor_bytes(n_traits, group_sizes)
        host = host_pinned_bytes(chunk_variants=chunk, depth=slots, n_traits=need,
                                 transfer_bytes_per_variant=transfer_bytes_per_variant, reduction_width=reduction_width)
        return device <= 0.85 * device_free_bytes and (not host_free_bytes or host * used <= 0.85 * host_free_bytes)

    # Double-buffered ring: a slot per reader being filled and as many filled
    # slots waiting for the GPU, checked at the largest candidate. With depth
    # equal to the readers, every slot is being decoded at once and the GPU
    # waits on the next fill; what matters is the variants in flight. H100,
    # hard-call store, full output K=512, 16 readers, scan seconds: chunk 2048
    # x depth 16 (32k variants) 38.8 / 60.8, 2048 x 32 and 4096 x 16 (65k)
    # 28.2 and 28.3 / 45.5, 4096 x 32 (131k) 32.4 / 81.8 (two rounds, load
    # 77-140). Readers first sized for the GPUs the CPU supply allowed also
    # deepen when fewer GPUs are used (full-scale full output, K=512, 1 GPU:
    # depth 8 from a 2-GPU plan held the readers at 8 where 15 fit).
    for slots in (2 * per_used, per_used):
        if slots <= depth:
            break
        if ring_fits(slots):
            depth = slots
            break
    result['reader_workers'] = used*min(depth, per_used)
    result['prefetch_chunks'] = depth
    result['devices_used'] = used
    result['why'] = why
    return result
