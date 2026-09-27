"""In-job user-CPU versus wall-time panel for chunk-size tuning.

Question: measured inside one job, does user CPU per variant rank chunk sizes
the way whole-job wall time does, with less repeat-to-repeat variation, and
does a short early window pick the faster size more reliably than early wall
time? The workload is the frozen H100 significant-pair job of the Sep-21
sustained panel: real PGEN prefix, synthetic null traits, empty output. It
runs the current source with explicit settings and blocking CUDA completion
events (so waiting threads do not accrue user time). No calculator or price
enters; these are empirical tuning observations, not component prices.

Each observation is a fresh process with verified cold private inputs.
A sampler thread records process user/sys CPU and faults every 0.25 s, and
per-thread CPU, runnable wait and GPU utilization every 1 s; chunk arrival
times come from the scan iterator. The sampler's own CPU is reported.
"""
import argparse, ctypes, hashlib, json, os, random, statistics, subprocess, sys, threading, time
from pathlib import Path
from unittest.mock import patch

ROOT = Path('results/usertime_panel_20260923')
DATA = Path('/data/zxie3/torchgwas_bench/significant_host_sustained_20260921')
OUT = Path('/data/zxie3/torchgwas_bench/usertime_panel_20260923')
REFERENCE = Path('/data484_4/zxie3/torchGWAS1.1/results/significant_host_sustained_20260921/reference.json')
INPUTS = ['input.pgen', 'input.pvar', 'input.psam', 'phenotype.npy', 'covariates.npy']
M, K, TILE = 1048576, 16385, 8193
ENVIRONMENT = dict(TORCHGWAS_SIGNIFICANCE_BACKEND='host', TORCHGWAS_PGEN_BACKEND='native',
                   TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_SCAN_PROFILE='0',
                   TORCHGWAS_BLOCKING_EVENTS='1', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1',
                   OPENBLAS_NUM_THREADS='1', OMP_WAIT_POLICY='PASSIVE',
                   # The frozen profile disables NumPy's huge-page advice. With THP
                   # defrag=madvise and fragmented memory, leaving it on made QC and
                   # tile materialization stall in direct compaction (smoke run r99:
                   # 1,139 s sys vs 70 s user, first chunk after 979 s).
                   NUMPY_MADVISE_HUGEPAGE='0')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=1, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def resident_pages(path):
    import mmap
    size = path.stat().st_size
    with path.open('rb') as f:
        mm = mmap.mmap(f.fileno(), size, flags=mmap.MAP_PRIVATE, prot=mmap.PROT_READ | mmap.PROT_WRITE)
        count = (size+mmap.PAGESIZE-1)//mmap.PAGESIZE
        vec = (ctypes.c_ubyte*count)(); buf = ctypes.c_char.from_buffer(mm)
        libc = ctypes.CDLL(None, use_errno=True)
        libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p]
        ret = libc.mincore(ctypes.addressof(buf), size, vec); error = ctypes.get_errno()
        del buf; mm.close()
    if ret:
        raise OSError(error, 'mincore failed')
    return sum(x & 1 for x in vec), count


def cold_inputs():
    records = []
    for name in INPUTS:
        path = DATA/name
        for attempt in range(5):
            with path.open('rb') as stream:
                os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
            resident, total = resident_pages(path)
            if not resident:
                break
            time.sleep(.05)
        if resident:
            raise RuntimeError('Cold precondition failed: '+name)
        records.append(dict(name=name, pages=total, resident_at_launch=resident, attempts=attempt+1))
    return records


class _Utilization(ctypes.Structure):
    _fields_ = [('gpu', ctypes.c_uint), ('memory', ctypes.c_uint)]


class Sampler(threading.Thread):
    """Process counters every `fast` s; threads, runnable wait and GPUs every `slow` s."""

    def __init__(self, uuids, fast=0.25, slow=1.0):
        super().__init__(daemon=True)
        self.fast_interval, self.slow_interval = fast, slow
        self.tick = os.sysconf('SC_CLK_TCK')
        self.stop_event = threading.Event()
        self.process, self.threads, self.gpus = [], [], []
        self.nvml = ctypes.CDLL('libnvidia-ml.so.1')
        assert self.nvml.nvmlInit_v2() == 0
        self.handles = {}
        for device, uuid in uuids.items():
            handle = ctypes.c_void_p()
            assert self.nvml.nvmlDeviceGetHandleByUUID(('GPU-'+uuid).encode(), ctypes.byref(handle)) == 0
            self.handles[device] = handle
        self.own_cpu_seconds = None

    @staticmethod
    def _stat(path):
        with open(path, 'rb') as stream:
            text = stream.read().decode()
        comm = text[text.index('(')+1:text.rindex(')')]
        fields = text[text.rindex(')')+2:].split()
        return comm, fields

    def _process(self):
        _, f = self._stat('/proc/self/stat')
        return (time.perf_counter(), int(f[11])/self.tick, int(f[12])/self.tick, int(f[7]), int(f[9]))

    def _threads(self):
        rows = {}
        for tid in os.listdir('/proc/self/task'):
            try:
                comm, f = self._stat(f'/proc/self/task/{tid}/stat')
                with open(f'/proc/self/task/{tid}/schedstat') as stream:
                    run_ns, wait_ns, _ = (int(x) for x in stream.read().split())
            except (FileNotFoundError, ProcessLookupError):
                continue  # thread exited between listing and reading
            rows[tid] = (comm, int(f[11])/self.tick, int(f[12])/self.tick, int(f[7]), run_ns, wait_ns)
        return rows

    def _gpu(self):
        out = {}
        for device, handle in self.handles.items():
            value = _Utilization()
            out[device] = (value.gpu if self.nvml.nvmlDeviceGetUtilizationRates(handle, ctypes.byref(value)) == 0
                           else None)
        return out

    def run(self):
        next_slow = 0.
        while not self.stop_event.is_set():
            now = time.perf_counter()
            self.process.append(self._process())
            if now >= next_slow:
                self.threads.append((now, self._threads()))
                self.gpus.append((now, self._gpu()))
                next_slow = now+self.slow_interval
            self.stop_event.wait(self.fast_interval)
        self.process.append(self._process())
        self.threads.append((time.perf_counter(), self._threads()))
        self.own_cpu_seconds = time.thread_time()
        self.nvml.nvmlShutdown()

    def stop(self):
        self.stop_event.set(); self.join()


def child(chunk, repeat, devices, cpus):
    os.sched_setaffinity(0, set(cpus))
    import numpy as np
    import torch
    torch.set_num_threads(4); torch.backends.cuda.matmul.allow_tf32 = False
    import torchgwas.api as api
    from torchgwas.api import run_linear_gwas
    from torchgwas.detailed_calibration import source_identity
    uuids = {d: str(torch.cuda.get_device_properties(d).uuid).lower().removeprefix('gpu-') for d in devices}
    for device in devices:
        torch.cuda.set_device(device)
        x = torch.ones((32, 32), device=device); (x@x).std(); torch.cuda.synchronize(device); del x
    reference = json.loads(REFERENCE.read_text())
    offsets = {d: i*TILE for i, d in enumerate(devices)}
    guard = threading.Lock(); checked = {}; events = []; original = api.linear_scan_streaming_chunks

    def scan(*args, **kwargs):
        device = kwargs['device']
        with guard:
            offset = offsets[device]; offsets[device] += len(devices)*TILE
        iterator, basis = original(*args, **kwargs)
        wanted = [(j, t-offset) for j, t in enumerate(reference['traits']) if offset <= t < offset+args[1].shape[1]]

        def inspect():
            try:
                for item in iterator:
                    with guard:
                        events.append((time.perf_counter(), device, int(item[0]), int(item[1])))
                    if len(item) == 6:  # dense beta, t, p, df from the host-selection path
                        for i, v in enumerate(reference['variants']):
                            if item[0] <= v < item[1]:
                                r = v-item[0]
                                for j, local in wanted:
                                    value = (float(item[2][r, local]), float(item[3][r, local]),
                                             float(np.asarray(item[5])[r, 0]))
                                    with guard:
                                        checked[i, j] = value
                    yield item
            finally:
                iterator.close()
        return inspect(), basis

    output = OUT/f'c{chunk}_r{repeat}'
    assert not output.exists()
    cold = cold_inputs()
    sampler = Sampler(uuids); sampler.start()
    started = time.perf_counter()
    kwargs = dict(genotype=str(DATA/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                  device=devices[0], trait_devices=list(devices), trait_block=TILE, chunk_size=chunk,
                  reader_workers=4, prefetch_chunks=4, compute_dtype='float32', reduce='significant',
                  significance_threshold=None, sumstats_format='binary', sumstats_fsync=True,
                  sumstats_queue_depth=1, sumstats_fields='beta+t')
    try:
        with patch('torchgwas.api.linear_scan_streaming_chunks', side_effect=scan):
            result = run_linear_gwas(phenotype=DATA/'phenotype.npy', covariates=DATA/'covariates.npy',
                                     output_dir=output, **kwargs)
    finally:
        api_seconds = time.perf_counter()-started
        sampler.stop()
    errors = []
    assert len(checked) == len(reference['variants'])*len(reference['traits']), len(checked)
    for (i, j), (beta, t, df) in checked.items():
        np.testing.assert_allclose([beta, t], [reference['beta'][i][j], reference['tstat'][i][j]],
                                   rtol=4e-4, atol=4e-5)
        assert df == reference['df'][i]; errors.append(abs(t-reference['tstat'][i][j]))
    timing = result.run_metadata['sumstats_write']
    record = dict(chunk=chunk, repeat=repeat, devices=devices, cpus=cpus, uuids=uuids,
                  environment={k: os.environ.get(k) for k in ENVIRONMENT}, api_kwargs=kwargs,
                  api_start=started, api_seconds=api_seconds,
                  executor_seconds=timing.get('setup_scan_and_write_seconds', timing['scan_and_write_seconds']),
                  scan_and_write_seconds=timing['scan_and_write_seconds'],
                  result_rows=result.run_metadata['n_result_rows'], checked_cells=len(checked),
                  max_t_abs_error=max(errors), cold=cold, chunk_events=events,
                  process_samples=sampler.process, thread_samples=sampler.threads, gpu_samples=sampler.gpus,
                  sampler_cpu_seconds=sampler.own_cpu_seconds,
                  peak_allocated={d: torch.cuda.max_memory_allocated(d) for d in devices},
                  source_sha256=source_identity(), harness_sha256=sha(__file__))
    save(ROOT/f'run_c{chunk}_r{repeat}.json', record)
    print(json.dumps({k: record[k] for k in ('chunk', 'repeat', 'executor_seconds', 'api_seconds',
                                              'result_rows', 'checked_cells', 'sampler_cpu_seconds')}), flush=True)


def observe(chunks, repeats, devices, cpus, seed):
    manifest = json.loads((DATA/'manifest.json').read_text())
    for name, row in manifest['files'].items():
        assert (DATA/name).stat().st_size == row['bytes'], name
        assert sha(DATA/name) == row['sha256'], name
    rng = random.Random(seed); schedule = []
    for repeat in range(repeats):
        order = list(chunks); rng.shuffle(order); schedule += [(c, repeat) for c in order]
    plan = dict(chunks=chunks, repeats=repeats, devices=devices, cpus=cpus, seed=seed,
                schedule=schedule, manifest=manifest, reference_sha256=sha(REFERENCE),
                harness_sha256=sha(__file__))
    plan_path = ROOT/'plan.json'
    if plan_path.exists():
        assert json.loads(plan_path.read_text()) == plan, 'Existing plan differs'
    else:
        save(plan_path, plan)
    for chunk, repeat in schedule:
        record, log = ROOT/f'run_c{chunk}_r{repeat}.json', ROOT/f'run_c{chunk}_r{repeat}.log'
        if record.exists():
            continue
        if log.exists():
            raise RuntimeError('Incomplete observation needs inspection: '+str(log))
        with log.open('x') as stream:
            subprocess.run([sys.executable, __file__, 'child', '--chunk', str(chunk), '--repeat', str(repeat),
                            '--devices', *devices, '--cpus', *map(str, cpus)],
                           stdout=stream, stderr=subprocess.STDOUT, check=True, env={**os.environ, **ENVIRONMENT})
        print(open(log).read().strip().splitlines()[-1], flush=True)
    report()


# ---------------------------------------------------------------- analysis

def _interp(samples, t, column):
    """Linear interpolation of a cumulative process counter at time t."""
    times = [s[0] for s in samples]
    if t <= times[0]:
        return samples[0][column]
    for a, b in zip(samples, samples[1:]):
        if a[0] <= t <= b[0]:
            w = 0. if b[0] == a[0] else (t-a[0])/(b[0]-a[0])
            return a[column]+w*(b[column]-a[column])
    return samples[-1][column]


def _wait(threads, t):
    """Summed runnable wait over threads alive at the nearest slow sample."""
    best = min(threads, key=lambda row: abs(row[0]-t))
    return sum(v[5] for v in best[1].values())/1e9


def window(record, low, high):
    """Counters between the arrival times of fractions `low` and `high` of all variants."""
    events = sorted(record['chunk_events'])
    total = sum(e-s for _, _, s, e in events)
    done, times = 0, {}
    for t, _, s, e in events:
        done += e-s
        for f in (low, high):
            if f not in times and done >= f*total:
                times[f] = t
    t0, t1 = times[low], times[high]
    p = record['process_samples']
    variants = (high-low)*total
    user = _interp(p, t1, 1)-_interp(p, t0, 1)
    sys_ = _interp(p, t1, 2)-_interp(p, t0, 2)
    return dict(wall=(t1-t0)/variants*1e6, user=user/variants*1e6, sys=sys_/variants*1e6,
                cpu=(user+sys_)/variants*1e6,
                minflt=(_interp(p, t1, 3)-_interp(p, t0, 3))/variants,
                runnable_wait=(_wait(record['thread_samples'], t1)-_wait(record['thread_samples'], t0))/variants*1e6,
                gpu_util=statistics.mean(v for t, row in record['gpu_samples'] if t0 <= t <= t1
                                         for v in row.values() if v is not None) if any(
                    t0 <= t <= t1 for t, _ in record['gpu_samples']) else None)


def _cv(values):
    return statistics.pstdev(values)/statistics.mean(values)


def report():
    plan = json.loads((ROOT/'plan.json').read_text())
    runs = {}
    for chunk, repeat in plan['schedule']:
        runs.setdefault(chunk, []).append(json.loads((ROOT/f'run_c{chunk}_r{repeat}.json').read_text()))
    metrics = ('wall', 'user', 'cpu', 'sys')
    table = {}
    for chunk, rows in sorted(runs.items()):
        steady = [window(r, .10, .90) for r in rows]
        early = [window(r, .02, .12) for r in rows]
        second = [window(r, .12, .22) for r in rows]
        scans = [max(e[0] for e in r['chunk_events'])-min(e[0] for e in r['chunk_events']) for r in rows]
        table[chunk] = dict(
            executor_seconds=[r['executor_seconds'] for r in rows],
            median_executor_seconds=statistics.median(r['executor_seconds'] for r in rows),
            cv_executor=_cv([r['executor_seconds'] for r in rows]),
            scan_seconds=scans, median_scan_seconds=statistics.median(scans), cv_scan=_cv(scans),
            pre_scan_seconds=[min(e[0] for e in r['chunk_events'])-r['api_start'] for r in rows],
            steady={m: dict(median=statistics.median(w[m] for w in steady), cv=_cv([w[m] for w in steady]))
                    for m in metrics},
            early={m: [w[m] for w in early] for m in metrics},
            steady_rows=steady, early_rows=early, second_rows=second,
            process_user_seconds=[r['process_samples'][-1][1]-r['process_samples'][0][1] for r in rows],
            process_sys_seconds=[r['process_samples'][-1][2]-r['process_samples'][0][2] for r in rows],
            sampler_cpu_seconds=[r['sampler_cpu_seconds'] for r in rows],
            max_t_abs_error=max(r['max_t_abs_error'] for r in rows),
            result_rows=sorted({r['result_rows'] for r in rows}))
    truth = sorted(table, key=lambda c: table[c]['median_executor_seconds'])
    truth_scan = sorted(table, key=lambda c: table[c]['median_scan_seconds'])
    rankings = {m: sorted(table, key=lambda c: table[c]['steady'][m]['median']) for m in metrics}
    # Early-window pairwise decisions across repeat pairings, as a JIT would
    # compare two sizes measured at different moments of the same job.
    chunks = sorted(table)
    accuracy = {}
    for rows_key in ('early_rows', 'second_rows'):
        for target in ('median_executor_seconds', 'median_scan_seconds'):
            decisions = {m: [0, 0] for m in metrics}
            for a in chunks:
                for b in chunks:
                    if a >= b:
                        continue
                    faster = a if table[a][target] < table[b][target] else b
                    for x in table[a][rows_key]:
                        for y in table[b][rows_key]:
                            for m in metrics:
                                guess = a if x[m] < y[m] else b
                                decisions[m][0] += guess == faster; decisions[m][1] += 1
            accuracy[rows_key+':'+target] = {m: v[0]/v[1] for m, v in decisions.items()}
    result = dict(chunks=chunks, true_order_by_median_executor=truth, true_order_by_median_scan=truth_scan,
                  steady_order=rankings, early_pairwise_accuracy=accuracy,
                  early_pairwise_decisions=sum(len(table[a]['early_rows'])*len(table[b]['early_rows'])
                                               for a in chunks for b in chunks if a < b), per_chunk=table,
                  units='wall/user/sys/cpu in microseconds per variant; minflt per variant',
                  scope='One H100 pair, one workload, empty significant output, blocking events, explicit settings; tuning observations, not component prices.')
    path = ROOT/'report.json'
    if path.exists():
        path.unlink()
    save(path, result)
    brief = {c: dict(median_executor=round(t['median_executor_seconds'], 2), cv_executor=round(t['cv_executor'], 3),
                     median_scan=round(t['median_scan_seconds'], 2), cv_scan=round(t['cv_scan'], 3),
                     **{f'{m}_us': round(t['steady'][m]['median'], 2) for m in metrics},
                     **{f'cv_{m}': round(t['steady'][m]['cv'], 3) for m in metrics}) for c, t in table.items()}
    print(json.dumps(dict(true_order=truth, true_order_scan=truth_scan, steady_order=rankings,
                          early_pairwise_accuracy=result['early_pairwise_accuracy'], per_chunk=brief), indent=1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['observe', 'child', 'report'])
    parser.add_argument('--chunk', type=int); parser.add_argument('--repeat', type=int)
    parser.add_argument('--chunks', type=int, nargs='+', default=[512, 1024, 4096])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--devices', nargs='+', default=['cuda:0', 'cuda:3'])
    parser.add_argument('--cpus', type=int, nargs='+', default=[0, 2, 4, 5, 7, 8, 9, 10])
    parser.add_argument('--seed', type=int, default=20260923)
    args = parser.parse_args()
    if args.stage == 'child':
        child(args.chunk, args.repeat, args.devices, args.cpus)
    elif args.stage == 'observe':
        observe(args.chunks, args.repeats, args.devices, args.cpus, args.seed)
    else:
        report()


if __name__ == '__main__':
    main()
