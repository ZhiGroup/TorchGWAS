"""Independent fixed API CPU with observed detached-thread intervals.

Run only with the supplied shared library in LD_AUDIT. Positive/negative
controls must pass before any Torch measurements are considered usable.
"""
import argparse, ast, ctypes, hashlib, json, os, random, statistics, threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

thread_environment = dict(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                          OMP_WAIT_POLICY='PASSIVE', GOMP_SPINCOUNT='0')
os.environ.update(thread_environment)
import torch


def worker(device, barrier, library, bank_kind='scan', checkpoint=None):
    torch.cuda.set_device(device)
    torch.backends.cuda.matmul.allow_tf32 = False
    held = ctypes.PyDLL(library)
    begin, end = held.gil_probe_begin_mode, held.gil_probe_end
    begin.argtypes = [ctypes.c_int]; begin.restype = None
    end.argtypes = [ctypes.POINTER(ctypes.c_double)]; end.restype = None
    released = ctypes.CDLL(library)
    native = released.gil_probe_native_work
    native.argtypes = [ctypes.c_uint]; native.restype = ctypes.c_double
    detached_native = held.gil_probe_detached_work
    detached_native.argtypes = [ctypes.c_uint]; detached_native.restype = ctypes.c_double
    nested = held.gil_probe_nested_work
    nested.argtypes = [ctypes.c_uint]; nested.restype = ctypes.c_double
    meter = held.gil_probe_meter_batch
    meter.argtypes = [ctypes.c_uint]; meter.restype = None
    state = (ctypes.c_double * 4)()
    previous = [None]

    def observe(fn, record=True):
        begin(int(record))
        previous[0] = fn()
        end(state)
        return dict(cpu_seconds=state[0], detached_cpu_seconds=state[1],
                    detached_intervals=int(state[2]), balance_errors=int(state[3]))

    def python_work():
        value = 0
        for i in range(100000): value += i
        return value

    # Resolve lazy PLT bindings and initialize the per-thread meter before
    # checking its steady operation. First-call setup is not native detached
    # CPU work; retain the same strict positive/negative thresholds afterward.
    for _ in range(3):
        observe(lambda:detached_native(1000000))
        observe(lambda:nested(300000))
        observe(python_work)
        observe(lambda:native(1000000))
    controls = dict(native=observe(lambda: detached_native(1000000)), nested=observe(lambda:nested(300000)), python=observe(python_work),
                    ctypes_glob_dat_unobserved=observe(lambda: native(1000000)))
    if checkpoint is not None:
        Path(checkpoint).write_text(json.dumps(dict(device=device,controls=controls,context_verified=False),indent=2))
    good = controls['native']['detached_intervals'] > 0 and controls['native']['balance_errors'] == 0
    good &= controls['native']['detached_cpu_seconds'] / controls['native']['cpu_seconds'] > .9
    good &= controls['python']['detached_intervals'] == controls['python']['balance_errors'] == 0
    good &= controls['nested']['detached_intervals'] == 3 and controls['nested']['balance_errors'] == 0
    good &= controls['nested']['detached_cpu_seconds'] / controls['nested']['cpu_seconds'] > .9
    if not good: raise RuntimeError('GIL audit positive/negative controls failed: ' + str(controls))
    meter_rows = []
    for repeat in range(7):
        for record in ([False,True] if repeat % 2 == 0 else [True,False]):
            measurement = observe(lambda: meter(20000), record)
            if measurement['balance_errors'] or measurement['detached_intervals'] != (20000 if record else 0):
                raise RuntimeError('Meter-only control did not conserve intervals: '+str(measurement))
            meter_rows.append(dict(measurement, repeat=repeat, recording=record, intervals=20000))

    if bank_kind=='jagwas':
        from direct_jagwas_host_primitives import build_bank
        bank=build_bank('cuda:'+str(device))
    elif bank_kind=='scan':
        source = Path('benchmarks/direct_calculator_host_primitives.py')
        wanted = {'a', 'b', 'mask', 'code', 'u8', 'i64', 'v', 'col', 'mask1', 'bank'}
        body = [node for node in ast.parse(source.read_text()).body if isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id in wanted for t in node.targets)]
        namespace = {'torch': torch}
        exec(compile(ast.Module(body=body, type_ignores=[]), str(source), 'exec'), namespace)
        bank = namespace['bank']
    else:raise ValueError('Unknown fixed primitive bank')
    names = list(bank)
    random.Random(920119).shuffle(names)
    for _ in range(20):
        for name in names: value = bank[name]()
    torch.cuda.synchronize()
    rows = []
    for repeat in range(5):
        barrier.wait(timeout=120)
        order = [False,True] if repeat % 2 == 0 else [True,False]
        for record in order:
            barrier.wait(timeout=120)
            samples = {name: [] for name in names}
            previous[0] = None
            empty = [observe(lambda: None,record) for _ in range(100)]
            for _ in range(100):
                for name in names: samples[name].append(observe(bank[name],record))
            torch.cuda.synchronize()
            overhead = statistics.fmean(row['cpu_seconds'] for row in empty)
            for name, values in samples.items():
                raw = statistics.fmean(row['cpu_seconds'] for row in values)
                detached = statistics.fmean(row['detached_cpu_seconds'] for row in values)
                cpu = max(0., raw - overhead)
                rows.append(dict(device=device, repeat=repeat, recording=record, primitive=name, sample_count=len(values),
                    raw_cpu_sum_seconds=sum(row['cpu_seconds'] for row in values), raw_cpu_seconds_per_call=raw,
                    empty_cpu_seconds_per_call=overhead, cpu_seconds_per_call=cpu,
                    detached_cpu_seconds_per_call=detached,
                    held_or_unknown_cpu_seconds_per_call=max(0., cpu-detached),
                    detached_intervals=sum(row['detached_intervals'] for row in values),
                    balance_errors=sum(row['balance_errors'] for row in values)))
    result=dict(device=device, controls=controls, meter_rows=meter_rows,
                context_verified=all(row['balance_errors']==0 for row in rows), rows=rows)
    if checkpoint is not None:Path(checkpoint).write_text(json.dumps(result,indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', required=True); parser.add_argument('--out', required=True)
    parser.add_argument('--bank',choices=['scan','jagwas'],default='scan')
    parser.add_argument('--devices',type=int,nargs='+',default=[1,2])
    args = parser.parse_args(); path = Path(args.out)
    if path.exists(): raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    library = str(Path(args.library).resolve())
    if os.environ.get('LD_AUDIT') != library: raise ValueError('Exact LD_AUDIT library required')
    if not args.devices or len(set(args.devices))!=len(args.devices):raise ValueError('Unique devices required')
    os.sched_setaffinity(0, range(12, 20)); torch.set_num_threads(4)
    results = []
    for mode, devices in [('single'+str(args.devices[0]),args.devices[:1]),('threads',args.devices)]:
        barrier = threading.Barrier(len(devices))
        with ThreadPoolExecutor(max_workers=len(devices)) as pool:
            futures = [pool.submit(worker, device, barrier, library,args.bank,
                str(path.with_name(path.stem+'.'+mode+'.'+str(device)+'.json'))) for device in devices]
            results.append(dict(mode=mode, workers=[future.result() for future in futures]))
    result = dict(results=results, context_verified=all(worker['context_verified'] for mode in results for worker in mode['workers']),
        torch_version=torch.__version__, affinity=sorted(os.sched_getaffinity(0)), thread_environment=thread_environment,
        primitive_bank=args.bank,fixed_shape=[32,32],host=os.uname().nodename,
        devices={str(d):dict(name=torch.cuda.get_device_name(d),capability=list(torch.cuda.get_device_capability(d))) for d in args.devices},
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),Path('benchmarks/direct_calculator_gil_audit.c'),
             Path('benchmarks/direct_jagwas_host_primitives.py' if args.bank=='jagwas' else 'benchmarks/direct_calculator_host_primitives.py'),Path(library)]},
        scope='Independent fixed32 API bank only. PLT-bound PyEval_SaveThread/RestoreThread and PyGILState_Ensure/Release intervals are observed with thread CPU clocks. The same explicit PLT routes are checked by native/nested controls. Recording on/off pairs expose clock-hook perturbation in the same audit process. GLOB_DAT relocations (_ctypes here), other detach APIs, static internal calls, hook overhead and reacquisition CPU remain held/unknown; not an exact GIL occupancy trace or a workload residual. No prices automatically applied to the model.')
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as stream: json.dump(result,stream,indent=2)
    print(json.dumps(dict(context_verified=result['context_verified'], controls=[worker['controls'] for mode in results for worker in mode['workers']]),indent=2))
