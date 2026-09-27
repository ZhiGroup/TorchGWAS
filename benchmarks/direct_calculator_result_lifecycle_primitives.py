"""Fixed 32-row exact-source finish and Event handoff primitives; no GWAS input.

The production finish body is extracted, not reimplemented. A no-op event
replaces CUDA synchronization; its separate CPU control permits replacement by
the independently measured ready-CUDA-event service. Shape-dependent status
work and allocator history beyond this tiny fixture are not measured here.
"""
import argparse, ast, ctypes, hashlib, inspect, json, os, platform
import statistics, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1',
                  NUMPY_MADVISE_HUGEPAGE='0')
import numpy as np
import torch
from torchgwas.native_scan import dosage_cuda_iterator


class ReadyEvent:
    def synchronize(self):
        pass


def finish_worker(index, barrier, library, reduction=None, return_beta=True):
    audit = ctypes.PyDLL(library)
    begin, end = audit.gil_probe_begin_mode, audit.gil_probe_end
    begin.argtypes = [ctypes.c_int]; begin.restype = None
    end.argtypes = [ctypes.POINTER(ctypes.c_double)]; end.restype = None
    native = audit.gil_probe_detached_work
    native.argtypes = [ctypes.c_uint]; native.restype = ctypes.c_double
    meter = audit.gil_probe_meter_batch
    meter.argtypes = [ctypes.c_uint]; meter.restype = None
    state = (ctypes.c_double * 4)()

    def observe(fn, record):
        begin(int(record)); value = fn(); del value; end(state)
        result = dict(cpu_seconds=state[0], detached_cpu_seconds=state[1],
                      detached_intervals=int(state[2]), balance_errors=int(state[3]))
        if result['balance_errors']:
            raise RuntimeError('Unbalanced finish GIL transition')
        return result

    controls = dict(native=observe(lambda: native(1000000), True),
                    python=observe(lambda: sum(range(10000)), True))
    assert controls['native']['detached_intervals'] == 1
    assert controls['native']['detached_cpu_seconds'] > .9 * controls['native']['cpu_seconds']
    assert controls['python']['detached_intervals'] == 0
    tree = ast.parse(inspect.getsource(dosage_cuda_iterator))
    bodies = [node for node in tree.body[0].body if isinstance(node, ast.FunctionDef) and node.name == 'finish']
    assert len(bodies) == 1
    event = ReadyEvent()
    if reduction=='jagwas':
        from torchgwas.reduce import JagwasReduction
        reducer=JagwasReduction()
        buffers=[reducer.host_buffers(32,1,pin_memory=False)]
        buffers[0][0].fill_(float('nan'));buffers[0][1].fill_(1.)
        buffers[0][2].zero_();buffers[0][3].zero_();buffers[0][4].fill_(20.)
    elif reduction is None:
        reducer=None
        buffers = [(torch.ones((32, 1)) if return_beta else None, torch.ones((32, 1)),
                    torch.zeros(32, dtype=torch.uint8), torch.full((32,), 20.))]
    else:raise ValueError('Unknown finish reduction')
    namespace = dict(np=np, result_done=[event], result_buffers=buffers,
                     reduction=reducer, compute_p_values=False, profiling=False,return_df=False)
    exec(compile(ast.Module(body=bodies, type_ignores=[]), 'exact_native_finish', 'exec'), namespace)
    finish = namespace['finish']
    meters, rows = [], []
    for repeat in range(7):
        for record in ([False, True] if repeat % 2 == 0 else [True, False]):
            result = observe(lambda: meter(20000), record)
            assert result['detached_intervals'] == (20000 if record else 0)
            meters.append(dict(result, repeat=repeat, recording=record, intervals=20000))
            for ownership in (['owned', 'borrowed'] if repeat % 2 == 0 else ['borrowed', 'owned']):
                namespace['borrow_results'] = ownership == 'borrowed'
                value = finish(0, 0, 32)
                assert value[0][3].flags.owndata == (ownership == 'owned')
                np.testing.assert_array_equal(value[0][3], np.ones((32, 1)))
                if reduction=='jagwas':
                    assert len(value[0])==6 and np.isnan(value[0][2]).all()
                    np.testing.assert_array_equal(value[0][5],np.zeros((32,1),dtype=np.int32))
                elif return_beta:np.testing.assert_array_equal(value[0][2], np.ones((32, 1)))
                else:assert value[0][2] is None
                assert value[1:3] == (0, 0)
                del value
                barrier.wait(timeout=120)
                empty = [observe(lambda: None, record) for _ in range(100)]
                dummy = [observe(event.synchronize, record) for _ in range(100)]
                samples = [observe(lambda: finish(0, 0, 32), record) for _ in range(300)]
                rows.append(dict(repeat=repeat, ownership=ownership, recording=record,
                    worker=index, samples=samples, empty_samples=empty,
                    dummy_event_samples=dummy))
    return dict(worker=index, controls=controls, meter_rows=meters, rows=rows,
                context_verified=True)


def acknowledgement_pair(index, barrier, loops=200):
    import queue
    commands, returned = queue.Queue(), queue.Queue()

    def receiver():
        while True:
            command = commands.get()
            if command is None:
                return
            event, stamp = command
            cpu = time.thread_time()
            assert event.wait(.05), 'Unexpected acknowledgement timeout'
            returned.put(dict(receiver_cpu_seconds=time.thread_time()-cpu,
                publication_to_return_seconds=time.perf_counter()-stamp['published']))

    thread = threading.Thread(target=receiver)
    thread.start(); rows = []
    try:
        barrier.wait(timeout=120)
        for _ in range(loops):
            cpu = time.thread_time(); event = threading.Event()
            create_cpu = time.thread_time()-cpu
            stamp = {}; commands.put((event, stamp))
            deadline = time.monotonic()+10
            while True:
                with event._cond:
                    waiting = bool(event._cond._waiters)
                if waiting:
                    break
                if time.monotonic() > deadline:
                    raise RuntimeError('Acknowledgement receiver did not block')
                time.sleep(0)
            cpu = time.thread_time(); stamp['published'] = time.perf_counter(); event.set()
            publish_cpu = time.thread_time()-cpu
            row = returned.get(timeout=10)
            rows.append(dict(row, create_cpu_seconds=create_cpu,
                publisher_cpu_seconds=publish_cpu, blocked_verified=True))
    finally:
        commands.put(None); thread.join(timeout=12)
        if thread.is_alive():
            raise RuntimeError('Acknowledgement receiver did not terminate')
    return dict(worker=index, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True); parser.add_argument('--library', required=True)
    parser.add_argument('--reduction',choices=['jagwas'])
    parser.add_argument('--omit-beta',action='store_true')
    parser.add_argument('--skip-acknowledgement',action='store_true')
    args = parser.parse_args()
    if args.omit_beta and args.reduction:parser.error('--omit-beta requires unreduced output')
    out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
    library = str(Path(args.library).resolve())
    assert os.environ.get('LD_AUDIT') == library
    os.sched_setaffinity(0, range(12, 20)); torch.set_num_threads(4)
    reports, handoffs = [], []
    for count in [1, 4]:
        barrier = threading.Barrier(count)
        with ThreadPoolExecutor(max_workers=count) as pool:
            futures = [pool.submit(finish_worker, index, barrier, library,args.reduction,not args.omit_beta) for index in range(count)]
            workers = [future.result() for future in futures]
        reports.append(dict(worker_count=count, workers=workers))
        (out/f'finish{count}.json').write_text(json.dumps(reports[-1], indent=2))
        print('FINISH_CONTEXT_COMPLETE', count, flush=True)
    for repeat in range(0 if args.skip_acknowledgement else 5):
        for pairs in [1, 2]:
            barrier = threading.Barrier(pairs)
            with ThreadPoolExecutor(max_workers=pairs) as pool:
                futures = [pool.submit(acknowledgement_pair, index, barrier) for index in range(pairs)]
                workers = [future.result() for future in futures]
            handoffs.append(dict(repeat=repeat, pairs=pairs, workers=workers))
    report = dict(results=reports, acknowledgement=handoffs, context_verified=True,return_df=False,
        numpy_version=np.__version__, torch_version=torch.__version__, python_version=sys.version,host=os.uname().nodename,
        libc=list(platform.libc_ver()), affinity=sorted(os.sched_getaffinity(0)),
        numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
        thread_environment={key:os.environ.get(key) for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
        source_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in
            [Path(__file__),Path('src/torchgwas/native_scan.py'),Path('src/torchgwas/reduce.py'),Path(library)]},
        timing_contract='Exact source finish(0,0,32), K=1, status clear, no reduction or p-values. Includes tensor slice/numpy conversion and destruction; dummy event synchronization is supplied separately. CPU prices must replace fixed finish plus slice/numpy control, not be added to them. Blocking CUDA wait is excluded.',
        scope='Independent tiny owned/borrowed finish and blocked Event acknowledgement only. No GWAS inputs, scan durations, timing fit or shape sweep. Larger status work, real CUDA synchronization, timeout retries and loaded-context latency remain separate. Raw signed control differences must be retained.')
    if args.reduction=='jagwas':
        from torchgwas.owned_result_work import owned_result_work
        report.update(reduction='jagwas',result_layout=owned_result_work(32,1,reduction='jagwas')['array_bytes'])
        report['timing_contract']=report['timing_contract'].replace('no reduction or p-values.','JAGWAS reduction, no p-values.')
    if args.omit_beta:
        from torchgwas.owned_result_work import owned_result_work
        report.update(return_beta=False,result_layout=owned_result_work(32,1,return_beta=False)['array_bytes'])
    (out/'primitives.json').write_text(json.dumps(report, indent=2))
    print('RESULT_PRIMITIVES_COMPLETE', flush=True)
