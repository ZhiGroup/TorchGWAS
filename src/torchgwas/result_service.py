"""Independent source finish and consumer acknowledgement service prices."""
import math
import statistics


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Finite nonnegative result-service observation required')
    return value


def sample(row, recording=True):
    cpu, detached = number(row['cpu_seconds']), number(row['detached_cpu_seconds'])
    intervals = row['detached_intervals']
    if type(intervals) is not int or intervals < 0 or row['balance_errors'] != 0 or detached > cpu:
        raise ValueError('Invalid result-service GIL observation')
    if not recording and (intervals or detached):
        raise ValueError('Disabled result meter reports detached work')
    return cpu, detached, intervals


def result_finish_prices(probe, *, workers, numpy_version, torch_version,
                         python_version, cpu_affinity, source_sha256, reduction=None, return_beta=True):
    if reduction not in (None,'jagwas'):
        raise ValueError('Unsupported finish reduction')
    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced finish layout')
    if probe.get('return_beta',True)!=return_beta:
        raise ValueError('Result finish beta layout mismatch')
    if not return_beta:
        from .owned_result_work import owned_result_work
        if (probe.get('return_df') is not False or
            probe.get('result_layout')!=owned_result_work(32,1,return_beta=False)['array_bytes']):
            raise ValueError('T-only finish requires the exact three-array layout without returned df')
    contract=('Exact source finish(0,0,32), K=1, status clear, no reduction or p-values.' if reduction is None else
              'Exact source finish(0,0,32), K=1, status clear, JAGWAS reduction, no p-values.')
    if reduction=='jagwas':
        from .owned_result_work import owned_result_work
        if (probe.get('return_df') is not False or
            probe.get('result_layout')!=owned_result_work(32,1,reduction='jagwas')['array_bytes']):
            raise ValueError('JAGWAS finish requires the exact five-array layout without returned df')
    if (probe.get('context_verified') is not True or probe['numpy_version'] != numpy_version
        or probe['torch_version'] != torch_version or probe['python_version'] != python_version
        or probe['affinity'] != list(cpu_affinity) or probe['numpy_madvise_hugepage'] is not False
        or probe['source_sha256'].get('src/torchgwas/native_scan.py') != source_sha256
        or probe.get('reduction') != reduction or not probe['timing_contract'].startswith(contract)):
        raise ValueError('Result finish source or runtime context mismatch')
    contexts = [row for row in probe['results'] if row['worker_count'] == workers]
    if type(workers) is not int or workers < 1 or len(contexts) != 1:
        raise ValueError('One complete finish worker context required')
    records = contexts[0]['workers']
    if sorted(row['worker'] for row in records) != list(range(workers)):
        raise ValueError('Unique complete finish workers required')
    grouped, coverage, meters = {}, [], []
    for worker in records:
        if worker.get('context_verified') is not True:
            raise ValueError('Unverified finish worker')
        for name, intervals in [('native', 1), ('python', 0)]:
            cpu, detached, count = sample(worker['controls'][name])
            if count != intervals or not cpu or (intervals and detached/cpu <= .9) or (not intervals and detached):
                raise ValueError('Result GIL positive/negative control failed')
        meter_pairs = {}
        for row in worker['meter_rows']:
            repeat, recording = row['repeat'], row['recording']
            if type(repeat) is not int or repeat < 0 or type(recording) is not bool:
                raise ValueError('Invalid result meter coordinates')
            cpu, detached, intervals = sample(row, recording)
            pair = meter_pairs.setdefault(repeat, {})
            if recording in pair or row['intervals'] != 20000 or intervals != (20000 if recording else 0):
                raise ValueError('Result meter interval accounting failed')
            pair[recording] = (cpu, detached)
        if len(meter_pairs) < 3 or any(set(pair) != {False, True} for pair in meter_pairs.values()):
            raise ValueError('Complete paired result meters required')
        overhead = statistics.median((p[True][0]-p[False][0])/20000 for p in meter_pairs.values())
        detached_overhead = statistics.median(p[True][1]/20000 for p in meter_pairs.values())
        if not 0 <= detached_overhead <= overhead:
            raise ValueError('Invalid result meter overhead partition')
        meters.append(dict(worker=worker['worker'], cpu_seconds=overhead, detached_seconds=detached_overhead))
        seen = set()
        for row in worker['rows']:
            repeat, ownership, recording = row['repeat'], row['ownership'], row['recording']
            key = (repeat, ownership, recording)
            if (type(repeat) is not int or repeat < 0 or ownership not in ('owned', 'borrowed')
                or type(recording) is not bool or row['worker'] != worker['worker'] or key in seen):
                raise ValueError('Invalid or duplicate result finish coordinates')
            seen.add(key)
            banks = {}
            for name in ['samples', 'empty_samples', 'dummy_event_samples']:
                if len(row[name]) < 3:
                    raise ValueError('Insufficient raw result-service observations')
                banks[name] = [sample(value, recording) for value in row[name]]
            if any(value[1] or value[2] for name in ['empty_samples', 'dummy_event_samples'] for value in banks[name]):
                raise ValueError('Dummy/empty result control unexpectedly detaches')
            raw, detached, intervals = [statistics.fmean(s[i] for s in banks['samples']) for i in range(3)]
            dummy = statistics.fmean(s[0] for s in banks['dummy_event_samples'])
            empty = statistics.fmean(s[0] for s in banks['empty_samples'])
            # Subtract the observed dummy event, not another tensor/control
            # estimate. Signed small-control differences remain in the record.
            cpu = raw-dummy-(intervals*overhead if recording else 0.)
            detached -= intervals*detached_overhead if recording else 0.
            if not 0 <= detached <= cpu:
                raise ValueError('Result meter correction exceeds observed service')
            grouped.setdefault(key, []).append(dict(cpu_seconds=cpu, serial_cpu_seconds=cpu-detached,
                raw_cpu_seconds=raw, empty_cpu_seconds=empty,
                dummy_less_empty_cpu_seconds=dummy-empty, intervals=intervals,
                sample_count=len(banks['samples'])))
        repeats = {key[0] for key in seen}
        if len(repeats) < 3 or seen != {(r, o, rec) for r in repeats for o in ['owned', 'borrowed'] for rec in [False, True]}:
            raise ValueError('Complete paired finish coverage required')
        coverage.append(seen)
    if any(keys != coverage[0] for keys in coverage):
        raise ValueError('Finish worker coverage differs')
    prices = {}
    for ownership in ['owned', 'borrowed']:
        repeats = []
        for repeat in sorted({key[0] for key in grouped}):
            pair = {rec: {name: statistics.fmean(r[name] for r in grouped[repeat, ownership, rec])
                          for name in ['cpu_seconds', 'serial_cpu_seconds', 'dummy_less_empty_cpu_seconds']}
                    for rec in [False, True]}
            repeats.append(dict(repeat=repeat, **pair[True], recording_disabled_cpu_seconds=pair[False]['cpu_seconds'],
                relative_recording_delta=pair[True]['cpu_seconds']/pair[False]['cpu_seconds']-1 if pair[False]['cpu_seconds'] else None))
        prices[ownership] = dict(cpu_seconds=statistics.median(r['cpu_seconds'] for r in repeats),
            serial_cpu_seconds=statistics.median(r['serial_cpu_seconds'] for r in repeats),
            baseline_copy_bytes=(544 if reduction=='jagwas' else 32*(9+4*int(return_beta))) if ownership == 'owned' else 0, repeat_means=repeats,
            replaces_fixed_finish_and_tensor_conversion=True, includes_ready_cuda_event=False)
        if reduction is not None:
            prices[ownership].update(reduction=reduction,result_arrays=5,baseline_rows=32,return_df=False)
        if not return_beta:
            prices[ownership].update(return_beta=False,result_arrays=3,baseline_rows=32,return_df=False)
    return dict(workers=workers, prices=prices, meter_overheads=meters,
        scope='Exact fixed 32-row finish, corrected per observed detached interval and dummy-event CPU. Replaces fixed NumPy finish plus tensor slice/numpy control service. Ready CUDA-event CPU remains separate. Median complete repeat means; larger status/QC work, shape transfer and loaded-context differences remain unpriced.')


def acknowledgement_prices(blocked, ready, *, pairs, python_version, cpu_affinity):
    if type(pairs) is not int or pairs not in (1, 2):
        raise ValueError('Measured acknowledgement context required')
    if any(p['python_version'] != python_version or p['affinity'] != list(cpu_affinity) for p in [blocked, ready]):
        raise ValueError('Acknowledgement runtime context mismatch')
    contexts = [r for r in ready['results'] if r['pairs'] == pairs]
    if len(contexts) != 1 or sorted(w['worker'] for w in contexts[0]['workers']) != list(range(pairs)):
        raise ValueError('Complete ready acknowledgement workers required')
    groups = {}; coverage = []
    for worker in contexts[0]['workers']:
        seen = set()
        for row in worker['rows']:
            repeat = row['repeat']
            if type(repeat) is not int or repeat < 0 or repeat in seen or type(row['calls']) is not int or row['calls'] < 1:
                raise ValueError('Invalid ready acknowledgement coverage')
            seen.add(repeat)
            loop = number(row['loop_cpu_seconds'])
            groups.setdefault(repeat, []).append([number(row[key])-loop for key in ['publish_cpu_seconds', 'receive_cpu_seconds']])
        if len(seen) < 3:
            raise ValueError('Insufficient ready acknowledgement repeats')
        coverage.append(seen)
    if any(s != coverage[0] for s in coverage):
        raise ValueError('Ready acknowledgement worker coverage differs')
    publish, receive = [statistics.median(statistics.fmean(r[i] for r in rows) for rows in groups.values()) for i in range(2)]
    if min(publish, receive) < 0:
        raise ValueError('Negative aggregate ready acknowledgement service')
    repeats = []; seen = set()
    for context in blocked['acknowledgement']:
        if context['pairs'] != pairs:
            continue
        repeat = context['repeat']
        if type(repeat) is not int or repeat < 0 or repeat in seen or sorted(w['worker'] for w in context['workers']) != list(range(pairs)):
            raise ValueError('Invalid blocked acknowledgement coverage')
        seen.add(repeat)
        rows = [r for w in context['workers'] for r in w['rows']]
        if not rows or any(r.get('blocked_verified') is not True for r in rows):
            raise ValueError('Acknowledgement blocking was not verified')
        means = {key:statistics.fmean(number(row[key]) for row in rows) for key in
                 ['create_cpu_seconds','publisher_cpu_seconds','receiver_cpu_seconds','publication_to_return_seconds']}
        repeats.append(dict(repeat=repeat, **means))
    if len(repeats) < 3:
        raise ValueError('Insufficient blocked acknowledgement repeats')
    extra = statistics.median(r['publication_to_return_seconds']-publish-receive for r in repeats)
    return dict(create_cpu_seconds=statistics.median(r['create_cpu_seconds'] for r in repeats),
        publish_cpu_seconds=publish, receive_cpu_seconds=receive, wakeup_seconds=max(0., extra),
        signed_extra_wakeup_seconds=extra, repeat_means=repeats,
        scope='One Event creation, set and ready wait per delivered chunk. Extra elapsed wakeup applies only when readiness follows the wait attempt. Parking/resume CPU placement, timeout retries and context transfer remain unpriced; no scan timings used.')
