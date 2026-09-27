"""Changing undecoded ranges must preserve samples, results and ownership."""
from contextlib import contextmanager
from dataclasses import FrozenInstanceError
import os
import threading
import time
from unittest.mock import patch

import numpy as np
import pytest
import torch

from torchgwas.adaptive_chunks import ChunkSizeControl, AlignedChunkSizeControl, _ChunkDelivery
from torchgwas.api import load_genotype
from torchgwas.linear import linear_scan, linear_scan_streaming_chunks, linear_scan_multigpu
from torchgwas.preprocess import residualize_and_standardize
from torchgwas.reduce import JagwasReduction, SignificantPairs, device_significant_pairs
from torchgwas.streaming import PinnedDosageLoader
from test_native_scan import Source
from test_pgen_native_reader import write_pgen_mixed
from test_device_significance import selected_arrays


class Direct(Source):
    allows_direct_native_fill = True
    validate_native_range = True
    native_dtype = np.int8
    decode_workers = 2

    def __init__(self, values):
        super().__init__(values)
        self.reads = []
        self.sessions = []
        self.lock = threading.Lock()

    @contextmanager
    def native_reader_session(self):
        with self.lock:
            session = {'closed': False}
            self.sessions.append(session)
        def fill(start, end, out):
            with self.lock:
                self.reads.append((start, end))
            np.copyto(out, self.values[:, start:end].T)
        try:
            yield fill
        finally:
            session['closed'] = True


@pytest.fixture(autouse=True)
def backend(monkeypatch):
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '0')
    monkeypatch.setenv('TORCHGWAS_PGEN_PACKED', '0')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.setenv('TORCHGWAS_SCAN_PROFILE', '0')


@pytest.fixture
def sample():
    rng = np.random.default_rng(391)
    calls = rng.integers(0, 3, (129, 137), dtype=np.int8)
    y = rng.normal(size=(129, 5)).astype(np.float32)
    cov = rng.normal(size=(129, 2)).astype(np.float32)
    return calls, y, cov


def require_cuda(count=1):
    if torch.cuda.device_count() < count:
        pytest.skip(f'{count} CUDA devices required')


def no_pipeline_threads():
    assert not any(t.name.startswith(('torchgwas-pinned-', 'torchgwas-copy-release',
        'torchgwas-result', 'torchgwas-shard-')) for t in threading.enumerate())


def check_ranges(ranges, first, last):
    cursor = first
    for start, end in sorted(ranges):
        assert start == cursor and end > start
        cursor = end
    assert cursor == last


def test_delivery_clock_boundary_excludes_observer_and_abandoned_chunks():
    observations = []
    delivery = _ChunkDelivery(3, 5, 4, 'cuda:0', (1., 2.), observations.append)
    payload = (3, 5, np.zeros((2, 1), np.float32), None)
    with patch('torchgwas.adaptive_chunks.time.perf_counter', side_effect=[9., 10., 13., 14., 17.]):
        yielded = delivery.deliver(payload)
        assert next(yielded) is payload and not observations
        with pytest.raises(StopIteration):
            next(yielded)
        delivery.complete()
    assert observations[0].consumer_seconds == 3.
    assert observations[0].first_result == 10. and observations[0].completed == 17.
    assert observations[0].result_bytes == 8 and observations[0].result_blocks == 1
    with pytest.raises(FrozenInstanceError):
        observations[0].end = 100
    abandoned = _ChunkDelivery(5, 7, 4, 'cuda:0', (1., 2.), observations.append)
    yielded = abandoned.deliver(payload)
    next(yielded)
    yielded.close()
    assert len(observations) == 1 and abandoned.consumer_seconds == 0.
    abandoned.complete()
    assert observations[1].consumer_cpu_seconds is None
    assert observations[1].consumer_runnable_wait_seconds is None
    assert observations[1].consumer_probe_wall_seconds is None


def test_delivery_retains_reader_cpu_and_wait_diagnostics():
    observations = []
    delivery = _ChunkDelivery(3, 5, 4, 'cuda:0',
        (1., 2., .6, .2, .0001), observations.append)
    delivery.complete()
    row = observations[0]
    assert row.read_cpu_seconds == .6
    assert row.read_runnable_wait_seconds == .2
    assert row.read_probe_wall_seconds == .0001


def test_delivery_measures_resumed_consumer_cpu_and_optional_wait():
    observations=[]
    delivery=_ChunkDelivery(3,5,4,'cuda:0',(1.,2.),observations.append)
    payload=(3,5,np.zeros((2,1),np.float32),None)
    with patch('torchgwas.streaming._thread_runnable_wait_sample',
               side_effect=[(10_000_000,True),(12_000_000,True)]),\
         patch('torchgwas.adaptive_chunks.time.thread_time',side_effect=[1.,1.25]),\
         patch('torchgwas.adaptive_chunks.time.perf_counter',side_effect=[9.,10.,13.,14.,17.]):
        yielded=delivery.deliver(payload)
        assert next(yielded) is payload
        with pytest.raises(StopIteration):next(yielded)
        delivery.complete()
    row=observations[0]
    assert row.consumer_seconds==3.
    assert row.consumer_cpu_seconds==.25
    assert row.consumer_runnable_wait_seconds==.002
    assert row.consumer_probe_wall_seconds==2.


def test_delivery_does_not_treat_ambiguous_wait_as_zero():
    observations=[]
    delivery=_ChunkDelivery(3,5,4,'cuda:0',(1.,2.),observations.append)
    payload=(3,5,np.zeros((2,1),np.float32),None)
    with patch('torchgwas.streaming._thread_runnable_wait_sample',
               side_effect=[(10_000_000,False),(10_000_000,False)]):
        list(delivery.deliver(payload))
        delivery.complete()
    assert observations[0].consumer_cpu_seconds is not None
    assert observations[0].consumer_runnable_wait_seconds is None
    assert observations[0].consumer_probe_wall_seconds is not None


def test_scheduler_wait_probe_is_optional():
    from io import StringIO
    from torchgwas.streaming import _thread_runnable_wait_sample
    with patch('builtins.open', side_effect=OSError('unavailable')):
        assert _thread_runnable_wait_sample() == (None, None)
    with patch('builtins.open', side_effect=[StringIO('0\n'), StringIO('100 0 2\n')]):
        assert _thread_runnable_wait_sample() == (0, False)
    with patch('builtins.open', side_effect=[StringIO('0\n'), StringIO('100 4321 2\n')]):
        assert _thread_runnable_wait_sample() == (4321, False)


@pytest.mark.parametrize('sizes,initial', [([], 1), ([2, 2], 2), ([0, 2], 2),
    ([True, 2], 2), ([2., 3], 3), ([2, 3], 4)])
def test_control_refuses_invalid_choices(sizes, initial):
    with pytest.raises(ValueError):
        ChunkSizeControl(sizes, initial=initial)


@pytest.mark.parametrize('fault', ['cpu', 'dtype', 'generic', 'packed', 'device_decode',
    'missing_capacity', 'selector', 'observer', 'small_ring'])
def test_unsupported_path_refused_before_preprocessing(sample, fault):
    calls, y, cov = sample
    source = Direct(calls)
    options = dict(device='cuda:0', compute_dtype='float32', chunk_size=17,
                   _chunk_size_selector=ChunkSizeControl([3, 17], initial=3))
    if fault == 'cpu': options['device'] = 'cpu'
    if fault == 'dtype': options['compute_dtype'] = 'float64'
    if fault == 'generic': source = Source(calls)
    if fault == 'packed': source.native_encoding = 'pgen_2bit'
    if fault == 'device_decode': source.iter_device_chunks = lambda **kw: None
    if fault == 'missing_capacity': options['chunk_size'] = None
    if fault == 'selector': options['_chunk_size_selector'] = 1
    if fault == 'observer': options['_chunk_observer'] = 1
    if fault == 'small_ring': options['chunk_size'] = 7
    with patch('torchgwas.linear.choose_device', side_effect=lambda device: torch.device(device)), \
         patch('torchgwas.linear.residualize_and_standardize', side_effect=AssertionError('preprocessed')):
        with pytest.raises(ValueError):
            linear_scan_streaming_chunks(source, y, cov, **options)
    assert not source.passes


@pytest.mark.parametrize('borrowed', [False, True])
def test_feedback_changes_future_ranges_and_preserves_fp64_results(sample, borrowed):
    require_cuda()
    calls, y, cov = sample
    source = Direct(calls)
    control = ChunkSizeControl([3, 7, 17], initial=3)
    observations = []; issued_before_change = []
    def observe(row):
        observations.append(row)
        if len(observations) == 1:
            with source.lock:
                issued_before_change.extend(source.reads)
            control.set_size(17)
    with patch('torchgwas.linear.residualize_and_standardize', wraps=residualize_and_standardize) as prep:
        chunks, _ = linear_scan_streaming_chunks(source, y, cov, chunk_size=17,
            device='cuda:0', prefetch_chunks=2, variant_range=(5, 132),
            borrow_results=borrowed, compute_p_values=False,
            _chunk_size_selector=control, _chunk_observer=observe)
        saved = []
        for row in chunks:
            # Borrowed arrays must be consumed before next(); owned arrays may
            # be retained. Both contracts must survive changing row counts.
            saved.append(tuple(a.copy() if borrowed and isinstance(a, np.ndarray) else a for a in row))
        assert prep.call_count == 1
    check_ranges([(r[0], r[1]) for r in saved], 5, 132)
    assert sorted(source.reads) == [(r[0], r[1]) for r in saved]
    assert len(set(source.reads)) == len(source.reads)
    assert [(r.start, r.end) for r in observations] == [(r[0], r[1]) for r in saved]
    assert any(b-a == 3 for a, b in source.reads) and any(b-a == 17 for a, b in source.reads)
    assert set(issued_before_change) <= set(source.reads)
    assert len(source.sessions) == 1 and source.sessions[0]['closed']
    for obs, result in zip(observations, saved):
        assert obs.read_started <= obs.read_finished <= obs.submitted <= obs.first_result <= obs.completed
        assert obs.capacity == 17 and obs.device == 'cuda:0' and obs.result_blocks == 1
        assert obs.result_bytes == sum(a.nbytes for a in result[2:] if a is not None)
    ref = linear_scan(calls[:, 5:132].astype(np.float64), y, cov, device='cpu', compute_dtype='float64')
    for field in (2, 3):
        np.testing.assert_allclose(np.concatenate([r[field] for r in saved]), ref[field-2], rtol=3e-4, atol=3e-5)
    no_pipeline_threads()


@pytest.mark.parametrize('failure', ['zero', 'overflow', 'float', 'raised', 'observer', 'close'])
def test_invalid_runtime_choice_or_consumer_failure_closes_leases(sample, failure):
    require_cuda()
    calls, y, cov = sample
    source = Direct(calls)
    observations = []
    def select(start, stop, capacity):
        if start >= 7 and failure in ('zero', 'overflow', 'float', 'raised'):
            if failure == 'raised': raise RuntimeError('injected selector failure')
            return {'zero': 0, 'overflow': 18, 'float': 2.5}[failure]
        return 7
    def observe(row):
        if failure == 'observer': raise RuntimeError('injected observer failure')
        observations.append(row)
    chunks, _ = linear_scan_streaming_chunks(source, y, cov, chunk_size=17,
        device='cuda:0', prefetch_chunks=2, compute_p_values=False,
        _chunk_size_selector=select, _chunk_observer=observe)
    if failure == 'close':
        next(chunks)
        chunks.close()
        assert not observations  # Last yield never resumed.
    else:
        with pytest.raises((ValueError, RuntimeError)):
            list(chunks)
    assert source.sessions and all(row['closed'] for row in source.sessions)
    no_pipeline_threads()


def test_native_pgen_ld_replay_with_arbitrary_chunk_boundaries(tmp_path, sample):
    require_cuda()
    calls, y, cov = sample
    path = tmp_path/'adaptive.pgen'
    forms = [0 if i % 5 == 0 else 2 if i % 2 else 3 for i in range(calls.shape[1])]
    write_pgen_mixed(path, calls.T.astype(np.uint8), forms)
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(len(calls))))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(calls.shape[1])))
    source = load_genotype(path, genotype_format='pgen', pgen_mode='hardcall', reader_workers=2)[0]
    observations = []
    def select(start, stop, capacity):
        return (1, 7, 17, 3)[start % 4]
    chunks, _ = linear_scan_streaming_chunks(source, y, cov, chunk_size=17,
        device='cuda:0', prefetch_chunks=2, variant_range=(3, 136),
        compute_p_values=False, _chunk_size_selector=select, _chunk_observer=observations.append)
    rows = list(chunks)
    check_ranges([(r[0], r[1]) for r in rows], 3, 136)
    ref = linear_scan(calls[:, 3:136].astype(np.float64), y, cov, device='cpu', compute_dtype='float64')
    for field in (2, 3):
        np.testing.assert_allclose(np.concatenate([r[field] for r in rows]), ref[field-2], rtol=3e-4, atol=3e-5)
    assert len(observations) == len(rows)
    no_pipeline_threads()


@pytest.mark.parametrize('threshold', [1., .2, 1e-30])
def test_device_selection_reports_once_per_source_chunk_and_preserves_all_pairs(sample, monkeypatch, threshold):
    require_cuda()
    monkeypatch.setenv('TORCHGWAS_SIGNIFICANCE_BACKEND', 'device')
    calls, y, cov = sample
    observations = []
    def small_blocks(*args, **kwargs):
        return device_significant_pairs(*args, **kwargs, max_cells=3)
    with patch('torchgwas.reduce.device_significant_pairs', side_effect=small_blocks):
        chunks, _ = linear_scan_streaming_chunks(Direct(calls), y, cov, chunk_size=17,
            device='cuda:0', prefetch_chunks=2, variant_range=(5, 38), compute_p_values=False,
            significance=SignificantPairs(threshold), significance_n_traits=y.shape[1],
            _chunk_size_selector=lambda start, stop, capacity: 7 if start % 2 else 3,
            _chunk_observer=observations.append)
        rows = list(chunks)
    check_ranges([(r.start, r.end) for r in observations], 5, 38)
    assert sum(r.result_blocks for r in observations) == len(rows)
    # Block count follows the shared selection geometry (3-row chunks pack into 5 blocks, not 6).
    from torchgwas.selection_geometry import device_selection_shape
    assert all(r.result_blocks == device_selection_shape(r.end-r.start, y.shape[1], 3)[2] for r in observations)
    assert sum(r.result_bytes for r in observations) == sum(sum(a.nbytes for a in row[2:]) for row in rows)
    ref = linear_scan(calls[:, 5:38].astype(np.float64), y, cov, device='cpu', compute_dtype='float64')
    indices = np.nonzero(np.isfinite(ref[1]) & (ref[2] <= threshold))
    got = selected_arrays(rows)
    np.testing.assert_array_equal(got[0], indices[0]+5)
    np.testing.assert_array_equal(got[1], indices[1])
    for value, want in zip(got[2:4], [ref[0][indices], ref[1][indices]]):
        np.testing.assert_allclose(value, want, rtol=3e-4, atol=3e-5)
    assert all('no durable-write guarantee' in r.boundary for r in observations)
    no_pipeline_threads()


def test_joint_two_gpu_adaptive_chunks_keep_full_factors_and_exact_coverage(sample):
    require_cuda(2)
    calls, y, cov = sample
    source = Direct(calls)
    control = ChunkSizeControl([3, 17], initial=3)
    observations = []; factors = []; lock = threading.Lock()
    def observe(row):
        with lock:
            observations.append(row)
        control.set_size(17)
    original = JagwasReduction.prepare
    def prepare(self, phenotype, device=None):
        result = original(self, phenotype, device=device)
        with lock:
            factors.append((str(self._inverse_cholesky.device), tuple(self._inverse_cholesky.shape)))
        return result
    with patch.object(JagwasReduction, 'prepare', prepare), \
         patch('torchgwas.linear.residualize_and_standardize', wraps=residualize_and_standardize) as prep:
        chunks, _ = linear_scan_multigpu(source, y, cov, devices=['cuda:0', 'cuda:1'],
            chunk_size=17, reader_workers=2, prefetch_chunks=2, compute_p_values=False,
            ordered=False, shared_queue_depth=1, reduction_factory=JagwasReduction,
            variant_range=(5, 132), _chunk_size_selector=control, _chunk_observer=observe)
        rows = sorted(list(chunks), key=lambda row: row[0])
        assert prep.call_count == 1
    check_ranges([(r[0], r[1]) for r in rows], 5, 132)
    assert sorted((r.start, r.end) for r in observations) == [(r[0], r[1]) for r in rows]
    assert sorted(source.reads) == [(r[0], r[1]) for r in rows]
    assert sorted(factors) == [('cuda:0', (5, 5)), ('cuda:1', (5, 5))]
    assert len(source.sessions) == 2 and all(row['closed'] for row in source.sessions)
    ref = linear_scan(calls[:, 5:132].astype(np.float64), y, cov, device='cpu', compute_dtype='float64')
    x = np.column_stack([np.ones(len(y)), cov.astype(np.float64)])
    yr = y - x @ np.linalg.lstsq(x, y, rcond=None)[0]
    correlation = np.corrcoef(yr.T)
    z = ref[1] / np.sqrt(1 + ref[1] ** 2 / (len(y) - np.linalg.matrix_rank(x) - 1))  # the score form
    expected = np.einsum('ij,ij->i', z, np.linalg.solve(correlation, z.T).T)
    np.testing.assert_allclose(np.concatenate([r[3][:, 0] for r in rows]), expected, rtol=3e-4, atol=3e-5)
    no_pipeline_threads()


def test_aligned_switches_use_only_declared_shapes_on_both_gpus(sample):
    require_cuda(2)
    from torchgwas.adaptive_chunks import aligned_chunk_shapes
    from torchgwas.linear import multigpu_variant_ranges
    calls,y,cov=sample;source=Direct(calls)
    control=AlignedChunkSizeControl([4,8,16],initial=4)
    observations=[];lock=threading.Lock()
    def observe(row):
        with lock:
            observations.append(row)
            control.set_size(16 if len(observations)>=2 else 8)
    chunks,_=linear_scan_multigpu(source,y,cov,devices=['cuda:0','cuda:1'],
        chunk_size=16,reader_workers=2,prefetch_chunks=2,compute_p_values=False,
        ordered=False,shared_queue_depth=1,reduction_factory=JagwasReduction,
        _chunk_size_selector=control,_chunk_observer=observe)
    rows=sorted(list(chunks),key=lambda row:row[0])
    check_ranges([(row[0],row[1]) for row in rows],0,137)
    spans=multigpu_variant_ranges(137,16,2)
    for index,(first,last) in enumerate(spans):
        seen=[row for row in observations if row.device=='cuda:'+str(index)]
        assert seen and all(row.capacity==16 for row in seen)
        shapes=aligned_chunk_shapes([4,8,16],last-first)
        assert all(row.end-row.start in shapes for row in seen)
        check_ranges([(row.start,row.end) for row in seen],first,last)
    assert {row.end-row.start for row in observations}>={4,16}
    ref=linear_scan(calls.astype(np.float64),y,cov,device='cpu',compute_dtype='float64')
    x=np.column_stack([np.ones(len(y)),cov.astype(np.float64)])
    yr=y-x@np.linalg.lstsq(x,y,rcond=None)[0];correlation=np.corrcoef(yr.T)
    z=ref[1]/np.sqrt(1+ref[1]**2/(len(y)-np.linalg.matrix_rank(x)-1))  # the score form
    expected=np.einsum('ij,ij->i',z,np.linalg.solve(correlation,z.T).T)
    np.testing.assert_allclose(np.concatenate([row[3][:,0] for row in rows]),expected,rtol=3e-4,atol=3e-5)
    assert all(row['closed'] for row in source.sessions)
    no_pipeline_threads()


@pytest.mark.parametrize('cuda_events', [False, True])
@pytest.mark.parametrize('mode', ['full', 'jagwas', 'significant'])
def test_initial_window_spreads_bounded_measurements_then_saves_for_later_job(
        sample, tmp_path, monkeypatch, cuda_events, mode):
    from torchgwas.adaptive_chunks import InitialChunkMeasurements
    from torchgwas.calibration_cache import CalibrationParameterCache
    require_cuda()
    calls, y, cov = sample
    source = Direct(calls)
    window = InitialChunkMeasurements(['cuda:0'], max_chunks_per_device=4,
        warmup_chunks=2, stride=3, max_window_seconds=100., cuda_events=cuda_events)
    options = {}
    if mode == 'jagwas': options['reduction'] = JagwasReduction()
    if mode == 'significant':
        monkeypatch.setenv('TORCHGWAS_SIGNIFICANCE_BACKEND', 'device')
        options['significance'] = SignificantPairs(.2)
    chunks, _ = linear_scan_streaming_chunks(source, y, cov, chunk_size=7,
        device='cuda:0', prefetch_chunks=2, compute_p_values=False,
        _chunk_size_selector=ChunkSizeControl([3, 7], initial=3),
        _chunk_observer=window, **options)
    rows = list(chunks)
    snapshot = window.snapshot()
    assert snapshot['reader_probe_protocol'] == 'thread_cpu_schedstat_v1'
    assert snapshot['consumer_probe_protocol'] == 'thread_cpu_schedstat_probe_wall_v2'
    assert snapshot['new_measurements_stopped']
    assert snapshot['reserved'] == {'cuda:0': 4} and snapshot['pending'] == []
    assert snapshot['source_chunks_seen'] == {'cuda:0': 46}
    expected = sorted(source.reads)[2:12:3]
    assert [(r['start'], r['end']) for r in snapshot['observations']] == expected
    for row in snapshot['observations']:
        assert row['read_cpu_seconds'] >= 0
        assert row['read_probe_wall_seconds'] >= 0
        assert (row['read_runnable_wait_seconds'] is None
                or row['read_runnable_wait_seconds'] >= 0)
        assert row['consumer_cpu_seconds'] >= 0
        assert (row['consumer_runnable_wait_seconds'] is None
                or row['consumer_runnable_wait_seconds'] >= 0)
        assert row['consumer_probe_wall_seconds'] >= 0
        if cuda_events:
            assert row['cuda']['h2d'] >= 0 and row['cuda']['conversion'] >= 0
            assert row['cuda']['statistics_and_reduction'] > 0
            if mode == 'significant': assert row['cuda']['result_transfer'] is None
            else: assert row['cuda']['result_transfer'] >= 0
        else:
            assert row['cuda'] is None
    # The scan includes every range after the measurement budget ends.
    check_ranges([(r[0], r[1]) for r in rows], 0, calls.shape[1])
    cache = CalibrationParameterCache(tmp_path/'cache')
    dependencies = dict(protocol='test-events-v1', device='cuda:0', mode=mode, source='test-direct')
    published = window.publish(cache, dependencies=dependencies,
        provenance=dict(job='first', numeric_checks='separate correctness tests'), max_age_seconds=60.)
    later = CalibrationParameterCache(tmp_path/'cache')
    hit = later.lookup('stage_observations', 'initial_chunks', dependencies=dependencies)
    assert hit['hit'] and hit['record_sha256'] == published['record_sha256']
    assert len(hit['record']['value']['observations']) == 4
    no_pipeline_threads()


def test_window_budgets_reservations_not_only_completed_measurements():
    from torchgwas.adaptive_chunks import InitialChunkMeasurements
    window = InitialChunkMeasurements(['cuda:0', 'cuda:1'], max_chunks_per_device=2,
        warmup_chunks=0, stride=1, max_window_seconds=3.)
    with patch('torchgwas.adaptive_chunks.time.perf_counter', return_value=10.):
        assert window.reserve_read(0, 3, 'cuda:0')
        assert window.reserve_read(3, 6, 'cuda:0')
        assert not window.reserve_read(6, 9, 'cuda:0')
        assert window.reserve_read(9, 12, 'cuda:1')
    with patch('torchgwas.adaptive_chunks.time.perf_counter', return_value=13.):
        assert not window.reserve_read(12, 15, 'cuda:1')
        report = window.snapshot()
        assert report['new_measurements_stopped'] and len(report['pending']) == 3
        assert report['reserved'] == {'cuda:0': 2, 'cuda:1': 1}
    with pytest.raises(ValueError, match='device'):
        window.reserve_read(12, 15, 'cuda:2')


def test_explicit_stop_preserves_inflight_measurement():
    from torchgwas.adaptive_chunks import InitialChunkMeasurements, ChunkObservation
    window = InitialChunkMeasurements(['cuda:0'], warmup_chunks=0, stride=1)
    assert window.reserve_read(0, 3, 'cuda:0')
    window.stop()
    assert not window.reserve_read(3, 6, 'cuda:0')
    row = ChunkObservation(0, 3, 7, 'cuda:0', 1., 2., 3., 4., 5., 1., 1, 60)
    window(row)
    assert window.snapshot()['pending'] == []
    assert len(window.snapshot()['observations']) == 1
    with pytest.raises(ValueError, match='reservation'):
        window(row)


def test_pinned_slots_are_sized_on_first_use_and_grow_once_to_the_capacity():
    require_cuda()  # page-locked host memory
    calls = np.random.default_rng(5).integers(0, 3, (17, 400), dtype=np.int8)
    control = AlignedChunkSizeControl([32, 64, 128], initial=32)
    loader = PinnedDosageLoader(Direct(calls), 128, depth=3, reader_workers=2, chunk_size_selector=control)
    assert loader.buffers == [None] * 3 and loader.slot_shape == (128, 17)
    seen = []
    try:
        for index, host, start, end in loader:
            np.testing.assert_array_equal(host.numpy(), calls[:, start:end].T)
            seen.append((end - start, loader.buffers[index].shape[0]))
            if end >= 96:
                control.set_size(64)
            loader.release(index)
    finally:
        loader.close()
    # The first three chunks (one per slot) are 32 rows, so every slot is
    # pinned at 32 and grows once, to the 128-row capacity, for a 64-row chunk.
    assert seen[:3] == [(32, 32)] * 3
    assert {rows for size, rows in seen if size == 64} == {128}
    assert sum(size for size, _ in seen) == 400
