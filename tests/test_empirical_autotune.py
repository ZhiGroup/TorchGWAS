"""Empirical autotune: segment trials pick the fastest chunk, layout follows memory and CPUs."""
from dataclasses import dataclass
import contextlib
import json
import random

import numpy as np
import pytest
import torch

from torchgwas.empirical_autotune import EmpiricalChunkTuner, plan_layout


@dataclass
class Obs:
    start: int
    end: int
    device: str
    completed: float


def simulate(tuner, cost, *, devices=('cuda:0',), rows=200_000, noise=0.0, seed=0, drift=0.0):
    """Devices pull chunks from their own ranges; per-row time depends on size."""
    rng = random.Random(seed)
    span = rows//len(devices)
    cursors = {d: i*span for i, d in enumerate(devices)}
    ends = {d: (i+1)*span for i, d in enumerate(devices)}
    clock = {d: 0. for d in devices}
    while any(cursors[d] < ends[d] for d in devices):
        device = min((d for d in devices if cursors[d] < ends[d]), key=clock.get)
        size = tuner.control(cursors[device], ends[device], tuner.capacity)
        start = cursors[device]
        cursors[device] += size
        load = 1+drift*clock[device]
        per_row = cost.get(size) or cost[min(k for k in cost if k >= size)]  # tails cost like the smallest fit
        clock[device] += size*per_row*load*(1+noise*rng.uniform(-1, 1))*1e-6
        tuner(Obs(start, start+size, device, clock[device]))
    return max(clock.values())


def test_trials_pick_fastest_size_and_commit_early():
    tuner = EmpiricalChunkTuner([512, 1024, 2048, 4096], total_rows=400_000, min_job_seconds=0)
    simulate(tuner, {512: 3.0, 1024: 2.0, 2048: 1.5, 4096: 2.5}, rows=400_000, noise=.05)
    audit = tuner.audit()
    assert audit['state'] == 'committed' and audit['choice'] == 2048
    assert audit['decided_fraction'] < .40
    assert [s['size'] for s in audit['segments']] == [512, 1024, 2048, 4096, 4096, 2048, 1024, 512]
    assert all(s['rows'] >= s['target'] for s in audit['segments'])


def test_linear_drift_cancels_with_forward_reverse_order():
    # Load doubles over the job; sizes cost the same. Forward-only order would
    # favour the first size; forward/reverse keeps the starting size.
    tuner = EmpiricalChunkTuner([512, 1024, 2048], total_rows=300_000, initial=1024, min_gain=.03, min_job_seconds=0)
    simulate(tuner, {512: 2.0, 1024: 2.0, 2048: 2.0}, rows=300_000, drift=2.0)
    assert tuner.audit()['choice'] == 1024


def test_small_gain_keeps_initial_and_short_job_skips():
    tuner = EmpiricalChunkTuner([512, 1024, 2048, 4096], total_rows=400_000, min_gain=.05, min_job_seconds=0)
    simulate(tuner, {512: 2.0, 1024: 2.0, 2048: 1.97, 4096: 2.0}, rows=400_000)
    assert tuner.audit()['choice'] == 1024 and tuner.reason == 'within_min_gain_of_initial'
    short = EmpiricalChunkTuner([512, 1024, 2048, 4096], total_rows=20_000)
    # Warmup throughput predicting a job under min_job_seconds also skips trials.
    brief = EmpiricalChunkTuner([512, 1024, 2048, 4096], total_rows=400_000, min_job_seconds=20)
    simulate(brief, {512: 2.0, 1024: 2.0, 2048: 1.0, 4096: 2.0}, rows=400_000)
    assert brief.state == 'skipped' and brief.reason == 'job_shorter_than_min_seconds' and not brief.segments
    assert short.state == 'skipped' and short.control(0, 20_000, 4096) == 1024
    # A medium job drops the largest sizes from its trials instead of skipping.
    medium = EmpiricalChunkTuner([512, 1024, 2048, 4096], total_rows=150_000, depth=4, concurrent=2)
    assert medium.state == 'warmup' and medium.trial_sizes[-1] < 4096 and 'dropped' in medium.reason
    assert medium.planned_rows <= 75_000


def test_deferred_trials_start_when_the_job_slows_down():
    # The first 10% runs 50x faster than the rest (a load change, as seen on
    # the H100): the warmup estimate is far too short, so trials are deferred,
    # then start once the measured rate shows enough remaining time.
    tuner = EmpiricalChunkTuner([512, 1024, 2048], total_rows=1_000_000, min_job_seconds=5)
    clock, cursor = 0., 0
    while cursor < 1_000_000:
        size = tuner.control(cursor, 1_000_000, tuner.capacity)
        per_row = {512: 3.0, 1024: 2.0, 2048: 1.0}.get(size, 3.0)*(1 if cursor < 100_000 else 50)
        clock += size*per_row*1e-6
        tuner(Obs(cursor, cursor+size, 'cuda:0', clock))
        cursor += size
    audit = tuner.audit()
    assert audit['estimated_job_seconds'] < 5 and audit['state'] == 'committed'
    assert audit['choice'] == 2048 and audit['reason'] == 'highest_measured_throughput'


def test_segments_wait_for_every_device_after_a_switch():
    tuner = EmpiricalChunkTuner([512, 1024, 2048], total_rows=600_000, min_job_seconds=0)
    simulate(tuner, {512: 3.0, 1024: 1.0, 2048: 2.0}, devices=('cuda:0', 'cuda:1'), rows=600_000)
    audit = tuner.audit()
    assert audit['devices'] == ['cuda:0', 'cuda:1'] and audit['choice'] == 1024
    assert all(s['other_size_rows'] == 0 for s in audit['segments'])


GIB = 1 << 30
BASE = dict(n_samples=35_365, covariate_rank=27, n_variants=1_048_576, capacity=4096, depth=4,
            transfer_bytes_per_variant=35_365., device_free_bytes=80*GIB, host_free_bytes=500*GIB)


def test_layout_rules():
    two = ['cuda:0', 'cuda:3']
    wide = plan_layout(mode='significant', n_traits=16_385, devices=two, cpus=16, **BASE)
    assert wide['trait_block'] == 8193 and wide['trait_devices'] == two
    # 4 readers per GPU, and a double-buffered ring of 8 slots per GPU when it fits.
    assert wide['reader_workers'] == 8 and wide['prefetch_chunks'] == 8
    busy = plan_layout(mode='significant', n_traits=16_385, devices=two, cpus=6, **BASE)
    # Without a decode measurement, 4 readers per GPU: a busy host is not a reader cap (time-slicing).
    assert busy['reader_workers'] == 8 and busy['prefetch_chunks'] == 8
    # A ring that fits only single-buffered keeps depth = readers.
    tight = plan_layout(mode='significant', n_traits=16_385, devices=two, cpus=16,
                        **dict(BASE, host_free_bytes=5*GIB))
    assert tight['reader_workers'] == 8 and tight['prefetch_chunks'] == 4
    many = [f'cuda:{i}' for i in range(8)]
    # Tiles reread the genotype: two tiles at K=8,192 even with 8 idle GPUs,
    # more only as the per-tile trait work grows.
    assert plan_layout(mode='significant', n_traits=8192, devices=many, cpus=40, **BASE)['trait_devices'] == many[:2]
    assert len(plan_layout(mode='significant', n_traits=100_000, devices=many, cpus=40, **BASE)['trait_devices']) == 4
    narrow = plan_layout(mode='significant', n_traits=300, devices=two, cpus=16, **BASE)
    assert narrow['trait_block'] is None and narrow['devices_used'] == 1
    jagwas = plan_layout(mode='jagwas', n_traits=5000, devices=two, cpus=16, **BASE)
    assert jagwas['variant_devices'] == two and jagwas['trait_block'] is None
    full = plan_layout(mode='full', n_traits=5000, devices=two, cpus=16, **BASE)
    assert full['variant_devices'] == two
    # Narrow panels: JAGWAS shards over at most 2 GPUs, dense output stays on 1.
    assert plan_layout(mode='jagwas', n_traits=300, devices=many, cpus=40, **BASE)['variant_devices'] == many[:2]
    assert plan_layout(mode='full', n_traits=64, devices=two, cpus=16, **BASE)['variant_devices'] is None
    # A saturated host trims GPUs, but never below two.
    cpu_bound = plan_layout(mode='jagwas', n_traits=5000, devices=['cuda:0', 'cuda:1', 'cuda:2'], cpus=1, **BASE)
    assert cpu_bound['variant_devices'] == ['cuda:0', 'cuda:1'] and cpu_bound['devices_used'] == 2
    tight = plan_layout(mode='significant', n_traits=2_000_000, devices=two, cpus=16,
                        **dict(BASE, device_free_bytes=40*GIB))
    assert tight['trait_block'] < 1_000_000 and len(tight['trait_devices']) == 2
    # Memory-limited tiles come in whole rounds: every GPU scans as many tiles.
    assert -(-2_000_000//tight['trait_block']) % 2 == 0
    four = plan_layout(mode='significant', n_traits=1_000_000, devices=many[:4], cpus=40,
                       **dict(BASE, device_free_bytes=24*GIB))
    assert -(-1_000_000//four['trait_block']) % 4 == 0 and len(four['trait_devices']) == 4
    # Chunk sizes whose rings do not fit are dropped before planning (float32
    # transfer at 600,000 samples on an 11 GiB card holds 512 but not 1,024).
    big = dict(BASE, n_samples=600_000, covariate_rank=10, transfer_bytes_per_variant=2.4e6,
               device_free_bytes=11*GIB, chunk_sizes=[512, 1024, 2048])
    small = plan_layout(mode='significant', n_traits=8192, devices=two, cpus=16, **big)
    assert small['chunk_sizes'] == [512] and small['chunk_sizes_dropped'] == [1024, 2048]
    with pytest.raises(ValueError, match='Not even chunk 512'):
        plan_layout(mode='significant', n_traits=8192, devices=two, cpus=16, **dict(big, device_free_bytes=2*GIB))
    with pytest.raises(ValueError, match='JAGWAS cannot tile'):
        plan_layout(mode='jagwas', n_traits=2_000_000, devices=two, cpus=16, **dict(BASE, device_free_bytes=4*GIB))
    # JAGWAS needs the full panel and its K x K factor on every GPU; no split.
    with pytest.raises(ValueError, match='JAGWAS factor'):
        plan_layout(mode='jagwas', n_traits=30_000, devices=many, cpus=40, **dict(BASE, device_free_bytes=8*GIB))
    with pytest.raises(ValueError, match='JAGWAS factor'):
        plan_layout(mode='jagwas', n_traits=30_000, devices=['cuda:0'], cpus=16, **dict(BASE, device_free_bytes=8*GIB))
    filtered = plan_layout(mode='full', n_traits=300, devices=two, cpus=16, allow_partitions=False, **BASE)
    assert filtered['variant_devices'] is None and filtered['trait_block'] is None


# ------------------------------------------------------------ real GPU runs

def _fixture(tmp_path, n=129, m=8192, k=40):
    from test_pgen_native_reader import write_pgen
    rng = np.random.default_rng(20260923)
    calls = rng.integers(0, 3, (n, m)).astype(np.uint8)
    y = rng.normal(size=(n, k)).astype(np.float32)
    y[:, :4] += .6*calls[:, :4].astype(np.float32)  # a few real signals
    c = rng.normal(size=(n, 2)).astype(np.float32)
    path = tmp_path/'input.pgen'
    write_pgen(path, calls.T.copy())
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    np.save(tmp_path/'y.npy', y); np.save(tmp_path/'c.npy', c)
    return path


@pytest.fixture
def native(monkeypatch):
    if not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    for key, value in dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0',
                           TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_SCAN_PROFILE='0').items():
        monkeypatch.setenv(key, value)


# Three sizes and a shallow ring keep the trial schedule, including the rows
# that drain after each switch, inside half of an 8,192-variant job.
OPTIONS = dict(chunk_sizes=[16, 32, 64], warmup_fraction=.02, trial_fraction=.2, min_job_seconds=0)
RING = dict(prefetch_chunks=2, reader_workers=4)
# cpus_per_device=1 in multi-GPU tests: they test partitioning, not how busy the host is.


def _run(tmp_path, name, path, **kwargs):
    from torchgwas.api import run_linear_gwas
    out = tmp_path/name
    result = run_linear_gwas(genotype=str(path), phenotype=tmp_path/'y.npy', covariates=tmp_path/'c.npy',
                             pgen_mode='hardcall', compute_dtype='float32', output_dir=out, **kwargs)
    return out, result, json.loads((out/'run.json').read_text())


def test_full_output_autotune_matches_fixed_chunks(tmp_path, native):
    from torchgwas.sumstats import open_binary_sumstats
    path = _fixture(tmp_path)
    _, _, _ = _run(tmp_path, 'fixed', path, device='cuda:0', chunk_size=32)
    _, _, run = _run(tmp_path, 'tuned', path, device='cuda:0', autotune=True, autotune_options=OPTIONS, **RING)
    tuned = run['autotune']
    assert tuned['method'] == 'empirical' and tuned['tuner_disabled_reason'] is None
    assert tuned['chunk']['state'] == 'committed' and tuned['chunk']['choice'] in (16, 32, 64)
    b0, t0, _logp, _ = open_binary_sumstats(tmp_path/'fixed'/'sumstats')
    b1, t1, _logp, _ = open_binary_sumstats(tmp_path/'tuned'/'sumstats')
    np.testing.assert_allclose(np.asarray(t1), np.asarray(t0), rtol=2e-5, atol=2e-5)
    np.testing.assert_allclose(np.asarray(b1), np.asarray(b0), rtol=2e-5, atol=2e-5)


def _pairs(directory, field='t_stat'):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    _, parts = open_indexed_sumstats(directory/'sumstats')
    rows = {}
    for part in parts:
        keys = zip(part['variant_index'], part['trait_index']) if 'trait_index' in part else (
            (v, None) for v in part['variant_index'])
        for key, value in zip(keys, part[field]):
            assert key not in rows
            rows[tuple(int(x) if x is not None else -1 for x in key)] = float(value)
    return rows


def test_significant_two_gpu_tiles_match_single_fixed(tmp_path, native):
    if torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    path = _fixture(tmp_path)
    common = dict(reduce='significant', significance_threshold=1e-3)
    _run(tmp_path, 'fixed', path, device='cuda:0', chunk_size=32, **common)
    _, _, run = _run(tmp_path, 'tuned', path, autotune=True,
                     autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8, cpus_per_device=1,
                                           split='traits'), **RING, **common)
    layout = run['autotune']['layout']
    assert layout['trait_block'] == 20 and layout['trait_devices'] == ['cuda:0', 'cuda:1']
    # Both tiles were fed by one decode pass.
    assert run['shared_decode']['enabled'] and run['shared_decode']['subscribers'] == 2
    assert run['shared_decode']['chunks'] > 0
    # gpu_fanout='auto': every GPU copies from pinned host memory (fan-out did
    # not pay end to end on lab-h100 or lab-a100); fan-out is opt-in.
    assert run['shared_decode']['transfer'] == 'pcie_per_gpu'
    assert run['autotune']['chunk']['state'] == 'committed'
    expected, actual = _pairs(tmp_path/'fixed'), _pairs(tmp_path/'tuned')
    assert expected and actual.keys() == expected.keys()
    np.testing.assert_allclose([actual[k] for k in expected], list(expected.values()), rtol=2e-4, atol=2e-4)


def test_significant_tiles_with_forced_gpu_fanout_match_single_fixed(tmp_path, native):
    # Fan-out copies each chunk once to the first tile GPU and the second tile
    # pulls it GPU to GPU (NVLink where present; PCIe peer or staged otherwise).
    if torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    path = _fixture(tmp_path)
    common = dict(reduce='significant', significance_threshold=1e-3)
    _run(tmp_path, 'fixed', path, device='cuda:0', chunk_size=32, **common)
    _, _, run = _run(tmp_path, 'fanout', path, autotune=True,
                     autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8,
                                           cpus_per_device=1, split='traits', gpu_fanout='peer'), **RING, **common)
    assert run['shared_decode']['transfer'] == 'pcie_to_root_then_gpu_peer'
    assert run['shared_decode']['fanout_device'] == 'cuda:0'
    # 'nvlink' fans out only where NVLink connects the tile GPUs.
    from torchgwas.shared_decode import nvlink_root
    _, _, nv = _run(tmp_path, 'nvlink', path, autotune=True,
                    autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8,
                                          cpus_per_device=1, split='traits', gpu_fanout='nvlink'), **RING, **common)
    linked = nvlink_root(['cuda:0', 'cuda:1']) is not None
    assert nv['shared_decode']['transfer'] == ('pcie_to_root_then_gpu_peer' if linked else 'pcie_per_gpu')
    # 'uplink' additionally needs the two GPUs to share a PCIe root port.
    from torchgwas.shared_decode import fanout_root
    _, _, up = _run(tmp_path, 'uplink', path, autotune=True,
                    autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8,
                                          cpus_per_device=1, split='traits', gpu_fanout='uplink'), **RING, **common)
    shares = fanout_root(['cuda:0', 'cuda:1']) is not None
    assert up['shared_decode']['transfer'] == ('pcie_to_root_then_gpu_peer' if shares else 'pcie_per_gpu')
    assert _pairs(tmp_path/'nvlink').keys() == _pairs(tmp_path/'fixed').keys()
    expected, actual = _pairs(tmp_path/'fixed'), _pairs(tmp_path/'fanout')
    assert expected and actual.keys() == expected.keys()
    np.testing.assert_allclose([actual[k] for k in expected], list(expected.values()), rtol=2e-4, atol=2e-4)


def test_jagwas_two_gpu_variant_shards_match_single_fixed(tmp_path, native):
    if torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    path = _fixture(tmp_path)
    _run(tmp_path, 'fixed', path, device='cuda:0', chunk_size=32, reduce='jagwas')
    # A larger probe share: upward probing with depth settling and the
    # start-size revisit costs more rows than 20% of this small job.
    _, _, run = _run(tmp_path, 'tuned', path, reduce='jagwas', autotune=True,
                     autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8, cpus_per_device=1,
                                           shard_setup_seconds=0, trial_fraction=0.9), **RING)
    assert run['autotune']['layout']['variant_devices'] == ['cuda:0', 'cuda:1']
    assert run['autotune']['chunk']['state'] == 'committed', (run['autotune']['chunk']['reason'], run['autotune']['chunk'].get('method'))
    expected, actual = _pairs(tmp_path/'fixed', 'chi2'), _pairs(tmp_path/'tuned', 'chi2')
    assert actual.keys() == expected.keys() and len(expected) == 8192
    np.testing.assert_allclose([actual[k] for k in expected], list(expected.values()), rtol=2e-4, atol=2e-4)


def test_autotune_argument_validation(tmp_path):
    from torchgwas.api import run_linear_gwas
    with pytest.raises(ValueError, match='requires output_dir'):
        run_linear_gwas(genotype='x.pgen', phenotype=np.zeros((2, 1)), autotune=True)
    with pytest.raises(ValueError, match='Unknown autotune_options'):
        run_linear_gwas(genotype='x.pgen', phenotype=np.zeros((2, 1)), autotune=True,
                        output_dir=tmp_path, autotune_options={'bogus': 1})
    with pytest.raises(ValueError, match='requires autotune'):
        run_linear_gwas(genotype='x.pgen', phenotype=np.zeros((2, 1)), autotune_options={})


def test_cpu_cap_from_measured_demand():
    from torchgwas.empirical_autotune import plan_layout
    seven = [f'cuda:{i}' for i in range(1, 8)]
    # JAGWAS K=8,192 on lab-2080ti: ~20 us decode CPU against ~175 us of GPU
    # GEMM per variant, so 5 spare cores feed all seven GPUs.
    heavy = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=6,
                        decode_cpu_per_variant=20e-6, gpu_seconds_per_variant=175e-6, **BASE)
    assert heavy['variant_devices'] == seven
    assert heavy['cpu_demand']['gpus_by_cpu'] >= 7 and heavy['cpu_demand']['cores_per_gpu'] < 0.5
    # Decode-bound (GPU time far below decode): the floor of two GPUs.
    light = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=6,
                        decode_cpu_per_variant=20e-6, gpu_seconds_per_variant=1e-6, **BASE)
    assert light['variant_devices'] == seven[:2]
    assert any('cores per GPU' in reason for reason in light['why'])
    # Balanced: (6 - 1) / (20/31.4 + 0.25) = 5 GPUs.
    middle = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=6,
                         decode_cpu_per_variant=20e-6, gpu_seconds_per_variant=31.4e-6, **BASE)
    assert middle['variant_devices'] == seven[:5]
    # An explicit cpus_per_device keeps the fixed rule; no measurement keeps the default of 3.
    fixed = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=6, cpus_per_device=3,
                        decode_cpu_per_variant=20e-6, gpu_seconds_per_variant=175e-6, **BASE)
    assert fixed['variant_devices'] == seven[:2] and 'cpu_demand' not in fixed
    assert plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=6, **BASE)['variant_devices'] == seven[:2]


def test_decode_probe_uses_native_fill():
    from torchgwas.empirical_autotune import decode_cpu_seconds_per_variant

    class Source:
        allows_direct_native_fill = True
        native_transfer_dtype = np.uint8
        native_row_width = 128
        shape = (500, 10_000)
        calls = []

        @contextlib.contextmanager
        def native_reader_session(self):
            def fill(start, end, out):
                assert out.shape == (end - start, 128) and out.flags.c_contiguous and out.dtype == np.uint8
                self.calls.append((start, end))
                out[:] = 1
                sum(range(20_000))  # some CPU per call
            yield fill

    source = Source()
    seconds = decode_cpu_seconds_per_variant(source, 100, 10_000, variants=512)
    assert seconds is not None and seconds >= 0
    assert source.calls == [(100, 356), (356, 612)]
    assert decode_cpu_seconds_per_variant(object(), 0, 100) is None
    assert decode_cpu_seconds_per_variant(source, 0, 10) is None  # too few variants to time


def test_framed_sources_get_frame_multiple_candidates_and_whole_frame_probes():
    from torchgwas.empirical_autotune import decode_cpu_seconds_per_variant, frame_aligned_sizes
    defaults = (512, 1024, 2048, 4096)
    assert frame_aligned_sizes(defaults, None) == list(defaults)
    assert frame_aligned_sizes(defaults, 2048) == [2048, 4096]   # hard-call store frames
    assert frame_aligned_sizes(defaults, 2500) == [2500, 5000]   # full-scale zstd store frames
    assert frame_aligned_sizes(defaults, 256) == list(defaults)

    class Framed:
        allows_direct_native_fill = True
        native_transfer_dtype = np.uint8
        native_row_width = 16
        shape = (16, 10_000)
        chunk_alignment_variants = 2500
        calls = []

        @contextlib.contextmanager
        def native_reader_session(self):
            def fill(start, end, out):
                self.calls.append((start, end))
            yield fill

    source = Framed()
    assert decode_cpu_seconds_per_variant(source, 100, 10_000) >= 0
    assert source.calls == [(2500, 5000), (5000, 7500)]  # whole frames, from the first boundary
    assert decode_cpu_seconds_per_variant(source, 100, 7000) is None

    class FramedPacked:
        allows_direct_native_fill = False
        _bytes_per_variant = 9
        shape = (33, 10_000)
        chunk_alignment_variants = 2048
        calls = []

        def read_packed_into(self, target, start, end):
            assert len(target) == (end - start) * 9
            self.calls.append((start, end))

    packed = FramedPacked()
    assert decode_cpu_seconds_per_variant(packed, 0, 10_000) >= 0
    assert packed.calls == [(0, 2048), (2048, 4096)]


def test_shard_count_balances_setup_against_work():
    from torchgwas.empirical_autotune import best_shard_count, plan_layout
    assert best_shard_count(50.0, 1.0, 7) == 7      # lab-2080ti JAGWAS K=8,192: long GPU work
    assert best_shard_count(1.6, 0.5, 7) == 2       # H100 JAGWAS: 1.6 s of work, 7 shards lost
    assert best_shard_count(0.1, 1.0, 7) == 1
    seven = [f'cuda:{i}' for i in range(1, 8)]
    base = dict(BASE, n_variants=200_000, capacity=2048)
    short = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=40,
                        gpu_seconds_per_variant=8e-6, shard_setup_seconds=0.5, **base)
    assert short['variant_devices'] == seven[:2] and short['shard_model']['best'] == 2
    assert any('shards minimize' in reason for reason in short['why'])
    long = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=40,
                       gpu_seconds_per_variant=250e-6, shard_setup_seconds=1.0, **base)
    assert long['variant_devices'] == seven
    one = plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=40,
                      gpu_seconds_per_variant=1e-7, shard_setup_seconds=1.0, **base)
    assert one['variant_devices'] is None and one['devices_used'] == 1
    # Without measurements the variant-count rule stands.
    assert plan_layout(mode='jagwas', n_traits=8192, devices=seven, cpus=40, **base)['variant_devices'] == seven


def test_cpu_supply_counts_fair_share_of_busy_cores():
    from torchgwas.empirical_autotune import plan_layout
    four = ['cuda:0', 'cuda:2', 'cuda:5', 'cuda:7']
    demand = dict(decode_cpu_per_variant=57.5e-6, gpu_seconds_per_variant=32.7e-6)  # A100 full scale, significant pairs
    # 96 cores at load 82 with 6 idle: our threads still get ~19 cores; 4 GPUs need ~9.
    busy = plan_layout(mode='jagwas', n_traits=8192, devices=four, cpus=6, cpu_cores=96, cpu_load=82., **demand, **BASE)
    assert busy['variant_devices'] == four and busy['cpu_demand']['supply'] == 'fair_share'
    assert busy['reader_workers'] == 16  # 4 per GPU (Erlang-C sizing within the attainable supply), not 6 // 4
    # Oversubscribed (load 97 on 48 cores): the floor of two.
    over = plan_layout(mode='jagwas', n_traits=8192, devices=four, cpus=1, cpu_cores=48, cpu_load=97., **demand, **BASE)
    assert over['variant_devices'] == four[:2]
    # Without the host load the idle-core supply applies: (6 - 1) / 2.01 -> 2.
    idle = plan_layout(mode='jagwas', n_traits=8192, devices=four, cpus=6, **demand, **BASE)
    assert idle['variant_devices'] == four[:2] and idle['cpu_demand']['supply'] == 'idle_cores'


def test_decode_probe_times_packed_reads_for_plink_sources():
    from torchgwas.empirical_autotune import decode_cpu_seconds_per_variant

    class Packed:
        allows_direct_native_fill = False
        _bytes_per_variant = 9
        shape = (33, 10_000)
        calls = []

        def read_packed_into(self, target, start, end):
            assert len(target) == (end - start) * 9
            self.calls.append((start, end))
            sum(range(10_000))

    source = Packed()
    assert decode_cpu_seconds_per_variant(source, 0, 10_000) >= 0
    assert source.calls == [(0, 256), (256, 512)]


def test_readers_follow_measured_decode_demand():
    from torchgwas.empirical_autotune import readers_per_gpu
    # Full output at K=512 on H100: GPU ~0.7 us per variant, packed decode ~8 us (offered load ~11 readers).
    assert readers_per_gpu(1, cpus=40, cpu_cores=48, cpu_load=6., decode_cpu_per_variant=8e-6,
                           gpu_seconds_per_variant=0.7e-6) == 16
    # Full-scale A100 significant pairs (57.5 us decode, 32.7 us GPU, load 82 on 96): 4, as measured
    # (4 readers 218-231 s, 2 readers 276 s).
    assert readers_per_gpu(4, cpus=6, cpu_cores=96, cpu_load=82., decode_cpu_per_variant=57.5e-6,
                           gpu_seconds_per_variant=32.7e-6) == 4
    # GEMM-bound (JAGWAS-like): few readers suffice.
    assert readers_per_gpu(4, cpus=40, cpu_cores=48, cpu_load=6., decode_cpu_per_variant=20e-6,
                           gpu_seconds_per_variant=250e-6) == 2
    # A crowded host is not a cap: each thread gets a share, so the offered load grows and so do the
    # readers (H100 at load 116-182: 2 capped readers lost 1.7-3.2x to 16).
    assert readers_per_gpu(1, cpus=2, cpu_cores=48, cpu_load=97., decode_cpu_per_variant=8e-6,
                           gpu_seconds_per_variant=0.7e-6) == 16
    assert readers_per_gpu(1, cpus=2, cpu_cores=48, cpu_load=150., decode_cpu_per_variant=42e-6,
                           gpu_seconds_per_variant=0.74e-6) == 16
    # No measurement: 4.
    assert readers_per_gpu(2, cpus=16) == 4 and readers_per_gpu(2, cpus=1) == 4


def test_readers_and_depth_follow_the_gpus_in_use():
    from torchgwas.empirical_autotune import plan_layout
    eight = [f'cuda:{i}' for i in range(8)]
    # Full output K=512 at full scale on H100: 8 idle GPUs, the CPU cap allows 2, narrow output uses 1;
    # the single GPU gets the readers the supply allows for one, and the ring is deep enough for them.
    layout = plan_layout(mode='full', n_traits=512, devices=eight, cpus=12, cpu_cores=48, cpu_load=34.3,
                         decode_cpu_per_variant=78.5e-6, gpu_seconds_per_variant=0.74e-6,
                         **dict(BASE, n_samples=35_365, n_variants=8_931_083, capacity=8192,
                                transfer_bytes_per_variant=8842.))
    assert layout['devices_used'] == 1
    assert layout['reader_workers'] > 8 and layout['prefetch_chunks'] >= layout['reader_workers']


def test_dense_output_shards_follow_the_measured_writer():
    # Full-scale H100 rates: one writer ~1.1 GB/s, four ~3 GB/s together; a
    # K=2,048 panel's GEMM is ~1.9 us per variant. One GPU would wait on its
    # writer for ~120 s, so shards pay; a disk that one writer saturates does not.
    four = ['cuda:0', 'cuda:1', 'cuda:2', 'cuda:3']
    base = dict(BASE, n_samples=22_250, covariate_rank=10, n_variants=8_086_101,
                transfer_bytes_per_variant=22_250., gpu_seconds_per_variant=1.9e-6, shard_setup_seconds=0.2)
    wide = dict(bytes_per_variant=2048 * 8 + 4, writer_bytes_per_second=1.1e9,
                aggregate_bytes_per_second=3.2e9, aggregate_writers=4)
    plan = plan_layout(mode='full', n_traits=2048, devices=four, cpus=40, output_rates=wide, **base)
    assert len(plan['variant_devices']) >= 3
    model = plan['shard_model']['dense_seconds']
    assert model['1'] > 2.5 * min(model.values())
    saturated = dict(wide, aggregate_bytes_per_second=1.1e9)  # one writer saturates the disk
    assert plan_layout(mode='full', n_traits=2048, devices=four, cpus=40, output_rates=saturated,
                       **base)['variant_devices'] is None
    # Without a write measurement the narrow-panel rule still holds.
    assert plan_layout(mode='full', n_traits=2048, devices=four, cpus=40, **base)['variant_devices'] is None


def test_dense_shard_model_takes_the_slower_of_gpu_and_writer():
    from torchgwas.empirical_autotune import dense_shard_seconds
    rates = dict(bytes_per_variant=4100, writer_bytes_per_second=1e9, aggregate_bytes_per_second=2e9,
                 aggregate_writers=4)
    one = dense_shard_seconds(1, n_variants=1_000_000, gpu_seconds_per_variant=1e-6, setup_seconds=0.5, **rates)
    assert one == 0.5 + 4.1  # writer-bound: 4.1 GB at 1 GB/s
    four = dense_shard_seconds(4, n_variants=1_000_000, gpu_seconds_per_variant=1e-6, setup_seconds=0.5, **rates)
    assert four == 2.0 + 4.1 / 2  # aggregate cap, not 4 GB/s
    gpu = dense_shard_seconds(2, n_variants=1_000_000, gpu_seconds_per_variant=1e-5, setup_seconds=0.5, **rates)
    assert gpu == 1.0 + 5.0  # GPU-bound: 10 s of GEMM over two shards


def test_grouped_jagwas_is_priced_by_its_groups(monkeypatch):
    import torchgwas.empirical_autotune as autotune
    from torchgwas.jagwas_blocks import projection_flops_per_variant
    groups = [1000] * 30
    # Factor memory: the largest group live, the others' retained factors.
    assert autotune.jagwas_factor_bytes(30_000) == 3 * 30_000 ** 2 * 8
    assert autotune.jagwas_factor_bytes(30_000, groups) == 3 * 1000 ** 2 * 8 + 8 * 29 * 1000 ** 2
    # A panel whose single factor cannot fit one GPU fits as 30 groups.
    many = [f'cuda:{i}' for i in range(4)]
    with pytest.raises(ValueError, match='JAGWAS factor'):
        plan_layout(mode='jagwas', n_traits=30_000, devices=many, cpus=40, **dict(BASE, device_free_bytes=8*GIB))
    plan_layout(mode='jagwas', n_traits=30_000, devices=many, cpus=40, group_sizes=groups,
                **dict(BASE, device_free_bytes=8*GIB))
    # GPU time: each group's projection, not one 30,000-wide one.
    monkeypatch.setattr(autotune, 'measured_gemm_rate', lambda device, dtype: 1e12)
    single = autotune.gpu_seconds_per_variant('cuda:0', mode='jagwas', n_samples=1000, n_traits=30_000)
    grouped = autotune.gpu_seconds_per_variant('cuda:0', mode='jagwas', n_samples=1000, n_traits=30_000,
                                              group_sizes=groups)
    assert grouped == pytest.approx((2 * 1000 * 30_000 + 30 * projection_flops_per_variant(1000)) / 1e12)
    # The FP32 scoring is the same either way. The projection is G/2 to G times
    # cheaper: a small group is one triangular block (~2 k^2), the wide panel many (~K^2).
    assert autotune.jagwas_projection_flops(30_000, groups) == 30 * projection_flops_per_variant(1000)
    assert autotune.jagwas_projection_flops(30_000) / autotune.jagwas_projection_flops(30_000, groups) > 10
    assert single / grouped > 5


def test_output_write_rates_uses_the_real_writer_and_cleans_up(tmp_path):
    from torchgwas.empirical_autotune import output_write_rates
    rates = output_write_rates(tmp_path, n_traits=64, writers=2, probe_bytes=1 << 20, chunk_rows=256)
    assert rates['bytes_per_variant'] == 64 * 12 + 4 and rates['aggregate_writers'] == 2
    assert rates['writer_bytes_per_second'] > 0 and rates['aggregate_bytes_per_second'] > 0
    assert rates['probe_rows'] % 256 == 0
    assert list(tmp_path.iterdir()) == []  # scratch stores removed
