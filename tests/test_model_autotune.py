"""Model tuner: per-chunk samples, load-adjusted comparison, re-planning on drift."""
from dataclasses import dataclass
import random

from torchgwas.model_autotune import ModelChunkTuner


@dataclass
class Obs:
    start: int
    end: int
    device: str
    completed: float
    read_cpu_seconds: float | None = None
    read_runnable_wait_seconds: float | None = None


def simulate(tuner, cost, *, devices=('cuda:0',), rows=400_000, noise=0.0, seed=0, load=None, change=None,
             startup=None):
    """Devices pull chunks of their own ranges; per-row cost (us) depends on size.

    load(rng) -> multiplier for the next chunk (reported as reader wait);
    change=(fraction, cost2) switches the size costs part-way through;
    startup=(fraction, factor) slows the first part of the job (not load).
    """
    rng = random.Random(seed)
    span = rows // len(devices)
    cursors = {d: i * span for i, d in enumerate(devices)}
    ends = {d: (i + 1) * span for i, d in enumerate(devices)}
    clock = {d: 0. for d in devices}
    factor = 1.0
    done = 0
    while any(cursors[d] < ends[d] for d in devices):
        device = min((d for d in devices if cursors[d] < ends[d]), key=clock.get)
        size = tuner.control(cursors[device], ends[device], tuner.capacity)
        start = cursors[device]
        cursors[device] += size
        table = change[1] if change and done >= change[0] * rows else cost
        per_row = table.get(size) or table[min(k for k in table if k >= size)]
        if load is not None:
            factor = load(rng, factor)
        ramp = startup[1] if startup and done < startup[0] * rows else 1.0
        seconds = size * per_row * factor * ramp * (1 + noise * rng.uniform(-1, 1)) * 1e-6
        clock[device] += seconds
        done += size
        cpu = seconds / factor
        tuner(Obs(start, start + size, device, clock[device], cpu, seconds - cpu))
    return max(clock.values())


def test_picks_fastest_size_early_from_per_chunk_samples():
    tuner = ModelChunkTuner([512, 1024, 2048], total_rows=400_000, min_job_seconds=0)
    simulate(tuner, {512: 3.0, 1024: 2.0, 2048: 1.5}, noise=.05)
    audit = tuner.audit()
    assert audit['state'] == 'committed' and audit['choice'] == 2048
    assert audit['decided_fraction'] < .15  # segment trials took ~40%
    # Only sizes above the start are probed; the start is revisited last.
    assert [(s['kind'], s['size']) for s in audit['segments']] == [('warmup', 1024), ('probe', 2048), ('revisit', 1024)]
    assert audit['model']['b_seconds_per_row'] > 0


def test_bursty_load_does_not_decide_the_size():
    # Load jumps between 1x and 4x at random, often enough to land inside
    # probes; 2048 is truly 10% faster.
    def bursts(rng, factor):
        return factor if rng.random() > .5 else rng.choice([1.0, 4.0])
    for seed in range(6):
        tuner = ModelChunkTuner([512, 1024, 2048], total_rows=600_000, min_job_seconds=0, seed=seed)
        simulate(tuner, {512: 2.2, 1024: 2.0, 2048: 1.8}, rows=600_000, seed=seed, load=bursts)
        audit = tuner.audit()
        assert audit['choice'] == 2048, (seed, audit['decisions'])
        assert audit['load_coefficient'] > 0.5  # the wait fraction explains the bursts


def test_small_gain_keeps_the_starting_size():
    tuner = ModelChunkTuner([512, 1024, 2048], total_rows=400_000, margin=.05, min_job_seconds=0)
    simulate(tuner, {512: 2.0, 1024: 2.0, 2048: 1.97})
    assert tuner.audit()['choice'] == 1024 and tuner.reason == 'within_margin_of_incumbent'


def test_per_chunk_cost_growth_reprobes_a_larger_size_and_switches():
    # 1024 is within the margin of 2048 until half-way; then a per-chunk cost
    # appears (per-row time a/c + b with a large a), and 2048 wins clearly.
    tuner = ModelChunkTuner([1024, 2048, 4096], total_rows=800_000, min_job_seconds=0, margin=.05)
    simulate(tuner, {1024: 2.0, 2048: 1.96, 4096: 1.95}, rows=800_000,
             change=(.5, {1024: 4.0, 2048: 2.6, 4096: 1.9}))
    audit = tuner.audit()
    assert [d['choice'] for d in audit['decisions']][:1] == [1024]
    assert audit['reprobes'] >= 1 and audit['choice'] in (2048, 4096)


def test_uniform_slowdown_at_the_largest_size_is_not_reprobed():
    # Host load doubles every size's per-row time half-way through. The
    # committed size is the largest and per-chunk cost is ~0, so no neighbour
    # can gain the margin: the drift is recorded, nothing is re-probed.
    tuner = ModelChunkTuner([1024, 2048, 4096], total_rows=800_000, min_job_seconds=0)
    simulate(tuner, {1024: 2.2, 2048: 1.6, 4096: 1.3}, rows=800_000,
             change=(.5, {1024: 4.4, 2048: 3.2, 4096: 2.6}))
    audit = tuner.audit()
    assert audit['choice'] == 4096 and audit['reprobes'] == 0
    assert audit['skipped_reprobes'] and audit['skipped_reprobes'][0]['reason'] == 'no neighbour can gain the margin'


def test_smaller_size_is_reprobed_only_when_the_fit_says_per_row_cost_grows():
    # Larger chunks are measurably slower per row (a < 0 beyond its error),
    # so after a drift the smaller neighbour may be re-probed.
    tuner = ModelChunkTuner([512, 1024, 2048], total_rows=900_000, min_job_seconds=0)
    simulate(tuner, {512: 1.0, 1024: 1.2, 2048: 3.0}, rows=900_000,
             change=(.5, {512: 1.0, 1024: 3.0, 2048: 3.0}))
    audit = tuner.audit()
    assert audit['decisions'][0]['choice'] == 1024
    assert audit['model']['a_seconds'] + audit['model']['a_error_seconds'] < 0
    assert audit['reprobes'] >= 1 and audit['choice'] == 512


def test_short_job_keeps_the_start_and_two_devices_settle():
    short = ModelChunkTuner([512, 1024, 2048], total_rows=20_000)
    assert short.control(0, 20_000, 2048) == 1024
    simulate(short, {512: 2.0, 1024: 2.0, 2048: 1.0}, rows=20_000)
    assert short.audit()['state'] in ('skipped', 'deferred', 'warmup') and short.audit()['choice'] in (None, 1024)
    pair = ModelChunkTuner([512, 1024, 2048], total_rows=600_000, concurrent=2, min_job_seconds=0)
    simulate(pair, {512: 3.0, 1024: 1.0, 2048: 2.0}, devices=('cuda:0', 'cuda:1'), rows=600_000)
    audit = pair.audit()
    assert audit['devices'] == ['cuda:0', 'cuda:1'] and audit['choice'] == 1024


def test_noisy_samples_do_not_switch_to_a_slower_size():
    # +-60% chunk-to-chunk noise (a loaded host); 512 is truly 5% slower than
    # the starting 1024, 2048 5% faster. Unresolved differences keep 1024.
    slower = 0
    for seed in range(8):
        tuner = ModelChunkTuner([512, 1024, 2048], total_rows=400_000, initial=1024, min_job_seconds=0, seed=seed)
        simulate(tuner, {512: 2.1, 1024: 2.0, 2048: 1.9}, noise=.6, seed=seed)
        first = tuner.audit()['decisions'][0]
        assert 'errors' in first and all(v >= 0 for v in first['errors'].values())
        slower += first['choice'] == 512
    assert slower <= 1


def test_noisy_samples_still_switch_on_a_large_gain():
    tuner = ModelChunkTuner([512, 1024, 2048], total_rows=400_000, initial=1024, min_job_seconds=0)
    simulate(tuner, {512: 3.0, 1024: 2.0, 2048: 1.2}, noise=.6)
    assert tuner.audit()['decisions'][0]['choice'] == 2048


def test_chunk_noise_alone_does_not_trigger_drift_reprobes():
    tuner = ModelChunkTuner([512, 1024, 2048], total_rows=800_000, initial=1024, min_job_seconds=0)
    simulate(tuner, {512: 2.5, 1024: 2.0, 2048: 2.2}, noise=.6, rows=800_000)
    assert tuner.audit()['reprobes'] == 0


def test_a_clearly_slower_size_ends_its_probe_early():
    tuner = ModelChunkTuner([512, 1024, 2048, 4096], total_rows=400_000, initial=1024, min_job_seconds=0,
                            probe_chunks=8)
    simulate(tuner, {512: 2.1, 1024: 2.0, 2048: 20.0, 4096: 1.9}, noise=.05)
    audit = tuner.audit()
    slow = next(s for s in audit['segments'] if s['size'] == 2048)
    assert slow['samples'] < 8  # stopped before the full visit
    assert audit['decisions'][0]['choice'] in (1024, 4096)
    assert all(s['size'] != 512 for s in audit['segments'])


def test_start_size_is_judged_by_its_revisit_not_its_start_up():
    # The warmup (first 2% of the job) runs 2x slower (pipeline start-up). Measured
    # only in warmup, 1024 would look worse than 2048; its revisit shows the
    # two within the margin, so the start size is kept.
    tuner = ModelChunkTuner([1024, 2048], total_rows=800_000, min_job_seconds=0, margin=.05)
    simulate(tuner, {1024: 2.0, 2048: 1.98}, rows=800_000,
             load=None, change=None, startup=(.02, 2.0))
    audit = tuner.audit()
    assert [s['kind'] for s in audit['segments']][:3] == ['warmup', 'probe', 'revisit']
    assert audit['decisions'][0]['choice'] == 1024
