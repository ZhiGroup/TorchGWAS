"""Bounded post-output steps must conserve the exact whole PGEN schedule."""
import pytest

from test_pgen_work_bounds import fixture
from torchgwas.incremental_pgen_schedule import IncrementalPgenSchedule
from torchgwas.pgen_work_bounds import (PgenHeaderWork,
    native_schedule_source_floor, paired_schedule_source_difference)


@pytest.mark.parametrize('span', [(0, 15), (1, 15), (2, 11)])
@pytest.mark.parametrize('chunk', [1, 2, 3, 8, 16])
def test_segments_conserve_exact_source_and_ld_restart_work(tmp_path,
                                                             span, chunk):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    direct = header.schedule_bounds(*span, chunk)
    staged = IncrementalPgenSchedule(header, *span, chunk,
        records_per_step=2 * chunk, max_chunks_per_step=2)
    assert staged.snapshot()['cursor'] == span[0]
    assert not staged.snapshot()['complete']
    with pytest.raises(ValueError, match='incomplete'):
        staged.finish()
    while not staged.snapshot()['complete']:
        before = staged.snapshot()['cursor']
        progress = staged.advance()
        assert before < progress['cursor'] <= min(before + 2 * chunk,
                                                  span[1])
    assert staged.advance()['completed_segments'] == staged.snapshot()[
        'completed_segments']
    merged = staged.finish()
    for key in direct:
        if key != 'scope':
            assert merged[key] == direct[key], key
    prices = {name: 1e-8 for name in merged['source_units']}
    profile = dict(decode_units=prices, cpu_fraction=.5,
                   depth=3, decode_workers=2,
                   cpu_available_cores=4.,
                   shared_dram_bytes_per_second=1e8,
                   read_bytes_per_second=1e7,
                   input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                              cpu_seconds_per_call=1e-5))
    capacities = dict(cpu=2., dram=5e7, input=5e6)
    actual = native_schedule_source_floor(merged, profile, capacities)
    expected = native_schedule_source_floor(direct, profile, capacities)
    assert actual['resource_work']['input_bytes'] == expected[
        'resource_work']['input_bytes']
    assert actual['resource_work']['cpu_seconds'] == pytest.approx(
        expected['resource_work']['cpu_seconds'])
    assert actual['source_stage_floor_seconds'] == pytest.approx(
        expected['source_stage_floor_seconds'])


def test_staged_schedules_preserve_paired_chunk_change(tmp_path):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 65)
    header = PgenHeaderWork(path)
    schedules = []
    for chunk in (1, 3):
        staged = IncrementalPgenSchedule(header, 1, 15, chunk,
            records_per_step=3 * chunk, max_chunks_per_step=3)
        while not staged.snapshot()['complete']:
            staged.advance()
        schedules.append(staged.finish())
    prices = {name: 1e-8 for row in schedules
              for name in row['source_units']}
    actual = paired_schedule_source_difference(*schedules, prices)
    expected = paired_schedule_source_difference(
        header.schedule_bounds(1, 15, 1),
        header.schedule_bounds(1, 15, 3), prices)
    assert actual == expected


@pytest.mark.parametrize('origin,cursors', [
    (0, [0, 1, 2, 3, 5, 8, 13]),
    (1, [1, 2, 4, 7, 12])])
def test_staged_primary_work_rebases_to_later_cursor_and_other_chunk_grid(
        tmp_path, origin, cursors):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    staged = IncrementalPgenSchedule(header, origin, 15, 3,
        records_per_step=6, max_chunks_per_step=2)
    while not staged.snapshot()['complete']:
        staged.advance()
    for cursor in cursors:
        for stop in sorted(set([min(15, cursor + 2),
                                min(15, cursor + 7), 15])):
            for chunk in (1, 2, 3, 8, 16):
                rebased = staged.rebase(cursor, chunk, stop=stop)
                direct = header.schedule_bounds(cursor, stop, chunk)
                for key in direct:
                    if key != 'scope':
                        assert rebased[key] == direct[key], \
                               (cursor, stop, chunk, key)
                if cursor > origin:
                    assert rebased['variant_range'][0] == cursor


def test_step_budget_and_changed_source_fail_closed(tmp_path):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    with pytest.raises(ValueError, match='segment budget'):
        IncrementalPgenSchedule(header, 0, 15, 1,
            records_per_step=2, max_segments=2)
    staged = IncrementalPgenSchedule(header, 0, 15, 3,
        records_per_step=3)
    staged.advance()
    original = path.read_bytes()
    path.write_bytes(original + b'changed')
    with pytest.raises(ValueError, match='changed'):
        staged.advance()


def test_no_unbounded_header_step_when_there_are_many_chunks(tmp_path,
                                                               monkeypatch):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    original = header.schedule_bounds
    calls = []

    def counted(start, stop, chunk, **kwargs):
        calls.append((start, stop, chunk))
        assert stop - start <= 4
        return original(start, stop, chunk, **kwargs)

    monkeypatch.setattr(header, 'schedule_bounds', counted)
    staged = IncrementalPgenSchedule(header, 0, 15, 2,
        records_per_step=4, max_chunks_per_step=2)
    while not staged.snapshot()['complete']:
        staged.advance()
    assert calls == [(0, 4, 2), (4, 8, 2), (8, 12, 2), (12, 15, 2)]
    assert staged.finish()['chunk_count'] == 8
