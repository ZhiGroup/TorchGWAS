"""First useful chunks can refute output assumptions without repricing them."""
import pytest

from torchgwas.productive_output_sample import ProductiveOutputSample
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition
from torchgwas.sumstats import DenseWriteProgress


PARTITIONS = [dict(id='a', device='cuda:1', variant_range=[100, 200], trait_range=[0, 4]),
              dict(id='b', device='cuda:2', variant_range=[100, 200], trait_range=[4, 7])]


def event(*, rows=0, start=100, end=104, device='cuda:1', traits=(0, 4),
          kind='significant', bound=True):
    partition = IndexedOutputPartition(device, (100, 200), traits) if bound else None
    return IndexedChunkWrite(start-100, end-100, kind, rows,
                             0 if rows == 0 else 100, None if rows == 0 else 'part.npz',
                             1., 1.25, rows > 0, partition,
                             (start, end) if bound else None)


def test_completed_output_refutes_only_incompatible_reduced_scenario():
    sample = ProductiveOutputSample(PARTITIONS, 'significant')
    sample.observe(event(rows=0))
    sample.check_scenarios({'low': 'empty'})
    with pytest.raises(ValueError, match='full-only'):
        sample.check_scenarios({'high': 'dense'})
    sample.observe(event(rows=2, start=104, end=108))
    with pytest.raises(ValueError, match='empty-only'):
        sample.check_scenarios({'low': 'empty'})
    with pytest.raises(ValueError, match='partial output'):
        sample.check_scenarios({'low': 'empty', 'high': 'dense'})
    sample.check_scenarios({'low': 'empty', 'middle': dict(
        retained_fraction=[1, 8], placement='spread'), 'high': 'dense'})
    row = sample.snapshot()
    assert row['counts']['empty'] == row['counts']['partial'] == 1
    assert row['events'][1]['writer_wall_seconds'] == pytest.approx(.25)
    assert row['events'][1]['possible'] == 16


def test_jagwas_uses_one_possible_result_per_variant_and_full_panel():
    sample = ProductiveOutputSample([dict(PARTITIONS[0], trait_range=[0, 7])], 'jagwas')
    sample.observe(event(rows=4, kind='jagwas', traits=(0, 7)))
    assert sample.snapshot()['counts']['full'] == 1
    sample.check_scenarios({'all': 'dense'})
    sample.observe(event(rows=1, start=104, end=108, kind='jagwas', traits=(0, 7)))
    with pytest.raises(ValueError, match='full-only'):
        sample.check_scenarios({'all': 'dense'})


def test_unbound_or_duplicate_event_disables_optional_planning():
    sample = ProductiveOutputSample(PARTITIONS, 'significant')
    sample.observe(event(bound=False))
    with pytest.raises(ValueError, match='cannot be bound'):
        sample.check_scenarios({'low': 'empty', 'high': 'dense'})
    other = ProductiveOutputSample(PARTITIONS, 'significant')
    other.observe(event())
    other.observe(event())
    assert other.snapshot()['counts']['invalid'] == 1
    with pytest.raises(ValueError, match='cannot be bound'):
        other.check_scenarios({'low': 'empty', 'high': 'dense'})


def test_sample_retention_is_bounded_while_counters_continue():
    sample = ProductiveOutputSample(PARTITIONS, 'significant', max_events=2)
    for i in range(5):
        sample.observe(event(start=100+4*i, end=104+4*i))
    result = sample.snapshot()
    assert len(result['events']) == 2
    assert result['counts']['indexed_events'] == result['counts']['empty'] == 5
    assert result['counts']['sampled'] == 2


def test_dense_progress_does_not_claim_per_chunk_occupancy():
    sample = ProductiveOutputSample(PARTITIONS, None)
    sample.observe(DenseWriteProgress(0, 4, (0, 4), 128, 4, 1., 'store'))
    assert sample.snapshot()['counts']['indexed_events'] == 0
