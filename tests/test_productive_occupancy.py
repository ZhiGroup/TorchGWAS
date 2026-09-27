"""Partial output has one immutable association ledger for both chunk sizes."""
import pytest

from torchgwas.productive_occupancy import (
    occupancy_kind, survivor_bins, validate_occupancy_scenarios)
from torchgwas.window_model import _coverage, _rebin_survivors


PARTITIONS = [dict(cursor=0, variant_range=[0, 16], trait_range=[0, 3])]


def windows(chunk):
    return [dict(trait_range=[0, 3], data=dict(encoded=dict(
        input_identity={'source': 'same'}, variant_range=[0, 8],
        chunk_ranges=[[start, start + chunk] for start in range(0, 8, chunk)])))]


@pytest.mark.parametrize('placement,counts', [
    ('spread', [2, 2, 2, 2]), ('clustered', [3, 2, 1, 3])])
def test_partial_bins_rebin_without_changing_retained_associations(placement, counts):
    scenario=dict(retained_fraction=[3, 10], placement=placement)
    bins=survivor_bins(PARTITIONS, 8, fine_chunk=2,
        reduction='significant', total_traits=3, scenario=scenario)
    earlier=survivor_bins(PARTITIONS, 4, fine_chunk=2,
        reduction='significant', total_traits=3, scenario=scenario)
    later=survivor_bins([dict(PARTITIONS[0],cursor=4)],4,fine_chunk=2,
        reduction='significant',total_traits=3,scenario=scenario)
    assert bins[:2]==earlier
    assert bins[2:]==later
    assert [row['retained'] for row in bins]==counts
    evidence=dict(input_identity={'source': 'same'}, reduction='significant',
        total_traits=3, significance_threshold=.01, bins=bins)
    assert _rebin_survivors(windows(2), evidence, _coverage(windows(2)),
        'significant', 3, .01, 256)==[counts]
    assert _rebin_survivors(windows(4), evidence, _coverage(windows(4)),
        'significant', 3, .01, 256)==[[sum(counts[:2]),sum(counts[2:])]]


def test_jagwas_sparse_bins_keep_complete_phenotype_panel():
    scenario=dict(retained_fraction=[3, 10], placement='spread')
    bins=survivor_bins(PARTITIONS, 8, fine_chunk=2,
        reduction='jagwas', total_traits=3, scenario=scenario)
    assert all(row['trait_range']==[0, 3] and 0<=row['retained']<=2 for row in bins)
    assert sum(row['retained'] for row in bins)==3
    with pytest.raises(ValueError,match='partition'):
        survivor_bins([dict(PARTITIONS[0],trait_range=[0,2])],8,fine_chunk=2,
            reduction='jagwas',total_traits=3,scenario=scenario)


def test_sparse_first_windows_and_whole_layout_share_global_tile_ledger():
    scenario=dict(retained_fraction=[3, 11], placement='spread')
    whole=[dict(cursor=4,variant_range=[0,20],trait_range=[0,5])]
    tiles=[dict(cursor=4,variant_range=[0,20],trait_range=trait)
           for trait in ([0,2],[2,5])]
    one=survivor_bins(whole,8,fine_chunk=2,reduction='significant',
        total_traits=5,scenario=scenario)
    split=survivor_bins(tiles,8,fine_chunk=2,reduction='significant',
        total_traits=5,scenario=scenario)
    expected=sum(
        ((position+1)*3+10)//11-(position*3+10)//11
        for variant in range(4,12) for trait in range(5)
        for position in [variant*5+trait])
    assert sum(row['retained'] for row in one)==expected
    assert sum(row['retained'] for row in split)==expected
    assert [one[i]['retained'] for i in range(4)]==[
        split[i]['retained']+split[4+i]['retained'] for i in range(4)]


@pytest.mark.parametrize('scenario', [
    {'retained_fraction':[0,10],'placement':'spread'},
    {'retained_fraction':[10,10],'placement':'spread'},
    {'retained_fraction':[1,1_000_000_000_001],'placement':'spread'},
    {'retained_fraction':[1,1_000_001],'placement':'clustered'},
    {'retained_fraction':[1.,10],'placement':'spread'},
    {'retained_fraction':[1,10],'placement':'unknown'},
    {'retained_fraction':[1,10]}])
def test_sparse_scenario_requires_bounded_exact_fraction(scenario):
    with pytest.raises(ValueError):occupancy_kind(scenario)


def test_scenario_and_expansion_budgets_are_explicit():
    validate_occupancy_scenarios({'empty':'empty','partial':dict(
        retained_fraction=[1,1000],placement='clustered'),'dense':'dense'})
    with pytest.raises(ValueError,match='budget'):
        survivor_bins(PARTITIONS,8,fine_chunk=2,reduction='significant',
            total_traits=3,scenario=dict(retained_fraction=[1,2],placement='spread'),
            max_bins=3)
