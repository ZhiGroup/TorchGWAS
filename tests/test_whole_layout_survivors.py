"""One sparse association pattern must survive candidate retile and sharding."""
import random
import pytest

from test_layout_output_floor import priced
from test_source_layout_floor import part, source
from torchgwas.productive_occupancy import whole_layout_survivors
from torchgwas.source_layout_floor import native_layout_source_floor


def layout(mode, traits, partitions):
    return dict(kind='torchgwas.pgen_layout_source_floor.v1',
                reduction=mode, total_traits=traits,
                input_identity={'test': 'conditional'},
                partitions=[dict(id=str(index), variant_range=list(variant),
                                 trait_range=list(trait))
                            for index, (variant, trait) in enumerate(partitions)])


def direct(lo, hi, traits, trait_range, numerator, denominator, placement):
    def selected(index):
        if placement == 'spread':
            return ((index + 1) * numerator + denominator - 1) // denominator > (
                index * numerator + denominator - 1) // denominator
        return index % denominator < numerator
    return sum(selected(v * traits + trait)
               for v in range(lo, hi)
               for trait in range(*trait_range))


@pytest.mark.parametrize('placement', ['spread', 'clustered'])
@pytest.mark.parametrize('fraction', [(1, 7), (3, 7), (5, 11)])
def test_significant_sparse_pattern_conserves_pairs_under_retile_and_sharding(
        placement, fraction):
    scenario=dict(retained_fraction=list(fraction), placement=placement)
    whole=whole_layout_survivors(
        layout('significant', 5, [((3, 12), (0, 5))]), scenario)
    retiled=whole_layout_survivors(
        layout('significant', 5, [
            ((3, 12), (0, 2)), ((3, 12), (2, 5))]), scenario)
    sharded=whole_layout_survivors(
        layout('significant', 5, [
            ((3, 7), (0, 5)), ((7, 12), (0, 5))]), scenario)
    expected=direct(3, 12, 5, (0, 5), *fraction, placement)
    assert whole['retained'] == retiled['retained'] == sharded['retained'] == expected
    assert retiled['retained_ranges']['0'][0] == direct(
        3, 12, 5, (0, 2), *fraction, placement)
    assert retiled['retained_ranges']['1'][0] == direct(
        3, 12, 5, (2, 5), *fraction, placement)


def test_sparse_rectangles_match_independent_pair_enumeration():
    rng = random.Random(5729)
    for _ in range(200):
        traits = rng.randint(1, 17)
        denominator = rng.randint(2, 31)
        numerator = rng.randint(1, denominator - 1)
        lo = rng.randint(0, 25)
        hi = lo + rng.randint(1, 24)
        a = rng.randrange(traits)
        b = rng.randint(a + 1, traits)
        placement = rng.choice(['spread', 'clustered'])
        scenario = dict(retained_fraction=[numerator, denominator],
                        placement=placement)
        report = whole_layout_survivors(
            layout('significant', traits, [((lo, hi), (a, b))]), scenario)
        assert report['retained'] == direct(
            lo, hi, traits, (a, b), numerator, denominator, placement)


@pytest.mark.parametrize('placement', ['spread', 'clustered'])
def test_jagwas_counts_variant_rows_without_phenotype_partition(placement):
    scenario=dict(retained_fraction=[3, 7], placement=placement)
    whole=whole_layout_survivors(
        layout('jagwas', 5, [((3, 12), (0, 5))]), scenario)
    sharded=whole_layout_survivors(
        layout('jagwas', 5, [
            ((3, 7), (0, 5)), ((7, 12), (0, 5))]), scenario)
    assert whole['retained'] == sharded['retained']
    assert whole['retained'] <= 9
    with pytest.raises(ValueError, match='complete phenotype'):
        whole_layout_survivors(
            layout('jagwas', 5, [((3, 12), (0, 3))]), scenario)


def test_large_spread_scenario_has_bounded_work_and_output_floor_accepts_counts(
        tmp_path):
    # The global pair ledger handles millions of markers without enumerating
    # their associations, then feeds the existing reduced-output calculator.
    huge=whole_layout_survivors(
        layout('significant', 600_000, [
            ((128, 8_086_101), (0, 300_000)),
            ((128, 8_086_101), (300_000, 600_000))]),
        dict(retained_fraction=[3, 1000], placement='spread'))
    total=(8_086_101 - 128) * 600_000
    assert huge['retained'] in (total * 3 // 1000, total * 3 // 1000 + 1)
    assert huge['clustered_variant_visits'] == 0

    floor=source(tmp_path)(0, 12, 3)
    fixed=native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    scenario=dict(retained_fraction=[3, 7], placement='spread')
    counts=whole_layout_survivors(fixed, scenario)
    output=priced(fixed, significant_backend='host',
                  retained_ranges=counts['retained_ranges'])
    assert sum(row['retained_rows'][0] for row in output['partitions']) == counts['retained']
    assert counts['retained_ranges']['low'][0] + counts['retained_ranges']['high'][0] == counts['retained']



def test_rare_spread_scenario_remains_exact_at_large_k():
    first, last, traits, denominator = 128, 8_086_101, 600_000, 100_000_000
    scenario = dict(retained_fraction=[1, denominator], placement='spread')
    full = whole_layout_survivors(layout('significant', traits,
        [((first, last), (0, traits))]), scenario)
    tiled = whole_layout_survivors(layout('significant', traits,
        [((first, last), (0, 100_000)),
         ((first, last), (100_000, traits))]), scenario)
    expected = ((last * traits + denominator - 1) // denominator -
                (first * traits + denominator - 1) // denominator)
    assert full['retained'] == tiled['retained'] == expected
    assert tiled['clustered_variant_visits'] == 0
def test_clustered_placement_has_explicit_work_budget():
    scenario=dict(retained_fraction=[1, 1_000_000], placement='clustered')
    with pytest.raises(ValueError, match='visit budget'):
        whole_layout_survivors(
            layout('significant', 3, [((0, 8_000_000), (0, 3))]),
            scenario, max_cluster_visits=4096)
    with pytest.raises(ValueError, match='Bounded reduced-output layout'):
        whole_layout_survivors(
            layout('significant', 3, [
                ((0, 2), (0, 1)), ((0, 2), (1, 3))]),
            'dense', max_partitions=1)
