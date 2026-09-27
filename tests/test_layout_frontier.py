"""JIT tile and GPU proposals conserve exactly the unissued associations."""
from copy import deepcopy
import time

import pytest

from test_source_layout_floor import part, source
from torchgwas.layout_frontier import unissued_frontier, bind_layout_to_frontier
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.source_layout_floor import native_layout_source_floor
from torchgwas.sumstats import DenseWriteProgress
from torchgwas.sumstats_indexed import IndexedChunkWrite


def productive_snapshot(traits=5, reduction=None):
    run = ProductiveTuningRun([
        dict(id='original', device='cuda:0', variant_range=[0, 12],
             trait_range=[0, traits])], chunk_sizes=[2, 4], initial=4)
    assert run.for_partition('original')(0, 12, run.capacity) == 4
    now = time.perf_counter()
    if reduction in ('jagwas', 'significant'):
        rows = 4 if reduction == 'jagwas' else 0
        event = IndexedChunkWrite(0, 4, reduction, rows, 100 if rows else 0,
            'part.npz' if rows else None, now, now, bool(rows))
    else:
        event = DenseWriteProgress(0, 4, (0, traits), 32 * traits, 4, now, 'cuda:0')
    run.output_written(event)
    snapshot = run.snapshot()
    run.finish(successful=False)
    return snapshot


def test_first_chunk_frontier_can_be_retile_without_duplicating_work(tmp_path):
    floor = source(tmp_path)
    source_floor = floor(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), source_floor),
        part('high', 'cuda:1', (2, 5), source_floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    frontier = unissued_frontier(productive_snapshot(),
        source_identity=source_floor['input_identity'],
        reduction=None, total_traits=5, job_variant_range=[0, 12])
    assert frontier['rectangles'][0]['variant_range'] == [4, 12]
    assert frontier['required_cells'] == 8 * 5
    coverage = bind_layout_to_frontier(frontier, layout)
    assert coverage['exact_coverage'] and coverage['required_cells'] == 40
    assert coverage['issued_revision'] == 1 and coverage['written_events'] == 1


def test_jagwas_unissued_full_panel_can_variant_shard(tmp_path):
    floor = source(tmp_path)
    left, right = floor(4, 8, 2), floor(8, 12, 2)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), left),
        part('right', 'cuda:1', (0, 5), right)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    frontier = unissued_frontier(productive_snapshot(reduction='jagwas'),
        source_identity=left['input_identity'],
        reduction='jagwas', total_traits=5, job_variant_range=[0, 12])
    assert bind_layout_to_frontier(frontier, layout)['exact_coverage']


def test_unequal_tile_cursors_preserve_each_issued_prefix(tmp_path):
    run = ProductiveTuningRun([
        dict(id='low', device='cuda:0', variant_range=[0, 12], trait_range=[0, 2]),
        dict(id='high', device='cuda:1', variant_range=[0, 12], trait_range=[2, 5])],
        chunk_sizes=[2, 4], initial=4)
    low, high = run.for_partition('low'), run.for_partition('high')
    assert low(0, 12, run.capacity) == 4
    assert low(4, 12, run.capacity) == 4
    assert high(0, 12, run.capacity) == 4
    now = time.perf_counter()
    run.output_written(DenseWriteProgress(0, 4, (2, 5), 96, 4, now, 'cuda:1'))
    snapshot = run.snapshot()
    run.finish(successful=False)
    floor = source(tmp_path)
    low_floor, high_floor = floor(8, 12, 2), floor(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), low_floor),
        part('high', 'cuda:1', (2, 5), high_floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    frontier = unissued_frontier(snapshot, source_identity=low_floor['input_identity'],
        reduction=None, total_traits=5, job_variant_range=[0, 12])
    assert frontier['required_cells'] == 4 * 2 + 8 * 3
    assert bind_layout_to_frontier(frontier, layout)['exact_coverage']
    rewound = deepcopy(layout)
    rewound['partitions'][0]['variant_range'] = [4, 12]
    with pytest.raises(ValueError, match='coverage'):
        bind_layout_to_frontier(frontier, rewound)


@pytest.mark.parametrize('damage', ['drop_trait', 'duplicate', 'issued_prefix',
                                    'extra_variant', 'source', 'mode', 'budget'])
def test_frontier_proof_rejects_changed_or_incomplete_candidate(tmp_path, damage):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    frontier = unissued_frontier(productive_snapshot(),
        source_identity=floor['input_identity'], reduction=None,
        total_traits=5, job_variant_range=[0, 12])
    layout = deepcopy(layout)
    if damage == 'drop_trait':
        layout['partitions'][1]['trait_range'] = [3, 5]
    elif damage == 'duplicate':
        layout['partitions'][1]['trait_range'] = [1, 5]
    elif damage == 'issued_prefix':
        layout['partitions'][0]['variant_range'] = [0, 12]
    elif damage == 'extra_variant':
        layout['partitions'][0]['variant_range'] = [4, 13]
    elif damage == 'source':
        layout['input_identity']['bytes'] += 1
    elif damage == 'mode':
        layout['reduction'] = 'significant'
    else:
        with pytest.raises(ValueError, match='budget'):
            bind_layout_to_frontier(frontier, layout, max_coverage_visits=1)
        return
    with pytest.raises(ValueError):
        bind_layout_to_frontier(frontier, layout)


@pytest.mark.parametrize('damage', ['before_output', 'gap', 'cursor', 'overlap',
                                    'whole_gap', 'whole_extra', 'jagwas_tile',
                                    'finished', 'stale'])
def test_frontier_rejects_invalid_productive_state(tmp_path, damage):
    floor = source(tmp_path)(0, 12, 2)
    snapshot = productive_snapshot()
    mode = None
    if damage == 'before_output':
        snapshot['first_written'] = None
    elif damage == 'gap':
        snapshot['partitions'][0]['ranges'] = [[1, 5]]
    elif damage == 'cursor':
        snapshot['partitions'][0]['cursor'] = 5
    elif damage == 'overlap':
        row = deepcopy(snapshot['partitions'][0])
        row['id'] = 'second'
        snapshot['partitions'].append(row)
    elif damage == 'whole_gap':
        snapshot['partitions'][0]['trait_range'] = [0, 4]
    elif damage == 'whole_extra':
        snapshot['partitions'][0]['variant_range'] = [0, 13]
    elif damage == 'jagwas_tile':
        mode = 'jagwas'
        snapshot['partitions'][0]['trait_range'] = [0, 4]
    elif damage == 'finished':
        snapshot['finished'] = {'successful': False}
    else:
        from pathlib import Path
        path = Path(floor['input_identity']['path'])
        path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(ValueError):
        unissued_frontier(snapshot, source_identity=floor['input_identity'],
                          reduction=mode, total_traits=5,
                          job_variant_range=[0, 12])
