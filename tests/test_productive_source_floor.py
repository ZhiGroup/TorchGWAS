"""Whole-source JIT work must cover precisely the post-output frontier."""
from copy import deepcopy
import time
import numpy as np
import pytest

from test_layout_frontier import productive_snapshot
from test_pgen_work_bounds import write_records
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.incremental_pgen_schedule import IncrementalPgenSchedule
from torchgwas.pgen_reader import pack_genovec
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.productive_source_stage import ProductiveSourceStage
from torchgwas.productive_source_floor import (
    productive_partial_floor, productive_source_floor)


def setup(tmp_path, reduction=None):
    path = tmp_path / 'plain.pgen'
    samples = 129
    records = [pack_genovec(np.full(samples, i % 3, dtype=np.uint8),
                             samples).tobytes() for i in range(12)]
    write_records(path, samples, [0] * 12, records)
    header = PgenHeaderWork(path)
    schedule = header.schedule_bounds(4, 12, 2)
    profile = dict(
        decode_units={name: 1e-8 for name in schedule['source_units']},
        cpu_fraction=.5, depth=2, decode_workers=2,
        cpu_available_cores=2.,
        shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7,
        input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                   cpu_seconds_per_call=1e-5))
    frontier = unissued_frontier(
        productive_snapshot(traits=5, reduction=reduction),
        source_identity=header.input_identity, reduction=reduction,
        total_traits=5, job_variant_range=[0, 12])
    return header, profile, frontier



def test_identical_trait_tile_source_prices_are_calculated_once(tmp_path, monkeypatch):
    import importlib

    module = importlib.import_module('torchgwas.productive_source_floor')
    header, profile, frontier = setup(tmp_path, reduction='significant')
    partitions = [dict(id='low', device='cuda:0', variant_range=[4, 12],
                       trait_range=[0, 2]),
                  dict(id='high', device='cuda:1', variant_range=[4, 12],
                       trait_range=[2, 5])]
    original = module.native_schedule_source_floor
    calls = []

    def counted(schedule, prices, capacities):
        calls.append(prices['decode_units']['copy_packed_byte'])
        return original(schedule, prices, capacities)

    monkeypatch.setattr(module, 'native_schedule_source_floor', counted)
    common = dict(chunk_markers=2,
                  shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
                  partition_axis='trait')
    equal = productive_source_floor(frontier, header, partitions,
        profiles={'low': profile, 'high': deepcopy(profile)}, **common)
    assert equal['unique_schedules'] == 1
    assert equal['unique_priced_source_floors'] == 1
    assert len(calls) == 1
    different = deepcopy(profile)
    different['decode_units']['copy_packed_byte'] *= 2
    split = productive_source_floor(frontier, header, partitions,
        profiles={'low': profile, 'high': different}, **common)
    assert split['unique_schedules'] == 1
    assert split['unique_priced_source_floors'] == 2
    assert len(calls) == 3
    assert (equal['source']['resource_work']['input_bytes'] ==
            split['source']['resource_work']['input_bytes'])
def test_post_output_trait_tiles_count_rereads_without_repeating_header_parse(
        tmp_path, monkeypatch):
    header, profile, frontier = setup(tmp_path)
    capacities = dict(cpu=1., dram=5e7, input=5e6)
    original = header.schedule_bounds
    calls = []

    def counted(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(header, 'schedule_bounds', counted)
    single = productive_source_floor(
        frontier, header,
        [dict(id='all', device='cuda:0', variant_range=[4, 12],
              trait_range=[0, 5])],
        chunk_markers=2, profiles={'all': profile},
        shared_capacities=capacities, partition_axis='trait')
    tiled = productive_source_floor(
        frontier, header, [
            dict(id='low', device='cuda:0', variant_range=[4, 12],
                 trait_range=[0, 2]),
            dict(id='high', device='cuda:1', variant_range=[4, 12],
                 trait_range=[2, 5])],
        chunk_markers=2, profiles={'low': profile, 'high': profile},
        shared_capacities=capacities, partition_axis='trait')
    assert calls == [(4, 12, 2), (4, 12, 2)]
    assert tiled['unique_schedules'] == 1
    assert tiled['unique_source_records'] == 8
    assert tiled['source']['resource_work']['input_bytes'] == 2 * single['source']['resource_work']['input_bytes']
    assert tiled['coverage']['exact_coverage']
    assert tiled['coverage']['required_cells'] == 40
    assert not tiled['selection_validated']


@pytest.mark.parametrize('reduction,axis', [
    (None, 'trait'), ('significant', 'trait'), ('jagwas', 'variant')])
def test_staged_primary_rebases_current_source_for_output_layouts(
        tmp_path, monkeypatch, reduction, axis):
    header, profile, frontier = setup(tmp_path, reduction=reduction)
    if axis == 'trait':
        partitions = [
            dict(id='low', device='cuda:0', variant_range=[4, 12],
                 trait_range=[0, 2]),
            dict(id='high', device='cuda:1', variant_range=[4, 12],
                 trait_range=[2, 5])]
    else:
        partitions = [
            dict(id='left', device='cuda:0', variant_range=[4, 8],
                 trait_range=[0, 5]),
            dict(id='right', device='cuda:1', variant_range=[8, 12],
                 trait_range=[0, 5])]
    options = dict(chunk_markers=3,
        profiles={row['id']: profile for row in partitions},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis=axis)
    direct = productive_source_floor(frontier, header, partitions, **options)
    staged = IncrementalPgenSchedule(header, 0, 12, 2,
                                     records_per_step=4)
    while not staged.snapshot()['complete']:
        staged.advance()
    with monkeypatch.context() as patch:
        patch.setattr(header, 'schedule_bounds',
            lambda *a, **k: pytest.fail('Staged rebase rescanned the source'))
        reused = productive_source_floor(frontier, header, partitions,
            staged_source=staged, **options)
    assert reused['source_schedule_method'] == 'staged_primary_rebase'
    assert reused['staged_segments_used'] == 3
    assert reused['prior_staged_calculation']['cpu_seconds'] > 0
    assert reused['coverage'] == direct['coverage']
    assert reused['source']['resource_work']['input_bytes'] == direct[
        'source']['resource_work']['input_bytes']
    assert reused['source']['resource_work']['cpu_seconds'] == pytest.approx(
        direct['source']['resource_work']['cpu_seconds'])
    assert reused['source']['source_stage_floor_seconds'] == pytest.approx(
        direct['source']['source_stage_floor_seconds'])


def test_completed_public_source_stage_charges_full_worker_before_source_rebase(
        tmp_path, monkeypatch):
    header, profile, frontier = setup(tmp_path)
    def header_factory():
        # Header preparation is real worker CPU even when an admitted PGEN
        # index is reused. Keep its cost distinct from schedule.advance().
        end = time.thread_time() + .005
        while time.thread_time() < end:
            pass
        return header
    config = dict(records_per_step=4, max_steps=3, max_cpu_seconds=10.,
                  max_window_seconds=10., max_retained_bytes=1 << 20,
                  extra_host_reserve_bytes=2 << 20)
    stage = ProductiveSourceStage(header_factory, header.input_identity, 2,
                                  config, on_step_observation=lambda: dict(valid=True))
    with pytest.raises(ValueError, match='Completed productive'):
        stage.completed_ledger()
    for _ in range(3):
        stage.output_written()
    deadline = time.monotonic() + 10.
    while not stage.snapshot()['complete'] and time.monotonic() < deadline:
        time.sleep(.01)
    audit = stage.finish()
    assert audit['complete'] and len(audit['steps']) == 3
    partitions = [dict(id='all', device='cuda:0', variant_range=[4, 12],
                       trait_range=[0, 5])]
    options = dict(chunk_markers=3, profiles={'all': profile},
                   shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
                   partition_axis='trait')
    direct = productive_source_floor(frontier, header, partitions, **options)
    with monkeypatch.context() as patch:
        patch.setattr(header, 'schedule_bounds',
                      lambda *a, **k: pytest.fail('Public staged rebase rescanned source'))
        reused = productive_source_floor(frontier, header, partitions,
                                         staged_source=stage, **options)
    assert reused['source_schedule_method'] == 'staged_primary_rebase'
    assert reused['coverage'] == direct['coverage']
    assert reused['source']['source_stage_floor_seconds'] == pytest.approx(
        direct['source']['source_stage_floor_seconds'])
    assert reused['prior_staged_calculation']['cpu_seconds'] == pytest.approx(
        audit['cpu_seconds'])
    assert reused['prior_staged_calculation']['wall_seconds'] == pytest.approx(
        audit['step_wall_seconds'])
    assert reused['prior_staged_calculation']['cpu_seconds'] > (
        stage.stage.snapshot()['calculation_cpu_seconds'] + .003)
    assert reused['productive_stage_worker_cost']['steps'] == 3


def test_jagwas_variant_shards_keep_complete_panel_and_budget_before_parse(
        tmp_path, monkeypatch):
    header, profile, frontier = setup(tmp_path, reduction='jagwas')
    rows = [
        dict(id='left', device='cuda:0', variant_range=[4, 8],
             trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[8, 12],
             trait_range=[0, 5])]
    options = dict(chunk_markers=2, profiles={'left': profile, 'right': profile},
                   shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
                   partition_axis='variant')
    result = productive_source_floor(frontier, header, rows, **options)
    assert result['coverage']['required_cells'] == 40
    assert result['unique_schedules'] == 2
    with pytest.raises(ValueError, match='complete phenotype'):
        productive_source_floor(frontier, header, [
            dict(rows[0], trait_range=[0, 3]), rows[1]], **options)
    with monkeypatch.context() as patch:
        patch.setattr(header, 'schedule_bounds',
                      lambda *a, **k: pytest.fail('Budget checked after parsing'))
        with pytest.raises(ValueError, match='record budget'):
            productive_source_floor(
                frontier, header, rows, max_unique_records=7, **options)
        with pytest.raises(ValueError, match='coverage'):
            productive_source_floor(
                frontier, header,
                [dict(rows[0], variant_range=[4, 9]), rows[1]], **options)


@pytest.mark.parametrize('reduction,axis', [
    ('significant', 'trait'), ('jagwas', 'variant')])
def test_one_global_sparse_scenario_composes_whole_source_compute_and_output(
        tmp_path, reduction, axis):
    header, profile, frontier = setup(tmp_path, reduction=reduction)
    common = dict(chunk_markers=2,
                  shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
                  partition_axis=axis)
    baseline_rows=[dict(id='whole', device='cuda:0',
                        variant_range=[4, 12], trait_range=[0, 5])]
    if reduction == 'significant':
        candidate_rows=[
            dict(id='low', device='cuda:0',
                 variant_range=[4, 12], trait_range=[0, 2]),
            dict(id='high', device='cuda:1',
                 variant_range=[4, 12], trait_range=[2, 5])]
    else:
        candidate_rows=[
            dict(id='left', device='cuda:0',
                 variant_range=[4, 8], trait_range=[0, 5]),
            dict(id='right', device='cuda:1',
                 variant_range=[8, 12], trait_range=[0, 5])]
    baseline=productive_source_floor(
        frontier, header, baseline_rows, profiles={'whole':profile}, **common)
    candidate=productive_source_floor(
        frontier, header, candidate_rows,
        profiles={row['id']:profile for row in candidate_rows}, **common)

    def partial(source):
        devices={row['device'] for row in source['source']['partitions']}
        compute=dict(covariate_rank=3,
            shared_h2d_bytes_per_second=1000.,
            per_device_h2d_bytes_per_second={d:800. for d in devices},
            peak_fp32_flops_per_second={d:100000. for d in devices})
        if reduction == 'jagwas':
            compute['peak_fp64_flops_per_second']={
                d:100000. for d in devices}
        output=dict(store_beta=True,
            shared_d2h_bytes_per_second=500.,
            per_device_d2h_bytes_per_second={d:400. for d in devices},
            output_bytes_per_second=200.)
        if reduction == 'significant':
            output['significant_backend']='host'
        return productive_partial_floor(
            frontier, source, compute_options=compute,
            output_options=output,
            occupancy_scenario=dict(
                retained_fraction=[3, 11], placement='spread'))

    old, new = partial(baseline), partial(candidate)
    assert old['survivors']['retained'] == new['survivors']['retained']
    assert old['envelope']['required_cells'] == new['envelope']['required_cells'] == 40
    assert old['envelope']['coverage']['exact_coverage']
    assert new['envelope']['coverage']['exact_coverage']
    assert new['envelope']['partial_floor_seconds'][0] > 0
    assert not new['prediction_complete'] and not new['selection_validated']


def test_dense_whole_source_composition_has_no_occupancy_assumption(tmp_path):
    header, profile, frontier = setup(tmp_path)
    source=productive_source_floor(
        frontier, header,
        [dict(id='dense', device='cuda:0', variant_range=[4, 12],
              trait_range=[0, 5])],
        chunk_markers=2, profiles={'dense':profile},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='trait')
    compute=dict(covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0':800.},
        peak_fp32_flops_per_second={'cuda:0':100000.})
    output=dict(store_beta=True,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0':400.},
        output_bytes_per_second=200.)
    result=productive_partial_floor(
        frontier, source, compute_options=compute, output_options=output)
    assert result['survivors'] is None
    assert result['output']['total_output_array_payload_bytes'][0] > 0
    assert result['envelope']['required_cells'] == 40
    boundary=dict(kind='torchgwas.productive_output_boundary.v1',
                  valid=True, reduction=None,
                  issued_revision=frontier['issued_revision'],
                  written_events=frontier['written_events'],
                  partitions=[dict(id='original',trait_range=[0,5],
                      issued_not_matrix_written_pairs=10,
                      issued_not_df_written_markers=2)])
    with_backlog=productive_partial_floor(
        frontier, source, compute_options=compute, output_options=output,
        output_boundary=boundary)
    assert with_backlog['issued_output_backlog']['array_payload_bytes_upper']==80
    assert with_backlog['array_payload_work_upper_bytes']==(
        result['output']['total_output_array_payload_bytes'][1]+80)
    stale_boundary=dict(boundary,issued_revision=boundary['issued_revision']+1)
    with pytest.raises(ValueError,match='output boundary'):
        productive_partial_floor(
            frontier,source,compute_options=compute,output_options=output,
            output_boundary=stale_boundary)
    with pytest.raises(ValueError, match='occupancy'):
        productive_partial_floor(
            frontier, source, compute_options=compute, output_options=output,
            occupancy_scenario='dense')
    stale=deepcopy(source)
    stale['coverage']['issued_revision']+=1
    with pytest.raises(ValueError, match='held productive frontier'):
        productive_partial_floor(
            frontier, stale, compute_options=compute, output_options=output)
