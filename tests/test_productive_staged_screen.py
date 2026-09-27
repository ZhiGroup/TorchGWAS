"""One post-output stage screens admitted layouts without rescanning PGEN."""
import time

import pytest

from test_layout_dense_writer_service import OPTIONS, profile as dense_writer_profile
from test_layout_jagwas_archive_floor import archive_price as jagwas_archive_price
from test_layout_jagwas_archive_floor import profile as jagwas_archive_profile
from test_layout_jagwas_selection_floor import prices as jagwas_selector_prices
from test_layout_significant_archive_floor import prices as significant_archive_prices
from test_layout_significant_archive_floor import profile as significant_archive_profile
from test_layout_significant_device_count_floor import price as count_transfer_price
from test_layout_significant_device_launch_floor import profile as launch_profile
from test_layout_significant_host_selection_floor import profile as significant_selector_profile
from test_productive_source_floor import setup
from test_significant_host_model import bank as significant_selector_prices
from torchgwas.productive_source_stage import ProductiveSourceStage
from torchgwas.productive_staged_screen import productive_staged_partial_screen


def staged(header):
    config = dict(records_per_step=4, max_steps=3, max_cpu_seconds=10.,
                  max_window_seconds=10., max_retained_bytes=1 << 20,
                  extra_host_reserve_bytes=2 << 20)
    stage = ProductiveSourceStage(lambda: header, header.input_identity, 2, config)
    for _ in range(3):
        stage.output_written()
    deadline = time.monotonic() + 10.
    while not stage.snapshot()['complete'] and time.monotonic() < deadline:
        time.sleep(.01)
    assert stage.finish()['complete']
    return stage


def candidate(mode, profile, chunk, split):
    axis = 'variant' if mode == 'jagwas' else 'trait'
    if split and mode == 'jagwas':
        parts = [dict(id='left', device='cuda:0', variant_range=[4, 8],
                      trait_range=[0, 5]),
                 dict(id='right', device='cuda:1', variant_range=[8, 12],
                      trait_range=[0, 5])]
    elif split:
        parts = [dict(id='low', device='cuda:0', variant_range=[4, 12],
                      trait_range=[0, 2]),
                 dict(id='high', device='cuda:1', variant_range=[4, 12],
                      trait_range=[2, 5])]
    else:
        parts = [dict(id='all', device='cuda:0', variant_range=[4, 12],
                      trait_range=[0, 5])]
    devices = {part['device'] for part in parts}
    compute = dict(covariate_rank=3, shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={device: 800. for device in devices},
        peak_fp32_flops_per_second={device: 100000. for device in devices})
    if mode == 'jagwas':
        compute['peak_fp64_flops_per_second'] = {
            device: 1000. for device in devices}
    output = dict(shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={device: 400. for device in devices},
        output_bytes_per_second=200.)
    if mode is None:
        output['dense_writer_options'] = OPTIONS
        services = dict(writer_profiles={
            device: dense_writer_profile() for device in devices})
    elif mode == 'jagwas':
        output['jagwas_writer_fsync'] = True
        services = dict(selection_prices=jagwas_selector_prices(),
            cpu_fraction_by_device={device: .5 for device in devices},
            archive_price=jagwas_archive_price(),
            archive_profiles={device: jagwas_archive_profile()
                              for device in devices})
    else:
        output.update(significant_backend='host',
                      significant_threshold_one=False,
                      significant_writer_fsync=True)
        services = dict(selection_prices=significant_selector_prices(),
            selection_profiles={device: significant_selector_profile()
                                for device in devices},
            archive_prices=significant_archive_prices(),
            archive_profiles={device: significant_archive_profile()
                              for device in devices})
    return dict(id='split' if split else 'baseline', chunk_markers=chunk,
        partitions=parts, partition_axis=axis,
        source_profiles={part['id']: profile for part in parts},
        compute_options=compute, output_options=output,
        mode_service_options=services)


@pytest.mark.parametrize('mode', [None, 'significant', 'jagwas'])
def test_staged_screen_joins_full_mode_service_and_exact_layout(
        tmp_path, monkeypatch, mode):
    header, profile, frontier = setup(tmp_path, reduction=mode)
    stage = staged(header)
    candidates = [candidate(mode, profile, 2, False),
                  candidate(mode, profile, 3, True)]
    scenario = None if mode is None else dict(
        retained_fraction=[1, 3], placement='spread')
    with monkeypatch.context() as patch:
        patch.setattr(header, 'schedule_bounds',
            lambda *a, **k: pytest.fail('Staged screen rescanned the source'))
        result = productive_staged_partial_screen(
            frontier, stage, candidates,
            shared_source_capacities=dict(cpu=1., dram=5e7, input=5e6),
            occupancy_scenario=scenario)
    assert result['stop_reason'] == 'complete'
    assert result['evaluated_candidates'] == 2
    assert result['prior_stage_cost_once']['cpu_seconds'] == pytest.approx(
        stage.snapshot()['cpu_seconds'])
    assert result['issued_revision'] == frontier['issued_revision']
    assert result['written_events'] == frontier['written_events']
    assert not result['prediction_complete'] and not result['selection_validated']
    baseline, split = [row['partial'] for row in result['candidates']]
    assert all(row['source_schedule_method'] == 'staged_primary_rebase'
               for row in result['candidates'])
    assert all(row['unique_priced_source_floors'] <=
               len(candidates[index]['partitions'])
               for index, row in enumerate(result['candidates']))
    assert all(partial['envelope']['coverage']['exact_coverage']
               for partial in (baseline, split))
    assert baseline['envelope']['required_cells'] == split['envelope']['required_cells'] == 40
    stages = split['envelope']['stage_floor_seconds']
    assert split['included_stage_services']
    if mode is None:
        assert stages['dense_writer_service'][0] > 0
        assert split['source']['resource_work']['input_bytes'] == (
            2 * baseline['source']['resource_work']['input_bytes'])
    elif mode == 'jagwas':
        assert stages['jagwas_selection_service'][0] > 0
        assert stages['jagwas_archive_service'][0] > 0
        assert stages['jagwas_single_consumer'][0] > 0
        assert baseline['survivors']['retained'] == split['survivors']['retained']
        assert all(row['trait_range'] == [0, 5] for row in split['source']['partitions'])
    else:
        assert stages['significant_host_selector_service'][0] > 0
        assert stages['significant_archive_service'][0] > 0
        assert stages['significant_single_consumer'][0] > 0
        assert baseline['survivors']['retained'] == split['survivors']['retained']


def test_screen_budget_stops_after_one_candidate(tmp_path):
    header, profile, frontier = setup(tmp_path)
    stage = staged(header)
    rows = [candidate(None, profile, 2, False),
            candidate(None, profile, 3, True)]
    result = productive_staged_partial_screen(
        frontier, stage, rows,
        shared_source_capacities=dict(cpu=1., dram=5e7, input=5e6),
        max_cpu_seconds=1e-12)
    assert result['stop_reason'] == 'screen_budget'
    assert result['evaluated_candidates'] == 1
    assert result['requested_candidates'] == 2


def test_screen_rejects_issued_prefix_and_jagwas_trait_split(tmp_path, monkeypatch):
    header, profile, frontier = setup(tmp_path, reduction='jagwas')
    stage = staged(header)
    row = candidate('jagwas', profile, 3, True)
    with monkeypatch.context() as patch:
        patch.setattr(header, 'schedule_bounds',
            lambda *a, **k: pytest.fail('Invalid layout parsed a source schedule'))
        row['partitions'][0]['variant_range'] = [3, 8]
        with pytest.raises(ValueError, match='coverage'):
            productive_staged_partial_screen(frontier, stage, [row],
                shared_source_capacities=dict(cpu=1., dram=5e7, input=5e6),
                occupancy_scenario='empty')
        row['partitions'][0]['variant_range'] = [4, 8]
        row['partitions'][0]['trait_range'] = [0, 4]
        with pytest.raises(ValueError, match='complete phenotype'):
            productive_staged_partial_screen(frontier, stage, [row],
                shared_source_capacities=dict(cpu=1., dram=5e7, input=5e6),
                occupancy_scenario='empty')


def test_device_significant_screen_charges_empty_block_barriers(tmp_path):
    header, profile, frontier = setup(tmp_path, reduction='significant')
    stage = staged(header)
    row = candidate('significant', profile, 3, True)
    row['output_options']['significant_backend'] = 'device'
    row['output_options']['device_selection_max_cells'] = 4
    row['mode_service_options'] = dict(
        count_transfer_prices={device: count_transfer_price()
                               for device in ('cuda:0', 'cuda:1')},
        launch_profiles={device: launch_profile()
                         for device in ('cuda:0', 'cuda:1')},
        archive_prices=significant_archive_prices(),
        archive_profiles={device: significant_archive_profile()
                          for device in ('cuda:0', 'cuda:1')})
    result = productive_staged_partial_screen(
        frontier, stage, [row],
        shared_source_capacities=dict(cpu=1., dram=5e7, input=5e6),
        occupancy_scenario='empty')
    partial = result['candidates'][0]['partial']
    assert partial['included_stage_services'] == [
        'significant_archive', 'significant_device_count',
        'significant_device_launch']
    assert partial['survivors']['retained'] == 0
    assert partial['envelope']['stage_floor_seconds'][
        'significant_device_count_barrier'][0] > 0
    assert partial['envelope']['stage_floor_seconds'][
        'significant_device_launch_count_serial'][0] > partial['envelope'][
        'stage_floor_seconds']['significant_device_count_barrier'][0]
    assert partial['envelope']['coverage']['exact_coverage']