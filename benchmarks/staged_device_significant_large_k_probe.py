"""Metadata-only 8M-by-600k device-significant staged-screen cost probe.

Eight indexed output events are simulated. No genotype scan, GPU selection or
NPZ write occurs. All service rates are synthetic, not calibration evidence.
"""
import argparse
import hashlib
import importlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time

import torch

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.productive_source_stage import ProductiveSourceStage
from torchgwas.productive_staged_screen import productive_staged_partial_screen
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_profile():
    units = ('copy_packed_byte', 'expand4_int8', 'expand_int8_tail_sample',
             'uleb1', 'uleb2', 'uleb3', 'uleb4', 'uleb5', 'set_category',
             'difflist_group_absolute_id', 'difflist_category_extract',
             'difflist_record_header', 'fill_packed_byte',
             'invert_packed_byte', 'native_onebit_tail_bit_test',
             'onebit8_native')
    return dict(decode_units={name: 1e-9 for name in units},
        cpu_fraction=.5, depth=2, decode_workers=2,
        cpu_available_cores=8., shared_dram_bytes_per_second=1e10,
        read_bytes_per_second=1e9)


def layout(name, chunk, tile, cursor, variants, traits, source):
    partitions = []
    for index, first in enumerate(range(0, traits, tile)):
        partitions.append(dict(id='tile-' + str(index),
            device='cuda:' + str(index % 2),
            variant_range=[cursor, variants],
            trait_range=[first, min(traits, first + tile)]))
    devices = {row['device'] for row in partitions}
    archive = dict(cpu_fraction=.5, shared_dram_bytes_per_second=1e10,
        fsync_seconds=.001,
        writeback_service=dict(pagecache_seconds_per_byte=1e-10,
                               storage_seconds_per_byte=1e-9),
        process_units=dict(numpy_copy_bytes=1e-10))
    price = dict(call_cpu_seconds=1e-6, byte_cpu_seconds=1e-10)
    return dict(id=name, chunk_markers=chunk, partitions=partitions,
        partition_axis='trait',
        source_profiles={row['id']: source for row in partitions},
        compute_options=dict(covariate_rank=3,
            shared_h2d_bytes_per_second=1e10,
            per_device_h2d_bytes_per_second={device: 1e10 for device in devices},
            peak_fp32_flops_per_second={device: 1e13 for device in devices}),
        output_options=dict(store_beta=True,
            significant_backend='device', device_selection_max_cells=1 << 20,
            significant_threshold_one=False, significant_writer_fsync=True,
            shared_d2h_bytes_per_second=1e10,
            per_device_d2h_bytes_per_second={device: 1e10 for device in devices},
            output_bytes_per_second=1e9),
        mode_service_options=dict(
            count_transfer_prices={device: dict(latency_seconds=1e-5,
                bytes_per_second=1e10, resources=['d2h']) for device in devices},
            launch_profiles={device: dict(torch_version=torch.__version__,
                cuda_runtime=torch.version.cuda,
                compute_capability=list(torch.cuda.get_device_capability(device)),
                kernel_launch_seconds=1e-6, gpu_fraction=1.) for device in devices},
            archive_prices={'True': dict(price), 'False': dict(price)},
            archive_profiles={device: archive for device in devices}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--variants', type=int, required=True)
    parser.add_argument('--traits', type=int, default=600000)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--matched-pricing-control', action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Probe output already exists')
    source = args.input.resolve(strict=True)
    mount = subprocess.check_output(['findmnt', '-T', str(source),
        '-o', 'TARGET,SOURCE,FSTYPE', '-n'], text=True).strip()
    if not mount.startswith('/data '):
        raise ValueError('Probe requires server-local /data PGEN')
    identity = input_identity(source)
    stage = ProductiveSourceStage(
        lambda: PgenHeaderWork(source, max_cached_signatures=32768,
                               max_cached_bounds=256),
        identity, 128,
        dict(records_per_step=1_048_576, max_steps=8,
            max_cpu_seconds=30., max_window_seconds=120.,
            max_retained_bytes=16 << 20, extra_host_reserve_bytes=528 << 20,
            max_cached_signatures=32768, max_cached_bounds=256))
    run = ProductiveTuningRun([dict(id='all', device='cuda:0',
        variant_range=[0, args.variants], trait_range=[0, args.traits])],
        chunk_sizes=[128, 512], initial=128,
        budget=IncrementalPlanningBudget(max_steps=2,
            max_cpu_seconds=60., max_window_seconds=120.))
    cursor = 0
    partition = IndexedOutputPartition('cuda:0', (0, args.variants),
                                       (0, args.traits))
    control = run.for_partition('all')
    for _ in range(8):
        width = control(cursor, args.variants, 512)
        end = cursor + width
        now = time.perf_counter()
        run.output_written(IndexedChunkWrite(cursor, end, 'significant',
            0, 0, None, now, now, False, partition,
            (cursor, end)))
        stage.output_written()
        cursor = end
    deadline = time.monotonic() + 120.
    while not stage.snapshot()['complete'] and time.monotonic() < deadline:
        observed = stage.snapshot()
        if observed['stop_reason'] is not None:
            raise RuntimeError('Stage stopped: ' + str(observed['stop_reason']))
        time.sleep(.05)
    stage_audit = stage.finish()
    if not stage_audit['complete'] or stage.stage.header._header.variant_ct != args.variants:
        raise ValueError('Exact large-source stage did not complete')
    frontier = unissued_frontier(run.snapshot(), source_identity=identity,
        reduction='significant', total_traits=args.traits,
        job_variant_range=[0, args.variants])
    profile = source_profile()
    candidates = [layout('six_tiles_128', 128, 100000, cursor,
                         args.variants, args.traits, profile),
                  layout('twelve_tiles_512', 512, 50000, cursor,
                         args.variants, args.traits, profile)]
    options = dict(shared_source_capacities=dict(cpu=8., dram=1e10, input=1e9),
        occupancy_scenario=dict(retained_fraction=[1, 100000000],
                                placement='spread'),
        max_cpu_seconds=60., max_wall_seconds=120.)
    matched = None
    if args.matched_pricing_control:
        module = importlib.import_module('torchgwas.productive_source_floor')
        original_digest = module._digest
        matched = []
        fingerprints = []
        try:
            for mode in ('cached', 'forced_no_reuse',
                         'forced_no_reuse', 'cached'):
                if mode == 'cached':
                    module._digest = original_digest
                else:
                    counter = itertools.count()
                    module._digest = lambda value: (
                        original_digest(value) + ':' + str(next(counter)))
                screen = productive_staged_partial_screen(
                    frontier, stage, candidates, **options)
                outcome = [dict(source=row['partial']['source']['resource_work'],
                    floor=row['partial']['envelope']['partial_floor_seconds'],
                    retained=row['partial']['survivors']['retained'])
                    for row in screen['candidates']]
                fingerprints.append(hashlib.sha256(json.dumps(
                    outcome, sort_keys=True).encode()).hexdigest())
                matched.append(dict(mode=mode,
                    screen_cpu_seconds=screen['screen_cpu_seconds'],
                    screen_wall_seconds=screen['screen_wall_seconds'],
                    source_cpu_seconds=[row['source_floor_cpu_seconds']
                                        for row in screen['candidates']],
                    priced_source_floors=[row['unique_priced_source_floors']
                                          for row in screen['candidates']],
                    model_fingerprint=fingerprints[-1]))
        finally:
            module._digest = original_digest
        if len(set(fingerprints)) != 1:
            raise ValueError('Matched pricing changed the conditional model')
    else:
        screen = productive_staged_partial_screen(
            frontier, stage, candidates, **options)
    if screen['evaluated_candidates'] != 2:
        raise ValueError('Both large-K candidate floors required')
    selected = [row['partial'] for row in screen['candidates']]
    if (selected[0]['survivors']['retained'] !=
            selected[1]['survivors']['retained'] or
            not all(row['envelope']['coverage']['exact_coverage'] for row in selected) or
            not all('significant_device_launch' in row['included_stage_services']
                    for row in selected)):
        raise ValueError('Large-K source/output screen invariants failed')
    result = dict(kind='torchgwas.staged_device_significant_large_k_probe.v1',
        input_identity=identity, mount=mount, traits=args.traits,
        stage=stage_audit, screen=screen, matched_source_pricing=matched,
        peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        script_sha256=sha(__file__),
        source_sha256={name: sha(Path('src/torchgwas') / name) for name in (
            'productive_staged_screen.py', 'productive_source_floor.py',
            'layout_significant_device_launch_floor.py',
            'layout_partial_envelope.py', 'selection_geometry.py')},
        scope='Metadata-only significant-pair post-output screen. Eight indexed completions are simulated; no genotype scan, GPU selection, or NPZ write. All service rates and occupancy are synthetic conditional scenarios, not calibrated runtime or a candidate winner.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    run.finish(successful=False)
    print(json.dumps(dict(stage_cpu_seconds=stage_audit['cpu_seconds'],
        screen_cpu_seconds=screen['screen_cpu_seconds'],
        screen_wall_seconds=screen['screen_wall_seconds'],
        evaluated=screen['evaluated_candidates'],
        retained=selected[0]['survivors']['retained'],
        peak_rss_kib=result['peak_process_rss_kib'])), flush=True)


if __name__ == '__main__':
    main()