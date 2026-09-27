"""Read-only large-PGEN metadata cost probe for the post-output JIT screen.

The output events are simulated, not a GWAS. All component rates below are
synthetic; only exact indexed source work and planner CPU/wall/RSS are measured.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.productive_source_stage import ProductiveSourceStage
from torchgwas.productive_staged_screen import productive_staged_partial_screen
from torchgwas.sumstats import DenseWriteProgress


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def candidate(name, chunk, partitions, source_profile):
    devices = {row['device'] for row in partitions}
    writer = dict(cpu_fraction=.5,
        writer_copy_service=dict(cpu_seconds_per_byte=1e-10,
                                 cpu_seconds_per_call=1e-6),
        process_units=dict(bytearray_zero_bytes=1e-10),
        executor_cpu_seconds=1e-6,
        writeback_service=dict(pagecache_seconds_per_byte=1e-10,
            storage_seconds_per_byte=1e-9, submit_seconds=1e-4,
            wait_seconds=1e-4, fadvise_seconds=1e-4,
            fadvise_eviction_seconds_per_byte=1e-10),
        fsync_seconds=.001)
    return dict(id=name, chunk_markers=chunk, partitions=partitions,
        partition_axis='trait',
        source_profiles={row['id']: source_profile for row in partitions},
        compute_options=dict(covariate_rank=3,
            shared_h2d_bytes_per_second=1e10,
            per_device_h2d_bytes_per_second={device: 1e10 for device in devices},
            peak_fp32_flops_per_second={device: 1e13 for device in devices}),
        output_options=dict(store_beta=True,
            shared_d2h_bytes_per_second=1e10,
            per_device_d2h_bytes_per_second={device: 1e10 for device in devices},
            output_bytes_per_second=1e9,
            dense_writer_options=dict(block_bytes=64 << 20, queue_depth=2,
                borrow_chunks=False, fsync=True, writeback_bytes=64 << 20,
                sync_file_range=True, store_variant_df=False)),
        mode_service_options=dict(writer_profiles={device: writer
                                                   for device in devices}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--variants', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Probe output already exists')
    source = args.input.resolve(strict=True)
    identity = input_identity(source)
    mount = subprocess.check_output(['findmnt', '-T', str(source),
        '-o', 'TARGET,SOURCE,FSTYPE', '-n'], text=True).strip()
    if not mount.startswith('/data '):
        raise ValueError('Probe requires server-local /data PGEN')
    initial = 128
    config = dict(records_per_step=1_048_576, max_steps=8,
        max_cpu_seconds=30., max_window_seconds=120.,
        max_retained_bytes=16 << 20, extra_host_reserve_bytes=528 << 20,
        max_cached_signatures=32768, max_cached_bounds=256)
    stage = ProductiveSourceStage(
        lambda: PgenHeaderWork(source, max_cached_signatures=32768,
                               max_cached_bounds=256),
        identity, initial, config)
    run = ProductiveTuningRun([dict(id='all', device='cuda:0',
        variant_range=[0, args.variants], trait_range=[0, 16])],
        chunk_sizes=[initial, 1024], initial=initial,
        budget=IncrementalPlanningBudget(max_steps=2,
            max_cpu_seconds=60., max_window_seconds=120.))
    control = run.for_partition('all')
    cursor = 0
    for _ in range(8):
        width = control(cursor, args.variants, 1024)
        end = cursor + width
        run.output_written(DenseWriteProgress(cursor, end, (0, 16),
            width * 16 * 8, end, time.perf_counter(), 'metadata-probe', 'cuda:0'))
        stage.output_written()
        cursor = end
    deadline = time.monotonic() + 120.
    while not stage.snapshot()['complete'] and time.monotonic() < deadline:
        observed = stage.snapshot()
        if observed['stop_reason'] is not None:
            raise RuntimeError('Stage stopped: ' + str(observed['stop_reason']) + ' ' +
                               str(observed['steps'][-1:] if observed['steps'] else ''))
        time.sleep(.05)
    stage_audit = stage.finish()
    if not stage_audit['complete'] or stage.stage.header._header.variant_ct != args.variants:
        raise ValueError('Exact large-source stage did not complete')
    snapshot = run.snapshot()
    frontier = unissued_frontier(snapshot, source_identity=identity,
        reduction=None, total_traits=16, job_variant_range=[0, args.variants])
    source_profile = dict(decode_units={name: 1e-9 for name in (
        'copy_packed_byte', 'expand4_int8', 'expand_int8_tail_sample',
        'uleb1', 'uleb2', 'uleb3', 'uleb4', 'uleb5', 'set_category',
        'difflist_group_absolute_id', 'difflist_category_extract',
        'difflist_record_header', 'fill_packed_byte',
        'invert_packed_byte', 'native_onebit_tail_bit_test',
        'onebit8_native')}, cpu_fraction=.5, depth=2, decode_workers=2,
        cpu_available_cores=8., shared_dram_bytes_per_second=1e10,
        read_bytes_per_second=1e9)
    remaining = [cursor, args.variants]
    baseline = candidate('baseline_128', 128,
        [dict(id='all', device='cuda:0', variant_range=remaining,
              trait_range=[0, 16])], source_profile)
    retiled = candidate('retiled_1024', 1024,
        [dict(id='low', device='cuda:0', variant_range=remaining,
              trait_range=[0, 8]),
         dict(id='high', device='cuda:1', variant_range=remaining,
              trait_range=[8, 16])], source_profile)
    screen = productive_staged_partial_screen(frontier, stage,
        [baseline, retiled],
        shared_source_capacities=dict(cpu=8., dram=1e10, input=1e9),
        max_cpu_seconds=60., max_wall_seconds=120.)
    result = dict(kind='torchgwas.staged_partial_screen_probe.v1',
        input_identity=identity, source_mount=mount, stage=stage_audit,
        issue_token=dict(issued_revision=frontier['issued_revision'],
                         written_events=frontier['written_events'],
                         required_cells=frontier['required_cells']),
        screen=screen, peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        script_sha256=sha(__file__),
        source_sha256={name: sha(Path('src/torchgwas') / name) for name in (
            'productive_staged_screen.py', 'productive_source_floor.py',
            'productive_source_stage.py', 'pgen_work_bounds.py')},
        scope='Metadata-only post-output JIT planner-cost probe. Events stand in for eight completed source chunks; no genotypes were scanned, no GPU work or output was written, and all model service rates are synthetic. Exact indexed PGEN work and planner CPU/wall/RSS only; no throughput, candidate winner or switch claim.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    run.finish(successful=False)
    print(json.dumps(dict(stage_cpu_seconds=stage_audit['cpu_seconds'],
        stage_steps=len(stage_audit['steps']),
        screen_cpu_seconds=screen['screen_cpu_seconds'],
        evaluated=screen['evaluated_candidates'],
        stop_reason=screen['stop_reason'],
        peak_rss_kib=result['peak_process_rss_kib'])), flush=True)


if __name__ == '__main__':
    main()
