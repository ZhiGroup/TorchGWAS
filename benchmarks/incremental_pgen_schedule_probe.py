"""Read-only matched cost probe for bounded whole-header PGEN accounting."""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.incremental_pgen_schedule import IncrementalPgenSchedule
from torchgwas.pgen_work_bounds import PgenHeaderWork


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _measure(build):
    wall = time.perf_counter()
    cpu = time.process_time()
    value = build()
    return value, dict(wall_seconds=time.perf_counter() - wall,
                       cpu_seconds=time.process_time() - cpu,
                       peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _staged(header, start, stop, chunk, records_per_step):
    plan = IncrementalPgenSchedule(header, start, stop, chunk,
        records_per_step=records_per_step,
        max_chunks_per_step=65536)
    steps = []
    while not plan.snapshot()['complete']:
        wall = time.perf_counter()
        cpu = time.process_time()
        progress = plan.advance()
        steps.append(dict(cursor=progress['cursor'],
                          wall_seconds=time.perf_counter() - wall,
                          cpu_seconds=time.process_time() - cpu))
    value = plan.finish()
    return value, steps, plan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int, required=True)
    parser.add_argument('--chunk-markers', type=int, required=True)
    parser.add_argument('--records-per-step', type=int, default=1048576)
    parser.add_argument('--rebase-cursor', type=int)
    parser.add_argument('--rebase-chunk-markers', type=int)
    parser.add_argument('--rebase-stop', type=int)
    parser.add_argument('--direct-first', action='store_true')
    parser.add_argument('--cached-signatures', type=int, default=1024)
    parser.add_argument('--cached-bounds', type=int, default=64)
    args = parser.parse_args()
    header, header_cost = _measure(lambda: PgenHeaderWork(args.input,
        max_cached_signatures=args.cached_signatures,
        max_cached_bounds=args.cached_bounds))
    direct_first = None
    if args.direct_first:
        direct_first = _measure(lambda: header.schedule_bounds(
            args.start, args.stop, args.chunk_markers,
            max_records=args.stop - args.start,
            max_chunks=200000))
    first, first_cost = _measure(lambda: _staged(
        header, args.start, args.stop, args.chunk_markers,
        args.records_per_step))
    direct, direct_cost = _measure(lambda: header.schedule_bounds(
        args.start, args.stop, args.chunk_markers,
        max_records=args.stop - args.start,
        max_chunks=200000))
    second, second_cost = _measure(lambda: _staged(
        header, args.start, args.stop, args.chunk_markers,
        args.records_per_step))
    def same(a, b):
        return {key: value for key, value in a.items() if key != 'scope'} == {
            key: value for key, value in b.items() if key != 'scope'}
    if (not same(first[0], direct) or not same(second[0], direct) or
            (direct_first is not None and not same(direct_first[0], direct))):
        raise ValueError('Incremental and direct source schedules differ')
    rebase = None
    if any(value is not None for value in (args.rebase_cursor,
            args.rebase_chunk_markers, args.rebase_stop)):
        if args.rebase_cursor is None or args.rebase_chunk_markers is None:
            raise ValueError('Both rebase cursor and chunk size required')
        rebase_stop = args.stop if args.rebase_stop is None else args.rebase_stop
        rebased, rebase_cost = _measure(lambda: second[2].rebase(
            args.rebase_cursor, args.rebase_chunk_markers,
            stop=rebase_stop))
        direct_rebased, direct_rebase_cost = _measure(lambda: header.schedule_bounds(
            args.rebase_cursor, rebase_stop, args.rebase_chunk_markers,
            max_records=rebase_stop - args.rebase_cursor,
            max_chunks=200000))
        if not same(rebased, direct_rebased):
            raise ValueError('Rebased and direct future schedules differ')
        rebase = dict(cursor=args.rebase_cursor, stop=rebase_stop,
                      chunk_markers=args.rebase_chunk_markers,
                      staged=rebase_cost, direct=direct_rebase_cost,
                      exact_work_match=True)
    report = dict(kind='torchgwas.incremental_pgen_schedule_probe.v1',
        input_identity=header.input_identity,
        input_path=str(Path(args.input).resolve()),
        variant_range=[args.start, args.stop],
        chunk_markers=args.chunk_markers,
        records_per_step=args.records_per_step,
        header_cache=dict(max_cached_signatures=args.cached_signatures,
                          max_cached_bounds=args.cached_bounds,
                          signature_work=header.cache_info(),
                          source_bounds=header.bounds_cache_info()),
        chunks=direct['chunk_count'],
        ld_replays=direct['ld_replay_count'],
        header=header_cost,
        staged_first=dict(first_cost, steps=first[1]),
        direct_first=(None if direct_first is None else direct_first[1]),
        direct=direct_cost,
        staged_second=dict(second_cost, steps=second[1]),
        rebase=rebase,
        exact_work_match=True,
        source_sha256={name: _sha(Path(__file__).parents[1] / name)
            for name in ('benchmarks/incremental_pgen_schedule_probe.py',
                         'src/torchgwas/incremental_pgen_schedule.py',
                         'src/torchgwas/pgen_work_bounds.py')},
        scope='Read-only indexed PGEN metadata accounting; no genotype payload, GPU, writer, cold process startup or association throughput measured. Matched same-header staged/direct/staged order is descriptive, not an end-to-end JIT timing guarantee.')
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in
        ('kind', 'input_path', 'chunks', 'ld_replays', 'exact_work_match',
         'header', 'direct')}, sort_keys=True))
    print(json.dumps({name: dict(wall_seconds=report[name]['wall_seconds'],
                                 cpu_seconds=report[name]['cpu_seconds'],
                                 max_step_wall_seconds=max(
                                     row['wall_seconds'] for row in report[name]['steps']),
                                 steps=len(report[name]['steps']))
                      for name in ('staged_first', 'staged_second')},
                     sort_keys=True))
    if direct_first is not None:
        print(json.dumps(dict(direct_first=direct_first[1]), sort_keys=True))
    if rebase is not None:
        print(json.dumps(rebase, sort_keys=True))


if __name__ == '__main__':
    main()
