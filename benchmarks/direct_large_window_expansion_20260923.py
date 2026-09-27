"""Read-only cost and conservation of full PGEN header chunk expansion."""
import gc
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork


SOURCE = Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
TARGET = Path('results/large_window_expansion_20260923')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    TARGET.mkdir(parents=True, exist_ok=False)
    before = input_identity(SOURCE)
    started, cpu = time.perf_counter(), time.process_time()
    header = PgenHeaderWork(SOURCE)
    parsed = dict(wall_seconds=time.perf_counter() - started,
                  cpu_seconds=time.process_time() - cpu)
    markers = int(header._header.variant_ct)
    rows = []
    for size in (1024, 4096):
        started, cpu = time.perf_counter(), time.process_time()
        expanded = header.window(0, markers, size, max_records=markers,
                                 max_chunks=10000, max_signatures=2048)
        expanded_time = dict(wall_seconds=time.perf_counter() - started,
                             cpu_seconds=time.process_time() - cpu)
        count = len(expanded['chunks'])
        read = sum(row['read_bytes'] for row in expanded['chunks'])
        decode = sum(row['decode_input_bytes'] for row in expanded['chunks'])
        units = {}
        for row in expanded['chunks']:
            for name, interval in row['source_units'].items():
                old = units.setdefault(name, [0, 0])
                old[0] += interval[0]
                old[1] += interval[1]
        del expanded
        gc.collect()
        started, cpu = time.perf_counter(), time.process_time()
        aggregate = header.schedule_bounds(0, markers, size, max_records=markers,
                                           max_signatures=65536, max_chunks=10000)
        aggregate_time = dict(wall_seconds=time.perf_counter() - started,
                              cpu_seconds=time.process_time() - cpu)
        assert count == aggregate['chunk_count']
        assert read == aggregate['read_bytes']
        assert decode == aggregate['decode_input_bytes']
        # Independently bounded difflist records can have wider intervals
        # when aggregated by signature; exact primitives must agree.
        mismatched_exact = {name: (value, aggregate['source_units'].get(name))
                            for name, value in units.items()
                            if value[0] == value[1] and
                            value != aggregate['source_units'].get(name)}
        assert not mismatched_exact, mismatched_exact
        rows.append(dict(chunk_markers=size, chunks=count, read_bytes=read,
                         decode_input_bytes=decode, expansion=expanded_time,
                         aggregate=aggregate_time, aggregate_cache=header.cache_info(),
                         peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         aggregate_encloses_chunk_intervals=all(
                             aggregate['source_units'].get(name, [0, 0])[0] <= value[0]
                             and value[1] <= aggregate['source_units'].get(name, [0, 0])[1]
                             for name, value in units.items())))
        print(json.dumps(rows[-1]), flush=True)
        assert rows[-1]['aggregate_encloses_chunk_intervals']
    assert input_identity(SOURCE) == before
    root = Path(__file__).parents[1]
    report = dict(source=before, header_parse=parsed, rows=rows,
                  source_model_sha256=sha(root / 'src/torchgwas/pgen_work_bounds.py'),
                  script_sha256=sha(__file__),
                  scope=__doc__ + ' No payload read or GWAS. Ordered shared-server '
                        'timings inform JIT planning cost, not runtime speedup.')
    with (TARGET / 'report.json').open('x') as out:
        json.dump(report, out, indent=2)


if __name__ == '__main__':
    main()
