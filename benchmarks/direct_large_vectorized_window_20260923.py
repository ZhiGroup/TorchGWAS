"""Read-only full PGEN vectorized chunk bounds versus compact aggregate."""
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork


SOURCE = Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
TARGET = Path('results/large_vectorized_window_20260923')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    TARGET.mkdir(parents=True, exist_ok=False)
    before = input_identity(SOURCE)
    header = PgenHeaderWork(SOURCE)
    markers = int(header._header.variant_ct)
    began, cpu = time.perf_counter(), time.process_time()
    window = header.vectorized_window(0, markers, 1024, max_records=markers,
                                      max_chunks=10000, max_signatures=2048,
                                      max_global_signatures=32768,
                                      max_signature_pairs=10000000)
    expanded = dict(wall_seconds=time.perf_counter() - began,
                    cpu_seconds=time.process_time() - cpu)
    rows = window['chunks']
    units = {}
    for row in rows:
        for name, pair in row['source_units'].items():
            result = units.setdefault(name, [0, 0])
            result[0] += pair[0]
            result[1] += pair[1]
    began, cpu = time.perf_counter(), time.process_time()
    aggregate = header.schedule_bounds(0, markers, 1024,
                                       max_records=markers, max_chunks=10000,
                                       max_signatures=65536)
    aggregate_time = dict(wall_seconds=time.perf_counter() - began,
                          cpu_seconds=time.process_time() - cpu)
    assert len(rows) == aggregate['chunk_count']
    assert sum(row['read_bytes'] for row in rows) == aggregate['read_bytes']
    assert sum(row['decode_input_bytes'] for row in rows) == aggregate['decode_input_bytes']
    assert sum(row['ld_replay'] is not None for row in rows) == aggregate['ld_replay_count']
    assert units == aggregate['source_units']
    assert input_identity(SOURCE) == before
    root = Path(__file__).parents[1]
    record = dict(source=before, vectorized=expanded, compact=aggregate_time,
                  chunks=len(rows), read_bytes=aggregate['read_bytes'],
                  ld_replays=aggregate['ld_replay_count'],
                  signature_cache=header.cache_info(),
                  peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  source_model_sha256=sha(root / 'src/torchgwas/pgen_work_bounds.py'),
                  script_sha256=sha(__file__),
                  scope=__doc__ + ' No payload read or GWAS. Conditional '
                        'source intervals only, not pipeline runtime.')
    with (TARGET / 'report.json').open('x') as out:
        json.dump(record, out, indent=2)
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
