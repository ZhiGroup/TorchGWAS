"""Read-only full PGEN header-window cost with a bounded signature cache."""
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork


SOURCE = Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
TARGET = Path('results/large_window_cache_20260923')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    TARGET.mkdir(parents=True, exist_ok=False)
    before = input_identity(SOURCE)
    began, cpu = time.perf_counter(), time.process_time()
    header = PgenHeaderWork(SOURCE, max_cached_signatures=16384)
    parse = dict(wall_seconds=time.perf_counter() - began,
                 cpu_seconds=time.process_time() - cpu)
    markers = int(header._header.variant_ct)
    began, cpu = time.perf_counter(), time.process_time()
    window = header.window(0, markers, 1024, max_records=markers,
                           max_chunks=10000, max_signatures=2048)
    expanded = dict(wall_seconds=time.perf_counter() - began,
                    cpu_seconds=time.process_time() - cpu)
    rows = window['chunks']
    counts = dict(chunks=len(rows), read_bytes=sum(row['read_bytes'] for row in rows),
                  decode_input_bytes=sum(row['decode_input_bytes'] for row in rows),
                  ld_replay_count=sum(row['ld_replay'] is not None for row in rows))
    assert input_identity(SOURCE) == before
    root = Path(__file__).parents[1]
    report = dict(source=before, header_parse=parse, expanded=expanded,
                  counts=counts, signature_cache=header.cache_info(),
                  peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  source_model_sha256=sha(root / 'src/torchgwas/pgen_work_bounds.py'),
                  script_sha256=sha(__file__),
                  scope=__doc__ + ' No payload read or GWAS; ordered shared-server '
                        'timing and process peak RSS only.')
    with (TARGET / 'report.json').open('x') as out:
        json.dump(report, out, indent=2)
    print(json.dumps(dict(expanded=expanded, counts=counts,
                          cache=report['signature_cache'],
                          peak_rss_kib=report['peak_rss_kib'])), flush=True)


if __name__ == '__main__':
    main()
