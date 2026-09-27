"""Audit header work bounds against payload counts; no GWAS timing calibration."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.decoder_work import decoder_work, native_read_layout
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.pgen_work_census import census


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]/'src'/'torchgwas'
    hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(source.rglob('*')) if p.suffix in ('.py', '.cpp')}
    identity = input_identity(args.input)
    wall = time.perf_counter(); cpu = time.process_time()
    index = PgenHeaderWork(args.input)
    setup = dict(wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu)
    n, m = int(index._header.sample_ct), int(index._header.variant_ct)
    rows = []
    for size in (128, 256, 512):
        for start in range(0, m, size):
            stop = min(start+size, m)
            wall = time.perf_counter(); cpu = time.process_time()
            bounds = index.bounds(start, stop)
            cost = dict(wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu)
            exact = census(args.input, size, variant_range=(start, stop))
            work = decoder_work(exact, 'torch_native_int8', restart_ld_bases=True)
            for key in bounds['source_units'].keys() | work['source_units'].keys():
                lo, hi = bounds['source_units'].get(key, [0, 0])
                assert lo <= work['source_units'].get(key, 0) <= hi, (start, size, key)
            for key in ('read_bytes', 'decode_input_bytes'):
                assert bounds[key] == native_read_layout(exact)[key]
            assert bounds['native_ld_base_update_bytes'] == work['native_ld_base_update_bytes']
            rows.append(dict(chunk_size=size, variant_range=[start, stop], calculation=cost,
                             work=bounds['structural_work'], all_exact_counts_enclosed=True))
    assert input_identity(args.input) == identity
    summary = {name: dict(median=statistics.median(row['calculation'][name] for row in rows),
                         maximum=max(row['calculation'][name] for row in rows))
               for name in ('wall_seconds', 'cpu_seconds')}
    native_source = source.parents[1]/'native'/'pgen_decode.c'
    native_library = Path(os.environ['TORCHGWAS_PGEN_LIBRARY']).resolve(strict=True)
    report = dict(samples=n, markers=m, input_identity=identity, index_setup=setup, steps=rows,
                  calculation_summary=summary, source_sha256=hashes,
                  native_source_sha256=hashlib.sha256(native_source.read_bytes()).hexdigest(),
                  native_library=dict(path=str(native_library), sha256=hashlib.sha256(native_library.read_bytes()).hexdigest()),
                  cpu_affinity=sorted(os.sched_getaffinity(0)),
                  benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_mount=subprocess.check_output(['findmnt', '-T', args.input, '-J'], text=True),
                  scope='All range intervals checked against independent exact payload censuses. Timings measure header arithmetic only; setup excluded from step cost. This audit deliberately reads payloads for verification and is not a production startup or a throughput benchmark.')
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(samples=n, markers=m, steps=len(rows), index_setup=setup,
                          calculation_summary=summary, source_files=len(hashes))))


if __name__ == '__main__':
    main()
