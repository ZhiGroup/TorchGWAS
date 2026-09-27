"""Read-only current-code timing of the large PGEN JIT memory admission core."""
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_memory_layout import (
    memory_layout, compact_rechunk_memory_layout, compact_shifted_reader_envelope)
from torchgwas.pgen_reader import read_header


SOURCE = Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
OUTPUT = Path('results/large_compact_admission_20260923/report.json')
SIZES = (128, 1024, 4096)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    identity = input_identity(SOURCE)
    if identity['bytes'] != 20_838_552_600:
        raise ValueError('Large source identity differs from the documented job')
    times = []

    def timed(name, call):
        wall, cpu = time.perf_counter(), time.process_time()
        value = call()
        times.append(dict(name=name, wall_seconds=time.perf_counter()-wall,
                          cpu_seconds=time.process_time()-cpu))
        return value

    header = timed('read_header', lambda: read_header(SOURCE))
    if header.variant_ct != 8_086_101 or header.sample_ct != 22_250:
        raise ValueError('Large source dimensions changed')
    layout = timed('compact_memory_layout', lambda: memory_layout(
        SOURCE, SIZES[0], header=header, compact=True, _defer_bases=True))
    rows = []
    for size in SIZES:
        fixed = timed(f'fixed_memory_{size}', lambda size=size:
                      compact_rechunk_memory_layout(layout, size))
        shifted = timed(f'shifted_memory_{size}', lambda size=size:
                        compact_shifted_reader_envelope(layout,
                            (0, header.variant_ct), size))
        rows.append(dict(chunk_markers=size, logical_chunks=fixed['logical_chunks'],
                         fixed_read_max_bytes=fixed['chunks'][0]['record_payload_bytes'],
                         shifted_read_max_bytes=shifted[0],
                         shifted_scratch_max_bytes=shifted[1]))
    if input_identity(SOURCE) != identity:
        raise ValueError('Large source changed during admission timing')
    root = Path(__file__).resolve().parents[1]
    report = dict(source=identity, samples=header.sample_ct,
                  variants=header.variant_ct, sizes=rows, times=times,
                  peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  source_sha256={name: digest(root/'src'/'torchgwas'/name)
                                 for name in ('pgen_reader.py', 'pgen_memory_layout.py')},
                  script_sha256=digest(Path(__file__)),
                  scope='Read-only server-local PGEN metadata, no genotype payload, GPU or GWAS. Ordered process CPU/wall times under uncontrolled page-cache and shared-server load; excludes complete public preparation, profile validation, tensor memory and output. No cold startup speedup claim.')
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open('x') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
