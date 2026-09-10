from __future__ import annotations

import argparse
import json
import resource
import tempfile
import time
from pathlib import Path

import numpy as np

from torchgwas.bed import PlinkBedGenotype, resolve_plink_triplet


def write_synthetic_bed(prefix: Path, n_samples: int, n_variants: int, seed: int) -> Path:
    bed, bim, fam = resolve_plink_triplet(prefix)
    fam.write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_samples)))
    bim.write_text("".join(f"1 rs{i} 0 {i + 1} A G\n" for i in range(n_variants)))
    bytes_per_variant = (n_samples + 3) // 4
    rng = np.random.default_rng(seed)
    with bed.open("wb", buffering=4 << 20) as handle:
        handle.write(b"\x6c\x1b\x01")
        for start in range(0, n_variants, 512):
            count = min(512, n_variants - start)
            handle.write(rng.integers(0, 256, size=(count, bytes_per_variant), dtype=np.uint8).tobytes())
    return bed


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark bounded parallel PLINK BED decoding")
    parser.add_argument("--n-samples", type=int, default=35_365)
    parser.add_argument("--n-variants", type=int, default=12_000)
    parser.add_argument("--chunk-size", type=int, default=1_000)
    parser.add_argument("--workers", default="1,2,4,8")
    parser.add_argument("--seed", type=int, default=20260910)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="torchgwas-bed-benchmark-") as tmpdir:
        bed = write_synthetic_bed(Path(tmpdir) / "synthetic", args.n_samples, args.n_variants, args.seed)
        for workers in [int(value) for value in args.workers.split(",")]:
            reader = PlinkBedGenotype(bed, reader_workers=workers, prefetch_chunks=max(2, workers))
            wall_start = time.perf_counter()
            cpu_start = time.process_time()
            emitted = 0
            checksum = 0.0
            for start, end, chunk in reader.iter_chunks(
                args.chunk_size,
                dtype=np.float32,
                reader_workers=workers,
                prefetch_chunks=max(2, workers),
            ):
                emitted += end - start
                checksum += float(np.nansum(chunk[:1, :]))
            wall = time.perf_counter() - wall_start
            cpu = time.process_time() - cpu_start
            if emitted != args.n_variants:
                raise RuntimeError(f"reader emitted {emitted} variants; expected {args.n_variants}")
            print(
                json.dumps(
                    {
                        "n_samples": args.n_samples,
                        "n_variants": args.n_variants,
                        "bed_bytes": bed.stat().st_size,
                        "chunk_size": args.chunk_size,
                        "reader_workers": workers,
                        "prefetch_chunks": max(2, workers),
                        "wall_seconds": wall,
                        "cpu_seconds": cpu,
                        "cpu_over_wall": cpu / wall,
                        "bed_mib_per_second": bed.stat().st_size / wall / (1 << 20),
                        "max_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                        "checksum": checksum,
                    },
                    sort_keys=True,
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
