"""Where a one-GPU significant-pairs scan spends host time between chunks.

Run under a sampling profiler that launches this process (py-spy record ...
-- python this.py), since the lab hosts allow ptrace only of descendants.
One GPU, K = 8,192 synthetic traits, p < 1e-5 with device selection, a prefix
of the full-scale hardcall PGEN. Prints the API wall time, result rows and the
process CPU time, so the profile can be read against the GPU floor.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time


def legacy_selector(beta, t, status, variant_df, critical, *, start=0, max_cells=1 << 20):
    """The per-block selector before 2026-09-26, for A/B runs: 1M-cell blocks,
    one nonzero and five blocking copies each (47 syncs per 1024 x 8192 chunk)."""
    import numpy as np
    import torch
    from torchgwas.selection_geometry import device_selection_shape
    rows, traits = t.shape
    width, height, _ = device_selection_shape(rows, traits, max_cells)
    for first in range(0, rows, height):
        last = min(rows, first + height)
        df = variant_df[first:last]
        limits = critical[df.to(torch.int64).clamp(0, len(critical) - 1)]
        valid = (status[first:last] == 0) & (df > 0) & torch.isfinite(df)
        for left in range(0, traits, width):
            right = min(traits, left + width)
            values = t[first:last, left:right]
            keep = torch.isfinite(values) & (values.abs() >= limits[:, None]) & valid[:, None]
            indices = keep.nonzero(as_tuple=False)
            if not indices.numel():
                yield (start + first, start + last, np.empty(0, np.int64), np.empty(0, np.int64),
                       np.empty(0, np.float32), np.empty(0, np.float32), np.empty(0, np.float32))
                continue
            ri, ti = indices[:, 0], indices[:, 1]
            yield (start + first, start + last, (ri + first + start).cpu().numpy(), (ti + left).cpu().numpy(),
                   beta[first:last, left:right][ri, ti].cpu().numpy(), values[ri, ti].cpu().numpy(),
                   df[ri].cpu().numpy())


def tree_digest(root):
    digest = hashlib.sha256()
    for path in sorted(Path(root).rglob('*')):
        if path.is_file() and path.suffix != '.json':
            digest.update(str(path.relative_to(root)).encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def content_digest(root):
    """Order-independent digest: every column after sorting by (variant, trait)."""
    import numpy as np
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    _, parts = open_indexed_sumstats(root)
    parts = list(parts)
    columns = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    trait_key = next(key for key in columns if key.startswith('trait'))
    order = np.lexsort((columns[trait_key], columns['variant_index']))
    digest = hashlib.sha256()
    for key in sorted(columns):
        digest.update(key.encode())
        digest.update(np.ascontiguousarray(columns[key][order]).tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, default=Path('/data/zxie3/torchgwas_bench/full_scale_20260925'))
    parser.add_argument('--variants', type=int, default=2_000_000)
    parser.add_argument('--chunk', type=int, default=1024)
    parser.add_argument('--readers', type=int, default=4)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--shards', type=int, default=1, help='variant shards over cuda:0..N-1 (readers scale by N)')
    parser.add_argument('--backend', choices=('host', 'device'), default='device')
    parser.add_argument('--selector', choices=('production', 'legacy'), default='production',
                        help='device selector: production (one nonzero per chunk) or the earlier per-block one')
    parser.add_argument('--cache-dir', type=Path, default=None, help='genotype metadata cache (skips the pvar parse)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    os.environ['TORCHGWAS_SIGNIFICANCE_BACKEND'] = args.backend
    from torchgwas.api import run_linear_gwas
    if args.selector == 'legacy':
        import torchgwas.reduce
        torchgwas.reduce.device_significant_pairs = legacy_selector
    manifest = json.loads((args.data / 'manifest.json').read_text())
    layout = (dict(device=args.device, reader_workers=args.readers) if args.shards == 1 else
              dict(variant_devices=[f'cuda:{i}' for i in range(args.shards)], reader_workers=args.readers * args.shards))
    with tempfile.TemporaryDirectory(dir=args.out.parent, prefix='significant-profile-') as scratch:
        started = time.perf_counter()
        result = run_linear_gwas(
            genotype=manifest.get('genotype', str(args.data / 'input.pgen')),
            genotype_format=manifest.get('genotype_format', 'auto'), pgen_mode='hardcall',
            phenotype=args.data / 'phenotype.npy', covariates=args.data / 'covariates.npy',
            compute_dtype='float32', output_dir=Path(scratch) / 'out', **layout,
            reduce='significant', significance_threshold=1e-5, chunk_size=args.chunk,
            prefetch_chunks=4, variant_range=(0, args.variants), genotype_cache_dir=args.cache_dir)
        elapsed = time.perf_counter() - started
        meta = result.run_metadata
        digest = tree_digest(Path(scratch) / 'out' / 'sumstats')
        content = content_digest(Path(scratch) / 'out' / 'sumstats')
    report = dict(variants=args.variants, chunk=args.chunk, readers=args.readers, device=args.device, backend=args.backend,
                  shards=args.shards, selector=args.selector, sumstats_sha256=digest, content_sha256=content,
                  api_seconds=elapsed, rows=meta['n_result_rows'],
                  write=meta.get('sumstats_write'), process_cpu=list(os.times()[:2]))
    args.out.write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps({key: report[key] for key in ('shards', 'selector', 'api_seconds', 'rows', 'process_cpu', 'content_sha256')}), (report['write'] or {}).get('scan_and_write_seconds'))


if __name__ == '__main__':
    main()
