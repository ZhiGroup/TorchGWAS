"""Full-scale dense store with -log10 P: correctness on sampled rows, and the write summary.

Runs one dense scan (fixed layout), then checks the stored neglog10p.f32
against the exact host tail at the store's df for evenly spaced variants,
and the stored t against its own -log10 P ordering. The output is deleted.

    python benchmarks/logp_store_check_20260927.py --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 \\
        --devices cuda:4 cuda:5 --out /dev/shm/zxie3_logp_check --record results/logp_store_check_20260927.json
"""
import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True)
    parser.add_argument('--rows', type=int, default=2000)
    args = parser.parse_args()
    os.environ.update(bench.ENVIRONMENT)
    kwargs = (dict(variant_devices=args.devices, reader_workers=8 * len(args.devices), prefetch_chunks=16)
              if len(args.devices) > 1 else dict(device=args.devices[0], reader_workers=16, prefetch_chunks=32))
    started = time.perf_counter()
    row = bench.child(args.data, args.out, dict(name='check', kwargs=dict(kwargs, chunk_size=4096)))
    wall = time.perf_counter() - started
    from torchgwas.sumstats import open_binary_df, open_binary_sumstats
    from torchgwas.tails import upper_tail_log10_from_t
    beta, t, logp, manifest = open_binary_sumstats(args.out / 'sumstats')
    df = open_binary_df(args.out / 'sumstats')
    n = t.shape[0]
    picks = np.linspace(0, n - 1, min(args.rows, n)).astype(np.int64)
    t_rows = np.stack([np.asarray(t[int(i)], dtype=np.float64) for i in picks])
    logp_rows = np.stack([np.asarray(logp[int(i)], dtype=np.float64) for i in picks])
    df_rows = np.broadcast_to(np.asarray(df), t.shape)
    df_rows = np.stack([np.asarray(df_rows[int(i)], dtype=np.float64) for i in picks])
    finite = np.isfinite(t_rows)
    reference = np.full_like(t_rows, np.nan)
    reference[finite] = upper_tail_log10_from_t(t_rows[finite], df_rows[finite])
    error = np.abs(logp_rows[finite] - reference[finite])
    relative = error / np.maximum(np.abs(reference[finite]), 1e-30)
    record = dict(data=str(args.data), devices=args.devices, wall_seconds=wall,
                  executor_seconds=row['executor_seconds'], api_seconds=row['api_seconds'],
                  arrays=sorted(manifest.get('arrays', [])) if isinstance(manifest.get('arrays'), list)
                  else sorted(manifest.get('arrays', {})), version=manifest.get('version'),
                  sampled_rows=int(len(picks)), sampled_cells=int(finite.sum()),
                  same_nan_pattern=bool(np.array_equal(np.isnan(t_rows), np.isnan(logp_rows))),
                  max_abs_error=float(error.max()), max_rel_error=float(relative.max()),
                  max_logp=float(np.nanmax(logp_rows)),
                  sumstats_write={key: value for key, value in row['sumstats_write'].items()
                                  if key in ('payload_bytes', 'write_seconds', 'cells', 'setup_scan_and_write_seconds')})
    print(json.dumps(record), flush=True)
    args.record.parent.mkdir(parents=True, exist_ok=True)
    args.record.write_text(json.dumps(record, indent=1) + '\n')
    shutil.rmtree(args.out)


if __name__ == '__main__':
    main()
