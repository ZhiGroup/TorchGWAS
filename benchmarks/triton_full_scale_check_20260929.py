"""The default statistics (Triton kernels, Triton tail, packed rows) against Torch, on the full cohort.

A correctness check, not a timing: the first `--variants` variants of the
full-scale store (22,250 samples), dense output at the store's own K and
min-p, each run with TORCHGWAS_STATS_BACKEND=torch (int8 rows, the Torch
statistics and the device tail) and with the defaults. Reports the largest
beta, t and -log10 P differences and, for min-p, how often the winner
agrees (float32 near-ties may flip; the value then agrees).

    python benchmarks/triton_full_scale_check_20260929.py --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --device cuda:0
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


def run(args, backend, mode, out):
    env = dict(os.environ, TORCHGWAS_PGEN_BACKEND='native', NUMPY_MADVISE_HUGEPAGE='0')
    for key in ('TORCHGWAS_NATIVE_STATS', 'TORCHGWAS_PGEN_PACKED', 'TORCHGWAS_STATS_BACKEND'):
        env.pop(key, None)
    if backend == 'torch':
        env.update(TORCHGWAS_STATS_BACKEND='torch', TORCHGWAS_PGEN_PACKED='0')
    code = f"""
from pathlib import Path
from torchgwas.api import run_linear_gwas
data = Path({str(args.data)!r})
run_linear_gwas(genotype=str(data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                genotype_cache_dir=str(data/'metadata_cache'), phenotype=data/'phenotype.npy',
                covariates=data/'covariates.npy', output_dir={str(out)!r}, compute_dtype='float32',
                chunk_size=4096, prefetch_chunks=4, reader_workers=8, device={args.device!r},
                variant_range=(0, {args.variants}), {'reduce="min-p",' if mode == 'min-p' else ''})
"""
    subprocess.run([sys.executable, '-c', code], check=True, env=env)


def dense(directory):
    from torchgwas.sumstats import open_binary_sumstats
    beta, t, logp, _ = open_binary_sumstats(directory/'sumstats')
    return np.asarray(beta), np.asarray(t), np.asarray(logp, dtype=np.float64)


def min_p(directory):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    _, parts = open_indexed_sumstats(directory/'sumstats')
    parts = list(parts)
    values = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    order = np.argsort(values['variant_index'], kind='stable')
    return {key: value[order] for key, value in values.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--variants', type=int, default=200_000)
    parser.add_argument('--out', type=Path, default=Path('/tmp/triton_check'))
    parser.add_argument('--modes', nargs='+', default=['full', 'min-p'])
    args = parser.parse_args()
    for mode in args.modes:
        outputs = {}
        for backend in ('torch', 'default'):
            outputs[backend] = args.out/f'{mode}_{backend}'
            shutil.rmtree(outputs[backend], ignore_errors=True)
            run(args, backend, mode, outputs[backend])
        if mode == 'full':
            (b0, t0, l0), (b1, t1, l1) = dense(outputs['torch']), dense(outputs['default'])
            finite = np.isfinite(t0) & np.isfinite(t1)
            report = dict(mode=mode, cells=int(t0.size), finite_agree=bool(np.array_equal(np.isfinite(t0), np.isfinite(t1))),
                          beta_abs=float(np.abs(b1 - b0)[finite].max()), t_abs=float(np.abs(t1 - t0)[finite].max()),
                          t_rel=float((np.abs(t1 - t0) / np.maximum(np.abs(t0), 1))[finite].max()),
                          logp_abs=float(np.abs(l1 - l0)[finite].max()),
                          logp_rel=float((np.abs(l1 - l0) / np.maximum(l0, 1))[finite].max()),
                          max_logp=float(l0[finite].max()))
        else:
            a, b = min_p(outputs['torch']), min_p(outputs['default'])
            same = a['trait_index'] == b['trait_index']
            report = dict(mode=mode, variants=int(len(a['variant_index'])),
                          rows_agree=bool(np.array_equal(a['variant_index'], b['variant_index'])),
                          winner_agreement=float(same.mean()),
                          logp_rel=float((np.abs(b['neg_log10_p'] - a['neg_log10_p'])
                                          / np.maximum(a['neg_log10_p'], 1)).max()),
                          t_abs_same_winner=float(np.abs(b['t_stat'] - a['t_stat'])[same].max()))
        print(json.dumps(report), flush=True)
        for directory in outputs.values():
            shutil.rmtree(directory, ignore_errors=True)


if __name__ == '__main__':
    main()
