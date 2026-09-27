"""Where does an autotuned job's time go: API phases versus the chunk flow.

Runs run_linear_gwas(autotune=True) once and prints the API phase clocks, the
tuner's first/last chunk completion (relative to scan start) and its decision.
Time before the first chunk is per-tile setup; time after the last is drain.
"""
import argparse, json, os, time
from pathlib import Path

ENVIRONMENT = dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
                   TORCHGWAS_SCAN_PROFILE='0', TORCHGWAS_BLOCKING_EVENTS='1', NUMPY_MADVISE_HUGEPAGE='0',
                   OMP_NUM_THREADS='4', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_WAIT_POLICY='PASSIVE')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reduce', default='significant')
    parser.add_argument('--threshold', type=float, default=1e-5)
    parser.add_argument('--options', default='{}')
    args = parser.parse_args()
    os.environ.update(ENVIRONMENT)
    from torchgwas.api import run_linear_gwas
    started = time.perf_counter()
    result = run_linear_gwas(genotype=str(args.data/'input.pgen'), phenotype=args.data/'phenotype.npy',
                             covariates=args.data/'covariates.npy', pgen_mode='hardcall', compute_dtype='float32',
                             output_dir=args.out, reduce=args.reduce,
                             significance_threshold=args.threshold if args.reduce == 'significant' else None,
                             autotune=True, autotune_options=json.loads(args.options))
    meta = result.run_metadata
    chunk = meta['autotune']['chunk'] or {}
    print(json.dumps(dict(api_seconds=time.perf_counter()-started, phases=meta['phase_seconds'],
                          layout=meta['autotune']['layout'].get('why'),
                          first_chunk_seconds=chunk.get('first_chunk_seconds'),
                          last_chunk_seconds=chunk.get('last_chunk_seconds'),
                          state=chunk.get('state'), reason=chunk.get('reason'), choice=chunk.get('choice'),
                          segments=[(s['size'], round(s['rows_per_second'])) for s in chunk.get('segments', [])]),
                     indent=1))


if __name__ == '__main__':
    main()
