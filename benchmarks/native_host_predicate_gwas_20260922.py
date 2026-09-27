"""Public GPU output equivalence; no inference of end-to-end speedup."""
import hashlib
import json
import os
from pathlib import Path
import threading
import time
from unittest.mock import patch
import numpy as np
import torch
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity, execution_context
from torchgwas.sumstats_indexed import open_indexed_sumstats
from torchgwas import native_host_predicate as native


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def output(path):
    manifest, parts = open_indexed_sumstats(path)
    parts = list(parts)
    if not parts:
        assert manifest['rows'] == 0
        return {}
    fields = ('variant_index', 'trait_index', 't_stat', 'df') + (('beta',) if 'beta' in parts[0] else ())
    values = {key:np.concatenate([p[key] for p in parts]) for key in fields}
    order = np.lexsort((values['trait_index'], values['variant_index']))
    return {key:array[order] for key,array in values.items()}


def main():
    fixture = Path('/data/zxie3/torchgwas_adaptive_candidate_fixture_v1_20260922')
    root = Path('results/native_host_predicate_gwas_20260922')
    target = Path('/data/zxie3/torchgwas_native_host_predicate_20260922')
    root.mkdir(parents=True, exist_ok=False); target.mkdir(parents=True, exist_ok=False)
    torch.set_num_threads(2); torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32 = False
    source = source_identity(); inputs = {n:sha(fixture/n) for n in
        ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    y = np.load(fixture/'phenotype.npy'); covariates = np.load(fixture/'covariates.npy')
    cases = [('single_sparse', .02, 'beta+t', {}),
        ('tiled_sparse', .02, 't', dict(trait_block=193, trait_devices=['cuda:1','cuda:2'])),
        ('tiled_dense', 1., 'beta+t', dict(trait_block=193, trait_devices=['cuda:1','cuda:2'])),
        ('tiled_empty', 1e-30, 't', dict(trait_block=193, trait_devices=['cuda:1','cuda:2']))]
    reports = []; original = native.fill_mask; lock = threading.Lock()
    for i,(label, threshold, fields, options) in enumerate(cases):
        reference = None
        for backend in (['numpy','native'] if i%2 == 0 else ['native','numpy']):
            os.environ['TORCHGWAS_HOST_PREDICATE'] = backend
            observed = []
            def fill(values, limits, mask):
                result = original(values, limits, mask)
                with lock: observed.append(dict(shape=list(values.shape), native=result))
                return result
            context = execution_context(['cuda:1','cuda:2'], input_path=fixture/'input.pgen', output_path=target)
            if backend == 'native': assert context['host_predicate'] == native.context()
            else: assert 'host_predicate' not in context
            directory = target/(label+'_'+backend)
            began = time.perf_counter()
            with patch.object(native, 'fill_mask', fill):
                result = run_linear_gwas(str(fixture/'input.pgen'), y, covariates,
                    genotype_format='pgen', device='cuda:1', chunk_size=128,
                    reader_workers=3, prefetch_chunks=3, reduce='significant',
                    significance_threshold=threshold, output_dir=directory,
                    sumstats_format='binary', sumstats_fields=fields,
                    sumstats_queue_depth=2, sumstats_fsync=True, **options)
            elapsed = time.perf_counter()-began
            data = output(directory/'sumstats')
            if reference is None: reference = data
            else:
                assert set(reference) == set(data)
                for name in reference: np.testing.assert_array_equal(data[name], reference[name])
            if backend == 'native':
                assert observed and all(r['native'] for r in observed)
                assert sum(np.prod(r['shape']) for r in observed) == 4097*512
            else: assert not observed
            row = dict(label=label, backend=backend, threshold=threshold, fields=fields,
                output_dir=str(directory), execution_context=context, predicate_calls=observed,
                api_seconds=elapsed, result_rows=result.run_metadata['n_result_rows'],
                writer=result.run_metadata['sumstats_write'])
            reports.append(row)
            print(json.dumps({k:row[k] for k in ('label','backend','api_seconds','result_rows')}), flush=True)
    assert source_identity() == source
    assert inputs == {n:sha(fixture/n) for n in inputs}
    (root/'report.json').write_text(json.dumps(dict(runs=reports, all_outputs_exact=True,
        source_sha256=source, inputs=inputs, harness_sha256=sha(__file__), scope=__doc__), indent=2)+'\n')


if __name__ == '__main__': main()
