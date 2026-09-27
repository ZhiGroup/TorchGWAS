"""Separate-process reuse audit; no new component measurement is performed."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import torch

from direct_jagwas_bounded_execution_20260922 import read_result
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import sha256_file,source_identity


def main():
    parser=argparse.ArgumentParser()
    for name in ('fixture','previous','reference','out','output-data'):
        parser.add_argument('--'+name,required=True)
    args=parser.parse_args();previous=Path(args.previous);out=Path(args.out)
    fixture=Path(args.fixture);destination=Path(args.output_data)
    out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    old=json.loads((previous/'report.json').read_text())
    profile=previous/'profile.json';config=json.loads((previous/'config.json').read_text())
    hashes={str(profile):sha256_file(profile),**old['profile']['component_artifacts']}
    y=np.load(fixture/'phenotype.npy');c=np.load(fixture/'covariates.npy')
    before=time.perf_counter()
    result=run_linear_gwas(fixture/'input.pgen',y,c,pgen_mode='hardcall',compute_dtype='float32',
        reduce='jagwas',output_dir=destination,sumstats_fields='t',sumstats_queue_depth=2,
        autotune_profile=profile,autotune_config=config)
    elapsed=time.perf_counter()-before
    _,actual=read_result(destination,4097);_,expected=read_result(Path(args.reference),4097)
    np.testing.assert_allclose(actual,expected,rtol=6e-5,atol=3e-4)
    audit=result.run_metadata['autotune'];state=audit['productive']
    assert state['finished']['successful'] and state['decisions'][0]['error'] is None
    assert all(p['cursor']==p['variant_range'][1] for p in state['partitions'])
    bound=state['forecast_attempts'][0]['scenarios'][0]['price_evidence']['bindings'][0]
    prior=old['results'][-1]['autotune']['productive']['forecast_attempts'][0]['scenarios'][0]['price_evidence']['bindings'][0]
    for field in ('record_sha256','observed_unix_seconds','created_unix_seconds','expires_unix_seconds'):
        assert bound[field]==prior[field]
    assert bound['age_seconds']>prior['age_seconds']
    assert hashes=={path:sha256_file(path) for path in hashes}
    assert source_identity()==old['source_sha256']
    (out/'report.json').write_text(json.dumps(dict(api_seconds=elapsed,
        first_written_after_api_seconds=state['first_written']-before,
        max_absolute_difference=float(np.max(np.abs(actual-expected))),autotune=audit,
        prior_age_seconds=prior['age_seconds'],age_seconds=bound['age_seconds'],
        preserved_artifact_sha256=hashes,source_sha256=source_identity(),
        benchmark_sha256=sha256_file(__file__),scope=__doc__),indent=2))
    print('COMPLETE separate-process reuse',elapsed,'age',bound['age_seconds'],flush=True)


if __name__=='__main__':main()
