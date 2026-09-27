import math,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'benchmarks'))
from direct_missing_marker_ties import branch_probabilities,apply_missingness

def test_iid_branches_partition_variants():
    b=branch_probabilities(8192,.001)
    assert b['complete']==pytest.approx((1-.001)**8192)
    assert sum(b[k] for k in ('gram','opener','sparse'))==pytest.approx(1)
    assert b['gram']>.9997 and b['sparse']<1e-7
    # Law of total expectation for observed sample count.
    assert b['complete']*8192+b['gram']*b['observed_given_missing']==pytest.approx(8192*.999)

def test_tiny_missing_rate_has_conditional_one_missing_call_limit():
    assert branch_probabilities(8192,1e-12)['observed_given_missing']==pytest.approx(8191,abs=1e-6)

def test_zero_missing_preserves_saved_profile_exactly():
    c={'torch':{'1':{'shape':[8192,1,8,2048,4,4]}},'plink':{'1':{'worker_elapsed_per_variant':.002}}}
    got,info=apply_missingness(c,0,[],0)
    assert got==c and got is not c and info['gram']==0

@pytest.mark.parametrize('rate',[-.01,1.,float('nan')])
def test_invalid_rate(rate):
    with pytest.raises(ValueError):branch_probabilities(8192,rate)


def test_missing_step_work_and_decoder_are_counted_once():
    c={'torch':{'1':{'shape':[100,1,8,2048,4,4]}},'plink':{'1':{'worker_elapsed_per_variant':8.}}}
    record={'samples':100,'predictors':10,'traits':1,'nm_samples':90,
      'timer':'CLOCK_THREAD_CPUTIME_ID per step','service_effective_cores':2.,
      'seconds_per_variant':{key:1. for key in ('nm_mask','expand','fill','syrk','xty','solve','vif','post')}}
    changed,info=apply_missingness(c,.01,[record],.4)
    b=branch_probabilities(100,.01)
    gram=(4+4*b['observed_given_missing']/90)/2+.1
    opener=(3+3*100/90)/2+.1
    expected=b['gram']*gram+b['opener']*opener+b['sparse']*2
    assert changed['plink']['1']['worker_elapsed_per_variant']==pytest.approx(4*expected)
    assert c['plink']['1']['worker_elapsed_per_variant']==8
