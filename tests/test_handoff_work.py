from copy import deepcopy
import pytest
from torchgwas.handoff_work import handoff_wakeup_prices


def fixture():
    return dict(python_version='py',affinity=[12],rows=[dict(kind=kind,pairs=1,ready=dict(producer_cpu_seconds=1.,receiver_cpu_seconds=2.),observations=[dict(blocked_verified=True,publication_to_return_seconds=10.)]) for kind in ['queue','future'] for _ in range(3)])


def test_elapsed_service_excludes_only_already_priced_ready_work():
    result=handoff_wakeup_prices(fixture(),pairs=1,python_version='py',cpu_affinity=[12])
    assert result['wakeup_seconds']==dict(queue=7.,future=8.)


@pytest.mark.parametrize('kind',['unverified','invalid','version','affinity','repeats'])
def test_handoff_prices_require_verified_matching_context(kind):
    probe=deepcopy(fixture())
    if kind=='unverified':probe['rows'][0]['observations'][0]['blocked_verified']=False
    elif kind=='invalid':probe['rows'][0]['observations'][0]['publication_to_return_seconds']=float('nan')
    elif kind=='version':probe['python_version']='different'
    elif kind=='affinity':probe['affinity']=[13]
    elif kind=='repeats':probe['rows'].pop()
    with pytest.raises(ValueError):handoff_wakeup_prices(probe,pairs=1,python_version='py',cpu_affinity=[12])
