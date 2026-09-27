import copy
import pytest
from torchgwas.cusolver_memory import jagwas_eigen_factor_workspace, jagwas_factor_workspace
from torchgwas.tensor_memory import eager_memory_plan


def evidence():
    census=dict(torch_version='2.5.1+cu124',cuda_version='12.4',preferred_linalg='_LinalgBackend.Default',
        compute_capability=[8,0],gpu='GPU',host='host',durations_recorded=False,
        library=dict(path='/runtime/libcusolver.so.11',version=[11,6,1],bytes=120220216,mtime_ns=1730487512825222275),
        rows=[dict(traits=7,device_workspace_bytes=160,host_workspace_bytes=0,matrix_allocated=False)])
    profile=dict(torch_version=census['torch_version'],cuda_version=census['cuda_version'],
        preferred_linalg=census['preferred_linalg'],compute_capability=census['compute_capability'],
        gpu_name=census['gpu'],host=census['host'],cusolver_library=copy.deepcopy(census['library']),
        sm_count=108,max_threads_per_sm=2048,cublas_workspace_config=None,cublas_handle_stream_pairs=1)
    return census,profile


def test_exact_query_requests_and_rounding():
    census,profile=evidence();result=jagwas_factor_workspace(census,7,profile)
    assert result['device_requested_bytes']==160
    assert result['device_rounded_bytes']==512
    assert result['host_requested_bytes']==0
    with pytest.raises(ValueError,match='not queried'):jagwas_factor_workspace(census,8,profile)


def test_workspace_increases_factor_preparation_without_counting_it_per_scan(monkeypatch):
    # The census is xpotrf's: the rounding cutoff's Cholesky factor.
    monkeypatch.setenv('TORCHGWAS_JAGWAS_RCOND', '0')
    census,profile=evidence();old=eager_memory_plan(257,13,7,2,2,profile,reduction='jagwas')
    census['rows'][0].update(device_workspace_bytes=2**30,host_workspace_bytes=32768)
    profile['jagwas_factor_workspace_census']=census
    new=eager_memory_plan(257,13,7,2,2,profile,reduction='jagwas')
    assert new['setup_bytes']==4*257*7+new['factor_setup']['distinct_temporary_bytes']+2**30+new['cublas']['total_bytes']
    assert new['scan']==old['scan']
    assert new['device_bytes']==new['setup_bytes']
    assert new['factor_host_workspace_bytes']==32768
    assert new['status']=='incomplete_memory_candidate'
    assert not new['prediction_complete']
    assert any('info/error' in term for term in new['unresolved_memory_terms'])


@pytest.mark.parametrize('change',['runtime','cuda','preferred','gpu','host','library','version','negative','float','duplicate','allocated','timing','shape'])
def test_unknown_or_mismatched_evidence_rejected(change):
    census,profile=evidence()
    if change=='runtime':profile['torch_version']='2.6.0'
    if change=='cuda':profile['cuda_version']='12.6'
    if change=='preferred':profile['preferred_linalg']='Magma'
    if change=='gpu':profile['gpu_name']='another GPU'
    if change=='host':profile['host']='another host'
    if change=='library':profile['cusolver_library']['mtime_ns']+=1
    if change=='version':
        census['library']['version']=[10,0,0];profile['cusolver_library']=copy.deepcopy(census['library'])
    if change=='negative':census['rows'][0]['device_workspace_bytes']=-1
    if change=='float':census['rows'][0]['host_workspace_bytes']=1.5
    if change=='duplicate':census['rows'].append(census['rows'][0])
    if change=='allocated':census['rows'][0]['matrix_allocated']=True
    if change=='timing':census['durations_recorded']=True
    if change=='shape':census['rows'][0]['traits']=True
    with pytest.raises(ValueError):jagwas_factor_workspace(census,7,profile)


def eigen_evidence():
    census,profile=evidence()
    census=dict(census,method='eigen',rows=[dict(traits=7,syevd_device_workspace_bytes=278512,syevd_host_workspace_bytes=0,
        geqrf_rows_queried=7,geqrf_max_device_workspace_bytes=1572992,geqrf_max_host_workspace_bytes=0,
        matrix_allocated=False)])
    return census,profile


def test_eigen_factor_needs_the_larger_of_its_two_workspaces():
    census,profile=eigen_evidence()
    result=jagwas_eigen_factor_workspace(census,7,profile)
    assert result['device_requested_bytes']==1572992
    assert result['device_rounded_bytes']==1573376
    census['rows'][0]['syevd_device_workspace_bytes']=2**30
    assert jagwas_eigen_factor_workspace(census,7,profile)['device_rounded_bytes']==2**30
    with pytest.raises(ValueError,match='not queried'):jagwas_eigen_factor_workspace(census,8,profile)


@pytest.mark.parametrize('change',['rows','allocated','method','negative'])
def test_eigen_census_must_cover_every_kept_row_count(change):
    census,profile=eigen_evidence()
    if change=='rows':census['rows'][0]['geqrf_rows_queried']=6
    if change=='allocated':census['rows'][0]['matrix_allocated']=True
    if change=='method':census.pop('method')
    if change=='negative':census['rows'][0]['syevd_host_workspace_bytes']=-1
    with pytest.raises(ValueError):jagwas_eigen_factor_workspace(census,7,profile)


def test_each_cutoff_reads_only_its_own_census():
    census,profile=eigen_evidence()
    with pytest.raises(ValueError,match='xpotrf'):jagwas_factor_workspace(census,7,profile)


def test_default_eigen_factor_adds_its_workspace_to_preparation_only(monkeypatch):
    monkeypatch.delenv('TORCHGWAS_JAGWAS_RCOND',raising=False)
    census,profile=eigen_evidence();old=eager_memory_plan(257,13,7,2,2,profile,reduction='jagwas')
    assert old['factor_setup']['method']=='eigen'
    assert any('syevd' in term for term in old['unresolved_memory_terms'])
    census['rows'][0]['syevd_device_workspace_bytes']=2**30
    profile['jagwas_eigen_workspace_census']=census
    new=eager_memory_plan(257,13,7,2,2,profile,reduction='jagwas')
    assert new['setup_bytes']==4*257*7+new['factor_setup']['distinct_temporary_bytes']+2**30+new['cublas']['total_bytes']
    assert new['scan']==old['scan']
    assert new['factor_workspace']['method']=='eigen'
    assert not any('syevd' in term for term in new['unresolved_memory_terms'])
    assert any('info tensors' in term for term in new['unresolved_memory_terms'])
