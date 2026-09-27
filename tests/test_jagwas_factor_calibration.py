"""Real independent observations feed preparation; remaining prices are controls."""
import copy
import json
from pathlib import Path
import pytest
from torchgwas.jagwas_preparation import attach_factor_calibration,factor_service
from torchgwas.jagwas_candidate import jagwas_candidate_runtime
from test_jagwas_preparation import input_path,choice,build
from test_jagwas_actual_candidate import writer_prices

FIXTURES=Path(__file__).parent/'fixtures'


def observations(index):
    return (json.loads((FIXTURES/f'jagwas_factor_phases_cuda{index}.json').read_text()),
            json.loads((FIXTURES/f'jagwas_fp64_capacity_cuda{index}.json').read_text()))


@pytest.mark.parametrize('count',[1,2])
def test_complete_real_components_attach_without_relabeling_synthetic_prices(input_path,count):
    c=choice(input_path,count);indices=[0,2][:count]
    c['devices']=['cuda:'+str(i) for i in indices]
    for tile,index in zip(c['tiles'],indices):
        tile['device']='cuda:'+str(index);bank,capacity=observations(index);before=copy.deepcopy(tile['profile'])
        tile['profile']=attach_factor_calibration(before,bank,capacity,device=tile['device'],factor_wait_cpu_fraction=1.)
        assert tile['profile']['gpu_resources']['fp64_flops_per_second']==capacity['resources']['fp64_flops_per_second']
        assert tile['profile']['gpu_resources']['fp32_flops_per_second']==before['gpu_resources']['fp32_flops_per_second']
        assert tile['profile']['factor_calibration_transfer_qualified'] is False
        assert len(tile['profile']['factor_observation_sha256'])==2
        scalar=factor_service(2049,512,tile['profile'],library_arithmetic='scalar')
        tensor=factor_service(2049,512,tile['profile'],library_arithmetic='tensor')
        assert scalar['seconds']>tensor['seconds']>0
    constructed=build(c)
    report=jagwas_candidate_runtime(c,writer_prices(),preparation=constructed['preparation'],occupancy='dense',host_serial_fraction=.5)
    assert report['retained_variants']==1025 and report['factor_preparation']['instances']==count
    assert not report['automatic_selection_ready']


@pytest.mark.parametrize('fault',['source','context','device','phase_missing','phase_duplicate','phase_clock','phase_summary',
    'capacity_missing','capacity_duplicate','work','rate','summary','flags','kernel','correctness'])
def test_corrupted_or_mixed_observations_refused(input_path,fault):
    bank,capacity=observations(0);device='cuda:0'
    p=choice(input_path,1)['tiles'][0]['profile']
    if fault=='source':bank['source_sha256']={}
    if fault=='context':capacity['affinity']=[0]
    if fault=='device':device='cuda:2'
    if fault=='phase_missing':bank['records'].pop()
    if fault=='phase_duplicate':bank['records'][0]=bank['records'][1]
    if fault=='phase_clock':bank['records'][0]['wall_seconds']+=1
    if fault=='phase_summary':bank['phases']['qr']['cpu_seconds']*=2
    if fault=='capacity_missing':capacity['rows'].pop()
    if fault=='capacity_duplicate':capacity['rows'][0]=capacity['rows'][1]
    if fault=='work':capacity['rows'][0]['dimension']=2048
    if fault=='rate':capacity['rows'][0]['flops_per_second']*=2
    if fault=='summary':capacity['resources']['fp64_flops_per_second']*=2
    if fault=='flags':capacity['selections']['scalar']['instruction_flags']=capacity['selections']['tensor']['instruction_flags']
    if fault=='kernel':capacity['kernels']['scalar']=capacity['kernels']['tensor']
    if fault=='correctness':capacity['max_abs_cpu_reference_error']['scalar']=1.
    with pytest.raises(ValueError):attach_factor_calibration(p,bank,capacity,device=device,factor_wait_cpu_fraction=1.)
