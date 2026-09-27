"""Captured statistics-tail geometry must be priced or explicitly rejected."""
import copy
import json
from pathlib import Path
import pytest
from torchgwas.tensor_service import DeviceService, gemm_work, host_primitive_name, tensor_stage_service
from torchgwas.tensor_work import eager_statistics_work

CAPTURE=json.loads((Path(__file__).parent/'fixtures'/'significant_singleton_n2049_b1_k1.json').read_text())


def test_narrow_singleton_counts_logical_work_and_actual_launches():
    row=CAPTURE['row']
    result=gemm_work(row['N'],row['B'],row['K']+row['C']+1,row['kernels'])
    assert result['useful_flops']==result['issued_flops']==2*2049*4
    assert result['arithmetic_accounting']=='logical_floor_not_verified_issued_work'
    assert result['workspace_write_bytes']==result['workspace_reduce_bytes']==128
    assert result['input_l2_bytes']==4*2049*5
    assert result['reduce_adds']==28 and result['kernel_count']==2
    assert result['launched_ctas']==result['grid_ctas']==8
    assert result['tile'] is None and result['unpriced_terms']
    assert not result['inner_k_tile_verified']
    # Exercise the complete tensor stage, including K=1 reductions/views.
    work=eager_statistics_work(2049,1,1,2,True)
    prices={host_primitive_name(call):1e-6 for call in work['host_calls']}
    resources=DeviceService(1e12,2e12,1e13,1e-6,2e-6,1<<20,108)
    stage=tensor_stage_service(work,resources,row['kernels'],host_primitives=prices)
    assert stage['estimated_span_seconds']>0
    assert stage['gemm']==result


@pytest.mark.parametrize('fault',['samples','markers','width','block','grid','split','specialization',
                                  'reduction_grid','reduction_block','missing_reduction','duplicate_main'])
def test_unobserved_narrow_singleton_contract_is_refused(fault):
    row=copy.deepcopy(CAPTURE['row']);kernels=row['kernels']
    main=next(k for k in kernels if 'gemvNSP_kernel<' in k['name'])
    reduction=next(k for k in kernels if 'splitKreduce_kernel<' in k['name'])
    n,b,width=2049,1,4
    if fault=='samples':n+=1
    if fault=='markers':b+=1
    if fault=='width':width+=1
    if fault=='block':main['geometry']['block']=[32,24,1]
    if fault=='grid':main['geometry']['grid'][0]+=1
    if fault=='split':main['geometry']['grid'][2]=4
    if fault=='specialization':main['name']=main['name'].replace('1, 32, 4, 1024','1, 16, 4, 1024')
    if fault=='reduction_grid':reduction['geometry']['grid'][1]=2
    if fault=='reduction_block':reduction['geometry']['block']=[16,32,1]
    if fault=='missing_reduction':kernels.remove(reduction)
    if fault=='duplicate_main':kernels.append(copy.deepcopy(main))
    with pytest.raises(ValueError):gemm_work(n,b,width,kernels)
