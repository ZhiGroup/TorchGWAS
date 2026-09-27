"""Full statistics and JAGWAS candidate checks with untimed launch censuses.

Component prices are synthetic controls; these are accounting, not calibration.
"""
import copy
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
import pytest
from torchgwas.jagwas_candidate import jagwas_candidate_runtime,detailed_jagwas_plan
from torchgwas.linear import multigpu_variant_ranges
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.layout_gpu_shape_service import native_layout_gpu_shape_service
from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_census import census
from torchgwas.tensor_work import eager_statistics_work
from torchgwas.tensor_service import host_primitive_name,gemm_work
from torchgwas.reduction_tensor_work import jagwas_projection_service
from test_jagwas_candidate import candidate,preparation
from test_jagwas_tensor_service import resources,prices as projection_prices
from test_jagwas_writer_service import bank,archive
from test_pgen_native_reader import write_pgen

CAPTURE=json.loads((Path(__file__).parent/'fixtures'/'jagwas_chunk_geometry.json').read_text())


@pytest.fixture
def input_path(tmp_path):
    path=tmp_path/'tail.pgen'
    write_pgen(path,(np.arange(1025*2049,dtype=np.uint32).reshape(1025,2049)%3).astype(np.uint8))
    return path


def actual_candidate(path,b,devices):
    result=candidate(path,devices)
    m=census(path,b)['markers']
    work=eager_statistics_work(2049,b,512,2,True)
    host={host_primitive_name(call):1e-6 for call in work['host_calls']}
    for tile,span in zip(result['tiles'],multigpu_variant_ranges(m,b,devices)):
        tile['variant_range']=list(span);tile['data']['markers']=span[1]-span[0]
        tile['data']['encoded']=census(path,b,variant_range=span,include_chunks=True)
        tile['profile'].update(chunk_markers=b,kernel_geometry=copy.deepcopy(CAPTURE['statistics']),
            host_primitives=host,joint_kernel_geometry=copy.deepcopy(CAPTURE['projection']))
    return result


def writer_prices():
    return dict(prices=bank(),archive=archive(),queue_cpu_seconds=dict(put=1e-6,get=1e-6))


@pytest.mark.parametrize('b',[128,256,512])
@pytest.mark.parametrize('devices',[1,2])
def test_actual_statistics_projection_and_tail_compose(input_path,b,devices):
    choice=actual_candidate(input_path,b,devices)
    reports=[torch_scan_work(tile['data'],tile['profile']) for tile in choice['tiles']]
    assert sum(len(work['blocks']) for work in reports)==(1025+b-1)//b
    for tile,work in zip(choice['tiles'],reports):
        assert sum(block['d2h_bytes'] for block in work['blocks'])==17*tile['data']['markers']
        for width,component in work['components'].items():
            stat=next(row for row in CAPTURE['statistics'] if row['B']==int(width))
            joint=next(row for row in CAPTURE['projection'] if row['B']==int(width))
            assert component['kernel_count']==len(stat['kernels'])+len(joint['kernels'])
            assert component['reduction_component']['gemm']['useful_flops']==2*int(width)*512**2
    for occupancy in ['empty','dense']:
        result=jagwas_candidate_runtime(choice,writer_prices(),preparation=preparation(choice),
            occupancy=occupancy,host_serial_fraction=.5)
        factor=result['factor_preparation']
        assert factor['instances']==devices
        assert factor['total_h2d_bytes']==devices*4*2049*512
        # FP64 Gram, then the eigen factor at k = K: eigh (2/3 K^3 + K^3) and QR (K^3 - K^3/3).
        assert factor['total_fp64_multiply_add_flops']==devices*(2*2049*512**2+(2*512**3)//3+512**3+512**3-512**3//3)
        assert result['retained_variants']==(1025 if occupancy=='dense' else 0)
        assert result['parts']==((1025+b-1)//b if occupancy=='dense' else 0)
        assert result['estimated_seconds']>0 and not result['automatic_selection_ready']


def test_compact_shape_service_matches_full_scan_components(input_path):
    choice = actual_candidate(input_path, 256, 2)
    tiles = choice['tiles']
    layout = dict(kind='torchgwas.pgen_layout_source_floor.v1',
                  input_identity=input_identity(input_path),
                  samples=2049, chunk_markers=256,
                  reduction='jagwas', total_traits=512,
                  partitions=[dict(id=str(index), device=tile['device'],
                                   variant_range=tile['variant_range'],
                                   trait_range=tile['trait_range'],
                                   chunks=(tile['data']['markers'] + 255) // 256)
                              for index, tile in enumerate(tiles)])
    compact = native_layout_gpu_shape_service(layout, covariate_rank=2,
        profiles={tile['device']: tile['profile'] for tile in tiles})
    for tile in tiles:
        expanded = torch_scan_work(tile['data'], tile['profile'])
        expected = sum(expanded['components'][block['markers']]
                       ['kernel_service_seconds']
                       for block in expanded['blocks'])
        actual = compact['per_device_service'][tile['device']]
        assert actual['kernel_service_seconds'] == pytest.approx(expected)
        assert actual['chunks'] == len(expanded['blocks'])
    assert compact['distinct_shapes'] <= 4
    assert not compact['prediction_complete']


def test_bound_optimizer_uses_real_geometry_for_all_chunks_and_devices(input_path):
    choices=[actual_candidate(input_path,b,devices) for b in [128,256,512] for devices in [1,2]]
    result=detailed_jagwas_plan(choices,writer_prices(),preparations=[preparation(c) for c in choices],
        occupancy_scenarios={'none':'empty','all':'dense'},host_scenarios={'half':dict(host_serial_fraction=.5)},
        cpu_workers=3,host_memory_bytes=1<<40,device_memory_bytes={'cuda:0':1<<40,'cuda:1':1<<40})
    scores=[max(s['estimate']['estimated_seconds'] for s in row['scenarios']) for row in result['candidates']]
    assert len(scores)==6 and scores==sorted(scores)
    assert not result['selection_validated']


def test_singleton_fp64_projection_uses_scalar_capacity():
    row=next(row for row in CAPTURE['projection'] if row['B']==1)
    def predict(device):
        return jagwas_projection_service(row['N'],1,row['K'],device,row['kernels'],
            host_primitives=projection_prices(row))
    fast=predict(replace(resources(),fp64_tensor_flops_per_second=None))
    slow=predict(replace(resources(),fp64_flops_per_second=1e10))
    op=lambda r:next(s for s in r['operations'] if s['op'] in ('aten.mm.default','aten.mm.out'))
    assert fast['gemm']['arithmetic_kind']=='fp64_scalar'
    assert fast['gemm']['input_l2_bytes']==8*512*513
    assert fast['gemm']['arithmetic_accounting']=='logical_floor_not_verified_issued_work'
    assert op(slow)['math_seconds']==pytest.approx(10*op(fast)['math_seconds'])
    assert fast['kernel_count']==len(row['kernels'])
    assert not fast['prediction_complete']
    with pytest.raises(ValueError,match='FP64'):predict(replace(resources(),fp64_flops_per_second=None))


@pytest.mark.parametrize('fault',['block','split','type','specialization','width'])
def test_singleton_fp64_unknown_geometry_refused(fault):
    row=copy.deepcopy(next(row for row in CAPTURE['projection'] if row['B']==1))
    kernel=next(k for k in row['kernels'] if 'internal::gemvx::kernel<' in k['name'])
    if fault=='block':kernel['geometry']['block']=[32,16,1]
    if fault=='split':kernel['geometry']['grid'][2]=2
    if fault=='type':kernel['name']=kernel['name'].replace('double','float')
    if fault=='specialization':kernel['name']=kernel['name'].replace('false, 8, false','false, 9, false')
    with pytest.raises(ValueError):gemm_work(512,1,513 if fault=='width' else 512,row['kernels'],dtype='float64')

@pytest.mark.parametrize('fault',['samples','block','grid','reduction','signature'])
def test_singleton_statistics_nsp_unknown_geometry_refused(fault):
    row=copy.deepcopy(next(row for row in CAPTURE['statistics'] if row['B']==1))
    kernel=next(k for k in row['kernels'] if 'gemvNSP_kernel<' in k['name'])
    if fault=='block':kernel['geometry']['block']=[32,16,1]
    if fault=='grid':kernel['geometry']['grid'][2]=4
    if fault=='reduction':row['kernels']=[k for k in row['kernels'] if 'splitKreduce' not in k['name']]
    if fault=='signature':kernel['name']=kernel['name'].replace('1024','512')
    with pytest.raises(ValueError):gemm_work(2050 if fault=='samples' else 2049,1,515,row['kernels'])
