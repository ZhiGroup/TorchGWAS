"""Actual duration-free FP64 geometry with independent synthetic resource controls."""
import copy
from dataclasses import replace
import json
from pathlib import Path
import pytest
from torchgwas.reduction_tensor_work import jagwas_tensor_work,jagwas_host_primitive_name,jagwas_projection_service
from torchgwas.jagwas_blocks import projection_flops_per_variant,triangular_blocks
from torchgwas.tensor_service import GEMM_OPS,DeviceService,gemm_work

CAPTURE=json.loads((Path(__file__).parent/'fixtures'/'jagwas_projection_geometry.json').read_text())


def resources():
    return DeviceService(hbm_bytes_per_second=1e30,l2_bytes_per_second=1e30,fp32_flops_per_second=1e12,
        kernel_launch_seconds=1e-6,host_dispatch_cpu_seconds=0.,available_l2_bytes=40<<20,sm_count=108,
        fp64_flops_per_second=1e11,fp64_tensor_flops_per_second=1e12)


def prices(row):
    work=jagwas_tensor_work(row['N'],row['B'],row['K'],phase='reduce',compute_dtype=row['compute_dtype'])
    return {jagwas_host_primitive_name(call):(i+1)*1e-7 for i,call in enumerate(work['host_calls'])}


def predict(row,**kwargs):
    opts=dict(compute_dtype=row['compute_dtype'],host_primitives=prices(row));opts.update(kwargs)
    return jagwas_projection_service(row['N'],row['B'],row['K'],opts.pop('resources',resources()),row['kernels'],**opts)


@pytest.mark.parametrize('row',CAPTURE['rows'],ids=lambda r:f"K{r['K']}-{r['compute_dtype']}")
def test_all_captured_projection_kernels_are_accounted_without_scan_rates(row):
    result=predict(row)
    work=jagwas_tensor_work(row['N'],row['B'],row['K'],phase='reduce',compute_dtype=row['compute_dtype'])
    bank=prices(row)
    assert result['kernel_count']==len(row['kernels'])
    # The block-triangular projection: one product per row block of L^-1.
    assert result['gemm']['useful_flops']==row['B']*projection_flops_per_variant(row['K'])
    assert result['gemm']['arithmetic_kind']=='fp64_tensor'
    assert result['gemm']['issued_flops']>=result['gemm']['useful_flops']
    assert result['host_dispatch_mode']=='source_grouped_typed'
    assert result['host_dispatch_cpu_seconds']==pytest.approx(sum(bank[jagwas_host_primitive_name(c)] for c in work['host_calls']))
    assert result['prediction_complete'] is False
    blocks=len(triangular_blocks(row['K']))
    assert sum(op['op'] in GEMM_OPS for op in result['operations'])==blocks
    if blocks>1:
        assert result['gemm']['products']==blocks and result['gemm']['split_k']==[1]*blocks
    else:assert result['gemm']['split_k']==1


def test_scalar_fp64_and_tensor_fp64_capacity_affect_separate_operations():
    row=next(r for r in CAPTURE['rows'] if r['K']==512 and r['compute_dtype']=='float32')
    fast=predict(row)
    slow_scalar=predict(row,resources=replace(resources(),fp64_flops_per_second=1e10))
    slow_tensor=predict(row,resources=replace(resources(),fp64_tensor_flops_per_second=1e11))
    def operation(result,name):return next(op for op in result['operations'] if op['op']==name)
    for name in ['aten.sum.dim_IntList']:
        assert operation(slow_scalar,name)['math_seconds']==pytest.approx(10*operation(fast,name)['math_seconds'])
        assert operation(slow_tensor,name)['math_seconds']==operation(fast,name)['math_seconds']
    # The finiteness guard and the score transform run on the FP32 t, before the FP64 cast.
    for name in ['aten.abs.default','aten.rsqrt_.default']:
        assert operation(slow_scalar,name)['math_seconds']==operation(fast,name)['math_seconds']
    assert operation(slow_tensor,'aten.mm.out')['math_seconds']==pytest.approx(10*operation(fast,'aten.mm.out')['math_seconds'])
    assert operation(slow_scalar,'aten.mm.out')['math_seconds']==operation(fast,'aten.mm.out')['math_seconds']


@pytest.mark.parametrize('field',['fp64_flops_per_second','fp64_tensor_flops_per_second'])
def test_missing_fp64_rate_is_not_substituted_with_fp32(field):
    row=CAPTURE['rows'][0]
    with pytest.raises(ValueError,match='FP64'):
        predict(row,resources=replace(resources(),**{field:None}))
    result=predict(row,resources=replace(resources(),**{field:0.}))
    assert result['estimated_span_seconds'] is None
    assert result['blocked_resources']==[field]


def test_missing_fixed_dispatch_price_and_unknown_compiled_family_are_refused():
    row=copy.deepcopy(CAPTURE['rows'][0])
    key='joint_isfinite_'+('fp32' if row['compute_dtype']=='float32' else 'fp64')
    bank=prices(row);bank.pop(key)
    with pytest.raises(ValueError,match=key):predict(row,host_primitives=bank)
    for kernel in row['kernels']:kernel['name']=kernel['name'].replace('d884gemm','h16816gemm')
    with pytest.raises(ValueError,match='GEMM launch census'):predict(row)


def test_fp64_operand_traffic_uses_eight_byte_values_and_swizzle_padding():
    main=dict(name='cutlass_80_tensorop_d884gemm_32x32_16x5_nn_align1',geometry=dict(grid=[32,2,1]))
    work=gemm_work(512,128,512,[main],dtype='float64')
    assert work['swizzle']==8 and work['grid_ctas']==64
    assert work['input_l2_bytes']==8*512*(128*16+512*4)
    with pytest.raises(ValueError,match='swizzled'):
        gemm_work(512,128,512,[dict(main,geometry=dict(grid=[31,2,1]))],dtype='float64')
