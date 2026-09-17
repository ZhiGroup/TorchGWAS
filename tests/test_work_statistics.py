import pytest
from torchgwas.work_statistics import native_pgen_work,form_service_cpu_seconds,CpuOperatingConditions

def test_ld_boundary_replay_counts_actual_native_work():
    work=native_pgen_work([0,2,2,0,3],[10,3,4,10,2],40,[(0,2),(2,4),(4,5)],[0,2,3,0,1])
    assert work['requested_records']==5
    assert work['decoded_records']==8
    assert work['replay_records']==3
    assert work['record_read_bytes']==52
    assert work['replay_read_bytes']==23
    assert work['packed_materialization_bytes']==80
    assert work['difference_patches']==8
    assert form_service_cpu_seconds(work,{0:2,2:3,3:4})==21

def test_single_range_avoids_repeated_base_prefix():
    work=native_pgen_work([0,2,2,0,3],[10,3,4,10,2],40,[(0,5)])
    assert work['record_read_bytes']==29 and work['replay_records']==0
    assert work['difference_patches'] is None
    with pytest.raises(ValueError,match='missing rate'):form_service_cpu_seconds(work,{0:1})

def test_cpu_capacity_is_explicit_operating_condition():
    assert CpuOperatingConditions(2.5,'independent CPU availability probe').capacity(4)==2.5
    assert CpuOperatingConditions(8,'scenario').capacity(4)==4
    with pytest.raises(ValueError):CpuOperatingConditions(0,'scenario').capacity(4)

def test_reject_unsupported_or_invalid_pgen_statistics():
    with pytest.raises(ValueError):native_pgen_work([2],[2],40,[(0,1)])
    with pytest.raises(ValueError):native_pgen_work([16],[10],40,[(0,1)])


def test_pipeline_explicit_work_and_capacity_replace_average():
    from dataclasses import replace
    from torchgwas.pipeline_model import Workload,InputProfile,Hardware,PipelinePlan,estimate
    w=Workload(8,40,1,0);plan=PipelinePlan(2,2,2,4,4)
    p=InputProfile('pgen-hardcall',100,40,1,cpu_decode_core_seconds_per_variant=99,
                  cpu_decode_core_seconds_total=12,genotype_record_read_bytes_total=160,
                  decode_work_shape=(8,40,2,2,2),decode_work_source='record index and independent form rates')
    h=Hardware(100,1e6,1e6,1e6,1e6,1e8,1e8,4,1e6)
    result=estimate(w,p,h,plan,operating_conditions=CpuOperatingConditions(2,'scenario'))
    assert result['resource_seconds']['cpu_decode']==6
    assert result['resource_seconds']['storage']==1.6
    with pytest.raises(ValueError,match='work statistics'):
        estimate(w,p,h,replace(plan,decode_variants=4))


def test_numpy_index_arrays_are_supported():
    import numpy as np
    work=native_pgen_work(np.array([0,2,2],dtype=np.uint8),np.array([10,3,4],dtype=np.uint32),40,[(2,3)])
    assert work['replay_records']==2 and work['record_read_bytes']==17
    import json
    json.dumps(work)

def test_wall_service_cannot_receive_a_second_load_correction():
    from torchgwas.pipeline_model import Workload,InputProfile,Hardware,PipelinePlan,ExecutionProfile,estimate
    w=Workload(8,40,1,0);plan=PipelinePlan(2,2,2,4,4)
    p=InputProfile('pgen-hardcall',100,40,1);h=Hardware(100,1e6,1e6,1e6,1e6,1e8,1e8,4,1e6)
    execution=ExecutionProfile(1,0,.1,(40,1,0,2,4,4),'probe',{})
    with pytest.raises(ValueError,match='already includes load'):
        estimate(w,p,h,plan,execution=execution,operating_conditions=CpuOperatingConditions(2,'scenario'))
