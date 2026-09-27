import copy,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'benchmarks'))
from direct_binary_ties import validate_profile


def test_profile_rejects_unsupported_geometry_and_incomplete_record_counts():
    c={'context':{'device_index':'0'},'torch':{str(k):{'shape':[8192,k,8,2048,4,4]} for k in range(1,33)}}
    s={'variants':65536,'record_forms':{'1':65536}}
    validate_profile(c,s)
    bad=copy.deepcopy(c);bad['torch']['2']['shape'][0]=10000
    with pytest.raises(ValueError,match='shape mismatch'):validate_profile(bad,s)
    with pytest.raises(ValueError,match='every variant'):validate_profile(c,{'variants':65536,'record_forms':{'1':65535}})
    with pytest.raises(ValueError,match='recalibrate'):validate_profile(c,{'variants':16384,'record_forms':{'1':16384}})
    bad=copy.deepcopy(c);bad['context']['device_index']='0,1'
    with pytest.raises(ValueError,match='one GPU'):validate_profile(bad,s)


def test_main_calculator_keeps_package_context_for_resource_profiles(tmp_path):
    import json,subprocess
    profile={
        'workload':{'variants':32,'samples':16,'traits':2,'covariates':1},
        'input':{'format':'pgen-hardcall','stored_bytes':128,'transfer_bytes_per_variant':16,'decoded_bytes_per_value':4},
        'hardware':{'disk_bytes_per_second':1e9,'h2d_bytes_per_second':1e10,'host_bytes_per_second':1e11,'gpu_flops_per_second':1e12,'gpu_bytes_per_second':1e11,'host_memory_bytes':2**30,'device_memory_bytes':2**30,'cpu_workers':1,'d2h_bytes_per_second':1e10},
        'plan':{'read_variants':8,'decode_variants':8,'chunk_variants':8,'workers':1,'depth':2}}
    path=tmp_path/'profile.json';path.write_text(json.dumps(profile))
    result=subprocess.run([sys.executable,str(Path(__file__).parents[1]/'runtime_calculator.py'),str(path)],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    output=json.loads(result.stdout)
    assert output['timing_device_count']==1 and output['memory_feasible']
    assert output['prediction_seconds'] is None  # no invented missing components
