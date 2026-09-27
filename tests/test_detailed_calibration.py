import copy
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from torchgwas.detailed_calibration import (SCHEMA,bind_detailed_profile,read_detailed_profile,
    sha256_file,source_identity,storage_identity,validate_detailed_profile,write_detailed_profile)


def fixture(tmp_path):
    artifact=tmp_path/'primitive.json';artifact.write_text('{"unit":"independent"}')
    contexts=[dict(name='one',devices=['cuda:1'],profiles={'cuda:1':{'primitive':1.}},shared_capacities={'cpu':4.})]
    execution=dict(devices={'cuda:1':dict(uuid='gpu-one',driver_version='d')},affinity=[1,2],torch_threads=4,
        environment={'OMP_NUM_THREADS':'4'},storage={'input':{'source':'disk-one'},'output':{'source':'disk-one'}},native_library_sha256='native')
    sources={'source.py':'source-hash'}
    profile=bind_detailed_profile(contexts,execution,component_artifacts={str(artifact):sha256_file(artifact)},
        limitations=['One independently measured context; no prediction guarantee.'],sources=sources)
    return profile,execution,sources,artifact


def test_binding_is_independent_of_workload_and_does_not_certify_accuracy(tmp_path):
    profile,execution,sources,_=fixture(tmp_path)
    assert profile['schema']==SCHEMA and profile['bound_at_utc']
    assert profile['selection_validated'] is profile['runtime_prediction_validated'] is False
    assert not {'workload','bounds','observations','selected'}&set(profile)
    assert validate_detailed_profile(profile,execution,sources=sources)['context_matches']
    execution['devices']['cuda:1']['uuid']='different'
    assert profile['execution_context']['devices']['cuda:1']['uuid']=='gpu-one'


def test_optional_shared_transfer_context_is_positive_and_immutable(tmp_path):
    profile,execution,sources,_=fixture(tmp_path)
    contexts=copy.deepcopy(profile['contexts'])
    contexts[0]['shared_transfer_capacities']=dict(h2d=1e10,d2h=8e9)
    bound=bind_detailed_profile(contexts,execution,
        component_artifacts=profile['component_artifacts'],limitations=[],sources=sources)
    assert validate_detailed_profile(bound,execution,sources=sources)['context_matches']
    assert bound['contexts'][0]['shared_transfer_capacities']==dict(h2d=1e10,d2h=8e9)
    for invalid in (dict(h2d=1e10),dict(h2d=0.,d2h=8e9),
                    dict(h2d=float('nan'),d2h=8e9)):
        contexts[0]['shared_transfer_capacities']=invalid
        with pytest.raises(ValueError,match='shared H2D/D2H'):
            bind_detailed_profile(contexts,execution,
                component_artifacts=profile['component_artifacts'],
                limitations=[],sources=sources)


@pytest.mark.parametrize('change,field',[
    (lambda e:e.update(affinity=[2,3]),'affinity'),
    (lambda e:e.update(torch_threads=8),'torch_threads'),
    (lambda e:e.update(torch_default_dtype='torch.float64'),'torch_default_dtype'),
    (lambda e:e['devices']['cuda:1'].update(uuid='same-name-different-card'),'uuid'),
    (lambda e:e['devices']['cuda:1'].update(driver_version='new'),'driver_version'),
    (lambda e:e['environment'].update(OMP_NUM_THREADS='8'),'OMP_NUM_THREADS'),
    (lambda e:e['environment'].update(CUDA_LAUNCH_BLOCKING='1'),'CUDA_LAUNCH_BLOCKING'),
    (lambda e:e['storage']['output'].update(source='other-disk'),'storage.output.source'),
    (lambda e:e.update(native_library_sha256='other-build'),'native_library_sha256'),
    (lambda e:e.update(unknown_setting=True),'unknown_setting')])
def test_runtime_changes_are_refused_with_exact_fields(tmp_path,change,field):
    profile,execution,sources,_=fixture(tmp_path);change(execution)
    with pytest.raises(ValueError,match=field):validate_detailed_profile(profile,execution,sources=sources)


def test_source_and_component_changes_fail_before_execution(tmp_path):
    profile,execution,sources,artifact=fixture(tmp_path)
    with pytest.raises(ValueError,match='source.source.py'):
        validate_detailed_profile(profile,execution,sources={'source.py':'changed'})
    artifact.write_text('changed')
    with pytest.raises(ValueError,match='Component artifact changed'):
        validate_detailed_profile(profile,execution,sources=sources)


def test_incomplete_device_contexts_and_unknown_schema_fail(tmp_path):
    profile,execution,sources,_=fixture(tmp_path)
    for change in [lambda p:p['contexts'][0]['profiles'].clear(),lambda p:p['execution_context']['devices'].clear(),
                   lambda p:p.update(schema='coarse'),lambda p:p['contexts'].append(copy.deepcopy(p['contexts'][0]))]:
        bad=copy.deepcopy(profile);change(bad)
        with pytest.raises(ValueError):validate_detailed_profile(bad,execution,sources=sources)


def test_atomic_profile_publication_does_not_replace_existing_file(tmp_path):
    profile,execution,sources,_=fixture(tmp_path);path=tmp_path/'profiles'/'calibration.json'
    assert write_detailed_profile(profile,path)==path
    assert read_detailed_profile(path)==profile
    before=path.read_bytes()
    with pytest.raises(FileExistsError):write_detailed_profile(dict(profile,scope='changed'),path)
    assert path.read_bytes()==before
    assert [p.name for p in path.parent.iterdir()]==['calibration.json']
    validate_detailed_profile(read_detailed_profile(path),execution,sources=sources)


def test_nonfinite_json_is_refused_and_no_partial_file_published(tmp_path):
    profile,_,_,_=fixture(tmp_path);path=tmp_path/'bad.json';profile['contexts'][0]['profiles']['cuda:1']['primitive']=float('nan')
    with pytest.raises(ValueError):write_detailed_profile(profile,path)
    assert not path.exists()
    path.write_text('{"schema":"'+SCHEMA+'","bad":NaN}')
    with pytest.raises(ValueError,match='Non-finite'):read_detailed_profile(path)
    path.write_text('{"schema":"'+SCHEMA+'","bad":1e999}')
    with pytest.raises(ValueError,match='Non-finite'):read_detailed_profile(path)
    path.write_text('{"schema":"'+SCHEMA+'","schema":"'+SCHEMA+'"}')
    with pytest.raises(ValueError,match='Duplicate profile key'):read_detailed_profile(path)


def test_binding_cannot_smuggle_an_old_plan_or_claim_validation(tmp_path):
    profile,execution,sources,_=fixture(tmp_path)
    for change in [lambda p:p.update(selected={'chunk_size':512}),lambda p:p.update(selection_validated=True),
                   lambda p:p['component_artifacts'].clear()]:
        bad=copy.deepcopy(profile);change(bad)
        with pytest.raises(ValueError):validate_detailed_profile(bad,execution,sources=sources)


def test_binding_refuses_relabeling_prices_for_different_hardware(tmp_path):
    profile,execution,sources,_=fixture(tmp_path)
    contexts=profile['contexts'];contexts[0]['profiles']['cuda:1']['gpu_resources']={'sm_count':132}
    execution['devices']['cuda:1']['sm_count']=108
    with pytest.raises(ValueError,match='Priced GPU properties'):
        bind_detailed_profile(contexts,execution,component_artifacts=profile['component_artifacts'],limitations=[],sources=sources)


def test_storage_check_uses_existing_ancestor_without_creating_output(tmp_path):
    path=tmp_path/'new'/'output';mount=dict(source='/dev/test',fstype='ext4',target=str(tmp_path))
    with patch('torchgwas.detailed_calibration.subprocess.check_output',return_value=json.dumps({'filesystems':[mount]})) as run:
        identity=storage_identity(path)
    assert identity['source']=='/dev/test' and not path.exists()
    assert run.call_args.args[0][3]==str(tmp_path.resolve())


def test_source_fingerprint_includes_actual_model_and_binder():
    identity=source_identity()
    assert all(len(identity[name])==64 for name in ['detailed_calibration.py','trait_candidate_space.py','trait_tiling_model.py','setup_work.py'])
