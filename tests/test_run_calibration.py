"""Public observation bindings keep evidence immutable and work-specific."""
import inspect
import json
import time
from pathlib import Path
from unittest.mock import patch
import pytest

from torchgwas.api import run_linear_gwas
from torchgwas.run_calibration import RunCalibration,calibration_options
from test_initial_calibration import feed


@pytest.fixture
def binding(tmp_path,monkeypatch):
    from torchgwas import detailed_calibration,adaptive_chunks
    context={'devices':['cuda:0','cuda:1'],'libraries':'test'}
    sources={'test.py':'original'}
    monkeypatch.setattr(detailed_calibration,'execution_context',lambda *a,**k:context)
    monkeypatch.setattr(detailed_calibration,'source_identity',lambda **kwargs:sources.copy())
    monkeypatch.setattr(adaptive_chunks,'validate_chunk_control',lambda *a:None)
    path=tmp_path/'input.pgen';path.write_bytes(b'fixture')
    config=dict(cache_dir=str(tmp_path/'cache'),validation_chunks_per_device=2,
        max_chunks_per_device=3,warmup_chunks=0,stride=1,max_window_seconds=10.)
    def make(**kwargs):
        obj=RunCalibration(config,inputs=kwargs.get('inputs',{'genotype':path}),
                           output_path=kwargs.get('output',tmp_path/'output'))
        obj.prepare(object(),devices=['cuda:0','cuda:1'],request=kwargs.get('request',{'capacity':10}))
        return obj
    return make,config,path,context,sources


def window(run,device='cuda:0',*,trait_range=(0,4),variant_range=(0,100)):
    return run.observer(device,variant_range=variant_range,trait_range=trait_range,
        reader_workers=2,capacity=10,depth=2)


def complete(run,device='cuda:0'):
    control=window(run,device)
    for i in range(3):feed(control,device,10*i)
    return run.finish(successful=True)


def publication(report,device='cuda:0'):
    return report['windows'][device]['publication']['publications'][device]


def test_public_reuse_preserves_original_evidence_and_age(binding):
    make,config,path,context,sources=binding
    first=complete(make());record=publication(first)
    old=Path(record['path']).read_bytes()
    original=json.loads(old)['observed_unix_seconds']
    second=complete(make(output=path.parent/'another-output'))
    report=second['windows']['cuda:0']
    decision=report['refresh_decisions']['cuda:0']
    assert decision['cache_hit'] and decision['state']=='cached_consistent'
    assert decision['previous_record_sha256']==record['record_sha256']
    assert report['reserved']=={'cuda:0':2}
    assert Path(record['path']).read_bytes()==old
    saved=json.loads(Path(publication(second)['path']).read_text())
    assert saved['name'].endswith(':validation')
    third=window(make()).snapshot()['refresh_decisions']['cuda:0']
    assert third['previous_record_sha256']==record['record_sha256']
    assert json.loads(Path(record['path']).read_text())['observed_unix_seconds']==original


@pytest.mark.parametrize('changed',['input','source','context','output-mode','tile','memory'])
def test_incompatible_evidence_is_not_reused(binding,changed):
    make,config,path,context,sources=binding
    complete(make())
    kwargs={}
    if changed=='input':path.write_bytes(b'changed fixture')
    if changed=='source':sources['test.py']='changed'
    if changed=='context':context['libraries']='changed'
    if changed=='output-mode':kwargs['request']={'capacity':10,'reduction':'significant'}
    if changed=='memory':kwargs['inputs']={'genotype':path,'phenotype':object()}
    control=window(make(**kwargs),trait_range=(4,8) if changed=='tile' else (0,4))
    assert not control.snapshot()['refresh_decisions']['cuda:0']['cache_hit']


def test_expiry_is_measured_from_original_observation(binding):
    make,config,*_=binding
    first=complete(make())
    saved=json.loads(Path(publication(first)['path']).read_text())
    with patch('torchgwas.calibration_cache.time.time',return_value=saved['observed_unix_seconds']+301.):
        state=window(make()).snapshot()['refresh_decisions']['cuda:0']
    assert not state['cache_hit'] and state['cache_reason']=='expired'


def test_only_first_tile_per_gpu_and_one_shared_wall_window(binding):
    make,*_=binding;run=make();a=window(run)
    assert window(run) is a
    assert window(run,trait_range=(4,8)) is None
    b=window(run,'cuda:1')
    with patch('torchgwas.run_calibration.time.perf_counter',return_value=10.):
        assert a.reserve_read(0,10,'cuda:0')
    with patch('torchgwas.run_calibration.time.perf_counter',return_value=19.):
        assert b.reserve_read(0,10,'cuda:1')
    with patch('torchgwas.run_calibration.time.perf_counter',return_value=20.):
        assert not a.reserve_read(10,20,'cuda:0')
        assert not b.reserve_read(10,20,'cuda:1')
    report=run.finish(successful=False)
    assert not report['unobserved_devices']
    assert not publication_if_any(report)


def publication_if_any(report):
    return [value for window in report['windows'].values()
            for value in window['publication']['publications'].values()]


def test_output_progress_is_bounded_audit_and_keeps_distinct_durability(binding):
    from torchgwas.sumstats import DenseWriteProgress
    from torchgwas.sumstats_indexed import IndexedChunkWrite
    make,*_=binding;run=make();base=run._created
    # Concurrent stores may deliver their timestamped notifications out of order.
    run.output_written(DenseWriteProgress(0,2,(0,4),64,0,base+2.,'/test/dense','cuda:0'))
    run.output_written(IndexedChunkWrite(0,2,'significant',0,0,None,base+.5,base+1.,False))
    run.output_written(IndexedChunkWrite(2,4,'significant',2,200,'part.npz',base+2.,base+3.,True))
    report=run.finish(successful=False)['output_progress']
    assert report['first_written_seconds']==1. and report['first_material_seconds']==2.
    assert report['first_fsynced_part_seconds']==3. and report['last_written_seconds']==3.
    assert report['dense_statistic_cells']==8 and report['dense_statistic_bytes']==64
    assert report['indexed_rows']==2 and report['indexed_part_bytes']==200
    assert report['events']==3 and report['dense_ranges']==1 and report['indexed_chunks']==2
    assert not list(run.cache.directory.glob('**/*.json'))


def output_request(**changes):
    return dict(capacity=10,variant_range=[10,110],phenotype_shape=[20,8],
        reduction='significant',significance_threshold=.01,**changes)


def output_event(*,device='cuda:0',start=10,traits=(0,4),rows=0):
    from torchgwas.sumstats_indexed import IndexedChunkWrite,IndexedOutputPartition
    now=time.perf_counter()
    return IndexedChunkWrite(start-10,start-8,'significant',rows,100 if rows else 0,
        'part.npz' if rows else None,now,now,bool(rows),
        IndexedOutputPartition(device,(10,110),traits),(start,start+2))


def output_sample(run):
    run.output_written(output_event())
    run.output_written(output_event(device='cuda:1',traits=(4,8),rows=2))
    return run.finish(successful=True)['output_occupancy']


def test_completed_output_counts_reuse_immutable_records_without_renewing_age(binding):
    make,*_=binding;first=output_sample(make(request=output_request()))
    assert first['status']=='published' and len(first['bins'])==2
    saved=Path(first['publication']['path']);contents=saved.read_bytes();record=json.loads(contents)
    assert record['kind']=='stage_observations' and record['observed_unix_seconds']<=record['created_unix_seconds']
    second=output_sample(make(request=output_request()))
    assert second['cache_hit'] and second['status']=='reused_without_renewal' and second['publication'] is None
    assert second['previous_record_sha256']==first['publication']['record_sha256']
    assert saved.read_bytes()==contents
    assert second['bins'][0]['retained']==0 and second['bins'][1]['trait_range']==[4,8]
    assert second['bins'][0]['variant_range']==[10,12]


def test_changed_survivors_publish_a_new_measurement_without_overwriting_old(binding):
    make,*_=binding;first=output_sample(make(request=output_request()))
    saved=Path(first['publication']['path']);old=saved.read_bytes()
    run=make(request=output_request());run.output_written(output_event(rows=1))
    report=run.finish(successful=True)['output_occupancy']
    assert report['cache_hit'] and report['status']=='published'
    assert report['publication']['record_sha256']!=first['publication']['record_sha256']
    assert saved.read_bytes()==old


def test_output_occupancy_expiry_and_changed_threshold_require_fresh_records(binding):
    make,*_=binding;first=output_sample(make(request=output_request()))
    old=json.loads(Path(first['publication']['path']).read_text())
    run=make(request=output_request())
    with patch('torchgwas.calibration_cache.time.time',return_value=old['observed_unix_seconds']+301.):
        run.output_written(output_event())
    expired=run.finish(successful=False)['output_occupancy']
    assert not expired['cache_hit'] and expired['cache_reason']=='expired'
    request=output_request();request['significance_threshold']=.02
    changed=output_sample(make(request=request))
    assert not changed['cache_hit'] and changed['status']=='published'


def test_output_count_collection_is_bounded_and_bad_overlap_cannot_be_published(binding):
    make,*_=binding;run=make(request=output_request())
    for start in range(10,100,2):run.output_written(output_event(start=start))
    report=run.finish(successful=True)['output_occupancy']
    assert len(report['bins'])==3 and report['max_bins_per_device']==3
    run=make(request=output_request());run.output_written(output_event());run.output_written(output_event())
    report=run.finish(successful=True)['output_occupancy']
    assert report['sample_status']=='invalid' and report['publication'] is None and report['bins']==[]


@pytest.mark.parametrize('fault',['failed','source_changed','write_error'])
def test_output_count_cache_errors_do_not_invalidate_finished_scientific_work(binding,fault):
    make,config,path,context,sources=binding
    run=make(request=output_request());run.output_written(output_event())
    if fault=='source_changed':sources['test.py']='changed'
    if fault=='write_error':
        with patch.object(run.cache,'store',side_effect=OSError('cache offline')):
            report=run.finish(successful=True)
        assert report['successful'] and report['output_occupancy']['status']=='cache_error'
    else:
        report=run.finish(successful=fault!='failed')
        assert report['output_occupancy']['publication'] is None
    assert not list(Path(config['cache_dir']).glob('**/*.json'))


@pytest.mark.parametrize('fault',['scan-failed','input-changed','source-changed'])
def test_failed_or_changed_run_does_not_publish(binding,fault):
    make,config,path,context,sources=binding
    run=make();control=window(run)
    for i in range(3):feed(control,'cuda:0',10*i)
    if fault=='input-changed':path.write_bytes(b'changed')
    if fault=='source-changed':sources['test.py']='changed'
    report=run.finish(successful=fault!='scan-failed')
    assert not publication_if_any(report)
    assert not list(Path(config['cache_dir']).glob('**/*.json'))


def test_optional_cache_failure_keeps_measurements_and_results(binding):
    make,*_=binding;run=make()
    with patch.object(run.cache,'lookup',side_effect=OSError('offline')):
        assert window(run) is None
    report=run.finish(successful=True)
    assert report['successful'] and report['cache_errors']['cuda:0']['stage']=='cache_lookup'


def test_strict_binding_option_and_optional_digest_publication_failure(binding):
    make,config,*_=binding;config['reuse_binding_digests']=False
    with patch('torchgwas.binding_digests.BindingDigestCache',side_effect=AssertionError('strict binding requested')):
        run=make();assert run.finish(successful=True)['binding_digests'] is None
    config['reuse_binding_digests']=True;run=make()
    with patch.object(run._binding_digests,'publish',side_effect=OSError('cache unavailable')):
        report=run.finish(successful=True)
    assert report['successful'] and report['binding_digests']['state']['closed']
    assert report['binding_digests']['publication']['status']=='cache_error'


def test_changed_binding_file_prevents_empirical_publication(binding):
    make,config,path,*_=binding;run=make();control=window(run)
    library=path.parent/'library.so';library.write_text('initial')
    run._binding_digests.digests([library]);library.write_text('updated')
    for i in range(3):feed(control,'cuda:0',10*i)
    report=run.finish(successful=True)
    assert report['successful'] and not report['publication_inputs_unchanged']
    assert not publication_if_any(report) and not list(Path(config['cache_dir']).glob('**/*.json'))


@pytest.mark.parametrize('config',[{},None,{'cache_dir':''},{'cache_dir':'x','unknown':1},
    {'cache_dir':'x','max_age_seconds':-1},{'cache_dir':'x','max_chunks_per_device':1},
    {'cache_dir':'x','max_window_seconds':0},{'cache_dir':'x','max_decision_cpu_seconds':0},
    {'cache_dir':'x','reuse_binding_digests':1}])
def test_invalid_window_settings_are_rejected_without_io(config):
    with pytest.raises(ValueError):calibration_options(config)


@pytest.mark.parametrize('field,value',[('chunk_size',None),('chunk_size',True),('genotype','x.bed'),
    ('pgen_mode','dosage'),('compute_dtype','float64'),('output_dir',None),
    ('sumstats_format','none'),('phenotype_table','x.tsv'),('pipeline_profile',{})])
def test_public_preflight_fails_before_loading_or_scanning(tmp_path,field,value):
    args={name:p.default for name,p in inspect.signature(run_linear_gwas).parameters.items()}
    args.update(genotype='x.pgen',phenotype='y.npy',pgen_mode='hardcall',chunk_size=128,
        output_dir=tmp_path/'output',initial_calibration={'cache_dir':str(tmp_path/'cache')})
    args[field]=value
    with patch('torchgwas.api.load_genotype',side_effect=AssertionError('must fail before load')):
        with pytest.raises(ValueError):run_linear_gwas(**args)
    assert not (tmp_path/'output').exists() and not (tmp_path/'cache').exists()
