"""Digest reuse detects ordinary file edits without caching runtime context."""
import json
import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pytest

from torchgwas import binding_digests as module
from torchgwas.binding_digests import BindingDigestCache
from torchgwas import detailed_calibration as calibration


@pytest.fixture
def files(tmp_path,monkeypatch):
    paths=[tmp_path/'a.so',tmp_path/'b.py']
    for path,text in zip(paths,['library','source']):path.write_text(text)
    monkeypatch.setattr(module,'_host_identity',lambda:'host-and-boot')
    # File-change checks use real stat metadata; only the one-second reuse
    # eligibility clock is advanced so unit tests do not sleep.
    monkeypatch.setattr(module.time,'time_ns',lambda:max(p.stat().st_ctime_ns for p in paths if p.exists())+2_000_000_000)
    return paths,tmp_path/'cache'


def artifacts(directory):return {str(p):p.read_bytes() for p in directory.glob('**/*.json')}


def test_cross_job_digest_reuse_preserves_records_and_avoids_byte_reads(files):
    paths,directory=files;cache=BindingDigestCache(directory)
    expected=[calibration.sha256_file(p) for p in paths]
    assert cache.digests(paths)==expected and not directory.exists()
    assert cache.snapshot()['hashed_files']==2
    saved=cache.publish(successful=True);before=artifacts(directory)
    assert len(saved['stored'])==1;cache.close()
    later=BindingDigestCache(directory)
    with patch.object(calibration,'sha256_file',side_effect=AssertionError('cached bytes reread')):
        assert later.digests(list(reversed(paths)))==list(reversed(expected))
        assert later.digests(paths)==expected
    assert later.snapshot()['disk_hits']==1 and later.snapshot()['memory_hits']==1
    assert later.snapshot()['hashed_files']==0
    assert later.publish(successful=True)['stored']==[] and artifacts(directory)==before
    record=json.loads(next(iter(before.values())))
    assert record['kind']=='source_work' and record['max_age_seconds'] is None
    assert 'observed_unix_seconds' not in record


@pytest.mark.parametrize('mutation',['same_size','restored_mtime','replacement','symlink','boot'])
def test_changed_file_or_host_boot_cannot_reuse_the_previous_digest(files,monkeypatch,mutation):
    paths,directory=files
    if mutation=='symlink':
        link=paths[0].parent/'link.so';link.symlink_to(paths[0]);selected=[link]
    else:selected=paths
    cache=BindingDigestCache(directory);old=cache.digests(selected);cache.publish(successful=True)
    before=artifacts(directory);metadata=paths[0].stat()
    if mutation=='boot':monkeypatch.setattr(module,'_host_identity',lambda:'another-boot')
    elif mutation=='symlink':selected[0].unlink();selected[0].symlink_to(paths[1])
    elif mutation=='replacement':
        replacement=paths[0].with_suffix('.new');replacement.write_text('changed');replacement.replace(paths[0])
    else:
        paths[0].write_text('changed')
        if mutation=='restored_mtime':os.utime(paths[0],ns=(metadata.st_atime_ns,metadata.st_mtime_ns))
    later=BindingDigestCache(directory);actual=later.digests(selected)
    assert actual==[calibration.sha256_file(p) for p in selected]
    assert later.snapshot()['disk_hits']==0 and later.snapshot()['hashed_files']==len(selected)
    if mutation!='boot':assert actual!=old
    assert all(artifacts(directory)[name]==body for name,body in before.items())


def test_recent_or_future_metadata_bypasses_reuse_and_publication(files,monkeypatch):
    paths,directory=files;monkeypatch.setattr(module,'_stable',lambda row:False)
    cache=BindingDigestCache(directory);cache.digests(paths);cache.digests(paths)
    assert cache.snapshot()['hashed_files']==4 and cache.snapshot()['bypasses']==2
    assert cache.publish(successful=True)['stored']==[] and not directory.exists()


def test_recent_same_tick_change_is_caught_by_content_recheck(files,monkeypatch):
    paths,directory=files;path=paths[0];identity=module.input_identity(path)
    monkeypatch.setattr(module,'_stable',lambda row:False)
    original=module.input_identity
    # Model coarse filesystem timestamps; bytes change but stat is identical.
    monkeypatch.setattr(module,'input_identity',lambda p:identity if Path(p)==path else original(p))
    cache=BindingDigestCache(directory);cache.digests([path]);path.write_text('changed')
    assert not cache.unchanged() and cache.publish(successful=True)['status']=='files_changed'
    assert not directory.exists()


@pytest.mark.parametrize('during',[True,False])
def test_symlink_retargeting_with_original_file_unchanged_is_detected(files,during):
    paths,directory=files;link=paths[0].parent/'alias';link.symlink_to(paths[0])
    cache=BindingDigestCache(directory);original=calibration.sha256_file
    def changing(path):
        value=original(path);link.unlink();link.symlink_to(paths[1]);return value
    if during:
        with patch.object(calibration,'sha256_file',side_effect=changing):
            with pytest.raises(ValueError,match='changed during'):cache.digests([link])
    else:
        cache.digests([link]);link.unlink();link.symlink_to(paths[1])
    assert not cache.unchanged() and cache.publish(successful=True)['status']=='files_changed'
    assert not directory.exists()


def test_changes_during_hashing_or_after_binding_refuse_publication(files):
    paths,directory=files;original=calibration.sha256_file
    def changing(path):
        value=original(path);Path(path).write_text('changed');return value
    cache=BindingDigestCache(directory)
    with patch.object(calibration,'sha256_file',side_effect=changing):
        with pytest.raises(ValueError,match='changed during'):cache.digests(paths)
    assert not cache.unchanged() and cache.publish(successful=True)['status']=='files_changed'
    later=BindingDigestCache(directory);later.digests(paths);paths[0].write_text('different again')
    assert not later.unchanged() and later.publish(successful=True)['stored']==[]
    assert not directory.exists()


def test_cache_corruption_failure_and_unsuccessful_jobs_fall_back_or_do_not_publish(files):
    paths,directory=files;expected=[calibration.sha256_file(p) for p in paths]
    failed=BindingDigestCache(directory);failed.digests(paths)
    assert failed.publish(successful=False)['stored']==[] and not directory.exists()
    cache=BindingDigestCache(directory)
    with patch.object(cache.cache,'lookup',side_effect=OSError('unavailable')):assert cache.digests(paths)==expected
    with patch.object(cache.cache,'store',side_effect=OSError('unavailable')):assert cache.publish(successful=True)['stored']==[]
    assert len(cache.snapshot()['errors'])==2
    saved=BindingDigestCache(directory);saved.digests(paths);record=saved.publish(successful=True)['stored'][0]
    Path(record['path']).write_text('invalid')
    later=BindingDigestCache(directory);assert later.digests(paths)==expected
    assert later.snapshot()['disk_hits']==0


def test_retention_is_bounded_and_closed_caches_cannot_be_reused(files):
    paths,directory=files;cache=BindingDigestCache(directory,max_groups=1,max_files=2)
    cache.digests([paths[0]]);cache.digests([paths[1]])
    assert cache.snapshot()['groups']==1 and cache.snapshot()['bypasses']==1
    extra=directory.parent/'extra';extra.write_text('third')
    with pytest.raises(ValueError,match='count exceeds'):cache.digests([extra])
    cache.close()
    assert not cache.unchanged()
    with pytest.raises(ValueError,match='closed'):cache.digests(paths)
    with pytest.raises(ValueError,match='closed'):cache.publish(successful=True)


def test_source_identity_default_still_hashes_and_cached_results_are_equal(tmp_path):
    cache=BindingDigestCache(tmp_path/'cache')
    strict=calibration.source_identity()
    assert calibration.source_identity(digest_cache=cache)==strict
    with patch.object(calibration,'sha256_file',wraps=calibration.sha256_file) as hasher:
        assert calibration.source_identity()==strict
        assert hasher.call_count==len(strict)


def test_runtime_settings_devices_and_mounts_remain_fresh_while_digests_hit(files,monkeypatch):
    import numpy as np
    import torch
    import threadpoolctl
    from torchgwas import gpu_identity,pgen_native
    paths,directory=files;cache=BindingDigestCache(directory)
    settings=dict(threads=2,driver='one',storage='/dev/first')
    uuid='12345678-1234-1234-1234-123456789abc';queries=[]
    monkeypatch.setattr(torch.cuda,'get_device_properties',lambda device:SimpleNamespace(
        uuid=uuid,name='fixture',total_memory=1<<30,multi_processor_count=1,
        max_threads_per_multi_processor=32,major=8,minor=0))
    def physical(uuids):
        queries.append(uuids)
        return {uuid:dict(pci_bus_id='00000000:01:00.0',driver_version=settings['driver'])}
    monkeypatch.setattr(gpu_identity,'physical_gpu_identity',physical)
    monkeypatch.setattr(pgen_native,'load_library',lambda:SimpleNamespace(_name=str(paths[0])))
    monkeypatch.setattr(threadpoolctl,'threadpool_info',lambda:[dict(internal_api='blas',prefix='test',
        version='1',num_threads=settings['threads'],filepath=str(paths[0]))])
    monkeypatch.setattr(torch,'get_num_threads',lambda:settings['threads'])
    monkeypatch.setattr(calibration,'storage_identity',lambda path:dict(source=settings['storage']))
    monkeypatch.setattr(calibration.os,'sched_getaffinity',lambda pid:{0})
    core=np._core._multiarray_umath
    features=dict(core.__cpu_features__)
    monkeypatch.setattr(core,'__cpu_features__',features)
    monkeypatch.delenv('NPY_DISABLE_CPU_FEATURES',raising=False)
    before=calibration.execution_context(['cuda:0'],input_path=paths[0],output_path=directory,digest_cache=cache)
    settings.update(threads=4,driver='two',storage='/dev/second')
    feature=next(iter(features));features[feature]=not features[feature]
    monkeypatch.setenv('NPY_DISABLE_CPU_FEATURES','test-declaration')
    with patch.object(calibration,'sha256_file',side_effect=AssertionError('digest bytes reread')):
        after=calibration.execution_context(['cuda:0'],input_path=paths[0],output_path=directory,digest_cache=cache)
    assert len(queries)==2 and cache.snapshot()['memory_hits']==4
    assert (before['torch_threads'],after['torch_threads'])==(2,4)
    assert after['cpu_pools'][0]['num_threads']==4
    assert after['devices']['cuda:0']['driver_version']=='two'
    assert after['storage']['input']['source']=='/dev/second'
    assert before['native_library_sha256']==after['native_library_sha256']
    assert before['numpy_core']['library_sha256']==after['numpy_core']['library_sha256']
    assert before['numpy_core']['cpu_features'][feature]!=after['numpy_core']['cpu_features'][feature]
    assert 'NPY_DISABLE_CPU_FEATURES' not in before['environment']
    assert after['environment']['NPY_DISABLE_CPU_FEATURES']=='test-declaration'
