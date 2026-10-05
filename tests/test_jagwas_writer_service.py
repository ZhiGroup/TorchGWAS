"""Source archive byte checks and resource conservation for joint output."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from torchgwas.reduced_output_work import jagwas_writer_work,jagwas_host_selection_service,jagwas_archive_service
from test_jagwas_schedule import shard,run

path=Path(__file__).parents[1]/'benchmarks'/'direct_jagwas_writer_primitives_20260921.py'


def writer_primitives():
    # A benchmark script outside the repository: only the test that needs it
    # skips, so the modules importing bank() and archive() still run.
    if not path.exists():
        pytest.skip(f'benchmarks/{path.name} is not in this repository')
    spec=importlib.util.spec_from_file_location('joint_writer_primitives',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def bank():
    return {name:dict(call_cpu_seconds=1e-6,unit_cpu_seconds=1e-9) for name in
        ['fp32_to_fp64_view','finite_fp64','flatnonzero_empty','flatnonzero_nonempty','index_add','fp64_gather']}


def profile():
    return dict(cpu_fraction=.5,shared_dram_bytes_per_second=1e8,fsync_seconds=.003,
        writeback_service=dict(pagecache_seconds_per_byte=1e-9,storage_seconds_per_byte=1e-8))


def archive():
    return dict(field_schema=[['variant_index','<i8'],['chi2','<f8']],call_cpu_seconds=1e-4,byte_cpu_seconds=1e-9)


@pytest.mark.parametrize('count',[1,37,65536])
def test_seekable_npz_header_rewrites_match_actual_numpy_submission(count):
    work=jagwas_writer_work(count,count);part=work['part']
    sink=writer_primitives().CountingFile()
    np.savez(sink,variant_index=np.arange(count,dtype=np.int64),chi2=np.ones(count,np.float64))
    assert sink.size==part['file_bytes']
    assert sink.written==part['file_bytes']+sum(row['local_header_bytes'] for row in part['arrays'])
    rows=jagwas_archive_service(work,archive(),profile(),host_serial_fraction=.25)
    cpu=sum(row['seconds']*row.get('resources',{}).get('cpu',0.) for row in rows)
    assert cpu==pytest.approx(archive()['call_cpu_seconds']+archive()['byte_cpu_seconds']*16*count+1e-9*sink.written)
    assert rows[1]['seconds']*rows[1]['resources']['output']==pytest.approx(sink.size)
    assert rows[2]['seconds']==profile()['fsync_seconds']


@pytest.mark.parametrize('retained',[0,11])
def test_empty_chunks_still_pay_selection_but_only_nonempty_parts_commit(retained):
    work=jagwas_writer_work(13,retained)
    rows=jagwas_host_selection_service(work,bank(),cpu_fraction=.5,dram_bytes_per_second=1e8,host_serial_fraction=.25)
    assert len(rows)==5 and all(row['seconds']>0 for row in rows)
    total=sum(row['seconds']*row['resources']['cpu'] for row in rows)
    assert total==pytest.approx(5e-6+1e-9*(3*13+2*retained))
    assert sum(row['seconds']*row['resources']['host_serial'] for row in rows)==pytest.approx(total/4)
    writers=jagwas_archive_service(work,archive(),profile(),host_serial_fraction=.25)
    assert bool(writers)==bool(retained)
    shards=[shard(chunks=2),shard('cuda:2',chunks=2)]
    for item in shards:
        item['outputs']=[[dict(cells=13,retained=retained,selection=rows,writer=writers)] for _ in range(2)]
    result=run(shards,shared_capacities=dict(cpu=2.,dram=1e8,host_serial=1.,output=1e8))
    assert result['retained_variants']==4*retained
    assert result['parts']==4*bool(retained)


def test_other_archive_schema_and_non_durable_boundary_are_rejected():
    wrong=archive();wrong['field_schema'].append(['df','<f4'])
    with pytest.raises(ValueError,match='two-array'):
        jagwas_archive_service(jagwas_writer_work(10,10),wrong,profile(),host_serial_fraction=1.)
    with pytest.raises(ValueError,match='durable'):
        jagwas_archive_service(jagwas_writer_work(10,10,fsync=False),archive(),profile(),host_serial_fraction=1.)
    bad=bank();bad['finite_fp64']['unit_cpu_seconds']=float('nan')
    with pytest.raises(ValueError,match='primitive price'):
        jagwas_host_selection_service(jagwas_writer_work(10,10),bad,cpu_fraction=1.,dram_bytes_per_second=1.,host_serial_fraction=1.)
    with pytest.raises(ValueError,match='CPU'):
        jagwas_archive_service(jagwas_writer_work(10,10),archive(),dict(profile(),cpu_fraction=0.),host_serial_fraction=1.)
