"""Shared-source floors retain physical rereads and full-panel JAGWAS."""
from copy import deepcopy

import numpy as np
import pytest

from torchgwas.pgen_reader import pack_genovec
from torchgwas.pgen_work_bounds import PgenHeaderWork,native_schedule_source_floor
from torchgwas.source_layout_floor import native_layout_source_floor
from test_pgen_work_bounds import fixture,write_records


def source(tmp_path, *, capacities=None):
    path=tmp_path/'plain.pgen';n=129
    records=[pack_genovec(np.full(n,i%3,dtype=np.uint8),n).tobytes()
             for i in range(12)]
    write_records(path,n,[0]*12,records)
    header=PgenHeaderWork(path)
    caps=(dict(cpu=1.,dram=5e7,input=5e6) if capacities is None
          else capacities)
    def floor(start,stop,chunk=3):
        schedule=header.schedule_bounds(start,stop,chunk)
        profile=dict(decode_units={name:1e-8 for name in schedule['source_units']},
            cpu_fraction=.5,depth=2,decode_workers=2,cpu_available_cores=2.,
            shared_dram_bytes_per_second=1e8,read_bytes_per_second=1e7,
            input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                       cpu_seconds_per_call=1e-5))
        return native_schedule_source_floor(schedule,profile,caps)
    return floor


def part(key,device,traits,floor):
    return dict(id=key,device=device,trait_range=list(traits),floor=floor)


def test_variant_shards_sum_shared_source_work_and_keep_jagwas_full_panel(tmp_path):
    floor=source(tmp_path);left=floor(0,6);right=floor(6,12);whole=floor(0,12)
    rows=[part('left','cuda:0',(0,128),left),
          part('right','cuda:1',(0,128),right)]
    combined=native_layout_source_floor(rows,total_traits=128,
        reduction='jagwas',partition_axis='variant')
    assert combined['resource_work']['input_bytes']==whole['resource_work']['input_bytes']
    assert combined['resource_work']['dram_bytes']==whole['resource_work']['dram_bytes']
    assert combined['resource_work']['cpu_seconds']==pytest.approx(
        [a+b for a,b in zip(left['resource_work']['cpu_seconds'],
                            right['resource_work']['cpu_seconds'])])
    assert set(combined['per_device_reader_floor_seconds'])=={'cuda:0','cuda:1'}
    assert combined['source_stage_floor_seconds'][0]>=combined['shared_resource_floor_seconds'][0]
    assert not combined['selection_validated']


def test_variant_shard_boundary_replays_the_same_ld_base_as_whole_scan(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path,129);header=PgenHeaderWork(path)
    schedules=[header.schedule_bounds(*span,3) for span in ((0,15),(0,6),(6,15))]
    def priced(schedule):
        profile=dict(decode_units={name:1e-8 for name in schedule['source_units']},
            cpu_fraction=.5,depth=2,decode_workers=2,cpu_available_cores=2.,
            shared_dram_bytes_per_second=1e8,read_bytes_per_second=1e7)
        return native_schedule_source_floor(schedule,profile,
            dict(cpu=1.,dram=5e7,input=5e6))
    full,left,right=map(priced,schedules)
    assert schedules[2]['ld_replay_count']>=1
    combined=native_layout_source_floor([
        part('left','cuda:0',(0,128),left),
        part('right','cuda:1',(0,128),right)],
        total_traits=128,reduction='jagwas',partition_axis='variant')
    assert combined['resource_work']['input_bytes']==full['resource_work']['input_bytes']
    assert combined['resource_work']['dram_bytes']==full['resource_work']['dram_bytes']


def test_trait_tiles_count_rereads_and_serial_per_device_workers(tmp_path):
    floor=source(tmp_path)(0,12)
    rows=[part('a','cuda:0',(0,64),floor),
          part('b','cuda:0',(64,128),floor)]
    serial=native_layout_source_floor(rows,total_traits=128,reduction=None,
                                      partition_axis='trait')
    assert serial['resource_work']['input_bytes']==2*floor['resource_work']['input_bytes']
    assert serial['per_device_reader_floor_seconds']['cuda:0']==pytest.approx(
        [2*x for x in floor['reader_worker_floor_seconds']])
    rows[1]['device']='cuda:1'
    parallel=native_layout_source_floor(rows,total_traits=128,reduction=None,
                                        partition_axis='trait')
    assert parallel['resource_work']==serial['resource_work']
    assert parallel['source_stage_floor_seconds'][0]<=serial['source_stage_floor_seconds'][0]
    assert parallel['shared_resource_floor_seconds']==serial['shared_resource_floor_seconds']


def test_large_trait_panel_uses_bounded_many_tile_layout(tmp_path):
    floor=source(tmp_path)(0,12)
    rows=[part(str(i),f'cuda:{i%4}',(i,i+1),floor) for i in range(96)]
    layout=native_layout_source_floor(rows,total_traits=96,reduction='significant',
                                      partition_axis='trait')
    assert len(layout['partitions'])==96
    assert layout['resource_work']['input_bytes']==96*floor['resource_work']['input_bytes']
    assert layout['per_device_reader_floor_seconds']['cuda:0']==pytest.approx(
        [24*x for x in floor['reader_worker_floor_seconds']])
    with pytest.raises(ValueError,match='Bounded'):
        native_layout_source_floor(rows,total_traits=96,reduction='significant',
                                   partition_axis='trait',max_partitions=95)


@pytest.mark.parametrize('damage',['jagwas_axis','jagwas_tile','overlap','source',
                                   'capacity','chunk','device','sample'])
def test_layout_source_floor_rejects_incompatible_partitions(tmp_path,damage):
    floor=source(tmp_path);left=floor(0,6);right=floor(6,12)
    rows=[part('a','cuda:0',(0,128),left),
          part('b','cuda:1',(0,128),right)]
    axis='variant';mode='jagwas'
    if damage=='jagwas_axis':axis='trait'
    elif damage=='jagwas_tile':rows[1]['trait_range']=[0,64]
    elif damage=='overlap':rows[1]['floor']=deepcopy(left)
    elif damage=='source':
        rows[1]['floor']=deepcopy(right)
        rows[1]['floor']['input_identity']['bytes']+=1
    elif damage=='capacity':
        rows[1]['floor']=deepcopy(right)
        rows[1]['floor']['capacities']['cpu']*=.5
    elif damage=='chunk':
        rows[1]['floor']=deepcopy(right)
        rows[1]['floor']['chunk_markers']=6
        rows[1]['floor']['chunk_count']=1
    elif damage=='sample':
        rows[1]['floor']=deepcopy(right)
        rows[1]['floor']['samples']+=1
    else:rows[1]['device']='cuda:0'
    with pytest.raises(ValueError):
        native_layout_source_floor(rows,total_traits=128,reduction=mode,
                                   partition_axis=axis)


def test_layout_source_floor_rejects_changed_file_after_partition_pricing(tmp_path):
    floor=source(tmp_path)(0,12)
    from pathlib import Path
    path=Path(floor['input_identity']['path'])
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='changed'):
        native_layout_source_floor([part('a','cuda:0',(0,128),floor)],
                                   total_traits=128,reduction='jagwas',
                                   partition_axis='variant')
