"""Adaptive ranges conserve independently parsed source and output work."""
from collections import Counter
from copy import deepcopy
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.pgen_work_census import census,rechunk_census,scheduled_census
from torchgwas.decoder_work import decoder_chunk_work,decoder_work,native_read_layout
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.binary_output_work import binary_output_work
from torchgwas.trait_tiling_model import trait_tiled_memory,trait_tiled_graph_chunks
from torchgwas.significant_host_work import significant_host_memory,indexed_part_work
from torchgwas.jagwas_candidate import jagwas_candidate_memory,jagwas_candidate_runtime
from torchgwas.reduced_output_work import jagwas_indexed_part_work
from test_native_ld_replay import fixture
from test_trait_tiling_model import candidate as dense_candidate,runtime as dense_runtime
from test_significant_host_model import runtime as significant_runtime
from test_jagwas_actual_candidate import actual_candidate,writer_prices
from test_jagwas_candidate import preparation
from test_pgen_native_reader import write_pgen
from test_mechanistic_shapes import component


@pytest.mark.parametrize('base_form',[0,1,4])
@pytest.mark.parametrize('span',[(0,45),(1,44),(3,45)])
def test_scheduled_counts_equal_independent_parser_per_range(tmp_path,base_form,span):
    path=tmp_path/'ld.pgen';fixture(path,base_form=base_form)
    source=census(path,2,span,include_chunks=True);before=deepcopy(source)
    ranges=[];cursor=span[0];sizes=[2,8,4,6];i=0
    while cursor<span[1]:
        end=min(cursor+sizes[i%len(sizes)],span[1]);ranges.append([cursor,end]);cursor=end;i+=1
    actual=scheduled_census(source,8,ranges)
    expected=[census(path,8,row,include_chunks=True)['chunks'][0] for row in ranges]
    assert actual['chunks']==expected and actual['chunk_ranges']==ranges
    assert source==before
    parts=decoder_chunk_work(actual,8)
    assert [row['census'] for row in parts]==expected
    units=Counter()
    for row in expected:units.update(decoder_work(row,'torch_native_int8',restart_ld_bases=True)['source_units'])
    assert dict(units)==decoder_work(actual,'torch_native_int8',restart_ld_bases=True)['source_units']
    assert native_read_layout(actual)['read_bytes']==sum(native_read_layout(row)['read_bytes'] for row in expected)
    assert actual['native_ld_replays']==[r for row in expected for r in row.get('native_ld_replays',[])]
    # Parent, child, caller ranges and original source do not alias.
    ranges[0][1]+=1
    actual['chunks'][1]['record_form_counts'][0]=999
    assert source==before and actual['chunk_ranges'][0][1]!=ranges[0][1]


def test_regular_schedule_matches_existing_view_and_remaining_subrange(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,2,include_chunks=True)
    ranges=[[lo,min(lo+8,45)] for lo in range(2,45,8)]
    actual=scheduled_census(source,8,ranges)
    assert actual.pop('chunk_ranges')==ranges
    assert actual==rechunk_census(source,8,(2,45))


@pytest.mark.parametrize('ranges',[
    [],[[False,4]],[[0]],[[0,0]],[[0,9]],[[0,4],[6,8]],[[0,4],[2,8]],
    [[0,4],[4,46]],[[0,3],[3,8]],[[0,4.0]],[[0,4],None]])
def test_malformed_or_unaligned_schedules_refused(tmp_path,ranges):
    path=tmp_path/'ld.pgen';fixture(path)
    with pytest.raises(ValueError):scheduled_census(census(path,2,include_chunks=True),8,ranges)


def test_budget_and_tampered_explicit_layout_refused(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,2,include_chunks=True)
    with pytest.raises(ValueError,match='budget'):scheduled_census(source,8,[[0,4]],max_source_chunks=2)
    actual=scheduled_census(source,8,[[0,4],[4,12]])
    for damage in ['missing_counts','range','capacity','child','replay']:
        broken=deepcopy(actual)
        if damage=='missing_counts':broken.pop('chunks')
        elif damage=='range':broken['chunk_ranges'][1][0]=3
        elif damage=='capacity':broken['chunk_ranges']=[[0,12]]
        elif damage=='child':broken['chunks'][1]['variant_range']=[8,12]
        else:broken['native_ld_replays'][0]['chunk_start']=5
        with pytest.raises(ValueError):decoder_chunk_work(broken,8)


def make_dense(tmp_path,count=2,block_bytes=None):
    path=tmp_path/'plain.pgen'
    write_pgen(path,np.arange(10*32,dtype=np.uint8).reshape(10,32)%3)
    choice=dense_candidate(path,count=count,block_bytes=block_bytes)
    encoded=scheduled_census(census(path,2,include_chunks=True),4,[[0,2],[2,6],[6,8],[8,10]])
    for tile in choice['tiles']:tile['data']['encoded']=deepcopy(encoded)
    return choice


@pytest.mark.parametrize('count',[1,2])
@pytest.mark.parametrize('block_bytes',[None,48])
def test_dense_and_significant_output_follow_actual_rows(tmp_path,count,block_bytes):
    choice=make_dense(tmp_path,count,block_bytes)
    before=deepcopy(choice)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(choice['tiles'][0]['data'],choice['tiles'][0]['profile'])
    assert [b['markers'] for b in work['blocks']]==[2,4,2,2]
    assert [b['variant_range'] for b in work['blocks']]==[[0,2],[2,6],[6,8],[8,10]]
    assert sum('reader_init_seconds' in b for b in work['blocks'])==work['workers']
    assert sum(b['h2d_bytes'] for b in work['blocks'])==32*10
    dense=dense_runtime(choice)
    assert [r['blocks'] for r in dense['tiles']]==[4,4,4]
    assert dense['binary_payload_bytes']==8*10*5+4*10*3
    # First and second equal-width tiles have different pinned-cache states;
    # the narrower third tile also has distinct work. All three are expanded.
    assert trait_tiled_graph_chunks(choice)==12
    assert choice==before
    with pytest.raises(ValueError,match='adaptive_candidate_memory'):trait_tiled_memory(choice)
    choice['output']['block_bytes']=None
    reduced=significant_runtime(choice)
    expected=sum(indexed_part_work(rows*tile['data']['traits_analyzed'])['file_bytes']
                 for tile in choice['tiles'] for rows in [2,4,2,2])
    assert reduced['parts']==12 and reduced['indexed_part_bytes']==expected
    assert significant_runtime(choice,'empty')['parts']==0
    with pytest.raises(ValueError,match='bounded graph'):significant_runtime(choice,max_source_chunks=11)
    with pytest.raises(ValueError,match='adaptive_candidate_memory'):significant_host_memory(choice)


@pytest.mark.parametrize('count',[1,2])
def test_joint_actual_geometry_full_panel_and_shifted_source_chunks(tmp_path,count):
    path=tmp_path/'joint.pgen';n,m=2049,1025
    write_pgen(path,(np.arange(n*m,dtype=np.uint32).reshape(m,n)%3).astype(np.uint8))
    choice=actual_candidate(path,512,count);source=census(path,128,include_chunks=True)
    all_rows=[]
    for tile in choice['tiles']:
        lo,hi=tile['variant_range'];ranges=[];cursor=lo
        for width in [128,256,512,128,1]:
            if cursor==hi:break
            # The last shard may contain only the file tail.
            end=min(cursor+width,hi);ranges.append([cursor,end]);cursor=end
        assert cursor==hi
        tile['data']['encoded']=scheduled_census(source,512,ranges)
        all_rows.extend(hi-lo for lo,hi in ranges)
        work=torch_scan_work(tile['data'],tile['profile'])
        assert sum(b['d2h_bytes'] for b in work['blocks'])==17*tile['data']['markers']
        assert sum(b['h2d_bytes'] for b in work['blocks'])==n*tile['data']['markers']
        assert [b['variant_range'] for b in work['blocks']]==ranges
    result=jagwas_candidate_runtime(choice,writer_prices(),preparation=preparation(choice),
        occupancy='dense',host_serial_fraction=.5)
    assert result['parts']==len(all_rows) and result['retained_variants']==m
    assert result['indexed_part_bytes']==sum(jagwas_indexed_part_work(rows)['file_bytes'] for rows in all_rows)
    assert result['factor_preparation']['instances']==count
    with pytest.raises(ValueError,match='bounded source-chunk'):
        jagwas_candidate_runtime(choice,writer_prices(),preparation=preparation(choice),
            occupancy='dense',host_serial_fraction=.5,max_source_chunks=len(all_rows)-1)
    with pytest.raises(ValueError,match='adaptive_candidate_memory'):jagwas_candidate_memory(choice)


def test_auto_writer_block_uses_first_actual_chunk_and_df_matches_layout():
    sizes=[128,512,256,1];k=4097
    work=binary_output_work(sum(sizes),k,512,borrow_chunks=False,store_variant_df=True,chunk_rows=sizes)
    assert work['block_bytes']==128*k*4
    assert work['binary_payload_bytes']==(8*k+4)*sum(sizes)
    assert work['staging_copy_bytes']==work['binary_payload_bytes']
    assert [r['markers'] for r in work['stream_work']['df']['chunks']]==sizes
    assert [r['markers'] for r in work['chunks']]==sizes
    assert work['allocated_staging_bytes']<binary_output_work(sum(sizes),k,512,store_variant_df=True)['allocated_staging_bytes']
    for bad in ([128,513,256],[0,512],[True,512],[128,512]):
        with pytest.raises(ValueError):binary_output_work(sum(sizes),k,512,chunk_rows=bad)


def test_actual_chunk_count_sets_active_reader_initializations(tmp_path):
    choice=make_dense(tmp_path,1);tile=choice['tiles'][0]
    profile=tile['profile'];profile.update(chunk_markers=16,depth=4,decode_workers=4)
    path=tile['data']['encoded']['path']
    tile['data']['encoded']=scheduled_census(census(path,2,include_chunks=True),16,
        [[0,2],[2,6],[6,8],[8,10]])
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(tile['data'],profile)
    assert work['workers']==4
    assert sum(b.get('reader_init_seconds',0.) for b in work['blocks'])==pytest.approx(
        4*10*profile['process_units']['pgen_index_records']/profile['cpu_fraction'])
