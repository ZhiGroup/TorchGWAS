"""Empty output and concurrent phenotype tiles retain exact producer identity."""
from dataclasses import FrozenInstanceError
from unittest.mock import patch
import numpy as np
import pytest

from torchgwas.api import _trait_blocked_significant_chunks
from torchgwas.sumstats_indexed import IndexedOutputPartition,PartitionedIndexedChunk,write_indexed_sumstats


def selected(start=10,end=12,rows=0):
    return (start,end,np.full(rows,start,dtype=np.int64),np.zeros(rows,dtype=np.int64),
            np.ones(rows,dtype=np.float32),np.ones(rows,dtype=np.float32),17)


@pytest.mark.parametrize('dtype',[np.int32,np.int64])
@pytest.mark.parametrize('layout',['contiguous','strided','readonly','empty'])
@pytest.mark.parametrize('devices',[['cuda:0'],['cuda:0','cuda:1']])
def test_trait_rebase_owns_outputs_without_mutating_borrowed_indices(dtype,layout,devices):
    base=np.zeros(8,dtype=dtype)
    incoming=base if layout=='contiguous' else base[::2]
    if layout=='readonly':incoming.flags.writeable=False
    if layout=='empty':incoming=incoming[:0]
    def scan(offset,width,device):
        for start in [10,12]:
            chunk=selected(start,start+2,len(incoming))
            yield (*chunk[:3],incoming,*chunk[4:])
    with patch('torchgwas.linear._significant_pairs_iterator',side_effect=lambda source,*a:source):
        chunks=list(_trait_blocked_significant_chunks(scan,None,5,2,17,devices=devices,queue_depth=1))
    assert len(chunks)==6
    np.testing.assert_array_equal(base,np.zeros_like(base))
    outputs=[chunk[3] for chunk in chunks]
    for array in outputs:
        assert array.dtype==np.int64 and array.flags.writeable
        assert not np.shares_memory(array,base)
    for i,array in enumerate(outputs):
        assert all(not np.shares_memory(array,other) for other in outputs[i+1:])
    if incoming.size:
        assert sorted(int(array[0]) for array in outputs)==[0,0,2,2,4,4]
        expected=[array.copy() for array in outputs]
        base.fill(99)
        for array,want in zip(outputs,expected):np.testing.assert_array_equal(array,want)


@pytest.mark.parametrize('devices',[['cuda:0'],['cuda:0','cuda:1']])
def test_producer_identity_survives_empty_chunks_and_shared_queue(tmp_path,devices):
    def scan(offset,width,device):yield selected(rows=0 if offset==0 else 1)
    def context(first,width,device):
        return IndexedOutputPartition(device or devices[0],(10,20),(first,first+width))
    with patch('torchgwas.linear._significant_pairs_iterator',side_effect=lambda source,*a:source):
        chunks=list(_trait_blocked_significant_chunks(scan,None,5,2,17,devices=devices,
            queue_depth=1,partition_context=context))
    assert len(chunks)==3 and all(isinstance(row,PartitionedIndexedChunk) for row in chunks)
    assert sorted((row.partition.trait_range,len(row[2])) for row in chunks)==[((0,2),0),((2,4),1),((4,5),1)]
    for chunk in chunks:
        assert chunk.partition.device==devices[(chunk.partition.trait_range[0]//2)%len(devices)]
        payload=(chunk[0]-10,chunk[1]-10,chunk[2]-10,*chunk[3:])
        translated=chunk.with_payload(payload)
        assert translated.partition is chunk.partition and translated[4] is chunk[4]
        events=[]
        directory=tmp_path/str(chunk.partition.trait_range[0])
        write_indexed_sumstats(directory,['a','b'],list('abcde'),20,[translated],kind='significant',df=17,
            on_chunk_written=events.append,variant_offset=10)
        event=events[0]
        assert (event.start,event.end)==(0,2) and event.source_variant_range==(10,12)
        assert event.partition==chunk.partition and event.rows==len(chunk[2])
        with pytest.raises(FrozenInstanceError):event.partition.device='cuda:9'
        if event.rows:
            with np.load(directory/event.part_file) as saved:
                assert saved['trait_index'][0]==chunk.partition.trait_range[0]
                assert saved['variant_index'][0]==0


def test_variant_resolver_uses_absolute_source_range_before_writing(tmp_path):
    calls=[];partition=IndexedOutputPartition('cuda:1',(10,12),(0,2));events=[]
    def resolve(a,b):calls.append((a,b));return partition
    write_indexed_sumstats(tmp_path,['a','b'],['x','y'],20,[(0,2,None,np.ones(2),None)],
        kind='jagwas',df=17,chi2_df=2,variant_offset=10,partition_for_range=resolve,on_chunk_written=events.append)
    assert calls==[(10,12)] and events[0].partition is partition


@pytest.mark.parametrize('fault',['conflict','outside','jagwas_tile','resolver_none'])
def test_invalid_identity_refuses_output_before_a_part_is_written(tmp_path,fault):
    partition=IndexedOutputPartition('cuda:0',(10,12),(0,2))
    resolved=(None if fault=='resolver_none' else IndexedOutputPartition('cuda:1',(10,12),(0,2)))
    if fault=='outside':partition=IndexedOutputPartition('cuda:0',(12,14),(0,2))
    if fault=='jagwas_tile':partition=IndexedOutputPartition('cuda:0',(10,12),(0,1))
    chunk=PartitionedIndexedChunk((0,2,None,np.ones(2),None),partition)
    with pytest.raises(ValueError):
        write_indexed_sumstats(tmp_path,['a','b'],['x','y'],20,[chunk],kind='jagwas',df=17,
            variant_offset=10,on_chunk_written=lambda event:pytest.fail('bad completion'),
            partition_for_range=(lambda a,b:resolved) if fault in ('conflict','resolver_none') else None)
    assert not list(tmp_path.glob('*.npz')) and not (tmp_path/'manifest.json').exists()


def test_partitioned_generator_closes_when_writer_fails(tmp_path):
    closed=[]
    def source():
        try:yield PartitionedIndexedChunk(selected(0,2,1),IndexedOutputPartition('cuda:0',(0,2),(0,1)))
        finally:closed.append(True)
    with patch('torchgwas.sumstats_indexed.os.fsync',side_effect=OSError('writer failure')):
        with pytest.raises(OSError):
            write_indexed_sumstats(tmp_path,['a','b'],['x'],20,source(),kind='significant',df=17,
                on_chunk_written=lambda row:pytest.fail('failed write has no event'))
    assert closed==[True]
