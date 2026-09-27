"""Header-only admission matches exact reader extents without payload access."""
import builtins
from unittest.mock import patch
import numpy as np
import pytest

from torchgwas.pgen_memory_layout import (memory_layout,rechunk_memory_layout,
    compact_rechunk_memory_layout,compact_shifted_reader_envelope)
from torchgwas.pgen_work_census import census
from torchgwas.pgen_reader import read_header
from torchgwas.decoder_work import native_read_layout,native_reader_workspace,decoder_work
from torchgwas.adaptive_candidate import _reader_envelope
from test_pgen_native_reader import write_pgen_mixed


def fixture(tmp_path):
    path=tmp_path/'memory.pgen';n,m=129,137
    calls=np.tile((np.arange(n)%3).astype(np.uint8),(m,1))
    forms=[0 if i%17==0 else 2+(i%2) for i in range(m)]
    for i in range(m):calls[i,i%n]=(calls[i,i%n]+1)%3
    write_pgen_mixed(path,calls,forms)
    return path


def test_admission_never_reads_or_maps_genotype_payload(tmp_path):
    path=fixture(tmp_path);header=read_header(path);first_payload=int(header.record_offsets[0])
    original=builtins.open;reads=[]
    class CheckedFile:
        def __init__(self,file):self.file=file
        def __enter__(self):self.file.__enter__();return self
        def __exit__(self,*args):return self.file.__exit__(*args)
        def read(self,size=-1):
            lo=self.file.tell()
            assert size>=0 and lo+size<=first_payload,'Genotype payload read during admission'
            value=self.file.read(size);reads.append((lo,len(value)));return value
        def __getattr__(self,key):return getattr(self.file,key)
    def checked(*args,**kwargs):return CheckedFile(original(*args,**kwargs))
    with patch('builtins.open',side_effect=checked),\
         patch('mmap.mmap',side_effect=AssertionError('Admission must not map genotype payloads')):
        result=memory_layout(path,4)
    assert reads and max(start+size for start,size in reads)<=first_payload
    assert result['index_bytes']==first_payload
    assert 'source_work' not in result and 'record_form_counts' not in result


@pytest.mark.parametrize('width',[1,4,7,16,128,137,200])
def test_vectorized_layout_preserves_each_ld_restart_extent(tmp_path,width):
    path=fixture(tmp_path);h=read_header(path)
    layout=memory_layout(path,width)
    bases=np.flatnonzero((h.vrtypes!=2)&(h.vrtypes!=3))
    final=int(h.record_offsets[-1])+int(h.record_lengths[-1])
    assert len(layout['chunks'])==(h.variant_ct+width-1)//width
    for row in layout['chunks']:
        lo,hi=row['variant_range']
        end=int(h.record_offsets[hi]) if hi<h.variant_ct else final
        assert row['record_payload_bytes']==end-int(h.record_offsets[lo])
        base=int(bases[np.searchsorted(bases,lo,side='right')-1])
        prefix=(int(h.record_offsets[lo])-int(h.record_offsets[base])
            if int(h.vrtypes[lo]) in (2,3) else 0)
        extra=((h.sample_ct+3)//4+8*(lo-base) if prefix else 0)
        assert (row['read_prefix_bytes'],row['reader_extra_workspace_bytes'])==(prefix,extra)


@pytest.mark.parametrize('fine,capacity',[(1,1),(1,7),(2,8),(4,4),(4,16),(4,32)])
def test_compact_memory_maxima_equal_full_grid_for_fixed_and_shifted_reads(tmp_path,fine,capacity):
    path=fixture(tmp_path);header=read_header(path)
    full=memory_layout(path,fine,header=header)
    compact=memory_layout(path,fine,header=header,compact=True)
    for lo,hi in ((0,137),(0,136),(4 if fine<=4 else fine,137)):
        if lo%fine or (hi%fine and hi!=137):continue
        fixed=rechunk_memory_layout(full,capacity,(lo,hi))
        summarized=compact_rechunk_memory_layout(compact,capacity,(lo,hi))
        assert summarized['record_payload_bytes']==fixed['record_payload_bytes']
        assert summarized['variant_range']==fixed['variant_range']
        assert native_reader_workspace(summarized['chunks'])==native_reader_workspace(fixed['chunks'])
        parts=full['chunks'][lo//fine:(hi+fine-1)//fine]
        assert compact_shifted_reader_envelope(compact,(lo,hi),capacity)==_reader_envelope(
            parts,capacity=capacity,fine=fine)


@pytest.mark.parametrize('capacity',[4,8,16,32])
def test_header_grouping_matches_exact_payload_and_ld_scratch_for_every_start(tmp_path,capacity):
    path=fixture(tmp_path);fine=memory_layout(path,4)
    for start in range(0,137,4):
        stop=min(start+capacity,137)
        grouped=rechunk_memory_layout(fine,capacity,(start,stop))
        exact=census(path,capacity,variant_range=(start,stop),include_chunks=True)
        assert grouped['record_payload_bytes']==exact['record_payload_bytes']
        for a,b in zip(grouped['chunks'],exact['chunks']):
            assert a['variant_range']==b['variant_range']
            known,measured=native_read_layout(a),native_read_layout(b)
            assert known['read_bytes']==measured['read_bytes']
            assert known['extra_workspace_bytes']==measured['extra_workspace_bytes']
            assert known['decode_input_bytes'] is None


def test_memory_extents_cannot_be_mistaken_for_decoder_work(tmp_path):
    path=fixture(tmp_path);layout=memory_layout(path,4)
    with pytest.raises(ValueError,match='cannot price decoder work'):
        decoder_work(layout,'torch_native_int8')
    from torchgwas.mechanistic_torch import torch_scan_work
    with pytest.raises(ValueError,match='cannot price scan work'):
        torch_scan_work({'encoded':layout},{})
    from torchgwas.adaptive_candidate import future_chunk_candidate
    with pytest.raises(ValueError,match='Memory-only'):
        future_chunk_candidate({},source_census=layout,chunk_sizes=[4,8],next_size=8,
            issued_chunks=[0],reduction=None)


@pytest.mark.parametrize('capacity,span',[(3,None),(8,(1,8)),(8,(0,7)),(8,(0,138)),(8,(8,4))])
def test_invalid_grouping_is_rejected(tmp_path,capacity,span):
    with pytest.raises(ValueError):rechunk_memory_layout(memory_layout(fixture(tmp_path),4),capacity,span)


def test_compact_admission_rejects_unsupported_record_type(tmp_path):
    from dataclasses import replace
    from torchgwas.pgen_reader import PgenFormatError
    path=fixture(tmp_path);header=read_header(path)
    types=header.vrtypes.copy();types[3]=5
    with pytest.raises(PgenFormatError,match='biallelic unphased hardcalls'):
        memory_layout(path,4,header=replace(header,vrtypes=types),compact=True)


@pytest.mark.parametrize('long_replay',[False,True])
def test_deferred_base_lookup_matches_full_compact_vectors(tmp_path,long_replay):
    from dataclasses import replace
    path=fixture(tmp_path);header=read_header(path)
    if long_replay:
        types=header.vrtypes.copy();types[:]=2;types[0]=0
        header=replace(header,vrtypes=types)
    for fine in (1,4,16):
        full=memory_layout(path,fine,header=header,compact=True)
        deferred=memory_layout(path,fine,header=header,compact=True,_defer_bases=True)
        for key in ('fine_payload_bytes','fine_prefix_bytes','fine_extra_workspace_bytes',
                    'payload_cumulative_bytes'):
            np.testing.assert_array_equal(deferred[key],full[key])
            assert not deferred[key].flags.writeable
        for capacity in (fine,4*fine):
            assert compact_rechunk_memory_layout(deferred,capacity)==compact_rechunk_memory_layout(full,capacity)
            assert compact_shifted_reader_envelope(deferred,(0,137),capacity)==compact_shifted_reader_envelope(full,(0,137),capacity)
