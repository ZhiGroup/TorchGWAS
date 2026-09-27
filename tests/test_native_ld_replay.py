"""LD restart work follows the decoder's last-non-LD-base semantics."""
from collections import Counter
from copy import deepcopy
from unittest.mock import patch

import numpy as np
import pytest

from test_pgen_native_reader import write_pgen_mixed
from torchgwas import pgen_native
from torchgwas.decoder_work import (decoder_work, decoder_chunk_work,
    native_ld_replays, native_read_layout, native_reader_workspace)
from torchgwas.pgen_native_reader import NativePgenReader
from torchgwas.pgen_reader import read_header, unpack_genovec
from torchgwas.pgen_work_census import census


def fixture(path, n=257, base_form=1):
    forms=[base_form]+[2,3]*20+[4,2,3,2]
    categories=np.zeros((len(forms),n),dtype=np.uint8)
    if base_form!=4:categories[0,::2]=2
    categories[0,3]=3
    base=categories[0]
    for i in range(1,len(forms)):
        if forms[i]==4:
            categories[i,3]=1;base=categories[i]
        else:
            categories[i]=base
            if forms[i]==3:
                categories[i,base==0]=2;categories[i,base==2]=0
            categories[i,(i*3)%n]=i%4
    write_pgen_mixed(path,categories,forms)
    return categories,forms


@pytest.mark.parametrize('n',[32,257,1025])
@pytest.mark.parametrize('base_form',[0,1,4])
@pytest.mark.skipif(not pgen_native.available(),reason='native PGEN decoder is not built')
def test_long_prefix_decodes_only_one_base_and_matches_pgenlib(tmp_path,n,base_form):
    pgenlib=pytest.importorskip('pgenlib')
    path=tmp_path/'ld.pgen'
    categories,forms=fixture(path,n,base_form)
    with pgenlib.PgenReader(str(path).encode()) as reference, NativePgenReader(path) as reader:
        for start in [40,1,42,17,0,41,44,3]:
            end=min(len(forms),start+3)
            expected=np.empty((end-start,n),np.int8)
            reference.read_range(start,end,expected)
            actual=np.empty_like(expected)
            with patch.object(pgen_native,'decode_range',wraps=pgen_native.decode_range) as decode:
                reader.read_range(start,end,actual)
            np.testing.assert_array_equal(actual,expected)
            wanted=categories[start:end].astype(np.int8)
            wanted[wanted==3]=-9
            np.testing.assert_array_equal(actual,wanted)
            if forms[start] in (2,3):
                assert len(decode.call_args_list)==2
                base=decode.call_args_list[0]
                assert len(base.args[1])==1 and base.args[6].shape==(1,(n+3)//4)
                assert base.kwargs['first_variant']==(0 if start<41 else 41)
            else:assert len(decode.call_args_list)==1
            packed=np.full((end-start,(n+3)//4+5),255,np.uint8)
            reader.read_packed_range_into(start,end,packed)
            for i,row in enumerate(categories[start:end]):
                # Form-3 inversion may set unused bits in the final byte.
                np.testing.assert_array_equal(unpack_genovec(packed[i,:(n+3)//4],n),row)
            assert not packed[:,(n+3)//4:].any()


@pytest.mark.parametrize('span',[(0,45),(1,40),(17,45)])
@pytest.mark.parametrize('chunk',[1,3,8,64])
def test_exact_base_counts_conserve_while_read_prefix_is_not_decoded(tmp_path,span,chunk):
    path=tmp_path/'ld.pgen';fixture(path)
    report=census(path,chunk,span,include_chunks=True)
    header=read_header(path)
    rows=native_ld_replays(report)
    for row in rows:
        at,base=row['chunk_start'],row['base_variant']
        assert row['base_record']==census(path,chunk,(base,base+1))
        assert row['read_prefix_bytes']==int(header.record_lengths[base:at].sum())
    plain=decoder_work(report,'torch_native_int8')
    restarted=decoder_work(report,'torch_native_int8',restart_ld_bases=True)
    expected=Counter(plain['source_units'])
    for row in rows:
        base=decoder_work(row['base_record'],'torch_native_int8')['source_units']
        expected.update({k:v for k,v in base.items() if not k.startswith('expand')})
    assert restarted['source_units']==dict(expected)
    assert restarted['logical_final_int8_bytes']==257*(span[1]-span[0])
    assert restarted['uncounted_mechanisms']==[]
    parts=decoder_chunk_work(report,chunk)
    units=Counter()
    for part in parts:
        units.update(part['decoder']['source_units'])
        assert part['census']==census(path,chunk,tuple(part['census']['variant_range']))
    assert dict(units)==dict(expected)
    for part in report['chunks']:
        layout=native_read_layout(part)
        lo,hi=part['variant_range']
        base=next((r['base_variant'] for r in rows if r['chunk_start']==lo),lo)
        assert layout['read_bytes']==int(header.record_lengths[base:hi].sum())
        expected_decode=int(header.record_lengths[lo:hi].sum())+(int(header.record_lengths[base]) if base<lo else 0)
        assert layout['decode_input_bytes']==expected_decode
        assert layout['extra_workspace_bytes']==((257+3)//4+8*(lo-base) if base<lo else 0)


@pytest.mark.parametrize('damage',['missing','count','bool_count','negative_count','prefix','base','order','span','geometry','form','file','sum','extra'])
def test_replay_ledger_refuses_incomplete_or_inconsistent_work(tmp_path,damage):
    path=tmp_path/'ld.pgen';fixture(path)
    report=census(path,8,include_chunks=True)
    row=report['native_ld_replays'][0]
    if damage=='missing':report.pop('native_ld_replays')
    elif damage=='count':report['ld_records_at_chunk_starts']+=1
    elif damage=='bool_count':report['ld_records_at_chunk_starts']=True
    elif damage=='negative_count':report['ld_records_at_chunk_starts']=-1
    elif damage=='prefix':row['read_prefix_bytes']=0
    elif damage=='base':row['base_variant']=row['chunk_start']
    elif damage=='order':report['native_ld_replays'].reverse()
    elif damage=='span':row['chunk_start']+=1
    elif damage=='geometry':report['chunk_markers']=0
    elif damage=='form':row['base_record']['record_form_counts']={2:1}
    elif damage=='file':row['base_record']['file_bytes']+=1
    elif damage=='sum':report['additional_base_record_bytes_if_every_chunk_restarts']+=1
    else:row['ignored']=123
    with pytest.raises(ValueError):native_read_layout(report)


def test_chunk_replay_read_bytes_must_conserve_parent(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    report=census(path,8,include_chunks=True)
    report['chunks'][1]['native_ld_replays'][0]['read_prefix_bytes']+=1
    with pytest.raises(ValueError,match='conserve native LD'):
        decoder_chunk_work(report,8)


def test_old_ld_census_remains_explicitly_unpriced(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    report=census(path,8)
    report.pop('native_ld_replays')
    work=decoder_work(report,'torch_native_int8',restart_ld_bases=True)
    assert work['uncounted_mechanisms']==['full source census of replayed LD-base records at reader restarts']


def test_tiled_memory_accounts_for_full_prefix_and_one_scratch_row(tmp_path):
    from test_trait_tiling_model import candidate
    from torchgwas.trait_tiling_model import trait_tiled_memory
    path=tmp_path/'ld.pgen';fixture(path,n=32)
    c=candidate(path,width=2,count=2,traits=5)
    b=c['tiles'][0]['profile']['chunk_markers']
    report=census(path,b,include_chunks=True)
    for tile in c['tiles']:
        tile['data'].update(samples=32,markers=45,encoded=report)
    actual=trait_tiled_memory(c)
    without=deepcopy(c)
    for tile in without['tiles']:
        for row in tile['data']['encoded']['chunks']:
            row.pop('native_ld_replays',None);row['ld_records_at_chunk_starts']=0
    baseline=trait_tiled_memory(without)
    payload,extra=native_reader_workspace(report['chunks'])
    old_payload=max(row['record_payload_bytes'] for row in report['chunks'])
    readers=sum(min(c['tiles'][i]['profile']['decode_workers'],c['tiles'][i]['profile']['depth'],(45+b-1)//b) for i in range(2))
    assert actual['host_bytes']-baseline['host_bytes']==readers*(payload+extra-old_payload)
    assert actual['device_bytes']==baseline['device_bytes']


def test_scan_charges_read_prefix_separately_from_decoded_base(tmp_path):
    from test_mechanistic_shapes import fixture as scan_fixture,component
    from torchgwas.mechanistic_torch import torch_scan_work
    path=tmp_path/'ld.pgen';fixture(path,n=32)
    data,profile=scan_fixture()
    report=census(path,8,(0,10),include_chunks=True)
    data['encoded']=report
    counts=decoder_work(report,'torch_native_int8',restart_ld_bases=True)
    profile['decode_units']={name:1e-9 for name in counts['source_units']}
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(data,profile)
    for block,part in zip(work['blocks'],report['chunks']):
        assert block['decode_read_seconds']==native_read_layout(part)['read_bytes']/profile['read_bytes_per_second']
    assert work['decoder']['native_ld_replay_records']==1
    assert work['host_workspace']['ld_replay_workspace_bytes_upper']==2*(8+8*8)
    assert any('LD-base restart' in term for term in work['allocator_unpriced_terms'])


def test_full_runtime_exposes_ld_restart_allocation_and_dispatch_limits(tmp_path):
    from test_mechanistic_shapes import fixture as scan_fixture,component
    from torchgwas.mechanistic_torch import torch_runtime
    path=tmp_path/'ld.pgen';fixture(path,n=32)
    data,profile=scan_fixture(k=1,c=8)
    data['encoded']=census(path,8,(0,10),include_chunks=True)
    data.update(traits_in_file=1,tables={'pvar':{'bytes':200,'field_characters':100},'psam':{'bytes':100}})
    counts=decoder_work(data['encoded'],'torch_native_int8',restart_ld_bases=True)
    profile['decode_units']={name:1e-9 for name in counts['source_units']}
    profile.update(write_bytes_per_second=1e8,npy_itemsize=4,initial_index_parses=2,
        timing_boundary='environment-ready',pin_cpu_seconds_per_page=1e-6,
        pin_driver_seconds_per_page=1e-6,tiny_setup_cpu_seconds=1e-5,
        tiny_setup_non_cpu_seconds=1e-5,first_use_seconds=.01,design_first_use_seconds=.1,fsync_seconds=.001)
    profile['gpu_resources']['gpu_fraction']=1.
    for name in ['pvar_rows','psam_rows','phenotype_qc_cells','covariate_qc_cells',
                 'covariate_basis_work','bytearray_zero_bytes','variant_id_rows']:
        profile['process_units'][name]=1e-8
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        result=torch_runtime(data,profile)
    assert result['estimated_seconds']>0
    assert any('LD-base restart' in term for term in result['unpriced_terms'])
    assert any('exact base records' in term for term in result['assumptions'])
