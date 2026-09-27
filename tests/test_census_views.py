"""Aligned source-count views must equal independently parsed PGEN censuses."""
from copy import deepcopy
from unittest.mock import patch

import pytest

from test_native_ld_replay import fixture
from test_trait_candidate_space import spec,prepare
from torchgwas.pgen_work_census import census,rechunk_census
from torchgwas.decoder_work import decoder_chunk_work


@pytest.mark.parametrize('fine,coarse',[(1,1),(1,7),(2,8),(3,9),(4,16),(8,24),(16,64)])
@pytest.mark.parametrize('span',['full','middle','tail'])
@pytest.mark.parametrize('base_form',[0,1,4])
def test_rechunk_and_range_are_exact_for_ld_and_sample_tails(tmp_path,fine,coarse,span,base_form):
    path=tmp_path/'ld.pgen';fixture(path,base_form=base_form)
    source=census(path,fine,include_chunks=True);before=deepcopy(source)
    extent=(0,45) if span=='full' else ((fine,45) if span=='tail' else (fine,min(2*fine,45)))
    actual=rechunk_census(source,coarse,extent)
    expected=census(path,coarse,extent,include_chunks=True)
    assert actual==expected
    assert source==before
    assert decoder_chunk_work(actual,coarse)==decoder_chunk_work(expected,coarse)


def test_shifted_source_range_and_nested_views_keep_correct_ld_bases(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,2,(1,44),include_chunks=True)
    coarse=rechunk_census(source,4,(3,44))
    assert coarse==census(path,4,(3,44),include_chunks=True)
    nested=rechunk_census(coarse,8,(7,44))
    assert nested==census(path,8,(7,44),include_chunks=True)


def test_view_does_not_alias_source_or_sibling_records(tmp_path):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,1,include_chunks=True);before=deepcopy(source)
    view=rechunk_census(source,8);sibling=rechunk_census(source,8)
    view['record_form_counts'][0]=123456
    view['chunks'][1]['native_ld_replays'][0]['base_record']['record_form_counts'][1]=123456
    assert source==before and sibling==census(path,8,include_chunks=True)
    assert view['native_ld_replays'][0]['base_record']['record_form_counts'][1]==1


@pytest.mark.parametrize('size,span',[(0,None),(True,None),(3,None),(4,(1,44)),(4,(0,43)),(4,(2,46)),(4,(4,4))])
def test_unaligned_or_invalid_view_is_refused(tmp_path,size,span):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,2,include_chunks=True)
    with pytest.raises(ValueError):rechunk_census(source,size,span)


@pytest.mark.parametrize('damage',['range','file','chunk_count','markers','missing_ld'])
def test_incomplete_source_chunks_are_refused(tmp_path,damage):
    path=tmp_path/'ld.pgen';fixture(path)
    source=census(path,2,include_chunks=True)
    if damage=='range':source['chunks'][1]['variant_range'][0]+=1
    elif damage=='file':source['chunks'][1]['file_bytes']+=1
    elif damage=='chunk_count':source['chunks'].pop()
    elif damage=='markers':source['markers']+=1
    else:source['chunks'][1].pop('native_ld_replays')
    with pytest.raises(ValueError):rechunk_census(source,8)


def test_unrelated_chunk_grids_keep_independent_full_passes_but_reuse_shards(tmp_path):
    value=spec(tmp_path)
    with patch('torchgwas.trait_candidate_space.census',wraps=census) as collect:
        result=prepare(value,chunks=[2,3,4],partition_axes=['trait','variant'])
    assert collect.call_count==2 and result['census_passes']==2
    assert result['census_views']==9
    for candidate in result['candidates']:
        for tile in candidate['tiles']:
            encoded=tile['data']['encoded'];span=tuple(encoded['variant_range'])
            expected=census(value['workload']['genotype'],tile['profile']['chunk_markers'],span,include_chunks=True)
            assert encoded==expected
