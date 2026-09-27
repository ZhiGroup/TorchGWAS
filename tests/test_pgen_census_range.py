import numpy as np
import pytest
from test_pgen_native_reader import write_pgen
from torchgwas.pgen_work_census import census


def test_shard_payload_and_full_index_are_distinct(tmp_path):
    path=tmp_path/'plain.pgen'
    write_pgen(path,np.arange(11*128,dtype=np.uint8).reshape(11,128)%4)
    full=census(path,chunk_markers=3)
    left=census(path,chunk_markers=3,variant_range=(0,6))
    right=census(path,chunk_markers=3,variant_range=(6,11))
    assert left['markers']==6 and right['markers']==5
    assert left['record_payload_bytes']+right['record_payload_bytes']==full['record_payload_bytes']
    assert right['file_markers']==11
    assert full['index_bytes']==left['index_bytes']==right['index_bytes']
    assert right['file_bytes']==full['file_bytes']
    assert right['record_form_counts']=={0:5}
    assert right['variant_range']==[6,11]
    for span in [(0,0),(-1,5),(0,12),(2,1),(True,2)]:
        with pytest.raises(ValueError,match='variant_range'):census(path,variant_range=span)


def test_range_start_finds_ld_base_before_shard_and_restarts_chunks(tmp_path):
    from test_pgen_native_reader import write_pgen_mixed
    from torchgwas.pgen_reader import read_header
    path=tmp_path/'mixed.pgen'
    categories=np.zeros((7,128),dtype=np.uint8)
    categories[1:,3]=1
    write_pgen_mixed(path,categories,[0,2,2,2,0,2,2])
    header=read_header(path)
    work=census(path,chunk_markers=2,variant_range=(1,7))
    assert work['record_form_counts']=={2:5,0:1}
    assert work['ld_records_at_chunk_starts']==3  # absolute1,3,5
    assert work['additional_base_record_bytes_if_every_chunk_restarts']==2*int(header.record_lengths[0])+int(header.record_lengths[4])
    assert work['record_payload_bytes']==int(header.record_lengths[1:].sum())
    whole=census(path,chunk_markers=2)
    prefix=census(path,chunk_markers=2,variant_range=(0,1))
    for key in ['record_payload_bytes','sample_delta_integer_count','total_difflist_groups']:
        assert whole[key]==prefix[key]+work[key]


@pytest.mark.parametrize('span', [None, (1, 10), (2, 9)])
def test_chunk_census_matches_independent_range_reads_and_conserves_work(tmp_path, span):
    from test_pgen_native_reader import write_pgen_mixed
    from torchgwas.decoder_work import decoder_chunk_work, decoder_work
    from collections import Counter
    path = tmp_path/'varied.pgen'
    n = 257
    forms = [0, 1, 4, 0, 1, 4, 4, 0, 1, 0, 4]
    categories = np.zeros((len(forms), n), dtype=np.uint8)
    for i, form in enumerate(forms):
        if form == 0:
            categories[i] = np.arange(n, dtype=np.uint16)%4
        elif form == 1:
            categories[i, ::2] = 2
            categories[i, 3:6] = 3
            categories[i, -1] = 2
        else:
            categories[i, i:i+2] = 1
    write_pgen_mixed(path, categories, forms)
    complete = census(path, 3, span, include_chunks=True)
    aggregate = {key:value for key,value in complete.items() if key not in {'chunks','chunk_census_scope'}}
    assert aggregate == census(path, 3, span)
    for part in complete['chunks']:
        assert part == census(path, 3, tuple(part['variant_range']))
    parts = decoder_chunk_work(complete, 3)
    units = Counter()
    for part in parts:
        units.update(part['decoder']['source_units'])
    assert dict(units) == decoder_work(complete, 'torch_native_int8', restart_ld_bases=True)['source_units']
    assert sum(part['census']['record_payload_bytes'] for part in parts) == complete['record_payload_bytes']
    assert sum(part['census']['native_onebit_tail_high_count'] for part in parts) == complete['native_onebit_tail_high_count']


def test_chunk_census_retains_ld_restart_obligations(tmp_path):
    from test_pgen_native_reader import write_pgen_mixed
    path = tmp_path/'ld.pgen'
    categories = np.zeros((7, 128), dtype=np.uint8)
    categories[1:, 3] = 1
    write_pgen_mixed(path, categories, [0, 2, 2, 2, 0, 2, 2])
    report = census(path, 2, (1, 7), include_chunks=True)
    assert sum(part['ld_records_at_chunk_starts'] for part in report['chunks']) == 3
    for part in report['chunks']:
        assert part == census(path, 2, tuple(part['variant_range']))


@pytest.mark.parametrize('field', ['payload', 'unit', 'range', 'file', 'geometry', 'missing_chunk'])
def test_chunk_decoder_rejects_inconsistent_supplied_work(tmp_path, field):
    from torchgwas.decoder_work import decoder_chunk_work
    path = tmp_path/'plain.pgen'
    write_pgen(path, np.zeros((11, 128), dtype=np.uint8))
    report = census(path, 3, include_chunks=True)
    if field == 'payload':report['chunks'][0]['record_payload_bytes'] += 1
    elif field == 'unit':report['chunks'][0]['record_form_counts'] = {4:3}
    elif field == 'range':report['chunks'][0]['variant_range'] = [1, 4]
    elif field == 'file':report['chunks'][0]['file_markers'] += 1
    elif field == 'geometry':report['chunk_markers'] += 1
    else:report['chunks'].pop()
    with pytest.raises(ValueError):decoder_chunk_work(report, 3)