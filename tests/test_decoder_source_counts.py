from torchgwas.decoder_work import decoder_work


def fixture(n=9):
    return dict(samples=n,markers=7,record_form_counts={str(f):1 for f in [0,1,2,3,4,6,7]},
                total_varint_lengths={'1':6},source_work=dict(varint_bytes=6,difflist_entries=0),
                total_difflist_groups=0,ld_records_at_chunk_starts=0,native_onebit_tail_high_count=1)


def test_native_base_refresh_includes_non_ld_records_only():
    work=decoder_work(fixture(),'torch_native_int8')
    # Seven logical three-byte rows: three input/base-read copies plus five
    # refreshes, since forms2/3 preserve the existing LD base.
    assert work['source_units']['copy_packed_byte']==24
    assert work['native_ld_base_update_bytes']==15
    assert work['logical_final_packed_bytes']==21


def test_existing_fill_and_onebit_tail_are_counted_once():
    work=decoder_work(fixture(),'torch_native_int8')
    assert work['source_units']['fill_packed_byte']==10  # three fills and one tail byte
    assert work['source_units']['set_category']==1
    assert work['source_units']['native_onebit_tail_bit_test']==1
    assert work['source_units']['expand_int8_tail_sample']==7


def test_native_refresh_does_not_change_other_reader_accounting():
    work=decoder_work(fixture(),'pgenlib_sse2')
    assert work['source_units']['copy_packed_byte']==9
    assert work['native_ld_base_update_bytes']==0


def test_form_zero_still_refreshes_base_without_any_ld_reference():
    census=dict(samples=8,markers=3,record_form_counts={'0':3},ld_records_at_chunk_starts=0)
    work=decoder_work(census,'torch_native_int8')
    assert work['source_units']['copy_packed_byte']==12
    assert work['native_ld_base_update_bytes']==6
