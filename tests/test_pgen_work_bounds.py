"""Header intervals enclose independent payload censuses and native outputs."""
import builtins
from collections import Counter
from copy import deepcopy
import json
import struct
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.decoder_work import decoder_work, native_read_layout, price_identified_units
from torchgwas.decoder_work import _native_decoder_memory_bytes
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.pgen_native_reader import NativePgenReader
from torchgwas.pgen_reader import read_header, _bytes_per_sample_id, pack_genovec
from torchgwas.pgen_work_bounds import (KIND,SCHEDULE_KIND,PgenHeaderWork,difflist_bounds,
    price_header_work,native_decoder_service_bounds,paired_schedule_source_difference,
    native_schedule_source_floor)
from torchgwas.pgen_work_census import census
from test_pgen_native_reader import _encode_difflist, _encode_onebit, _uleb128


def write_records(path, n, forms, records):
    m = len(records); data_start = 20 + 5*m
    path.write_bytes(b'\x6c\x1b\x10' + struct.pack('<II', m, n) + b'\x07'
        + struct.pack('<Q', data_start) + bytes(forms)
        + b''.join(struct.pack('<I', len(record)) for record in records) + b''.join(records))


def fixture(path, n):
    forms = [0, 2, 3, 1, 2, 4, 3, 6, 2, 7, 3, 1, 2, 0, 3]
    calls = []; records = []; base = None
    rng = np.random.default_rng(915)
    for i, form in enumerate(forms):
        if form in (4, 6, 7):
            background = {4: 0, 6: 2, 7: 3}[form]
            row = np.full(n, background, np.uint8)
            row[np.arange(i % 3, n, 17)] = (background+1) % 4
            record = _encode_difflist([(int(s), int(row[s])) for s in np.flatnonzero(row != background)],
                                     _bytes_per_sample_id(n))
        elif form == 1:
            row = rng.integers(0, 3, n, dtype=np.uint8)
            record = _encode_onebit(row, n, _bytes_per_sample_id(n), 0, 1)
        elif form == 0:
            row = rng.integers(0, 4, n, dtype=np.uint8)
            record = pack_genovec(row, n).tobytes()
        else:
            row = base.copy(); row[i % n] = (row[i % n]+1) % 4
            wanted = row.copy()
            if form == 3:
                wanted[row == 0] = 2; wanted[row == 2] = 0
            record = _encode_difflist([(int(s), int(wanted[s])) for s in np.flatnonzero(wanted != base)],
                                     _bytes_per_sample_id(n))
        records.append(record); calls.append(row)
        if form not in (2, 3):
            base = row
    write_records(path, n, forms, records)
    return np.stack(calls)


def assert_encloses(bounds, exact):
    work = decoder_work(exact, 'torch_native_int8', restart_ld_bases=True)
    for name in set(bounds['source_units']) | set(work['source_units']):
        lo, hi = bounds['source_units'].get(name, [0, 0])
        assert lo <= work['source_units'].get(name, 0) <= hi, (name, lo, hi, work['source_units'])
    for key in ('native_ld_base_update_bytes', 'logical_final_packed_bytes', 'logical_final_int8_bytes'):
        assert bounds[key] == work[key]
    layout = native_read_layout(exact)
    assert bounds['read_bytes'] == layout['read_bytes']
    assert bounds['decode_input_bytes'] == layout['decode_input_bytes']
    assert bounds['native_ld_replay_packed_bytes'] == work.get('native_ld_replay_packed_bytes', 0)
    return work


@pytest.mark.parametrize('n', [1, 3, 4, 7, 8, 9, 63, 64, 65, 128, 129, 256, 257, 2049])
def test_all_forms_every_chunk_start_and_tail_enclose_exact_work(tmp_path, n):
    path = tmp_path/'input.pgen'; fixture(path, n); index = PgenHeaderWork(path)
    uncached = PgenHeaderWork(path, max_cached_signatures=0)
    for start in range(15):
        for size in (1, 3, 8, 15):
            stop = min(15, start+size)
            bounds = index.bounds(start, stop)
            assert bounds == uncached.bounds(start, stop)
            exact = census(path, size, variant_range=(start, stop))
            work = assert_encloses(bounds, exact)
            prices = {name: (i+1)*1e-8 for i, name in enumerate(sorted(bounds['source_units']))}
            priced = price_header_work(bounds, prices)
            exact_price = price_identified_units(work, prices)['identified_cpu_seconds']
            lo, hi = priced['cpu_seconds']
            assert lo <= exact_price+1e-12 and exact_price <= hi+1e-12


@pytest.mark.parametrize('n',[1,65,129,2049])
@pytest.mark.parametrize('span',[(0,15),(1,15),(2,11)])
@pytest.mark.parametrize('chunk',[1,2,3,8,16])
def test_whole_header_schedule_encloses_every_ld_restart(tmp_path,n,span,chunk):
    path=tmp_path/'input.pgen';fixture(path,n)
    header=PgenHeaderWork(path)
    bound=header.schedule_bounds(*span,chunk)
    exact=census(path,chunk,variant_range=span)
    work=assert_encloses(bound,exact)
    assert bound['kind']==SCHEDULE_KIND
    assert bound['chunk_count']==(span[1]-span[0]+chunk-1)//chunk
    assert bound['ld_replay_count']==exact['ld_records_at_chunk_starts']
    assert bound['structural_work']['payload_bytes_read']==0
    prices={name:(i+1)*1e-8 for i,name in enumerate(sorted(bound['source_units']))}
    lo,hi=price_header_work(bound,prices)['cpu_seconds']
    exact_price=price_identified_units(work,prices)['identified_cpu_seconds']
    assert lo<=exact_price+1e-12 and exact_price<=hi+1e-12
    with pytest.raises(ValueError,match='intervals'):
        decoder_work(bound,'torch_native_int8')


def test_whole_header_schedule_budgets_and_unpriced_units(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129)
    header=PgenHeaderWork(path)
    for kwargs in ({'max_records':14},{'max_chunks':2},{'max_signatures':1}):
        with pytest.raises(ValueError,match='budget'):
            header.schedule_bounds(0,15,3,**kwargs)
    bound=header.schedule_bounds(0,15,3)
    possible=next(name for name,interval in bound['source_units'].items() if interval[1])
    prices={name:1e-8 for name in bound['source_units'] if name!=possible}
    priced=price_header_work(bound,prices)
    assert priced['cpu_seconds'] is None and possible in priced['unpriced_source_units']
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='changed'):
        header.schedule_bounds(0,15,3)


@pytest.mark.parametrize('chunk',[1,3,8])
@pytest.mark.parametrize('buffered_cpu',[False,True])
def test_whole_schedule_source_floor_conserves_exact_stage_work(tmp_path,chunk,buffered_cpu):
    path=tmp_path/'input.pgen';fixture(path,129)
    span=(1,15);schedule=PgenHeaderWork(path).schedule_bounds(*span,chunk)
    prices={name:(i+1)*1e-8 for i,name in enumerate(sorted(schedule['source_units']))}
    profile=dict(decode_units=prices,cpu_fraction=.5,depth=3,decode_workers=2,
        cpu_available_cores=4.,shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7)
    if buffered_cpu:
        profile['input_read_cpu_prices']=dict(cpu_seconds_per_byte=1e-8,
                                               cpu_seconds_per_call=1e-5)
    caps=dict(cpu=2.,dram=5e7,input=5e6)
    report=native_schedule_source_floor(schedule,profile,caps)
    assert report['samples']==129
    exact=census(path,chunk,variant_range=span)
    decoder=decoder_work(exact,'torch_native_int8',restart_ld_bases=True)
    exact_cpu=price_identified_units(decoder,prices)['identified_cpu_seconds']
    read_cpu=(schedule['read_bytes']*1e-8+schedule['chunk_count']*1e-5
              if buffered_cpu else 0.)
    low,high=report['resource_work']['cpu_seconds']
    assert low-1e-12<=exact_cpu+read_cpu<=high+1e-12
    assert report['resource_work']['read_cpu_seconds']==pytest.approx(read_cpu)
    assert report['resource_work']['input_bytes']==native_read_layout(exact)['read_bytes']
    expected_dram=_native_decoder_memory_bytes(129,span[1]-span[0],
        schedule['read_bytes'],schedule['decode_input_bytes'],
        schedule['native_ld_base_update_bytes'],schedule['native_ld_replay_packed_bytes'],
        separate_buffered_read=buffered_cpu)
    if buffered_cpu:expected_dram+=2*schedule['read_bytes']
    assert report['resource_work']['dram_bytes']==expected_dram
    assert report['chunk_count']==(span[1]-span[0]+chunk-1)//chunk
    assert report['source_stage_floor_seconds'][0]<=report['source_stage_floor_seconds'][1]
    assert all(a<=b for a,b in zip(report['resource_floor_seconds'],
                                    report['source_stage_floor_seconds']))
    assert not report['prediction_complete'] and not report['selection_validated']


def test_whole_schedule_source_floor_rejects_missing_price_capacity_and_staleness(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129)
    schedule=PgenHeaderWork(path).schedule_bounds(0,15,3)
    prices={name:1e-8 for name in schedule['source_units']}
    profile=dict(decode_units=prices,cpu_fraction=.5,depth=3,decode_workers=2,
        cpu_available_cores=4.,shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7)
    caps=dict(cpu=2.,dram=5e7,input=5e6)
    missing=next(name for name,interval in schedule['source_units'].items() if interval[1])
    with pytest.raises(ValueError,match='Unpriced'):
        native_schedule_source_floor(schedule,dict(profile,decode_units={
            key:value for key,value in prices.items() if key!=missing}),caps)
    with pytest.raises(ValueError,match='exceeds'):
        native_schedule_source_floor(schedule,profile,dict(caps,cpu=5.))
    altered=deepcopy(schedule);altered['read_bytes']+=1
    with pytest.raises(ValueError,match='conserve'):
        native_schedule_source_floor(altered,profile,caps)
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='changed'):
        native_schedule_source_floor(schedule,profile,caps)


@pytest.mark.parametrize('buffered_cpu',[False,True])
def test_whole_schedule_source_floor_is_below_finite_worker_graph(tmp_path,buffered_cpu):
    path=tmp_path/'plain.pgen';n=129;m=9
    records=[pack_genovec(np.full(n,i%3,dtype=np.uint8),n).tobytes()
             for i in range(m)]
    write_records(path,n,[0]*m,records)
    header=PgenHeaderWork(path);schedule=header.schedule_bounds(0,m,3)
    prices={name:1e-8 for name in schedule['source_units']}
    profile=dict(decode_units=prices,cpu_fraction=.5,depth=2,decode_workers=2,
        cpu_available_cores=2.,shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7)
    if buffered_cpu:
        profile['input_read_cpu_prices']=dict(cpu_seconds_per_byte=1e-8,
                                               cpu_seconds_per_call=1e-5)
    caps=dict(cpu=1.,dram=5e7,input=5e6)
    floor=native_schedule_source_floor(schedule,profile,caps)
    graph=ExecutionGraph();graph.capacities=caps
    for i,row in enumerate(header.window(0,m,3)['chunks']):
        service=native_decoder_service_bounds(row,profile)
        read=service['read'];read_name=f'read:{i}';decode_name=f'decode:{i}'
        prior=[f'decode:{i-2}'] if i>=2 else []
        graph.add(read_name,read['seconds'],prior,read['resources'])
        cpu=service['decode_resource_work']['cpu'][0]
        dram=service['decode_resource_work']['dram'][0]
        seconds=service['decode_seconds'][0]
        graph.add(decode_name,seconds,[read_name],
                  {'cpu':cpu/seconds,'dram':dram/seconds})
    assert floor['source_stage_floor_seconds'][0]<=graph.solve()['seconds']+1e-12


@pytest.mark.parametrize('n',[1,65,129,2049])
@pytest.mark.parametrize('span',[(0,15),(1,15),(2,11)])
@pytest.mark.parametrize('chunk',[1,2,3,8,16])
def test_vectorized_window_matches_independent_chunk_bounds(tmp_path,n,span,chunk):
    path=tmp_path/'input.pgen';fixture(path,n)
    header=PgenHeaderWork(path,max_cached_signatures=1)
    ordinary=header.window(*span,chunk,max_chunks=16,max_signatures=32)
    bulk=header.vectorized_window(*span,chunk,max_chunks=16,max_signatures=32)
    assert bulk['chunk_ranges']==ordinary['chunk_ranges']
    assert bulk['chunks']==ordinary['chunks']
    assert bulk['structural_work']==ordinary['structural_work']
    schedule=header.schedule_bounds(*span,chunk)
    assert sum(row['read_bytes'] for row in bulk['chunks'])==schedule['read_bytes']
    assert sum(row['decode_input_bytes'] for row in bulk['chunks'])==schedule['decode_input_bytes']


def test_vectorized_window_limits_and_stale_source(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129)
    header=PgenHeaderWork(path)
    for kwargs in ({'max_records':14},{'max_chunks':2},{'max_signatures':1},
                   {'max_global_signatures':1},{'max_signature_pairs':1}):
        with pytest.raises(ValueError,match='budget'):
            header.vectorized_window(0,15,3,**kwargs)
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='changed'):
        header.vectorized_window(0,15,3)


@pytest.mark.parametrize('span',[(0,15),(1,15),(2,11)])
@pytest.mark.parametrize('sizes',[(1,3),(2,8),(3,16)])
def test_paired_schedule_cancels_common_source_and_encloses_exact_delta(tmp_path,span,sizes):
    path=tmp_path/'input.pgen';fixture(path,129)
    header=PgenHeaderWork(path)
    first,second=(header.schedule_bounds(*span,size) for size in sizes)
    prices={name:(i+1)*1e-8 for i,name in enumerate(sorted(set(first['source_units'])|set(second['source_units'])))}
    delta=paired_schedule_source_difference(first,second,prices)
    same=paired_schedule_source_difference(first,first,prices)
    assert same['source_unit_delta_intervals']=={}
    assert same['candidate_minus_baseline_decoder_cpu_seconds']==[0.,0.]
    assert same['decoder_cpu_delta_lower_bound']==same['decoder_cpu_delta_upper_bound']==0.
    assert same['read_bytes_delta']==0
    exact=[]
    for size in sizes:
        census_row=census(path,size,variant_range=span)
        exact.append((decoder_work(census_row,'torch_native_int8',restart_ld_bases=True),
                      native_read_layout(census_row)))
    cpu=[price_identified_units(row[0],prices)['identified_cpu_seconds'] for row in exact]
    assert delta['candidate_minus_baseline_decoder_cpu_seconds'][0]-1e-12<=cpu[1]-cpu[0]
    assert cpu[1]-cpu[0]<=delta['candidate_minus_baseline_decoder_cpu_seconds'][1]+1e-12
    assert delta['read_bytes_delta']==exact[1][1]['read_bytes']-exact[0][1]['read_bytes']
    assert delta['decode_input_bytes_delta']==exact[1][1]['decode_input_bytes']-exact[0][1]['decode_input_bytes']
    assert delta['chunk_count_delta']==second['chunk_count']-first['chunk_count']
    if delta['source_unit_delta_intervals']:
        missing=next(iter(delta['source_unit_delta_intervals']))
        partial=paired_schedule_source_difference(first,second,{k:v for k,v in prices.items() if k!=missing})
        assert partial['candidate_minus_baseline_decoder_cpu_seconds'] is None
        assert missing in partial['unpriced_source_units']
        assert partial['decoder_cpu_delta_upper_bound'] is not None
        assert cpu[1]-cpu[0]<=partial['decoder_cpu_delta_upper_bound']+1e-12
        reverse=paired_schedule_source_difference(second,first,{k:v for k,v in prices.items() if k!=missing})
        assert reverse['decoder_cpu_delta_lower_bound'] is not None
        assert cpu[0]-cpu[1]>=reverse['decoder_cpu_delta_lower_bound']-1e-12
    altered=deepcopy(second);altered['input_identity']['bytes']+=1
    with pytest.raises(ValueError):paired_schedule_source_difference(first,altered,prices)
    altered=deepcopy(second);altered['source_units'].clear()
    with pytest.raises(ValueError,match='conserve'):
        paired_schedule_source_difference(first,altered,prices)


def test_paired_schedule_rejects_inconsistent_shared_restart(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129)
    header=PgenHeaderWork(path)
    fine=header.schedule_bounds(0,15,1)
    coarse=header.schedule_bounds(0,15,2)
    assert coarse['additional_replay_work']['entries']
    altered=deepcopy(coarse)
    altered['additional_replay_work']['entries'][0][2]+=1
    with pytest.raises(ValueError,match='conserve'):
        paired_schedule_source_difference(fine,altered,{})


@pytest.mark.parametrize('limit', [0, 1, 3, 1024])
def test_signature_cache_is_bounded_and_does_not_expose_mutable_rows(tmp_path, limit):
    path = tmp_path/'input.pgen'; fixture(path, 129)
    index = PgenHeaderWork(path, max_cached_signatures=limit,max_cached_bounds=0)
    expected = PgenHeaderWork(path, max_cached_signatures=0,max_cached_bounds=0).window(0, 15, 3)
    for _ in range(2):
        result = index.window(0, 15, 3)
        assert result == expected
        for row in result['chunks']:
            for value in row['source_units'].values(): value[:] = [-100, -100]
            row['input_identity'].clear()
            row['varint_constraints'].clear()
        assert index.cache_info()['currsize'] <= limit
    assert index.window(0, 15, 3) == expected
    info = index.cache_info()
    assert info['maxsize'] == limit and info['misses'] > 0
    if limit == 0: assert info['hits'] == 0
    if limit == 1024: assert info['hits'] > 0


def test_overlapping_window_bounds_cache_is_bounded_and_file_bound(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129)
    index=PgenHeaderWork(path,max_cached_bounds=2)
    expected=PgenHeaderWork(path,max_cached_bounds=0).bounds(0,3)
    first=index.bounds(0,3)
    first['source_units'].clear();first['input_identity'].clear()
    assert index.bounds(0,3)==expected
    assert index.bounds_cache_info()==dict(hits=1,misses=1,entries=1,max_entries=2)
    index.bounds(3,6);index.bounds(6,9)
    assert index.bounds_cache_info()['entries']==2
    assert index.bounds(0,3)==expected
    assert index.bounds_cache_info()['misses']==4
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='PGEN input changed'):
        index.bounds(0,3)


def test_prepared_header_reuses_bound_locators_without_reopening_index(tmp_path):
    import torchgwas.pgen_work_bounds as module
    from torchgwas.analytical_plan_cache import input_identity
    path=tmp_path/'input.pgen';fixture(path,129)
    fresh=PgenHeaderWork(path)
    assert fresh._bases is None
    expected=fresh.window(0,15,3)
    prepared=(input_identity(path),read_header(path))
    with patch.object(module,'read_header',side_effect=AssertionError('index reread')):
        reused=PgenHeaderWork(path,_prepared_header=prepared)
    assert reused.window(0,15,3)==expected
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='Prepared PGEN header differs'):
        PgenHeaderWork(path,_prepared_header=prepared)


def test_admitted_index_skips_repeated_whole_file_validation(tmp_path):
    import torchgwas.pgen_work_bounds as module
    from torchgwas.analytical_plan_cache import input_identity
    from torchgwas.pgen_memory_layout import memory_layout
    path=tmp_path/'input.pgen';fixture(path,129)
    expected=PgenHeaderWork(path).window(0,15,3)
    identity=input_identity(path);header=read_header(path);bases=[]
    memory_layout(path,3,header=header,_validated_bases_receiver=bases.append)
    assert len(bases)==1
    compact_bases=bases[0].astype('uint32');compact_bases.flags.writeable=False
    prepared=(identity,header);certificate=(identity,header,compact_bases)
    with (patch.object(module,'read_header',side_effect=AssertionError('index reread')),
          patch.object(module.np,'unique',side_effect=AssertionError('whole-index validation'))):
        reused=PgenHeaderWork(path,_prepared_header=prepared,_prepared_index=certificate)
    assert reused.window(0,15,3)==expected
    original_search=np.searchsorted
    def checked_search(array,value,**kwargs):
        if array is certificate[2]:
            assert isinstance(value,np.uint32)
        return original_search(array,value,**kwargs)
    with patch.object(module.np,'searchsorted',side_effect=checked_search):
        assert reused.bounds(1,2)['ld_replay']['base_variant']==0
    with pytest.raises(ValueError,match='Prepared PGEN index differs'):
        PgenHeaderWork(path,_prepared_header=prepared,
            _prepared_index=(identity,read_header(path),compact_bases))
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='Prepared PGEN header differs'):
        PgenHeaderWork(path,_prepared_header=prepared,_prepared_index=certificate)


def test_signature_cache_distinguishes_ld_base_replay_from_expansion(tmp_path):
    import torchgwas.pgen_work_bounds as module
    path = tmp_path/'input.pgen'; fixture(path, 129)
    with patch.object(module, '_record_units', wraps=module._record_units) as build:
        index = PgenHeaderWork(path)
        first = index.bounds(0, 1)
        replay = index.bounds(1, 2)
        assert first['logical_final_int8_bytes'] == replay['logical_final_int8_bytes'] == 129
        modes = [call.kwargs['expand'] for call in build.call_args_list if call.args[1] == 0]
        assert modes == [True, False]
        calls = build.call_count
        assert index.bounds(0, 1) == first and index.bounds(1, 2) == replay
        assert build.call_count == calls


@pytest.mark.parametrize('limit', [-1, True, 1.5, None])
def test_invalid_signature_cache_limit(tmp_path, limit):
    with pytest.raises(ValueError, match='max_cached_signatures'):
        PgenHeaderWork(tmp_path/'not_read.pgen', max_cached_signatures=limit)
    with pytest.raises(ValueError, match='max_cached_bounds'):
        PgenHeaderWork(tmp_path/'not_read.pgen', max_cached_bounds=limit)
    with pytest.raises(ValueError, match='max_cached_bounds'):
        PgenHeaderWork(tmp_path/'not_read.pgen', max_cached_bounds=257)


@pytest.mark.parametrize('n', [9, 129, 257])
def test_fixture_decodes_natively_with_ld_restart(tmp_path, n):
    path = tmp_path/'input.pgen'; calls = fixture(path, n)
    # Check real encoded records, including constant-category forms 6 and 7.
    with NativePgenReader(path) as reader:
        for start in (0, 2, 4, 6, 8, 10, 12, 14):
            stop = min(start+3, 15)
            got = np.empty((stop-start, n), dtype=np.int8)
            reader.read_range(start, stop, got)
            expected = calls[start:stop].astype(np.int8)
            expected[expected == 3] = -9
            np.testing.assert_array_equal(got, expected)


@pytest.mark.parametrize('n', [1, 64, 65, 127, 128, 256, 257, 65536, 65537, 2**24+1, 2**32-1])
def test_length_inequalities_cover_varint_and_group_boundaries(n):
    rng = np.random.default_rng(n)
    for e in sorted({0, 1, n, *[min(n, v) for v in (2, 63, 64, 65, 127, 128, 129, 16383, 16384)]}):
        g = (e+63)//64; d = e-g; width = _bytes_per_sample_id(n)
        for header_len in range(max(1, (e.bit_length()+6)//7), 6):
            for delta_len in (1, 2, 3, 4, 5):
                length = header_len + g*width + max(g-1, 0) + (e+3)//4 + delta_len*d
                result = difflist_bounds(n, length)
                assert result['entries'][0] <= e <= result['entries'][1]
                counts = Counter({header_len: 1}); counts[delta_len] += d
                for k in range(1, 6):
                    lo, hi = result['source_units'][f'uleb{k}']
                    assert lo <= counts[k] <= hi
        # Mixed integer lengths rather than only one delta length per record.
        if d <= 20000:
            lengths = rng.integers(1, 6, size=d).tolist() + [max(1, (e.bit_length()+6)//7)]
            result = difflist_bounds(n, sum(lengths)+g*width+max(g-1,0)+(e+3)//4)
            for k, count in Counter(lengths).items():
                assert result['source_units'][f'uleb{k}'][0] <= count <= result['source_units'][f'uleb{k}'][1]


@pytest.mark.parametrize('length', [1, 2, 3, 4, 5])
def test_actual_sparse_record_with_every_integer_width(tmp_path, length):
    n = max(2, (1 << (7*(length-1)))+1)
    record = _encode_difflist([(0, 1), (n-1, 2)], _bytes_per_sample_id(n))
    path = tmp_path/'sparse.pgen'; write_records(path, n, [4], [record])
    work = assert_encloses(PgenHeaderWork(path).bounds(0, 1), census(path, 1))
    assert work['source_units'][f'uleb{length}'] >= 1


def test_noncanonical_five_byte_count_is_not_assumed_one_byte(tmp_path):
    path = tmp_path/'overlong.pgen'
    write_records(path, 3, [4], [b'\x80\x80\x80\x80\x00'])
    work = assert_encloses(PgenHeaderWork(path).bounds(0, 1), census(path, 1))
    assert work['source_units']['uleb5'] == 1


def test_no_genotype_payload_is_read_or_mapped(tmp_path):
    path = tmp_path/'input.pgen'; fixture(path, 129)
    first = int(read_header(path).record_offsets[0]); original = builtins.open; reads = []
    class Checked:
        def __init__(self, file): self.file = file
        def __enter__(self): self.file.__enter__(); return self
        def __exit__(self, *args): return self.file.__exit__(*args)
        def __getattr__(self, key): return getattr(self.file, key)
        def read(self, size=-1):
            lo = self.file.tell()
            assert size >= 0 and lo+size <= first
            value = self.file.read(size); reads.append((lo, len(value))); return value
    with patch('builtins.open', side_effect=lambda *a, **kw: Checked(original(*a, **kw))), \
         patch('mmap.mmap', side_effect=AssertionError('Payload mapping forbidden')):
        index = PgenHeaderWork(path)
        for start in range(15): index.bounds(start, min(15, start+4))
        index.schedule_bounds(0,15,3)
    assert reads and max(lo+size for lo, size in reads) <= first


def test_structural_bounds_cache_preserves_original_record_and_dependencies(tmp_path):
    path = tmp_path/'input.pgen'; fixture(path, 129); index = PgenHeaderWork(path)
    bounds = index.bounds(2, 7)
    deps = dict(input=bounds['input_identity'], schema=KIND, decoder_source='test-fixture', span=[2, 7])
    cache = CalibrationParameterCache(tmp_path/'cache')
    with patch('torchgwas.calibration_cache.time.time', return_value=100):
        stored = cache.store('source_work', 'header-bounds', bounds, dependencies=deps,
                             provenance={'method': 'header_only'})
    from pathlib import Path
    artifact = Path(stored['path']); original = artifact.read_bytes()
    with patch('torchgwas.calibration_cache.time.time', return_value=1e9):
        hit = cache.lookup('source_work', 'header-bounds', dependencies=deps)
        assert hit['hit'] and hit['age_seconds'] == 1e9-100
        assert hit['record']['value']['source_units'] == bounds['source_units']
        assert artifact.read_bytes() == original
        changed = deepcopy(deps); changed['decoder_source'] = 'changed'
        assert not cache.lookup('source_work', 'header-bounds', dependencies=changed)['hit']
    # A filesystem identity change makes this already-open header unusable.
    with path.open('ab') as stream: stream.write(b'\0')
    with pytest.raises(ValueError, match='changed'): index.bounds(2, 7)
    assert artifact.read_bytes() == original


def test_budgets_and_unknown_prices_fail_explicitly(tmp_path):
    path = tmp_path/'input.pgen'; fixture(path, 129); index = PgenHeaderWork(path)
    for kw in ({'max_records': 4}, {'max_signatures': 1}):
        with pytest.raises(ValueError, match='budget'): index.bounds(0, 15, **kw)
    bounds = index.bounds(3, 4)
    prices = {key: 1e-8 for key in bounds['source_units']}
    maybe = next(key for key, (lo, hi) in bounds['source_units'].items() if lo == 0 < hi)
    del prices[maybe]
    result = price_header_work(bounds, prices)
    assert result['cpu_seconds'] is None and maybe in result['unpriced_source_units']
    with pytest.raises(ValueError, match='intervals'): decoder_work(bounds, 'torch_native_int8')
    bad = deepcopy(bounds); bad['source_units']['set_category'] = [2, 1]
    with pytest.raises(ValueError, match='interval'): price_header_work(bad, prices)


@pytest.mark.parametrize('samples,length', [(0, 1), (True, 1), (2**32, 1), (3, 0), (1, 100)])
def test_invalid_diff_extents(samples, length):
    with pytest.raises(ValueError): difflist_bounds(samples, length)


def test_plain_records_collapse_to_exact_counts(tmp_path):
    path = tmp_path/'plain.pgen'; n = 9
    write_records(path, n, [0]*3, [b'\0'*3]*3)
    bounds = PgenHeaderWork(path).bounds(0, 3)
    work = assert_encloses(bounds, census(path, 3))
    assert bounds['source_units'] == {key: [value, value] for key, value in work['source_units'].items()}
    # Intervals survive serialization without relying on integer dict keys.
    decoded = json.loads(json.dumps(bounds))
    assert price_header_work(decoded, {key: 1 for key in bounds['source_units']})['cpu_seconds'][0] == sum(work['source_units'].values())


@pytest.mark.parametrize('cpu_fraction', [.25, 1.])
@pytest.mark.parametrize('dram', [1e5, 1e12])
@pytest.mark.parametrize('buffered', [False, True])
def test_service_bounds_enclose_existing_scan_resource_accounting(tmp_path, cpu_fraction, dram, buffered):
    from test_mechanistic_shapes import fixture as profile_fixture, component
    from torchgwas.mechanistic_torch import torch_scan_work
    path = tmp_path/'input.pgen'; fixture(path, 129); index = PgenHeaderWork(path)
    data, profile = profile_fixture(k=2, c=8)
    data.update(samples=129, markers=15, encoded=census(path, 4, include_chunks=True))
    profile.update(chunk_markers=4, cpu_fraction=cpu_fraction, shared_dram_bytes_per_second=dram,
                   kernel_geometry=[dict(N=129, B=b, K=2, C=8, kernels=[]) for b in (4, 3)])
    rows = [index.bounds(start, min(start+4, 15)) for start in range(0, 15, 4)]
    names = sorted({name for row in rows for name in row['source_units']})
    profile['decode_units'] = {name: (i+1)*1e-8 for i, name in enumerate(names)}
    if buffered:
        profile['input_read_cpu_prices'] = dict(cpu_seconds_per_byte=1e-9, cpu_seconds_per_call=1e-6)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service', side_effect=component):
        scan = torch_scan_work(data, profile)
    for row, block in zip(rows, scan['blocks']):
        bounded = native_decoder_service_bounds(row, profile)
        lo, hi = bounded['decode_seconds']
        assert lo-1e-12 <= block['decode_seconds'] <= hi+1e-12
        for name, (lo, hi) in bounded['decode_resource_work'].items():
            amount = block['decode_seconds']*block['decode_resources'][name]
            assert lo-1e-9 <= amount <= hi+1e-9
        assert bounded['read']['seconds'] == block['decode_read_seconds']
        assert bounded['read']['resources'] == block['read_resources']


def test_joint_integer_count_tightens_marginal_histogram_prices(tmp_path):
    path = tmp_path/'input.pgen'; fixture(path, 129)
    bounds = PgenHeaderWork(path).bounds(3, 4)
    prices = {key: 1 for key in bounds['source_units']}
    constrained = price_header_work(bounds, prices)['cpu_seconds']
    marginal = deepcopy(bounds); marginal.pop('varint_constraints')
    separate = price_header_work(marginal, prices)['cpu_seconds']
    assert constrained[0] >= separate[0] and constrained[1] < separate[1]
    assert_encloses(bounds, census(path, 1, variant_range=(3, 4)))
