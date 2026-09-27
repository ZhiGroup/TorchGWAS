"""Price the exact source ranges retained by a completed adaptive GPU run.

This is a structural replay audit. Synthetic prices exercise the graph; no
association duration or cached stage span is converted to a service rate.
"""
import argparse
import copy
import json
from pathlib import Path

from torchgwas.decoder_work import decoder_work,decoder_chunk_work,native_read_layout
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.jagwas_candidate import jagwas_candidate_runtime
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.pgen_work_census import census,scheduled_census
from test_jagwas_actual_candidate import actual_candidate,writer_prices
from test_jagwas_candidate import preparation


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--execution',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();root=Path(args.execution);out=Path(args.out)
    out.mkdir(parents=True,exist_ok=False)
    report_path=root/'report.json';previous=json.loads(report_path.read_text())
    source=source_identity()
    for name,digest in previous['inputs'].items():assert sha256_file(name)==digest
    path=Path(next(name for name in previous['inputs'] if name.endswith('.pgen')))
    fine=census(path,128,include_chunks=True);choice=actual_candidate(path,512,2)
    fixed=copy.deepcopy(choice);rows=[];exact_reads=regular_reads=0
    for tile,old in zip(choice['tiles'],fixed['tiles']):
        lo,hi=tile['variant_range']
        ranges=[row for row in previous['source_ranges'] if lo<=row[0]<row[1]<=hi]
        actual=scheduled_census(fine,512,ranges)
        # Reopen and parse each actual range independently of the regrouping.
        expected=[census(path,512,row,include_chunks=True)['chunks'][0] for row in ranges]
        assert actual['chunks']==expected
        decoder_chunk_work(actual,512)
        tile['data']['encoded']=actual
        needed=set(decoder_work(actual,'torch_native_int8',restart_ld_bases=True)['source_units'])
        tile['profile']['decode_units']={key:1e-9 for key in needed}
        old['profile']['decode_units']=copy.deepcopy(tile['profile']['decode_units'])
        work=torch_scan_work(tile['data'],tile['profile'])
        assert [b['variant_range'] for b in work['blocks']]==ranges
        read=native_read_layout(actual)['read_bytes']
        regular=native_read_layout(old['data']['encoded'])['read_bytes']
        exact_reads+=read;regular_reads+=regular
        rows.append(dict(device=tile['device'],variant_range=[lo,hi],chunk_ranges=ranges,
            chunks=len(ranges),read_bytes=read,fixed_grid_read_bytes=regular,
            ld_restarts=actual['ld_records_at_chunk_starts'],
            fixed_grid_ld_restarts=old['data']['encoded']['ld_records_at_chunk_starts'],
            h2d_bytes=sum(b['h2d_bytes'] for b in work['blocks']),
            d2h_bytes=sum(b['d2h_bytes'] for b in work['blocks']),
            reader_initializations=sum('reader_init_seconds' in b for b in work['blocks'])))
    result=jagwas_candidate_runtime(choice,writer_prices(),preparation=preparation(choice),
        occupancy='dense',host_serial_fraction=.5)
    actual_parts=list((root/'adaptive'/'sumstats').glob('part_*.npz'))
    assert result['parts']==len(actual_parts)==len(previous['source_ranges'])
    assert result['indexed_part_bytes']==sum(p.stat().st_size for p in actual_parts)
    assert result['retained_variants']==previous['retained_variants']
    assert exact_reads>regular_reads
    assert sum(row['h2d_bytes'] for row in rows)==previous['dimensions']['N']*previous['dimensions']['M']
    assert sum(row['d2h_bytes'] for row in rows)==17*previous['dimensions']['M']
    assert source_identity()==source
    write_record(out/'report.json',dict(source_sha256=source,benchmark_sha256=sha256_file(__file__),
        execution_report_sha256=sha256_file(report_path),execution_source_sha256=previous['source_sha256'],
        inputs=previous['inputs'],dimensions=previous['dimensions'],shards=rows,
        exact_read_bytes=exact_reads,regular_grid_read_bytes=regular_reads,
        extra_read_bytes=exact_reads-regular_reads,parts=result['parts'],
        indexed_part_bytes=result['indexed_part_bytes'],retained_variants=result['retained_variants'],
        source_ranges=previous['source_ranges'],synthetic_graph_seconds=result['estimated_seconds'],
        scope='Exact source-count and output-size replay of a previously completed real two-GPU adaptive run. Independent range parses and physical indexed-part sizes match. Component prices and setup are synthetic controls; this does not validate runtime prediction, live continuation or an automatic switching policy.'))
    print(json.dumps(dict(parts=result['parts'],retained_variants=result['retained_variants'],
        extra_read_bytes=exact_reads-regular_reads,indexed_part_bytes=result['indexed_part_bytes'])),flush=True)


if __name__=='__main__':main()
