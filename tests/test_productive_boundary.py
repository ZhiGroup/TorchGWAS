"""Issued work and written output remain distinct at a productive frontier."""
import time
from dataclasses import replace

from torchgwas.productive_boundary import ProductiveBoundaryProgress
from torchgwas.productive_output_backlog import productive_output_backlog
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.sumstats import DenseWriteProgress
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition


def dense_event(start, end, traits, device, df):
    return DenseWriteProgress(start, end, traits,
        (end-start)*(traits[1]-traits[0])*8, df,
        time.perf_counter(), 'sumstats', device)


def test_dense_boundaries_use_source_coordinates_and_keep_df_separate():
    partitions=[
        dict(id='low', device='cuda:0', variant_range=[100, 108],
             trait_range=[0, 3]),
        dict(id='high', device='cuda:1', variant_range=[100, 108],
             trait_range=[3, 5])]
    run=ProductiveTuningRun(partitions, chunk_sizes=[2, 4], initial=4)
    progress=ProductiveBoundaryProgress(partitions, None)
    assert run.for_partition('low')(100,108,4)==4
    assert run.for_partition('low')(104,108,4)==4
    assert run.for_partition('high')(100,108,4)==4
    for event in (dense_event(0,2,(0,3),'cuda:0',1),
                  dense_event(2,4,(0,3),'cuda:0',3),
                  dense_event(0,4,(3,5),'cuda:1',2)):
        run.output_written(event)
        progress.observe(event)
    boundary=progress.bind(run.snapshot())
    assert boundary['valid'] and boundary['written_events']==3
    low,high=boundary['partitions']
    assert low['matrix_written_to']==104 and low['df_written_to']==103
    assert low['issued_not_matrix_written']==[104,108]
    assert low['issued_not_matrix_written_pairs']==12
    assert low['issued_not_df_written_markers']==5
    assert high['matrix_written_to']==104 and high['df_written_to']==102
    assert high['issued_not_matrix_written']==[104,104]
    assert boundary['issued_not_written_pairs']==12
    backlog=productive_output_backlog(boundary,total_traits=5)
    assert backlog['array_payload_bytes_upper']==96
    assert backlog['partitions'][0]['pending_pairs']==12
    with_df=productive_output_backlog(boundary,total_traits=5,
                                      store_variant_df=True)
    assert with_df['array_payload_bytes_upper']==124
    run.finish(successful=False)


def test_indexed_completed_ranges_exclude_only_matching_reserved_chunks():
    partitions=[dict(id='j', device='cuda:0',
        variant_range=[20,32], trait_range=[0,5])]
    run=ProductiveTuningRun(partitions, chunk_sizes=[2,4], initial=4)
    progress=ProductiveBoundaryProgress(partitions,'jagwas')
    control=run.for_partition('j')
    assert [control(start,32,4) for start in (20,24,28)]==[4,4,4]
    owner=IndexedOutputPartition('cuda:0',(20,32),(0,5))
    for lo,hi,rows in ((20,24,0),(28,32,2)):
        event=IndexedChunkWrite(lo-20,hi-20,'jagwas',rows,
            0 if rows==0 else 128,None if rows==0 else 'part.npz',
            time.perf_counter(),time.perf_counter(),bool(rows),owner,(lo,hi))
        run.output_written(event)
        progress.observe(event)
    boundary=progress.bind(run.snapshot())
    row=boundary['partitions'][0]
    assert boundary['valid'] and row['issued_not_indexed_written']==[[24,28]]
    assert row['indexed_written_ranges']==[[20,24],[28,32]]
    assert row['indexed_rows']==2 and row['fsynced_parts']==1
    assert row['issued_not_indexed_written_pairs']==20
    assert boundary['issued_not_written_pairs']==20
    backlog=productive_output_backlog(boundary,total_traits=5,
                                      occupancy_scenario='dense')
    assert backlog['partitions'][0]['selected_rows']==4
    assert backlog['array_payload_bytes_upper']==64
    run.finish(successful=False)


def test_invalid_optional_output_progress_refuses_binding_without_raising():
    partitions=[dict(id='x', device='cuda:0',
        variant_range=[0,8], trait_range=[0,3])]
    run=ProductiveTuningRun(partitions, chunk_sizes=[2,4], initial=4)
    progress=ProductiveBoundaryProgress(partitions,None)
    assert run.for_partition('x')(0,8,4)==4
    event=dense_event(0,2,(0,3),'cuda:1',1)
    run.output_written(event)
    progress.observe(event)
    assert progress.bind(run.snapshot())['valid'] is False
    assert progress.snapshot()['invalid_events']==1
    run.finish(successful=False)


def test_native_dense_writer_prefix_binds_only_with_unique_source_owner():
    partitions=[dict(id='a',device='cuda:0',variant_range=[0,4],trait_range=[0,3]),
                dict(id='b',device='cuda:1',variant_range=[4,8],trait_range=[0,3])]
    run=ProductiveTuningRun(partitions,chunk_sizes=[2,4],initial=4)
    progress=ProductiveBoundaryProgress(partitions,None)
    assert run.for_partition('a')(0,4,4)==4
    assert run.for_partition('b')(4,8,4)==4
    event=dense_event(0,4,(0,3),None,4)
    run.output_written(event);progress.observe(event)
    bound=progress.bind(run.snapshot())
    assert bound['valid'] and bound['partitions'][0]['matrix_written_to']==4
    assert bound['partitions'][1]['matrix_written_to']==4
    crossing=dense_event(4,8,(0,3),None,8)
    run.output_written(crossing);progress.observe(crossing)
    assert progress.bind(run.snapshot())['valid']
    run.finish(successful=False)


def test_native_dense_writer_cross_shard_prefix_stays_unbound():
    partitions=[dict(id='a',device='cuda:0',variant_range=[0,4],trait_range=[0,3]),
                dict(id='b',device='cuda:1',variant_range=[4,8],trait_range=[0,3])]
    run=ProductiveTuningRun(partitions,chunk_sizes=[2,4],initial=4)
    progress=ProductiveBoundaryProgress(partitions,None)
    assert run.for_partition('a')(0,4,4)==4
    assert run.for_partition('b')(4,8,4)==4
    event=dense_event(0,8,(0,3),None,8)
    run.output_written(event);progress.observe(event)
    assert not progress.bind(run.snapshot())['valid']
    assert progress.snapshot()['invalid_events']==1
    run.finish(successful=False)


def test_dense_queue_observations_keep_separate_writer_owners():
    partitions=[dict(id='a',device='cuda:0',variant_range=[0,4],trait_range=[0,3]),
                dict(id='b',device='cuda:1',variant_range=[4,8],trait_range=[0,3])]
    run=ProductiveTuningRun(partitions,chunk_sizes=[2,4],initial=4)
    progress=ProductiveBoundaryProgress(partitions,None)
    assert run.for_partition('a')(0,4,4)==4
    assert run.for_partition('b')(4,8,4)==4
    queue=dict(kind='torchgwas.dense_writer_queue_observation.v1',valid=True)
    for lo,hi,directory in ((0,4,'writer-a'),(4,8,'writer-b')):
        event=replace(dense_event(lo,hi,(0,3),None,hi),directory=directory,
                      writer_queue=queue)
        run.output_written(event);progress.observe(event)
    bound=progress.bind(run.snapshot())
    assert bound['valid'] and set(bound['writer_queues'])=={'writer-a','writer-b'}
    assert bound['writer_queues']['writer-a']['event_partition_id']=='a'
    assert bound['writer_queues']['writer-b']['event_partition_id']=='b'
    run.finish(successful=False)


def test_significant_pending_scenario_uses_global_coordinates_across_chunks():
    scenario=dict(retained_fraction=[2,7],placement='clustered')
    base=dict(kind='torchgwas.productive_output_boundary.v1',valid=True,
              issued_revision=9,written_events=2,reduction='significant')
    def boundary(spans):
        return dict(base,partitions=[dict(id='tile',trait_range=[2,5],
                                          issued_not_indexed_written=spans)])
    whole=productive_output_backlog(boundary([[20,30]]),total_traits=6,
                                    occupancy_scenario=scenario)
    split=productive_output_backlog(boundary([[20,24],[24,30]]),total_traits=6,
                                    occupancy_scenario=scenario)
    assert whole['array_payload_bytes_upper']==split['array_payload_bytes_upper']
    assert whole['partitions'][0]['selected_rows']==split['partitions'][0]['selected_rows']
