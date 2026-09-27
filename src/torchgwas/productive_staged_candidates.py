"""Construct post-output chunk candidates from one active measured context.

This does no PGEN work and does not choose or apply a new layout. The current
unissued partitions stay on their original GPU and phenotype/variant ranges;
only future chunk width varies. Tile and device moves need separate memory
and output-ownership admission before becoming runtime candidates.
"""
from copy import deepcopy
import math

from .layout_frontier import KIND as FRONTIER_KIND


_SOURCE_FIELDS=('decode_units','cpu_fraction','depth','decode_workers',
                'cpu_available_cores','shared_dram_bytes_per_second',
                'read_bytes_per_second')


def _positive(value,name):
    if (isinstance(value,bool) or not isinstance(value,(int,float)) or
            not math.isfinite(value) or value<=0):
        raise ValueError('Positive measured '+name+' required')
    return value


def fixed_partition_chunk_candidates(frontier,context,*,chunk_sizes,
                                     partition_axis,covariate_rank,output,
                                     writer_options=None,reduction_prices=None,
                                     significance_threshold=None):
    """Build at most four exact-frontier, independently priced chunk layouts.

    Transfer and output ceilings come from the active context, never from GPU
    identity or a fitted job runtime. The caller must still audit original
    measurement ages, source/work price targets, memory and final completion.
    """
    if (not isinstance(frontier,dict) or frontier.get('kind')!=FRONTIER_KIND or
            not isinstance(context,dict) or
            not isinstance(chunk_sizes,list) or not 0<len(chunk_sizes)<=4 or
            any(type(size) is not int or size<1 for size in chunk_sizes) or
            len(set(chunk_sizes))!=len(chunk_sizes) or
            partition_axis not in ('trait','variant') or
            type(covariate_rank) is not int or covariate_rank<0 or
            not isinstance(output,dict) or
            set(output)!={'block_bytes','queue_depth','store_beta','fsync'}):
        raise ValueError('Bounded active staged chunk construction required')
    mode=frontier.get('reduction')
    if (mode not in (None,'significant','jagwas') or
            (mode=='jagwas' and partition_axis!='variant') or
            (mode=='significant' and partition_axis!='trait')):
        raise ValueError('Staged chunk axis differs from output mode')
    rectangles=frontier.get('rectangles')
    if not isinstance(rectangles,list) or not 0<len(rectangles)<=16:
        raise ValueError('Bounded unissued partitions required')
    ids=[row.get('id') for row in rectangles if isinstance(row,dict)]
    devices={row.get('device') for row in rectangles if isinstance(row,dict)}
    profiles=context.get('profiles')
    if (len(ids)!=len(rectangles) or len(set(ids))!=len(ids) or
            not all(isinstance(key,str) and key for key in ids) or
            not isinstance(profiles,dict) or not devices<=set(profiles)):
        raise ValueError('Every unissued partition needs an active device profile')
    shared=context.get('shared_capacities')
    transfer=context.get('shared_transfer_capacities')
    if (not isinstance(shared,dict) or not isinstance(transfer,dict) or
            set(transfer)!={'h2d','d2h'}):
        raise ValueError('Measured shared source/output and transfer capacities required')
    for name in ('cpu','dram','input','output'):
        _positive(shared.get(name),'shared '+name)
    for name in ('h2d','d2h'):
        _positive(transfer[name],'shared '+name)
    if (type(output['store_beta']) is not bool or
            type(output['fsync']) is not bool):
        raise ValueError('Explicit output fields and durability required')
    if mode is None:
        if (not isinstance(writer_options,dict) or
                set(writer_options)!={'block_bytes','queue_depth',
                                      'borrow_chunks','fsync','writeback_bytes',
                                      'sync_file_range','store_variant_df'} or
                writer_options['fsync']!=output['fsync'] or
                writer_options['block_bytes']!=output['block_bytes'] or
                writer_options['queue_depth']!=output['queue_depth']):
            raise ValueError('Actual dense writer settings required')
    elif writer_options is not None or not isinstance(reduction_prices,dict):
        raise ValueError('Current reduced-output primitive record required')
    if mode=='significant':
        if (isinstance(significance_threshold,bool) or
                not isinstance(significance_threshold,(int,float)) or
                not math.isfinite(significance_threshold) or
                not 0<significance_threshold<=1):
            raise ValueError('Bound significant-pair threshold required')
    elif significance_threshold is not None:
        raise ValueError('Significance threshold applies only to selected pairs')
    source_profiles={}
    for row in rectangles:
        active=profiles[row['device']]
        if not isinstance(active,dict):
            raise ValueError('Active device profile required')
        if any(field not in active for field in _SOURCE_FIELDS):
            raise ValueError('Measured source prices required on every GPU')
        source_profiles[row['id']]={field:active[field] for field in _SOURCE_FIELDS}
        if active.get('input_read_cpu_prices') is not None:
            source_profiles[row['id']]['input_read_cpu_prices']=active[
                'input_read_cpu_prices']
    device_profiles={device:profiles[device] for device in devices}
    if any(not isinstance(profile.get('gpu_resources'),dict)
           for profile in device_profiles.values()):
        raise ValueError('Measured GPU resources required on every device')
    h2d={device:_positive(profile.get('h2d_bytes_per_second'),
                          device+' H2D') for device,profile in device_profiles.items()}
    d2h={device:_positive(profile.get('d2h_bytes_per_second'),
                          device+' D2H') for device,profile in device_profiles.items()}
    fp32={device:_positive(profile.get('gpu_resources',{}).get(
        'fp32_flops_per_second'),device+' FP32')
        for device,profile in device_profiles.items()}
    links=deepcopy(context.get('shared_links',()))
    candidates=[]
    for size in chunk_sizes:
        compute=dict(covariate_rank=covariate_rank,
            shared_h2d_bytes_per_second=transfer['h2d'],
            per_device_h2d_bytes_per_second=h2d,
            peak_fp32_flops_per_second=fp32,shared_links=links)
        result=dict(store_beta=output['store_beta'],
            shared_d2h_bytes_per_second=transfer['d2h'],
            per_device_d2h_bytes_per_second=d2h,
            output_bytes_per_second=shared['output'],shared_links=links)
        if mode is None:
            result['dense_writer_options']=deepcopy(writer_options)
            services=dict(writer_profiles=device_profiles)
        elif mode=='jagwas':
            writer=reduction_prices.get('writer_prices')
            if not isinstance(writer,dict) or not {'prices','archive'}<=set(writer):
                raise ValueError('Complete JAGWAS writer primitive bank required')
            compute['peak_fp64_flops_per_second']={device:_positive(
                profile.get('gpu_resources',{}).get('fp64_flops_per_second'),
                device+' FP64') for device,profile in device_profiles.items()}
            result['jagwas_writer_fsync']=output['fsync']
            services=dict(selection_prices=writer['prices'],
                cpu_fraction_by_device={device:profile['cpu_fraction']
                                        for device,profile in device_profiles.items()},
                archive_price=writer['archive'],archive_profiles=device_profiles)
        else:
            if 'archive' not in reduction_prices:
                raise ValueError('Complete significant-pair primitive bank required')
            result.update(significant_backend='host',
                significant_threshold_one=(significance_threshold==1),
                significant_writer_fsync=output['fsync'])
            services=dict(selection_prices=reduction_prices,
                selection_profiles=device_profiles,
                archive_prices=reduction_prices['archive'],
                archive_profiles=device_profiles)
        candidates.append(dict(id='chunk:'+str(size),chunk_markers=size,
            partitions=deepcopy(rectangles),partition_axis=partition_axis,
            source_profiles=source_profiles,compute_options=compute,
            output_options=result,mode_service_options=services,
            shape_profiles={device:dict(profile,chunk_markers=size)
                            for device,profile in device_profiles.items()}))
    return candidates
