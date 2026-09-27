"""Duration-free reduced-output work for the shared analytical calculator.

Native FP32 complete-phenotype execution only. These source counts do not yet
price the selection/writer graph or authorize detailed autotune for reductions.
"""
import math
import io
import numpy as np
from .mechanistic_plan import _integer
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS, device_selection_shape


def significant_output_work(samples,markers,traits,chunk_markers,*,backend,
                            retained_per_block=None,store_beta=True,fsync=True,
                            max_selection_cells=DEVICE_SELECTION_MAX_CELLS,include_blocks=True,max_blocks=100000,significance_threshold=None,critical_table_reused=False):
    for name,value in [('samples',samples),('markers',markers),('traits',traits),
                       ('chunk_markers',chunk_markers),('max_selection_cells',max_selection_cells)]:
        _integer(name,value)
    if significance_threshold is not None and (not math.isfinite(significance_threshold) or not 0<significance_threshold<=1):
        raise ValueError('Significance threshold must be in (0,1]')
    if backend not in ('host','device'):raise ValueError('Explicit host/device selection required')
    if type(critical_table_reused) is not bool:raise ValueError('Boolean critical_table_reused required')
    if type(store_beta) is not bool or type(fsync) is not bool:raise ValueError('Boolean writer settings required')
    _integer('max_blocks',max_blocks)
    if type(include_blocks) is not bool:raise ValueError('Boolean include_blocks required')
    if not include_blocks and retained_per_block is not None:
        raise ValueError('Explicit retained counts require explicit blocks')
    chunks=(markers+chunk_markers-1)//chunk_markers
    full,tail=divmod(markers,chunk_markers)
    if backend=='host':
        block_count=chunks
        maximum_block_cells=min(markers,chunk_markers)*traits
    else:
        full_shape=(device_selection_shape(chunk_markers,traits,max_selection_cells)
                    if full else None)
        tail_shape=(device_selection_shape(tail,traits,max_selection_cells)
                    if tail else None)
        block_count=(full*(full_shape[2] if full_shape else 0)+
                     (tail_shape[2] if tail_shape else 0))
        maximum_block_cells=max(
            (shape[0]*shape[1] for shape in (full_shape,tail_shape)
             if shape is not None), default=0)
    if include_blocks and block_count>max_blocks:raise ValueError('Selection census exceeds max_blocks; use aggregate bounds')
    blocks=[]
    if include_blocks:
        for start in range(0,markers,chunk_markers):
            stop=min(markers,start+chunk_markers)
            width,height,_=(
                (traits,chunk_markers,1) if backend=='host' else
                device_selection_shape(stop-start,traits,max_selection_cells))
            for first in range(start,stop,height):
                last=min(stop,first+height)
                for left in range(0,traits,width):
                    right=min(traits,left+width)
                    blocks.append(dict(variant_range=[first,last],trait_range=[left,right],cells=(last-first)*(right-left)))
    if retained_per_block is not None:
        if not isinstance(retained_per_block,(list,tuple)) or len(retained_per_block)!=len(blocks):
            raise ValueError('One exact retained count per selection block required')
        for block,count in zip(blocks,retained_per_block):
            _integer('retained count',count,0)
            if count>block['cells']:raise ValueError('Retained count exceeds block cells')
            block['retained']=count
        retained=sum(retained_per_block);parts=sum(v>0 for v in retained_per_block)
    else:retained=parts=None
    payload_per_pair=28 if store_beta else 24  # int64 variant/trait, FP32 beta/t/df
    maximum=markers*traits
    return dict(mode='significant',backend=backend,source_chunks=chunks,blocks=blocks if include_blocks else None,
        selection_cells=maximum,selection_blocks=block_count,retained_pairs=retained,
        selection_count_d2h_bytes=4*block_count if backend=='device' else 0,
        retained_pair_bounds=[0,maximum],maximum_selection_block_cells=maximum_block_cells,
        critical_value_evaluations=0 if significance_threshold==1 or critical_table_reused else samples if backend=='device' else markers,
        critical_table_reused=critical_table_reused,
        critical_lookup_entries=markers if backend=='device' or critical_table_reused else 0,
        significance_threshold=significance_threshold,
        critical_h2d_bytes=4*(samples+1) if backend=='device' else 0,
        dense_result_ring=backend=='host',
        result_payload_d2h_bytes=(4*(1+int(store_beta))*maximum+5*markers if backend=='host' else None if retained is None else markers+28*retained),
        result_payload_d2h_bounds=([4*(1+int(store_beta))*maximum+5*markers]*2 if backend=='host' else [markers,markers+28*maximum]),
        indexed_array_payload_bytes=None if retained is None else payload_per_pair*retained,
        indexed_array_payload_bounds=[0,payload_per_pair*maximum],nonempty_parts=parts,
        part_fsync_calls=None if parts is None else parts*fsync,
        unpriced_terms=['GPU selection kernels and nonzero synchronization; dynamic-count transport is outside payload bytes',
            'Host comparisons, nonzero, gathers, threshold evaluation and Python dispatch',
            'Variable output staging, bounded multi-device result queue and single consumer',
            'NPY/ZIP headers, metadata publication, filesystem page-cache work and durable commit',
            'Selection workspace and allocator retention must be admitted before detailed execution'],
        scope='Exact native FP32 source payload/count ledger conditional on per-block survivor counts. No timing coefficient, observed selectivity estimate or autotune authorization.')


def jagwas_output_work(samples,markers,traits,chunk_markers,*,covariate_rank=0,retained_variants=None,fsync=True):
    for name,value in [('samples',samples),('markers',markers),('traits',traits),('chunk_markers',chunk_markers)]:_integer(name,value)
    _integer('covariate_rank',covariate_rank,0)
    if covariate_rank>=samples-1:raise ValueError('No residual phenotype rank')
    if retained_variants is not None:
        _integer('retained_variants',retained_variants,0)
        if retained_variants>markers:raise ValueError('Too many retained variants')
    if type(fsync) is not bool:raise ValueError('Boolean fsync required')
    chunks=(markers+chunk_markers-1)//chunk_markers
    return dict(mode='jagwas',trait_separable=False,maximum_full_rank_traits=samples-covariate_rank-1,
        full_rank_possible=traits<=samples-covariate_rank-1,
        correlation_fp32_matmul_flops=2*samples*traits*traits,
        setup_fp64_cholesky_order=traits,setup_fp64_triangular_solve_rhs=traits,
        persistent_inverse_cholesky_bytes=8*traits*traits,
        setup_explicit_matrix_live_bytes=28*traits*traits,
        projection_fp64_matmul_flops=2*markers*traits*traits,
        projected_square_cells=markers*traits,projection_sum_additions=markers*(traits-1),
        native_result_payload_d2h_bytes=17*markers,
        indexed_array_payload_bytes=None if retained_variants is None else 16*retained_variants,
        indexed_array_payload_bounds=[0,16*markers],source_chunks=chunks,
        part_fsync_call_bounds=[0,chunks if fsync else 0],
        unpriced_terms=['FP64 Cholesky and triangular solve service and vendor workspace',
            'Compiled projection/finite-check/nan replacement/square/reduction geometry',
            'Indexed writer serialization, metadata and durable storage commit',
            'Replicated setup and factors for any supported variant partition'],
        scope='Source operation/payload ledger. Explicit live setup matrices exclude phenotype, allocator and library workspace; this is not a complete peak-memory bound or a timing prediction.')


def significant_execution_layout(traits, trait_block, devices, reader_workers, *, queue_depth=None):
    """Shared executor/calculator contract for bounded significant trait workers.

    Readers are a global budget. Each device processes strided whole tiles;
    all active workers feed one bounded queue and one synchronous writer.
    The caller still owes separate complete host/device memory admission.
    """
    for name, value in [('traits', traits), ('trait_block', trait_block), ('reader_workers', reader_workers)]:
        _integer(name, value)
    if (not isinstance(devices, (list, tuple)) or not devices
            or any(not isinstance(device, str) or not device for device in devices)
            or len(set(devices)) != len(devices)):
        raise ValueError('Unique explicit significant tile devices required')
    tiles = (traits + trait_block - 1) // trait_block
    devices = list(devices[:tiles])
    if reader_workers < len(devices):
        raise ValueError('reader_workers must provide at least one reader per active device')
    if queue_depth is None:
        queue_depth = 4 * len(devices)
    _integer('queue_depth', queue_depth)
    per_device, remainder = divmod(reader_workers, len(devices))
    return dict(tiles=tiles, devices=devices,
        readers_per_device=[per_device + (index < remainder) for index in range(len(devices))],
        queue_depth=queue_depth if len(devices) > 1 else 0,
        writer_workers=1, producer_pending_slots=len(devices) if len(devices) > 1 else 0,
        consumer_slots=1, covariate_basis_evaluations=1)

def _indexed_array_part_work(rows, fields):
    """Uncompressed one-dimensional NPZ arrays, including NPY/ZIP headers.

    NumPy's force_zip64 local headers are 20 bytes longer than ordinary local
    headers. This bounded contract excludes members/offsets requiring ZIP64
    central-directory records; the caller must reject larger parts explicitly.
    """
    _integer('rows', rows, 0)
    if not rows:
        return dict(rows=0, arrays=[], array_payload_bytes=0, file_bytes=0, fsync_calls=0)
    arrays = []
    offset = 0
    central = 0
    for field, dtype in fields:
        header = io.BytesIO()
        np.lib.format.write_array_header_1_0(header, dict(descr=dtype, fortran_order=False, shape=(rows,)))
        payload = rows * np.dtype(dtype).itemsize
        member = len(header.getvalue()) + payload
        if member >= (1 << 31) - 1 or offset >= (1 << 31) - 1:
            raise ValueError('Indexed part exceeds bounded ZIP member/offset contract')
        filename_bytes = len((field + '.npy').encode('ascii'))
        local = 30 + filename_bytes + 20
        arrays.append(dict(field=field, dtype=dtype, payload_bytes=payload,
            npy_header_bytes=len(header.getvalue()), local_header_bytes=local,
            serialization_buffer_bytes=min(payload, 16 << 20)))
        offset += local + member
        central += 46 + filename_bytes
    if offset + central >= (1 << 31) - 1:
        raise ValueError('Indexed part exceeds bounded ZIP directory contract')
    return dict(rows=rows, arrays=arrays, array_payload_bytes=sum(a['payload_bytes'] for a in arrays),
        file_bytes=offset + central + 22, fsync_calls=1)


def jagwas_indexed_part_work(rows):
    """The actual writer stores absolute int64 indices and FP64 chi-square."""
    return _indexed_array_part_work(rows,[('variant_index','<i8'),('chi2','<f8')])


def jagwas_writer_work(markers, retained, *, fsync=True):
    """Host selection and one indexed part for a reduced native FP32 chunk.

    The scan transports 17 bytes per variant but releases status/df before
    handing the 12-byte beta/stat/index tuple to the writer. Unknown validity
    must be assessed over [0, markers], never silently assumed empty.
    """
    _integer('markers',markers);_integer('retained',retained,0)
    if retained>markers:raise ValueError('Retained variants exceed chunk width')
    if type(fsync) is not bool:raise ValueError('Boolean fsync required')
    part=jagwas_indexed_part_work(retained)
    return dict(markers=markers,retained=retained,float32_to_float64_cells=markers,
        finite_predicate_cells=markers,nonzero_cells=markers,
        nonzero_indices=retained,index_add_elements=retained,statistic_gather_elements=retained,
        native_result_d2h_bytes=17*markers,owned_chunk_array_bytes=12*markers,
        indexed_array_payload_bytes=part['array_payload_bytes'],part=part,
        part_fsync_calls=int(fsync and retained>0),
        # Sum all distinct writer arrays and both serialization buffers. This
        # bounds explicit arrays, not allocator RSS. The previous iteration's
        # converted values and keep indices can survive assignment; reserve
        # both at full chunk width even when this chunk retains no variants.
        writer_array_bytes_upper=25*markers+24*retained+2*min(8*retained,16<<20),
        scope='Source-counted native FP32 JAGWAS writer work for explicit survivor count. No timing, metadata, allocator or filesystem-memory bound.')


def _joint_cpu_step(cpu,traffic,cpu_fraction,dram_bytes_per_second,host_serial_fraction):
    for name,value in [('CPU work',cpu),('DRAM work',traffic),('CPU fraction',cpu_fraction),
                       ('DRAM capacity',dram_bytes_per_second),('host serial fraction',host_serial_fraction)]:
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
            raise ValueError('Invalid '+name)
    if not 0<cpu_fraction<=1 or not dram_bytes_per_second or host_serial_fraction>1:
        raise ValueError('Positive CPU/DRAM and bounded host serialization required')
    seconds=max(cpu/cpu_fraction,traffic/dram_bytes_per_second)
    return dict(seconds=seconds,resources=dict(cpu=cpu/seconds if seconds else 0.,
        host_serial=cpu*host_serial_fraction/seconds if seconds else 0.,dram=traffic/seconds if seconds else 0.))


def jagwas_host_selection_service(work,prices,*,cpu_fraction,dram_bytes_per_second,host_serial_fraction):
    """Price the exact FP32-to-FP64, validity and gather sequence after dequeue.

    Rates describe independently measured fixed generic operations. Retained
    counts are supplied scenarios; neither a threshold nor prior GWAS timings
    determine them. DRAM traffic is logical source work, not a physical bound.
    """
    b,h=work['markers'],work['retained']
    terms=[('fp32_to_fp64_view',b,12*b),('finite_fp64',b,9*b),
        ('flatnonzero_nonempty' if h else 'flatnonzero_empty',b,b+8*h),
        ('index_add',h,16*h),('fp64_gather',h,24*h)]
    steps=[]
    for name,count,traffic in terms:
        price=prices[name]
        if set(price)!={'call_cpu_seconds','unit_cpu_seconds'}:
            raise ValueError('Joint primitive requires fixed and per-unit CPU prices')
        if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in price.values()):
            raise ValueError('Invalid joint primitive price: '+name)
        cpu=price['call_cpu_seconds']+count*price['unit_cpu_seconds']
        steps.append(_joint_cpu_step(cpu,traffic,cpu_fraction,dram_bytes_per_second,host_serial_fraction))
    return steps


def jagwas_archive_service(work,price,profile,*,host_serial_fraction):
    """Two-array uncompressed NPZ serialization, page-cache copy and commit.

    The field signature prevents reuse of significant-pairs archive costs.
    CUDA transport, host selection, metadata publication and final directory
    commit are separate stages. Non-durable buffering is not modeled here.
    """
    part=work['part']
    if not part['rows']:return []
    schema=[[row['field'],row['dtype']] for row in part['arrays']]
    if price.get('field_schema')!=schema:
        raise ValueError('Independent archive price must match the two-array JAGWAS schema')
    if work['part_fsync_calls']!=1:
        raise ValueError('Joint archive service requires the durable part boundary')
    names=['call_cpu_seconds','byte_cpu_seconds']
    values=[price[name] for name in names]
    writeback=profile['writeback_service']
    values += [writeback['pagecache_seconds_per_byte'],writeback['storage_seconds_per_byte'],profile['fsync_seconds']]
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in values):
        raise ValueError('Invalid independent archive/storage price')
    fixed,bulk,page,storage,fsync=values
    if not storage:raise ValueError('Positive storage service required')
    # Seekable ZIP output rewrites each local header at close. Final storage
    # extent counts it once; CPU page-cache copying sees both submissions.
    submitted=part['file_bytes']+sum(row['local_header_bytes'] for row in part['arrays'])
    cpu=fixed+bulk*part['array_payload_bytes']+page*submitted
    traffic=5*part['array_payload_bytes']+2*(submitted-part['array_payload_bytes'])
    return [_joint_cpu_step(cpu,traffic,profile['cpu_fraction'],profile['shared_dram_bytes_per_second'],host_serial_fraction),
        dict(seconds=part['file_bytes']*storage,resources=dict(output=1/storage)),
        dict(seconds=fsync)]
