"""CUDA selector phases from installed launches and source/resource equations.

The launch partition is exact for the checked PyTorch 2.5.1/CUB 2.3 policy.
Traffic, cache attainment and arithmetic issue are explicitly modeled work,
not recovered physical counters or a fitted selector-duration table.
"""
import math
from .mechanistic_plan import _integer
from .first_principles import positive
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS

SOURCES={
    'nonzero':'https://github.com/pytorch/pytorch/blob/v2.5.1/aten/src/ATen/native/cuda/Nonzero.cu',
    'reduce':'https://github.com/NVIDIA/cccl/blob/v2.3.1/cub/cub/device/dispatch/dispatch_reduce.cuh',
    'select_policy':'https://github.com/NVIDIA/cccl/blob/v2.3.1/cub/cub/device/dispatch/tuning/tuning_select_if.cuh',
    'select_agent':'https://github.com/NVIDIA/cccl/blob/v2.3.1/cub/cub/agent/agent_select_if.cuh',
    'prefix':'https://github.com/NVIDIA/cccl/blob/v2.3.1/cub/cub/agent/single_pass_scan_operators.cuh'}


def selection_gpu_census(work,census,current_context):
    """Exact source shape, occupancy and installed library/device identity."""
    if census.get('durations_recorded') is not False:raise ValueError('Untimed selector census required')
    context=census.get('context',{})
    fields=['torch_version','cuda_runtime','device_uuid','compute_capability','sm_count','library_sha256']
    if any(key not in current_context or context.get(key)!=current_context[key] for key in fields):
        raise ValueError('Selector kernel census differs from current installed context')
    if context['cuda_runtime']!='12.4':raise ValueError('Unverified CUB toolkit source policy')
    if any(census['source_sha256'].get(key)!=value for key,value in work['source_sha256'].items()):
        raise ValueError('Selector source changed since kernel census')
    counts=[b['retained'] for b in work['blocks']]
    matches=[r for r in census['rows'] if [r.get(k) for k in ['N','B','K','max_cells']]==
        [work[k] for k in ['samples','markers','traits','max_cells']] and [b['retained'] for b in r['blocks']]==counts]
    if len(matches)!=1:raise ValueError('One exact selector geometry/occupancy census required')
    row=matches[0]
    if row.get('durations_recorded') is not False or row.get('independent_cpu_arrays_equal') is not True:
        raise ValueError('Selector geometry requires untimed independent array checks')
    return row['kernels']


def _geometry(row):
    if not isinstance(row,dict) or set(row)!={'name','geometry'}:
        raise ValueError('Duration-free kernel name and geometry only')
    if not isinstance(row['name'],str) or not row['name']:raise ValueError('Compiled kernel name required')
    geometry=row['geometry']
    if set(geometry)-{'grid','block','registers per thread','shared memory'}:
        raise ValueError('Unexpected timed or unsupported kernel geometry')
    for name in ['grid','block']:
        shape=geometry.get(name)
        if not isinstance(shape,list) or len(shape)!=3 or any(type(v) is not int or v<1 for v in shape):
            raise ValueError('Positive three-dimensional compiled '+name+' required')
    for name in ['registers per thread','shared memory']:
        if name in geometry:_integer(name,geometry[name],0)
    return geometry


def selection_gpu_work(work,kernels,*,compute_capability,lookback_windows=1):
    """Associate every recorded launch with a source operation or CUB phase.

    lookback_windows declares predecessor windows per non-first CUB tile;
    retry spinning remains unpriced. No launch count is inferred from a time.
    Captured input policy and grids must match the supported source policy.
    """
    if str(work['torch_version']).split('+')[0]!='2.5.1':raise ValueError('PyTorch 2.5.1 source policy required')
    policy={(8,0):(384,6),(9,0):(384,11)}.get(tuple(compute_capability))
    if policy is None:raise ValueError('Unverified CUB architecture policy')
    _integer('lookback_windows',lookback_windows)
    if not isinstance(kernels,list):raise ValueError('Explicit compiled launch census required')
    for kernel in kernels:_geometry(kernel)
    index=0;rows=[];operations={};groups=[]
    blocks_by_nonzero=iter(work['blocks'])
    def take(family):
        nonlocal index
        if index>=len(kernels) or family not in kernels[index]['name']:
            raise ValueError('Compiled launch does not match source phase '+family+' at '+str(index))
        value=kernels[index];index+=1;return value
    def add(kernel,phase,read,write,ops,**extra):
        row=dict(kernel=kernel,phase=phase,read_bytes=read,write_bytes=write,
            logical_bytes=read+write,modeled_scalar_ops=ops,**extra)
        rows.append(row);return len(rows)-1
    for step_index,step in enumerate(work['steps']):
        if step['alias_only'] or step['allocation_only'] or step['host_copy']:continue
        if step['op']!='aten.nonzero.default':
            # The packed payload's stack is one batched concatenation copy.
            kernel=take('CatArrayBatchedCopy' if step['op']=='aten.stack.default' else 'elementwise_kernel')
            if not kernel['name'].startswith('void at::native::'):
                raise ValueError('Unknown eager CUDA kernel family')
            cells=max((math.prod(t['shape']) for t in step['outputs']),default=0)
            if not cells:raise ValueError('Zero-extent eager operation unexpectedly launches a kernel')
            scalar=cells*(2 if step['op']=='aten.clamp.default' else 1)
            operations[str(step_index)]=[add(kernel,'eager',step['read_bytes'],step['write_bytes'],scalar,
                op=step['op'],step_index=step_index)]
            continue
        block=next(blocks_by_nonzero);cells,retained=block['cells'],block['retained']
        # PyTorch 2.5.1 nonzero is one count reduce and one flagged select at
        # any extent below INT_MAX; above 1M cells the count grid saturates at
        # CUB's occupancy bound, which the partial check below admits.
        if cells>DEVICE_SELECTION_MAX_CELLS:raise ValueError('Selector block exceeds the CUDA nonzero limit')
        phases={};count=[]
        first=take('DeviceReduce')
        if 'Policy600' not in first['name'] or 'NonZeroOp<bool>' not in first['name']:
            raise ValueError('Unknown bool count-reduction policy')
        if first['geometry']['block']!=[256,1,1]:raise ValueError('Unexpected count-reduction block')
        if 'DeviceReduceSingleTileKernel<' in first['name']:
            if cells>4096 or first['geometry']['grid']!=[1,1,1]:raise ValueError('Invalid single-tile count geometry')
            count.append(add(first,'count',cells,4,2*cells-1))
        elif 'DeviceReduceKernel<' in first['name']:
            partials=math.prod(first['geometry']['grid'])
            if cells<=4096 or not 1<=partials<=math.ceil(cells/4096):raise ValueError('Invalid partial count geometry')
            count.append(add(first,'count',cells,4*partials,2*cells-partials))
            second=take('DeviceReduceSingleTileKernel<')
            if 'Policy600' not in second['name'] or 'NonZeroOp<bool>' in second['name'] or second['geometry']['grid']!=[1,1,1] or second['geometry']['block']!=[256,1,1]:
                raise ValueError('Invalid final partial-count reduction')
            count.append(add(second,'count',4*partials,4,partials-1))
        else:raise ValueError('Unrecognized count-reduction kernel')
        phases['count']=count
        init=take('DeviceCompactInitKernel<');sweep=take('DeviceSelectSweepKernel<')
        threads,items=policy;tiles=math.ceil(cells/(threads*items))
        if (init['geometry']['block']!=[128,1,1] or init['geometry']['grid']!=[math.ceil(tiles/128),1,1]
                or sweep['geometry']['block']!=[threads,1,1] or sweep['geometry']['grid']!=[tiles,1,1]):
            raise ValueError('Compiled select geometry differs from the architecture policy')
        if not all(text in sweep['name'] for text in ['CountingInputIterator<long','NonZeroOp<bool>','ScanTileState<int, true>','Policy900']):
            raise ValueError('Unknown flagged-selection types or policy')
        init_row=add(init,'flagged_select',0,8*(tiles+32)+4,tiles+33)
        # One 8-byte descriptor for the first tile (unless also last); every
        # subsequent tile publishes partial and inclusive descriptors.
        writes=8*(1+2*(tiles-1)) if tiles>1 else 0
        lookback=32*8*lookback_windows*(tiles-1)
        # Register reduction + scan and five warp-scan rounds. This explicit
        # source algorithm work proxy excludes compiler instruction expansion.
        padded=tiles*threads*items
        scan_ops=2*padded+5*threads*tiles
        sweep_row=add(sweep,'flagged_select',cells+lookback,8*retained+4+writes,
            cells+scan_ops+32*lookback_windows*(tiles-1),cub_tiles=tiles,
            lookback_bytes=lookback,lookback_windows=lookback_windows,
            shared_scatter_bytes_upper=16*retained,items_per_thread=items)
        phases['flagged_select']=[init_row,sweep_row]
        if retained:
            scatter=take('write_indices<')
            if scatter['geometry']['grid']!=[math.ceil(retained/256),1,1] or scatter['geometry']['block']!=[256,1,1]:
                raise ValueError('Compiled coordinate scatter differs from source extent')
            phases['coordinate_scatter']=[add(scatter,'coordinate_scatter',8*retained,16*retained,2*retained,
                int64_divmod_pairs=2*retained)]
        groups.append(phases)
    if index!=len(kernels):raise ValueError('Unassigned compiled CUDA launches')
    if len(groups)!=len(work['blocks']):raise ValueError('Missing source nonzero phase group')
    return dict(kernels=rows,operations=operations,nonzero_groups=groups,kernel_count=len(rows),
        compute_capability=list(compute_capability),logical_bytes=sum(r['logical_bytes'] for r in rows),
        source_sha256=work['source_sha256'],sources=SOURCES,prediction_complete=False,
        unpriced_terms=['CUB prefix retry spinning and scheduler-dependent predecessor progress',
            'Issued instruction expansion, integer addressing and intra-CTA scan/barrier latency',
            'Survivor placement within CUB tiles changes shared-scatter work',
            'Physical cache transactions, write allocation and gather coalescing'],
        scope='Exact source launch association with explicit CUB policy, partial buffers and predecessor-window traffic. Arithmetic and traffic are analytical proxies, not measured physical counters.')


def selection_gpu_service(ledger,resources,*,traffic_mode,int64_divmod_per_second=None):
    """Use the existing independent GPU resource bank; never a timing grid.

    logical_hbm and logical_l2 are cache-attainment scenarios, not guaranteed
    bounds. Div/mod capacity must be independent when coordinates are emitted.
    CPU API calls and all count/payload transfers remain in the caller graph.
    """
    if traffic_mode not in ('logical_hbm','logical_l2'):raise ValueError('Explicit cache-attainment scenario required')
    fraction=positive('GPU fraction',resources.gpu_fraction)
    memory_rate=positive('memory capacity',getattr(resources,'hbm_bytes_per_second' if traffic_mode=='logical_hbm' else 'l2_bytes_per_second'))
    arithmetic=positive('scalar issue capacity',resources.fp32_flops_per_second)
    launch=positive('kernel launch',resources.kernel_launch_seconds,True)
    reports=[]
    for row in ledger['kernels']:
        active=min(1.,math.prod(row['kernel']['geometry']['grid'])/resources.sm_count)
        scalar=row['modeled_scalar_ops']/(arithmetic*active*fraction)
        if row.get('int64_divmod_pairs'):
            rate=positive('independent int64 div/mod capacity',int64_divmod_per_second) if int64_divmod_per_second is not None else None
            if rate is None:raise ValueError('Independent int64 div/mod capacity required for nonempty coordinates')
            scalar+=row['int64_divmod_pairs']/(rate*active*fraction)
        memory=row['logical_bytes']/(memory_rate*fraction)
        seconds=launch/fraction+max(memory,scalar)
        reports.append(dict(seconds=seconds,memory_seconds=memory,arithmetic_seconds=scalar,
            phase=row['phase'],kernel_name=row['kernel']['name']))
    operations={key:dict(op=ledger['kernels'][indices[0]]['op'],seconds=sum(reports[i]['seconds'] for i in indices)) for key,indices in ledger['operations'].items()}
    nonzero=[{phase:dict(seconds=sum(reports[i]['seconds'] for i in indices)) for phase,indices in group.items()} for group in ledger['nonzero_groups']]
    return dict(operation_services=operations,nonzero_services=nonzero,kernel_reports=reports,
        gpu_service_seconds=sum(r['seconds'] for r in reports),traffic_mode=traffic_mode,
        unpriced_terms=ledger['unpriced_terms']+['Scalar/integer proxy shares the existing FP32-equivalent issue resource',
            'CUB internal shared-memory traffic and CTA synchronization lack an independent attainment service'],
        prediction_complete=False,scope='Conditional sum of GPU kernel services from exact installed launches and independent capacities. Host barriers, copies and dispatch are priced separately.')
