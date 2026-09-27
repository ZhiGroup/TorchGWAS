"""Installed CUDA nonzero allocation requests, without duration-based fitting."""
from .mechanistic_plan import _integer
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS, device_selection_shape


def device_nonzero_workspace(census, cells, context):
    """Verify exact requested allocations at this C and installed context.

    No size interpolation: preparation must capture each distinct mask extent.
    The CUB count allocation remains alive while selection scratch is allocated.
    Coordinate output is priced separately from these private workspaces.
    """
    _integer('cells',cells)
    if census.get('duration_fields_recorded') is not False:
        raise ValueError('Untimed allocation census required')
    recorded=census['context']
    for key in ['torch_version','cuda_runtime','device_uuid','compute_capability','sm_count','library_sha256']:
        if key not in context or context[key]!=recorded.get(key):
            raise ValueError('Nonzero allocation context differs: '+key)
    if context['torch_version'].split('+')[0]!='2.5.1':
        raise ValueError('Nonzero allocation source contract requires PyTorch 2.5.1')
    rows=[row for row in census['rows'] if row['cells']==cells]
    if sorted(row['retained'] for row in rows)!=sorted({0,1,cells}):
        raise ValueError('Exact empty, single and full occupancy allocation controls required')
    round_up=lambda n:((n+511)//512)*512
    scratch=[];output_excess=[]
    for row in rows:
        if row['shape']!=[1,cells] or row['dtype']!='torch.bool' or row['duration_fields_recorded'] is not False:
            raise ValueError('Contiguous two-dimensional bool source contract required')
        if any('time_us' in event or 'duration' in event for event in row['trace']):
            raise ValueError('Discard allocator timestamps before supplying the memory census')
        count=row['retained']
        if row['output_bytes']!=16*count:
            raise ValueError('Unexpected nonzero coordinate extent')
        allocations=[e for e in row['trace'] if e['action']=='alloc']
        frees=[e for e in row['trace'] if e['action']=='free_requested']
        if len(allocations)!=(4 if count else 3) or len(frees)!=3:
            raise ValueError('Installed nonzero allocation sequence differs from source')
        if allocations[0]['size']!=4:
            raise ValueError('Missing device count allocation')
        if count:
            coordinates=allocations[2]
            if coordinates['addr']!=row['output_ptr'] or coordinates['size']!=16*count:
                raise ValueError('Coordinate output does not match its allocation')
        private=[allocations[0],allocations[1],allocations[-1]]
        for event in private: _integer('allocation request',event['size'])
        if [(e['addr'],e['size']) for e in frees]!=[(e['addr'],e['size']) for e in [private[1],private[2],private[0]]]:
            raise ValueError('Nonzero private allocation release order changed')
        allocated=row.get('output_allocation_bytes')
        if not isinstance(allocated,int) or isinstance(allocated,bool) or allocated<round_up(16*count) or (not count and allocated):
            raise ValueError('Live coordinate allocation block must be captured separately from its request')
        output_excess.append(allocated-round_up(16*count))
        if row['extra_allocated_peak_bytes']!=sum(round_up(e['size']) for e in private)+allocated:
            raise ValueError('Nonzero allocation peak does not reconcile to source overlap')
        scratch.append([e['size'] for e in private])
    if any(row!=scratch[0] for row in scratch):
        raise ValueError('Private nonzero workspace unexpectedly depends on survivor count')
    count_bytes,reduce_bytes,select_bytes=scratch[0]
    return dict(cells=cells,count_requested_bytes=count_bytes,
        reduce_requested_bytes=reduce_bytes,select_requested_bytes=select_bytes,
        private_requested_bytes=sum(scratch[0]),private_rounded_bytes=sum(round_up(n) for n in scratch[0]),
        observed_coordinate_block_excess_bytes=output_excess,
        source='PyTorch 2.5.1 Nonzero.cu count allocation and overlapping count/select scratch',
        context=recorded,
        scope='Exact installed private allocation requests at this mask extent, with empty/single/dense controls. Coordinate output, mask, pinned count, allocator reservations and driver memory are separate.')


def device_selection_memory(samples, markers, traits, *, census, context, max_cells=DEVICE_SELECTION_MAX_CELLS):
    """Bound selector-local storages by two consecutive worst-case blocks.

    Statistics inputs/resident critical values belong to the native scan ledger.
    No occupancy assumption relaxes this bound. CUB requests must be captured
    at every distinct full/tail mask extent, not interpolated from a timing grid.
    """
    from .device_significance_work import device_significant_tensor_work
    for name,value in [('samples',samples),('markers',markers),('traits',traits),('max_cells',max_cells)]:
        _integer(name,value)
    if max_cells>DEVICE_SELECTION_MAX_CELLS:raise ValueError('Selection limit exceeds the production bound')
    width,height,_=device_selection_shape(markers,traits,max_cells)
    heights={min(markers,height)};widths={width}
    if markers%height:heights.add(markers%height)
    if traits%width:widths.add(traits%width)
    blocks=[];source={}
    for h in sorted(heights):
        for w in sorted(widths):
            cells=h*w
            work=device_significant_tensor_work(samples,h,w,[cells],max_cells=max_cells)
            workspace=device_nonzero_workspace(census,cells,context)
            blocks.append(dict(rows=h,traits=w,cells=cells,
                tensor_temporary_bytes=work['distinct_temporary_bytes'],workspace=workspace,
                private_rounded_bytes=workspace['private_rounded_bytes']))
            source.update(work['source_sha256'])
    maximum=max(row['cells'] for row in blocks)
    tensor=max(row['tensor_temporary_bytes'] for row in blocks)
    private=max(row['private_rounded_bytes'] for row in blocks)
    return dict(selection_gpu_bytes=2*(tensor+private),
        single_block_tensor_bytes=tensor,single_block_private_bytes=private,
        maximum_selection_cells=maximum,block_shapes=blocks,source_sha256=source,
        selected_payload_bytes_per_block=20*maximum,
        pinned_count_requested_bytes=4,pinned_count_allocator_bytes=4,
        scope='Conservative sum of every distinct source intermediate in a dense block, plus exact installed CUB requests; doubled for previous/current Python references and stream work. Native dense statistics, persistent critical values, host output queue, allocator reservations and driver memory are separate.')
