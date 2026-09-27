"""Typed independent primitive names for the device selector's Python calls."""
from collections import Counter


def device_selection_host_primitives(work):
    """Resolve API semantics rather than a workload size or measured duration."""
    rows = []
    for call in work['host_calls']:
        steps = [work['steps'][i] for i in call['step_indices']]
        ops = tuple(row['op'] for row in steps)
        name = call['name']
        primitive = None
        if name == 'cpu_numpy':
            if ops != ('aten._to_copy.default',): raise ValueError('Unknown host-copy trace')
            primitive = {'torch.float32':'copy_cpu_numpy_fp32', 'torch.int64':'copy_cpu_numpy_int64',
                         'torch.int32':'copy_cpu_numpy_int32'}[call['input_dtypes'][0]]
        elif ops == ('aten._to_copy.default',):
            # Casts differ only by their output type: df -> int64 table index,
            # coordinates -> int32 packed payload.
            primitive = {'torch.int64':'df_cast_int64', 'torch.int32':'coordinate_cast_int32'}.get(
                steps[0]['outputs'][0]['dtype'])
        elif name == 'nonzero':
            if ops != ('aten.nonzero.default',): raise ValueError('Unknown nonzero trace')
            primitive = 'nonzero_nonempty' if steps[0]['outputs'][0]['shape'][0] else 'nonzero_empty'
        elif name == '__getitem__':
            if ops == ('aten.index.Tensor',):
                primitive = 'gather_matrix' if len(call['input_shapes'][0]) == 2 else 'gather_vector'
            else:
                primitive = {
                    ('aten.alias.default',):'view_alias',
                    ('aten.slice.Tensor',):'view_slice',
                    ('aten.slice.Tensor','aten.slice.Tensor'):'view_two_slices',
                    ('aten.slice.Tensor','aten.unsqueeze.default'):'view_column',
                    ('aten.slice.Tensor','aten.select.int'):'view_coordinate',
                    ('aten.unsqueeze.default',):'view_unsqueeze',
                }.get(ops)
        else:
            primitive = {
                ('aten.clamp.default',):'clamp_int64',
                ('aten.eq.Scalar',):'status_equal_zero',
                ('aten.gt.Scalar',):'df_greater_zero',
                ('aten.bitwise_and.Tensor',):'and_bool',
                ('aten.abs.default','aten.ne.Scalar','aten.eq.Tensor','aten.mul.Tensor'):'isfinite_fp32',
                ('aten.abs.default',):'abs_fp32',
                ('aten.ge.Tensor',):'compare_cutoff',
                ('aten.add.Tensor',):'index_add_int64',
                ('aten.scalar_tensor.default','aten.where.self'):'where_limit',
                ('aten.new_empty.default',):'empty_mask',
                ('aten.ge.Tensor_out',):'compare_cutoff_out',
                ('aten.lt.Scalar',):'finite_upper',
                ('aten.bitwise_and_.Tensor',):'and_bool_inplace',
                ('aten.unbind.int',):'view_unbind',
                ('aten.view.dtype','aten.detach.default'):'view_dtype',
                ('aten.stack.default',):'stack_int32',
            }.get(ops)
        if primitive is None:
            raise ValueError('Unpriced device-selection Python call: '+str((name,ops)))
        rows.append(dict(call_id=call['id'], primitive=primitive, step_indices=call['step_indices']))
    return dict(calls=rows, counts=dict(Counter(row['primitive'] for row in rows)),
        scope='Source API grouping for independent fixed primitive measurements. Dynamic nonzero and blocking copy prices must preserve their synchronization boundaries.')


def device_selection_graph(work, *, host_prices, operation_services, nonzero_services,
                           transfer_prices, yield_cpu_seconds, cpu_fraction=1.,
                           host_serial_fraction=1., wait_cpu_fraction=0.,
                           capacities=None):
    """Source-ordered CPU/CUDA graph with explicit count and payload barriers.

    Services are independently supplied components, not complete-selector API
    timings. Blocking host prices must separate pre-barrier and post-barrier
    CPU work; a whole blocking API CPU/wall measurement cannot be used twice.
    Return yield/resume ports so the existing indexed writer can pause Python
    without draining pending CUDA work from an empty nonzero result.
    """
    import math
    from .execution_graph import ExecutionGraph
    from .first_principles import positive

    for name, value in [('cpu_fraction',cpu_fraction),('host_serial_fraction',host_serial_fraction),
                        ('wait_cpu_fraction',wait_cpu_fraction)]:
        if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('Invalid '+name)
    if not cpu_fraction: raise ValueError('Positive CPU capacity required')
    positive('yield CPU',yield_cpu_seconds,True)
    resolved = device_selection_host_primitives(work)
    by_call = {row['call_id']:row for row in resolved['calls']}
    calls = {call['id']:call for call in work['host_calls']}
    expected_ops = {str(i) for i,step in enumerate(work['steps'])
        if not (step['alias_only'] or step['allocation_only'] or step['host_copy'] or step['op']=='aten.nonzero.default')}
    if set(operation_services) != expected_ops:
        raise ValueError('One independently priced service per active non-copy tensor step required')
    if len(nonzero_services) != len(work['blocks']):
        raise ValueError('One count/select service group per source block required')
    if set(transfer_prices) != {'count','payload'}:
        raise ValueError('Separate pinned-count and pageable-payload transfer prices required')
    for kind, price in transfer_prices.items():
        if set(price) != {'latency_seconds','bytes_per_second','resources'}:
            raise ValueError('Transfer service requires explicit latency, capacity and resource names')
        positive(kind+' latency',price['latency_seconds'],True)
        positive(kind+' capacity',price['bytes_per_second'])
        if not isinstance(price['resources'],(list,tuple)) or len(set(price['resources']))!=len(price['resources']):
            raise ValueError('Unique transfer resource names required')
    graph = ExecutionGraph(); graph.capacities = dict(capacities or {})
    graph.capacities.setdefault('cpu',cpu_fraction)
    graph.capacities.setdefault('host_serial',1.)
    host = gpu = graph.add('begin')
    yielded = []; resumed = []; gpu_tails = []; used_calls = set()
    cpu_total = 0.; payload_bytes = 0; count_bytes = 0

    def cpu_node(name, seconds, after):
        nonlocal cpu_total
        positive(name,seconds,True); cpu_total += seconds
        return graph.add(name,seconds/cpu_fraction,[after],
            dict(cpu=cpu_fraction,host_serial=cpu_fraction*host_serial_fraction))

    def gpu_node(name, row, after):
        if not isinstance(row,dict) or set(row)-{'seconds','resources','op','phase'} or 'seconds' not in row:
            raise ValueError('Explicit independently priced GPU service required')
        return graph.add(name,row['seconds'],after,row.get('resources'))

    def transfer(name, size, kind, after):
        price = transfer_prices[kind]
        seconds = price['latency_seconds']+size/price['bytes_per_second']
        resources = {key:size/seconds for key in price['resources']}
        return graph.add(name,seconds,after,resources)

    def wait(name, submitted, completed):
        if wait_cpu_fraction:
            graph.resource_waits[name+':active'] = dict(after=submitted,until=completed,
                resources=dict(cpu=cpu_fraction*wait_cpu_fraction))
        return graph.add(name,after=[submitted,completed])

    for index, block in enumerate(work['blocks']):
        lo, hi = block['step_range']
        call_ids = list(dict.fromkeys(step['host_call_id'] for step in work['steps'][lo:hi]))
        for call_id in call_ids:
            if call_id in used_calls: raise ValueError('Python API call crosses a source yield')
            used_calls.add(call_id)
            call = calls[call_id]; primitive = by_call[call_id]['primitive']
            price = host_prices.get(primitive)
            if not isinstance(price,dict): raise ValueError('Missing independent host primitive '+primitive)
            prefix = 'call:'+str(call_id)
            blocking = primitive.startswith(('copy_cpu_numpy_','nonzero_'))
            expected = {'before_cpu_seconds','after_cpu_seconds'} if blocking else {'cpu_seconds'}
            if set(price) != expected:
                raise ValueError('Blocking CPU work must be partitioned around the barrier: '+primitive)
            host = cpu_node(prefix+':before',price['before_cpu_seconds'] if blocking else price['cpu_seconds'],host)
            if primitive.startswith('nonzero_'):
                phases = nonzero_services[index]
                names = ['count','flagged_select']+(['coordinate_scatter'] if block['retained'] else [])
                if set(phases)!=set(names): raise ValueError('Nonzero services differ from count/select/scatter source phases')
                gpu = gpu_node(prefix+':count',phases['count'],[gpu,host])
                gpu = transfer(prefix+':count_d2h',4,'count',[gpu,host]); count_bytes += 4
                host = wait(prefix+':count_ready',host,gpu)
                host = cpu_node(prefix+':after',price['after_cpu_seconds'],host)
                for phase in names[1:]: gpu = gpu_node(prefix+':'+phase,phases[phase],[gpu,host])
            elif primitive.startswith('copy_cpu_numpy_'):
                size = work['steps'][call['step_indices'][0]]['inputs'][0]['bytes']
                gpu = transfer(prefix+':payload_d2h',size,'payload',[gpu,host]); payload_bytes += size
                host = wait(prefix+':payload_ready',host,gpu)
                host = cpu_node(prefix+':after',price['after_cpu_seconds'],host)
            else:
                for step_index in call['step_indices']:
                    step = work['steps'][step_index]
                    if step['alias_only'] or step['allocation_only']: continue
                    service = operation_services[str(step_index)]
                    if service.get('op') != step['op']: raise ValueError('GPU service source-operation mismatch')
                    gpu = gpu_node(prefix+':step:'+str(step_index),service,[gpu,host])
        host = cpu_node(f'block:{index}:yield',yield_cpu_seconds,host)
        yielded.append(host); gpu_tails.append(gpu)
        host = graph.add(f'block:{index}:resume',after=[host]); resumed.append(host)
    if used_calls != set(calls): raise ValueError('Unassigned source Python API calls')
    if payload_bytes != work['selected_payload_d2h_bytes'] or count_bytes != work['nonzero_count_d2h_bytes']:
        raise ValueError('Source transfer conservation failed')
    graph.add('drained',after=[host,gpu])
    return dict(graph=graph,yield_nodes=yielded,resume_nodes=resumed,gpu_tail=gpu,
        block_gpu_tails=gpu_tails,retained_per_block=[b['retained'] for b in work['blocks']],
        selected_payload_d2h_bytes=payload_bytes,nonzero_count_d2h_bytes=count_bytes,
        host_cpu_seconds=cpu_total,source_sha256=work['source_sha256'],prediction_complete=False,
        scope='Conditional independent-service source graph. Empty yields need not wait for flagged selection; nonempty payloads do. Resume ports belong to the queue/writer. Count and payload CPU waits consume only the explicit wait scenario.')
