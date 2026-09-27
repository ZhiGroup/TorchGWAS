"""Source-derived device selection work with explicit dynamic survivor shapes."""
import ast
import hashlib
import inspect
from pathlib import Path
import textwrap

from .reduced_output_work import significant_output_work
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS, device_selection_shape

NONZERO_SOURCE = 'https://github.com/pytorch/pytorch/blob/v2.5.1/aten/src/ATen/native/cuda/Nonzero.cu'


def cuda_nonzero_work(cells, retained, *, torch_version):
    """Explicit bool-mask traffic around PyTorch 2.5.1's CUB calls.

    The count reduction precedes a blocking four-byte D2H. Flagged selection
    still scans the mask when the count is zero. Coordinate scatter follows
    selection for nonempty two-dimensional masks. CUB internal passes and
    workspace are not inferred from these wrapper-level accesses.
    """
    from .mechanistic_plan import _integer
    _integer('cells', cells); _integer('retained', retained, 0)
    if retained > cells: raise ValueError('Retained count exceeds mask cells')
    if cells >= (1 << 31)-1: raise ValueError('CUDA nonzero requires fewer than INT_MAX cells')
    if str(torch_version).split('+')[0] != '2.5.1':
        raise ValueError('CUDA nonzero source accounting requires PyTorch 2.5.1')
    phases = [
        dict(phase='count', gpu_read_bytes=cells, gpu_write_bytes=4,
             predicate_evaluations=cells, sum_additions=cells-1),
        dict(phase='count_to_host', gpu_read_bytes=4, gpu_write_bytes=0, d2h_bytes=4,
             blocking=True),
        dict(phase='flagged_select', gpu_read_bytes=cells, gpu_write_bytes=8*retained+4,
             predicate_evaluations=cells, counting_iterator_materialized_bytes=0),
    ]
    if retained:
        phases.append(dict(phase='coordinate_scatter', gpu_read_bytes=8*retained,
            gpu_write_bytes=16*retained, integer_divisions=2*retained,
            integer_remainders=2*retained, integer_multiplications=2*retained,
            grid=[(retained+255)//256, 1, 1], block=[256, 1, 1]))
    return dict(cells=cells, retained=retained, phases=phases, count_d2h_bytes=4,
        output_bytes=16*retained, device_count_allocation_bytes=4, pinned_count_allocation_bytes=4,
        flat_indices_alias_output=True,
        gpu_logical_bytes=sum(p['gpu_read_bytes']+p['gpu_write_bytes'] for p in phases),
        synchronization='Count transfer is blocking; selection and coordinate scatter follow that barrier.',
        source=NONZERO_SOURCE, source_version='2.5.1',
        unpriced_terms=['CUB count/select internal traffic, prefix arithmetic, launch policy and workspace',
            'Allocator rounding, cache reuse, pinned-count lifecycle and CPU dispatch',
            'Integer source operators are not a count of issued GPU instructions'],
        scope='Source-level accesses for one contiguous two-dimensional bool mask. '
              'No observed timing, calibrated service or complete memory bound.')


def device_significant_tensor_work(samples, markers, traits, retained_per_block,
        *, max_cells=DEVICE_SELECTION_MAX_CELLS, max_blocks=1024):
    """Trace the production selector on metadata, with caller-supplied counts.

    Two boundaries cannot execute on meta tensors: dynamic nonzero extents and
    blocking CPU/NumPy copies. Substitute declared nonzero shapes and record
    copy requests at their original source positions. All comparisons, gathers,
    slicing, loop bounds and index arithmetic execute from the source function.
    The GPU census must separately verify these substitutions and the compiled
    kernel work. No real phenotype/statistics allocation or timing is used.
    """
    import torch
    from torch.utils._python_dispatch import TorchDispatchMode
    from torch.overrides import TorchFunctionMode
    from torch.utils._pytree import tree_flatten
    from .reduce import device_significant_pairs

    ledger = significant_output_work(samples, markers, traits, markers,
        backend='device', retained_per_block=retained_per_block,
        max_selection_cells=max_cells, max_blocks=max_blocks)
    if ledger['retained_pairs'] is None:
        raise ValueError('Explicit per-block survivor counts required for meta selection')
    source = Path(inspect.getsourcefile(device_significant_pairs))
    geometry_source = Path(inspect.getsourcefile(device_selection_shape.__wrapped__))
    syntax = ast.parse(textwrap.dedent(inspect.getsource(device_significant_pairs)))
    copies = []
    steps = []
    host_calls = []
    active_call = []
    storage_ids = {}
    keep = []
    block_index = 0
    nonzero_index = 0

    def describe(value):
        storage = value.untyped_storage(); key = storage._cdata
        if key not in storage_ids: storage_ids[key] = len(storage_ids)
        return dict(storage=storage_ids[key], shape=list(value.shape), stride=list(value.stride()),
            dtype=str(value.dtype), device_type='cuda', bytes=value.numel()*value.element_size(),
            storage_bytes=storage.nbytes(), offset_bytes=value.storage_offset()*value.element_size())

    class HostRecord(TorchFunctionMode):
        def __torch_function__(self, function, types, args=(), kwargs=None):
            kwargs = kwargs or {}
            inputs = [v for v in tree_flatten((args, kwargs))[0] if isinstance(v, torch.Tensor)]
            call_id = len(host_calls)
            host_calls.append(dict(id=call_id, name=getattr(function, '__name__', str(function)),
                input_dtypes=[str(v.dtype) for v in inputs], input_shapes=[list(v.shape) for v in inputs],
                kwargs={key:str(value) for key,value in kwargs.items() if not isinstance(value, torch.Tensor)},
                step_indices=[]))
            active_call.append(call_id)
            try:
                return function(*args, **kwargs)
            finally:
                active_call.pop()

    def host_copy(value):
        before = describe(value)
        call_id = len(host_calls)
        host_calls.append(dict(id=call_id, name='cpu_numpy', input_dtypes=[before['dtype']],
            input_shapes=[before['shape']], kwargs={}, step_indices=[]))
        after = dict(before, device_type='cpu', storage=None, stride=[1], offset_bytes=0,
                     storage_bytes=before['bytes'])
        row = dict(op='aten._to_copy.default', inputs=[before], outputs=[after],
            alias_only=False, allocation_only=False, host_copy=True, block=block_index, host_call_id=call_id,
            read_bytes=before['bytes'], write_bytes=after['bytes'],
            logical_bytes=before['bytes']+after['bytes'], gpu_logical_bytes=before['bytes'])
        steps.append(row); copies.append(dict(block=block_index, bytes=before['bytes'], dtype=before['dtype']))
        keep.append(value)
        return value  # Only yielded; no later tensor computation consumes this host result.

    class ReplaceCopies(ast.NodeTransformer):
        def __init__(self): self.count = 0
        def visit_Call(self, node):
            node = self.generic_visit(node)
            if (isinstance(node.func, ast.Attribute) and node.func.attr == 'numpy'
                    and not node.args and not node.keywords
                    and isinstance(node.func.value, ast.Call)):
                call = node.func.value
                if (isinstance(call.func, ast.Attribute) and call.func.attr == 'cpu'
                        and not call.args and not call.keywords):
                    self.count += 1
                    return ast.copy_location(ast.Call(func=ast.Name(id='_trace_host_copy', ctx=ast.Load()),
                        args=[call.func.value], keywords=[]), node)
            return node

    replacement = ReplaceCopies(); syntax = replacement.visit(syntax)
    if replacement.count != 1:
        raise ValueError('Device selector host-copy source contract changed')
    ast.fix_missing_locations(syntax)
    # The packed copy's host split is NumPy work after the recorded copy. On
    # meta tensors it yields five same-length placeholders and no tensor ops.
    namespace = dict(device_significant_pairs.__globals__, _trace_host_copy=host_copy,
                     _owned_pairs=lambda packed, first_variant, first_trait: (range(packed.shape[1]),)*5)
    exec(compile(syntax, str(source), 'exec'), namespace)
    selector = namespace[device_significant_pairs.__name__]

    class Record(TorchDispatchMode):
        def __torch_dispatch__(self, func, types, args=(), kwargs=None):
            nonlocal nonzero_index
            kwargs = kwargs or {}
            inputs = [v for v in tree_flatten((args, kwargs))[0] if isinstance(v, torch.Tensor)]
            before = [describe(v) for v in inputs]
            name = str(func)
            if name == 'aten.nonzero.default':
                if nonzero_index >= len(ledger['blocks']):
                    raise ValueError('Source emits more nonzero calls than its block contract')
                count = ledger['blocks'][nonzero_index]['retained']; nonzero_index += 1
                result = torch.empty_strided((count, 2), (1, max(1, count)), dtype=torch.int64, device='meta')
            else:
                result = func(*args, **kwargs)
            outputs = [v for v in tree_flatten(result)[0] if isinstance(v, torch.Tensor)]
            after = [describe(v) for v in outputs]
            keep.extend(inputs); keep.extend(outputs)
            alias = bool(after) and all(v['storage'] in {r['storage'] for r in before} for v in after) and not func._schema.is_mutable
            allocate = name.startswith(('aten.empty', 'aten.new_empty'))
            reads = 0 if alias or allocate else sum(v['bytes'] for v in {r['storage']:r for r in before}.values())
            writes = 0 if alias or allocate else sum(v['bytes'] for v in {r['storage']:r for r in after}.values())
            if name == 'aten.index.Tensor':
                # Gather reads selected source values, not the entire base array.
                reads = sum(v['bytes'] for v in after) + sum(v['bytes'] for v in before[1:])
            steps.append(dict(op=name, inputs=before, outputs=after, alias_only=alias,
                allocation_only=allocate, host_copy=False, block=block_index,
                host_call_id=active_call[-1] if active_call else None,
                read_bytes=reads, write_bytes=writes, logical_bytes=reads+writes, gpu_logical_bytes=reads+writes))
            return result

    beta = torch.empty((markers, traits), dtype=torch.float32, device='meta')
    values = torch.empty_like(beta)
    status = torch.empty(markers, dtype=torch.uint8, device='meta')
    df = torch.empty(markers, dtype=torch.float32, device='meta')
    critical = torch.empty(samples+1, dtype=torch.float32, device='meta')
    initial = [describe(v) for v in [beta, values, status, df, critical]]
    blocks = []
    with HostRecord(), Record():
        generator = selector(beta, values, status, df, critical, max_cells=max_cells)
        for block_index, block in enumerate(ledger['blocks']):
            first = len(steps)
            row = next(generator)
            if list(row[:2]) != block['variant_range'] or any(len(v) != block['retained'] for v in row[2:]):
                raise ValueError('Meta output differs from declared source block')
            backend = cuda_nonzero_work(block['cells'], block['retained'], torch_version=torch.__version__)
            blocks.append(dict(block, step_range=[first, len(steps)], nonzero_backend_work=backend))
        if next(generator, None) is not None:
            raise ValueError('Source emits additional selection blocks')
    if nonzero_index != len(blocks): raise ValueError('Incomplete source nonzero census')
    for index, step in enumerate(steps):
        if step['host_call_id'] is None:
            raise ValueError('Source tensor operation lacks its Python API call')
        host_calls[step['host_call_id']]['step_indices'].append(index)
    active_host_calls = [call for call in host_calls if call['step_indices']]
    initial_ids = {v['storage'] for v in initial}
    temporary = {}
    for step in steps:
        for value in step['inputs'] + step['outputs']:
            key = value['storage']
            if key is not None and key not in initial_ids:
                temporary[key] = max(temporary.get(key, 0), value['storage_bytes'])
    return dict(samples=samples, markers=markers, traits=traits, max_cells=max_cells,
        blocks=blocks, steps=steps, initial_storages=initial, host_calls=active_host_calls,
        distinct_temporary_bytes=sum(((v+511)//512)*512 for v in temporary.values()),
        source_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (source, geometry_source)},
        torch_version=torch.__version__, selection_blocks=len(blocks),
        nonzero_calls=nonzero_index, selected_copy_calls=len(copies), selected_copies=copies,
        selected_payload_d2h_bytes=sum(row['bytes'] for row in copies),
        gpu_logical_bytes=sum(row['gpu_logical_bytes'] for row in steps),
        nonzero_count_d2h_bytes=sum(block['nonzero_backend_work']['count_d2h_bytes'] for block in blocks),
        selector_d2h_bytes=sum(row['bytes'] for row in copies)+4*len(blocks),
        gpu_logical_bytes_with_nonzero_source=sum(row['gpu_logical_bytes'] for row in steps if row['op'] != 'aten.nonzero.default')
            +sum(block['nonzero_backend_work']['gpu_logical_bytes'] for block in blocks),
        retained_pairs=ledger['retained_pairs'], prediction_complete=False,
        unpriced_terms=['CUB internal count/select passes, arithmetic, workspace and timing across the known count barrier',
            'Compiled kernel geometry, instruction throughput and physical HBM/cache traffic',
            'CPU dispatch, allocation, NumPy wrapping and blocking pageable-copy service',
            'Critical table preparation/upload and the native source status transfer',
            'Queue, output serialization, durable storage, cleanup and allocator retention'],
        scope='Exact source loops with declared dynamic survivor shapes and recorded blocking host-copy requests. '
              'Logical gather accesses count selected values; nonzero backend internals remain separate. '
              'All-intermediate storage sum is conservative and excludes vendor workspace; not a peak-memory admission or timing qualification.')
