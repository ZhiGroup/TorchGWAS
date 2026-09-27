"""Untimed operation and transfer census of bounded device pair selection."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
import torch
from torch.utils._python_dispatch import TorchDispatchMode
from torch.utils._pytree import tree_flatten

from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import kernel_census, write_record
from torchgwas.reduce import device_significant_pairs
from torchgwas.selection_geometry import DEVICE_SELECTION_MAX_CELLS

parser = argparse.ArgumentParser()
parser.add_argument('--out', required=True)
parser.add_argument('--device', default='cuda:0')
# B,K,max_cells; max_cells 'max' is the production whole-chunk block.
parser.add_argument('--shapes', nargs='+', default=['13,7,17', '5,37,11', '257,4093,1048576'])
parser.add_argument('--memory-limit-mib', type=int, default=512)
args = parser.parse_args()
shapes = [tuple(int(v) if v != 'max' else DEVICE_SELECTION_MAX_CELLS for v in shape.split(',')) for shape in args.shapes]
root = Path(args.out); root.mkdir(parents=True, exist_ok=False)
device = args.device
limit = args.memory_limit_mib << 20
torch.cuda.set_device(device)
torch.cuda.set_per_process_memory_fraction(limit / torch.cuda.get_device_properties(device).total_memory, device)
if torch.cuda.mem_get_info(device)[0] < limit:
    raise ValueError('Insufficient free memory for bounded synthetic census')
torch.set_num_threads(2); torch.set_num_interop_threads(1)
source = source_identity()
prop=torch.cuda.get_device_properties(device)
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device=device,
    device_uuid=str(prop.uuid),name=prop.name,compute_capability=[prop.major,prop.minor],
    sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
(root / 'harness.py').write_bytes(Path(__file__).read_bytes())


def describe(value):
    return dict(shape=list(value.shape), stride=list(value.stride()),
        dtype=str(value.dtype), device_type=value.device.type,
        bytes=value.numel()*value.element_size())


class Observe(TorchDispatchMode):
    def __init__(self): self.steps = []
    def __torch_dispatch__(self, func, types, args=(), kwargs=None):
        kwargs = kwargs or {}
        before = [describe(v) for v in tree_flatten((args, kwargs))[0] if isinstance(v, torch.Tensor)]
        result = func(*args, **kwargs)
        after = [describe(v) for v in tree_flatten(result)[0] if isinstance(v, torch.Tensor)]
        self.steps.append(dict(op=str(func), inputs=before, outputs=after))
        return result


rows = []
for b, k, cells in shapes:
    for mode in ['empty', 'sparse', 'dense', 'invalid']:
        ids = np.arange(b*k, dtype=np.int64).reshape(b, k)
        beta = (ids % 1009).astype(np.float32)
        values = np.zeros((b, k), np.float32)
        if mode == 'sparse': values.flat[::17] = 2.
        elif mode in ('dense', 'invalid'): values.fill(2.)
        status = np.zeros(b, np.uint8)
        df = np.full(b, 38., np.float32)
        if mode == 'invalid':
            status[1::5] = 1; status[2::5] = 2
            df[::7] = 0.
            values.flat[::13] = np.nan
        table = np.ones(41, np.float32); table[0] = np.inf
        expected_indices = np.nonzero(np.isfinite(values) & (np.abs(values) >= table[df.astype(np.int64), None])
            & (status[:, None] == 0) & (df[:, None] > 0))
        expected = [expected_indices[0]+19, expected_indices[1], beta[expected_indices],
                    values[expected_indices], df[expected_indices[0]]]
        tensors = [torch.as_tensor(v, device=device) for v in [beta, values, status, df, table]]
        def run():
            return list(device_significant_pairs(*tensors, start=19, max_cells=cells))
        outputs = run(); torch.cuda.synchronize(device)
        actual = [np.concatenate([part[i] for part in outputs]) for i in range(2, 7)]
        order = np.lexsort((actual[1], actual[0]))
        for observed, target in zip(actual, expected):
            np.testing.assert_array_equal(observed[order], target)
        blocks = [dict(variant_range=list(part[:2]), retained=len(part[2]),
            payload_bytes=sum(v.nbytes for v in part[2:])) for part in outputs]
        del outputs, actual, order, expected
        observer = Observe()
        with observer: outputs = run()
        torch.cuda.synchronize(device); del outputs
        gc.collect(); torch.cuda.synchronize(device)
        baseline = torch.cuda.memory_allocated(device)
        torch.cuda.reset_peak_memory_stats(device)
        with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,
                torch.profiler.ProfilerActivity.CUDA]) as profiler:
            outputs = run(); torch.cuda.synchronize(device)
        peak = torch.cuda.max_memory_allocated(device)
        del outputs
        with tempfile.TemporaryDirectory(prefix='device-selection-census-') as temporary:
            trace = Path(temporary) / 'trace.json'; profiler.export_chrome_trace(str(trace))
            events = json.loads(trace.read_text())['traceEvents']
            kernels = kernel_census(events)
            runtime = [e['name'] for e in events if e.get('cat') == 'cuda_runtime']
            transfers = [dict(name=e['name'], bytes=e.get('args', {}).get('bytes'))
                for e in events if e.get('cat') in ('gpu_memcpy', 'gpu_memset')]
        row = dict(N=40, B=b, K=k, max_cells=cells, mode=mode, blocks=blocks,
            observed_tensor_steps=observer.steps, kernels=kernels, cuda_runtime_calls=runtime,
            transfer_events=transfers, extra_allocated_bytes=peak-baseline, durations_recorded=False,
            independent_cpu_arrays_equal=True)
        rows.append(row)
        write_record(root / f'{b}_{k}_{cells}_{mode}.json', row)
        print(json.dumps(dict(B=b, K=k, max_cells=cells, mode=mode, blocks=len(blocks),
            retained=sum(v['retained'] for v in blocks), tensor_operations=len(observer.steps),
            kernels=len(kernels), transfer_events=len(transfers))), flush=True)
        del tensors, profiler, observer, row, events, kernels, transfers, runtime
        gc.collect(); torch.cuda.empty_cache()
if source != source_identity(): raise ValueError('Source changed during device census')
write_record(root / 'census.json', dict(source_sha256=source,
    harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    torch_version=torch.__version__, device=torch.cuda.get_device_properties(device).name,context=context,
    memory_probe_limit_bytes=limit, durations_recorded=False, rows=rows,
    scope='Bounded synthetic selector operations, dynamic output shapes, CUDA launch and transfer census. '
          'Independent complete CPU reference comparisons; no timing prices, GWAS durations or selection qualification.'))
