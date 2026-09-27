"""Independent fixed API primitives; retain blocking and loaded-host boundaries."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import threading
import time

import torch
from torch.utils._python_dispatch import TorchDispatchMode
from torch.utils._pytree import tree_flatten
from direct_device_significance_primitives import build_bank
from torchgwas.detailed_calibration import source_identity, sha256_file
from torchgwas.geometry_collection import write_record


def telemetry():
    command = ['nvidia-smi', '--query-gpu=index,uuid,utilization.gpu,memory.used,power.draw,clocks.sm,clocks.mem', '--format=csv,noheader,nounits']
    result = subprocess.run(command, text=True, capture_output=True, timeout=20)
    return dict(utc=datetime.now(timezone.utc).isoformat(), command=command,
        returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)


def bulk_bank(device):
    # Fixed 1M elements, not selected from a GWAS candidate or timing result.
    fp = torch.ones(1 << 20, device=device, dtype=torch.float32)
    integer = torch.arange(1 << 20, device=device, dtype=torch.int64)
    status = torch.zeros(1 << 20, device=device, dtype=torch.uint8)
    empty = torch.zeros((1024,1024), device=device, dtype=torch.bool)
    sparse = empty.clone(); sparse[0,0] = True
    dense = torch.ones_like(empty)
    return {
        'loop_control': lambda: None,
        'copy_cpu_numpy_fp32': lambda: fp.cpu().numpy(),
        'copy_cpu_numpy_int64': lambda: integer.cpu().numpy(),
        'copy_cpu_numpy_uint8': lambda: status.cpu().numpy(),
        'nonzero_empty': lambda: empty.nonzero(as_tuple=False),
        'nonzero_single': lambda: sparse.nonzero(as_tuple=False),
        'nonzero_dense': lambda: dense.nonzero(as_tuple=False),
    }


class Signature(TorchDispatchMode):
    def __init__(self): self.steps = []
    def __torch_dispatch__(self, function, types, args=(), kwargs=None):
        def tensors(values):
            return [dict(shape=list(v.shape), stride=list(v.stride()), dtype=str(v.dtype), device=v.device.type)
                for v in tree_flatten(values)[0] if isinstance(v, torch.Tensor)]
        kwargs = kwargs or {}
        before = tensors((args,kwargs))
        output = function(*args,**kwargs)
        self.steps.append(dict(op=str(function), inputs=before, outputs=tensors(output)))
        return output


def observe(function, count, device):
    torch.cuda.synchronize(device)
    wall = time.perf_counter_ns(); thread = time.thread_time_ns(); process = time.process_time_ns()
    for _ in range(count):
        result = function()
        del result
    process_elapsed = time.process_time_ns()-process
    thread_elapsed = time.thread_time_ns()-thread
    wall_elapsed = time.perf_counter_ns()-wall
    drain_wall = time.perf_counter_ns(); drain_cpu = time.thread_time_ns()
    torch.cuda.synchronize(device)
    drain_cpu = time.thread_time_ns()-drain_cpu
    drain_wall = time.perf_counter_ns()-drain_wall
    return dict(calls=count, api_return_wall_seconds=wall_elapsed*1e-9/count,
        caller_thread_cpu_seconds=thread_elapsed*1e-9/count,
        diagnostic_process_cpu_seconds=process_elapsed*1e-9/count,
        external_drain_wall_seconds=drain_wall*1e-9/count,
        external_drain_thread_cpu_seconds=drain_cpu*1e-9/count)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    parser.add_argument('--devices', nargs='+', default=['cuda:0','cuda:2'])
    args = parser.parse_args()
    if len(args.devices) != 2 or len(set(args.devices)) != 2:
        raise ValueError('This protocol declares exactly two distinct devices')
    root = Path(args.out); root.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2); torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.set_float32_matmul_precision('highest')
    source = source_identity()
    paths = [Path(__file__),Path(__file__).with_name('direct_device_significance_primitives.py')]
    for path in paths: (root/path.name).write_bytes(path.read_bytes())
    hashes = {p.name:sha256_file(p) for p in paths}
    limit = 512 << 20
    banks = {}; signatures = {}; properties = {}
    for device in args.devices:
        torch.cuda.set_device(device)
        prop = torch.cuda.get_device_properties(device)
        torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,device)
        if torch.cuda.mem_get_info(device)[0] < limit: raise ValueError('Insufficient free GPU memory')
        properties[device] = dict(name=prop.name, uuid=str(prop.uuid), total_memory=prop.total_memory,
            multiprocessors=prop.multi_processor_count, capability=[prop.major,prop.minor])
        banks[device] = {'tiny':build_bank(device), 'bulk':bulk_bank(device)}
        signatures[device] = {}
        for scale, bank in banks[device].items():
            for name,function in bank.items():
                for _ in range(8): function()
                torch.cuda.synchronize(device)
                signature = Signature()
                with signature: function()
                torch.cuda.synchronize(device)
                signatures[device][scale+':'+name] = signature.steps
    groups = [[args.devices[0]],[args.devices[1]],args.devices]
    plan = []
    rng = random.Random(9211829)
    for index, devices in enumerate(groups):
        rounds = []
        for repeat in range(9):
            orders = {}
            for device in devices:
                order = [(scale,name,128 if scale=='tiny' else 8) for scale,bank in banks[device].items() for name in bank]
                rng.shuffle(order); orders[device] = order
            rounds.append(dict(repeat=repeat,orders=orders))
        plan.append(dict(group=index,devices=devices,rounds=rounds))
    context = dict(torch_version=torch.__version__, cuda_runtime=torch.version.cuda,
        devices=properties, affinity=sorted(os.sched_getaffinity(0)),
        torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),
        python=os.sys.version, memory_limit_per_device_bytes=limit,
        environment={k:v for k,v in sorted(os.environ.items()) if k.startswith(('CUDA_','TORCH','OMP_','MKL_','OPENBLAS_','GOMP_','PYTORCH_'))})
    protocol = dict(source_sha256=source,harness_sha256=hashes,context=context,plan=plan,signatures=signatures,
        boundary='Warm API calls including immediate output destruction. Every batch starts drained. External final GPU drain is recorded separately. Blocking copy and nonzero count barriers remain inside the API boundary.',
        cpu_boundary='Caller thread CPU measured directly. Process CPU is diagnostic and overlaps across concurrent workers. Neither CPU nor wall differences are clipped or converted into GIL/kernel/transfer prices.',
        primitive_extents='Tiny tensors use 32 elements or 32x32. Bulk tensors always use 1048576 elements. No candidate-shaped timings or association outputs are used.',
        limitations=['Warm allocator and library state', 'Pending asynchronous GPU work can delay a later API call',
            'Thread CPU includes busy driver waiting; it is not pure dispatch', 'Python iteration and timers are included; loop control is reported without subtraction',
            'Other GPU/CPU work on the shared host is not controlled', 'These measurements do not partition GIL service or CUB kernel stages'],
        transfer_qualified=False)
    write_record(root/'protocol.json',protocol)
    observations = []; states = []; gc.collect(); gc.disable()
    try:
        for group in plan:
            devices = group['devices']
            for round_ in group['rounds']:
                states.append(dict(group=group['group'],repeat=round_['repeat'],phase='before',telemetry=telemetry()))
                barrier = threading.Barrier(len(devices))
                def worker(device):
                    torch.cuda.set_device(device)
                    torch.cuda.synchronize(device)
                    barrier.wait(timeout=30)
                    rows = []
                    for scale,name,count in round_['orders'][device]:
                        row = observe(banks[device][scale][name],count,device)
                        row.update(group=group['group'],active_devices=devices,repeat=round_['repeat'],device=device,scale=scale,primitive=name)
                        rows.append(row)
                    return rows
                with ThreadPoolExecutor(max_workers=len(devices)) as pool:
                    futures = [pool.submit(worker,device) for device in devices]
                    rows = [row for future in futures for row in future.result()]
                observations.extend(rows)
                write_record(root/f"group_{group['group']}_repeat_{round_['repeat']}.json",dict(observations=rows))
                states.append(dict(group=group['group'],repeat=round_['repeat'],phase='after',telemetry=telemetry()))
                print(json.dumps(dict(group=group['group'],devices=devices,repeat=round_['repeat'],observations=len(rows))),flush=True)
    finally:
        gc.enable()
        write_record(root/'telemetry.json',states)
    if source != source_identity() or hashes != {p.name:sha256_file(p) for p in paths}:
        raise ValueError('Source changed during primitive collection')
    summaries = []
    keys = sorted({(r['group'],r['device'],r['scale'],r['primitive']) for r in observations})
    metrics = ['api_return_wall_seconds','caller_thread_cpu_seconds','external_drain_wall_seconds','external_drain_thread_cpu_seconds']
    for key in keys:
        rows = [r for r in observations if (r['group'],r['device'],r['scale'],r['primitive']) == key]
        values = {metric:dict(median=statistics.median(r[metric] for r in rows),minimum=min(r[metric] for r in rows),maximum=max(r[metric] for r in rows)) for metric in metrics}
        summaries.append(dict(group=key[0],device=key[1],scale=key[2],primitive=key[3],repetitions=len(rows),metrics=values))
    write_record(root/'report.json',dict(protocol_sha256=sha256_file(root/'protocol.json'),observations=observations,summaries=summaries,
        transfer_qualified=False,scope='Independent raw fixed-primitive observations. No full scan measurements, fitted candidate rankings, or automatic service decomposition.'))
    print(json.dumps(dict(complete=True,observations=len(observations),primitive_contexts=len(summaries))),flush=True)


if __name__ == '__main__': main()
