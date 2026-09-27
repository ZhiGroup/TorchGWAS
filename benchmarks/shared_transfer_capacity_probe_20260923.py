"""Diagnostic pinned H2D/D2H single- and multi-GPU transfer observations.

This writes raw wall and CUDA-event spans, not a calibrated capacity profile.
Selected GPUs must look idle at launch; other server work is recorded rather
than silently treated as an uncontended hardware ceiling.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import threading
import time

import torch


def gpu_snapshot():
    fields='index,uuid,pci.bus_id,utilization.gpu,memory.used'
    output=subprocess.check_output(['nvidia-smi','--query-gpu='+fields,
        '--format=csv,noheader,nounits'],text=True)
    rows={}
    for line in output.splitlines():
        index,uuid,bus,util,memory=(field.strip() for field in line.split(','))
        rows[int(index)]=dict(uuid=uuid.lower(),pci_bus_id=bus.lower(),
            utilization_percent=int(util),memory_used_mib=int(memory))
    return rows


def _uuid(value):
    value=str(value).lower()
    return value.removeprefix('gpu-')


def _buffers(devices,size):
    result={}
    for device in devices:
        with torch.cuda.device(device):
            host=torch.empty(size,dtype=torch.uint8,pin_memory=True)
            host.fill_(int(device)%251)
            gpu=torch.empty(size,dtype=torch.uint8,device=device)
            gpu.fill_(int(device)%251)
            stream=torch.cuda.Stream(device=device)
            torch.cuda.synchronize(device)
        result[device]=(host,gpu,stream)
    return result


def observe(group,buffers,*,direction,size,copies,pool):
    gate=threading.Barrier(len(group)+1)
    def copy(device):
        host,gpu,stream=buffers[device]
        src,dst=(host,gpu) if direction=='h2d' else (gpu,host)
        with torch.cuda.device(device),torch.cuda.stream(stream):
            begin=torch.cuda.Event(enable_timing=True)
            end=torch.cuda.Event(enable_timing=True)
            gate.wait()
            wall_begin=time.perf_counter()
            begin.record(stream)
            for _ in range(copies):dst.copy_(src,non_blocking=True)
            end.record(stream)
            end.synchronize()
            wall_end=time.perf_counter()
            return dict(device=device,wall_start=wall_begin,wall_end=wall_end,
                        event_seconds=begin.elapsed_time(end)/1000.)
    pending=[pool.submit(copy,device) for device in group]
    gate.wait()
    spans=[future.result() for future in pending]
    wall=max(row['wall_end'] for row in spans)-min(row['wall_start'] for row in spans)
    if wall<=0:raise ValueError('Positive simultaneous copy span required')
    return dict(direction=direction,devices=list(group),bytes_per_copy=size,
        copies_per_device=copies,total_bytes=size*copies*len(group),
        wall_seconds=wall,aggregate_bytes_per_second=size*copies*len(group)/wall,
        per_device_event_seconds={str(row['device']):row['event_seconds']
                                  for row in spans},
        start_skew_seconds=max(row['wall_start'] for row in spans)-
                           min(row['wall_start'] for row in spans))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--devices',type=int,nargs='+',required=True)
    parser.add_argument('--size-mib',type=int,default=32)
    parser.add_argument('--copies',type=int,default=8)
    parser.add_argument('--warmups',type=int,default=3)
    parser.add_argument('--repeats',type=int,default=12)
    parser.add_argument('--idle-utilization-limit',type=int,default=20)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    devices=list(args.devices)
    if (not devices or len(set(devices))!=len(devices) or
            any(device<0 for device in devices) or
            not 1<=args.size_mib<=256 or not 1<=args.copies<=64 or
            not 0<=args.warmups<=20 or not 1<=args.repeats<=100 or
            not 0<=args.idle_utilization_limit<=100):
        raise ValueError('Bounded explicit transfer probe arguments required')
    before=gpu_snapshot()
    if not set(devices)<=set(before):
        raise ValueError('Requested GPU is unavailable')
    cuda_uuids={device:_uuid(torch.cuda.get_device_properties(device).uuid)
                for device in devices}
    if any(cuda_uuids[device]!=_uuid(before[device]['uuid']) for device in devices):
        raise ValueError('CUDA device indices differ from physical nvidia-smi indices')
    size=args.size_mib<<20
    report=dict(kind='torchgwas.shared_transfer_diagnostic.v1',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        torch_version=torch.__version__,cuda_version=torch.version.cuda,
        cpu_affinity=(sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity')
                      else None),arguments=dict(devices=devices,size_bytes=size,
        copies=args.copies,warmups=args.warmups,repeats=args.repeats,
        idle_utilization_limit=args.idle_utilization_limit),
        gpu_before=before,cuda_uuids=cuda_uuids,observations=[],
        scope='Pinned host tensors and explicit copy streams. Raw service observations under recorded server load; no guaranteed capacity ceiling or immutable profile price.')
    busy=[device for device in devices if
          before[device]['utilization_percent']>args.idle_utilization_limit]
    if busy:
        report.update(status='skipped_busy_selected_gpu',busy_devices=busy,
                      gpu_after=gpu_snapshot())
    else:
        buffers=_buffers(devices,size)
        groups=[[device] for device in devices]
        if len(devices)>1:groups.append(devices)
        for direction in ('h2d','d2h'):
            for group in groups:
                with ThreadPoolExecutor(max_workers=len(group)) as pool:
                    for _ in range(args.warmups):
                        observe(group,buffers,direction=direction,size=size,
                                copies=args.copies,pool=pool)
                    rows=[observe(group,buffers,direction=direction,size=size,
                                  copies=args.copies,pool=pool)
                          for _ in range(args.repeats)]
                report['observations'].append(dict(direction=direction,
                    devices=list(group),rows=rows,
                    median_aggregate_bytes_per_second=statistics.median(
                        row['aggregate_bytes_per_second'] for row in rows),
                    median_wall_seconds=statistics.median(
                        row['wall_seconds'] for row in rows)))
        report.update(status='diagnostic_complete',gpu_after=gpu_snapshot())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')
        stream.flush();os.fsync(stream.fileno())
    print(json.dumps(dict(status=report['status'],output=str(args.output),
        groups=[dict(direction=row['direction'],devices=row['devices'],
                     median_aggregate_bytes_per_second=row[
                         'median_aggregate_bytes_per_second'])
                for row in report['observations']])))


if __name__=='__main__':main()
