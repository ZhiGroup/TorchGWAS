"""Host-to-GPU delivery rate: per-GPU PCIe copies versus one PCIe copy plus peer fan-out.

Usage: python fanout_transfer_probe_20260923.py cuda:0 cuda:1 [--mib 256 --chunks 24]

Delivers the same pinned chunks to every listed GPU two ways, pipelined like
SharedDecodeHub: 'pcie' (each GPU copies each chunk from pinned host memory on
its own stream) and 'fanout' (each chunk crosses PCIe once into a two-slot
ring on the first GPU; the others pull it from there on dedicated source-side
streams while the next chunk arrives). Reports GB/s delivered per GPU.
Transfer only: no decode, no compute.
"""
import argparse
import json
import time

import torch


def run(devices, mib, chunks, mode):
    host = torch.empty(mib << 20, dtype=torch.uint8).pin_memory()
    dst = {d: [torch.empty_like(host, device=d) for _ in range(2)] for d in devices}
    copy = {d: torch.cuda.Stream(device=d) for d in devices}
    root = devices[0]
    peer = {d: torch.cuda.Stream(device=root) for d in devices[1:]}
    arrived = [torch.cuda.Event() for _ in range(2)]
    done = {d: [torch.cuda.Event() for _ in range(2)] for d in devices}
    for d in devices:
        torch.cuda.synchronize(d)
    started = time.perf_counter()
    for index in range(chunks):
        slot = index % 2
        if mode == 'pcie':
            for d in devices:
                with torch.cuda.stream(copy[d]):
                    if index >= 2:
                        copy[d].wait_event(done[d][slot])
                    dst[d][slot].copy_(host, non_blocking=True)
                    done[d][slot].record(copy[d])
        else:
            with torch.cuda.stream(copy[root]):
                if index >= 2:  # every peer finished reading this root slot
                    for d in devices[1:]:
                        copy[root].wait_event(done[d][slot])
                dst[root][slot].copy_(host, non_blocking=True)
                arrived[slot].record(copy[root])
            for d in devices[1:]:
                with torch.cuda.stream(peer[d]), torch.cuda.stream(copy[d]):
                    copy[d].wait_event(arrived[slot])
                    dst[d][slot].copy_(dst[root][slot], non_blocking=True)
                    done[d][slot].record(copy[d])
    for d in devices:
        torch.cuda.synchronize(d)
    seconds = time.perf_counter()-started
    return chunks*host.numel()/seconds/1e9


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('devices', nargs='+')
    parser.add_argument('--mib', type=int, default=256)
    parser.add_argument('--chunks', type=int, default=24)
    parser.add_argument('--rounds', type=int, default=5)
    args = parser.parse_args()
    devices = [str(torch.device(d)) for d in args.devices]
    peer_ok = all(torch.cuda.can_device_access_peer(torch.device(d).index, torch.device(devices[0]).index)
                  for d in devices[1:])
    rows = {'alone': [], 'pcie': [], 'fanout': []}
    run(devices, args.mib, 4, 'pcie')  # warm contexts and peer mappings
    run(devices, args.mib, 4, 'fanout')
    for _ in range(args.rounds):
        rows['alone'].append(run(devices[:1], args.mib, args.chunks, 'pcie'))
        rows['pcie'].append(run(devices, args.mib, args.chunks, 'pcie'))
        rows['fanout'].append(run(devices, args.mib, args.chunks, 'fanout'))
    summary = {k: sorted(round(x, 2) for x in v) for k, v in rows.items()}
    print(json.dumps(dict(devices=devices, peer_access=peer_ok, mib=args.mib, chunks=args.chunks,
                          gb_per_s_per_gpu=summary)))


if __name__ == '__main__':
    main()
