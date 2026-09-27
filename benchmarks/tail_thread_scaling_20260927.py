"""Does the device tail serialize across shard threads? Calls per second, 1 vs N threads.

Each thread owns one GPU and calls tails.neg_log10_p_device on a (4096, 1)
column (a min-p chunk's winners) and on a (4096, 512) chunk (dense K = 512),
`--calls` times, synchronizing at the end. If the host side of a call holds
the GIL, N threads take about N times as long as one.

    python benchmarks/tail_thread_scaling_20260927.py cuda:4 cuda:5 cuda:6 cuda:7
"""
import json
import sys
import threading
import time

import torch

from torchgwas import tails


def run(device, shape, calls, out):
    t = torch.randn(shape, device=device)
    df = torch.full((shape[0], 1), 22_238.0, device=device)
    result = torch.empty(shape, dtype=torch.float32, device=device)
    tails.neg_log10_p_device(t, df, out=result)
    torch.cuda.synchronize(device)
    started = time.perf_counter()
    for _ in range(calls):
        tails.neg_log10_p_device(t, df, out=result)
    torch.cuda.synchronize(device)
    out[str(device)] = time.perf_counter() - started


def main():
    devices = [torch.device(name) for name in sys.argv[1:]]
    for device in devices:
        tails.prepare_device_tail(device)
    for shape, calls in (((4096, 1), 200), ((4096, 512), 100)):
        for count in (1, len(devices)):
            out = {}
            threads = [threading.Thread(target=run, args=(device, shape, calls, out)) for device in devices[:count]]
            started = time.perf_counter()
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            wall = time.perf_counter() - started
            print(json.dumps(dict(shape=shape, threads=count, calls_per_thread=calls, wall_seconds=round(wall, 3),
                                  ms_per_call=round(1e3 * wall / calls, 3),
                                  kind=tails._KINDS.get(devices[0]))), flush=True)


if __name__ == '__main__':
    main()
