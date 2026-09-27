"""Host (GIL) seconds to issue one chunk's device work, against its GPU seconds.

The same synthetic chunk as empirical_autotune.scan_device_seconds_per_variant,
timed on the host (wall time to issue the ops, no synchronization) and on the
device (CUDA events), for min-p and dense output, chunk 1024 and 4096,
K = 512 and 8,192, one thread and one per GPU on `devices` at once. A
full-scale min-p run at K = 512 on four GPUs behaved as if every chunk cost
~3.5 ms of shared host time (1,975 chunks in 7.1 s at 4096; the tuner
estimated 27.9 s at 1024).

    python benchmarks/chunk_host_dispatch_20260927.py cuda:4 cuda:5 cuda:6 cuda:7
"""
import json
import sys
import threading
import time

import torch


def measure(device, mode, width, chunk, n_samples, repeats, out):
    from torchgwas.linear import _dosage_statistics
    from torchgwas.min_p import MinPReduction
    from torchgwas.tails import neg_log10_p_device
    torch.cuda.set_device(device)
    codes = torch.randint(0, 3, (chunk, n_samples), device=device, dtype=torch.uint8)
    design = torch.randn(n_samples, width + 1, device=device) / n_samples ** 0.5
    phenotype_ss = torch.ones(width, device=device)
    reduction = MinPReduction()

    def once():
        beta, t, status, variant_df = _dosage_statistics(codes.float(), design, phenotype_ss, width,
                                                         n_samples - 3, False, covariate_rank=1)
        if mode == 'min-p':
            reduction.reduce(beta, t, status, variant_df, 1, log10_p=(None, None, torch.float32))
        else:
            neg_log10_p_device(t, variant_df.reshape(-1, 1))

    once()
    torch.cuda.synchronize(device)
    host, gpu = [], []
    for _ in range(repeats):
        start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
        started = time.perf_counter()
        start.record()
        once()
        stop.record()
        host.append(time.perf_counter() - started)
        stop.synchronize()
        gpu.append(start.elapsed_time(stop) / 1e3)
    out[str(device)] = (min(host), min(gpu))


def main():
    from torchgwas import tails
    devices = [torch.device(name) for name in sys.argv[1:]]
    for device in devices:
        tails.prepare_device_tail(device)
    for mode in ('min-p', 'dense'):
        for width in (512, 8192):
            for chunk in (1024, 4096):
                for count in (1, len(devices)):
                    out = {}
                    threads = [threading.Thread(target=measure, args=(d, mode, width, chunk, 22_250, 5, out))
                               for d in devices[:count]]
                    for thread in threads:
                        thread.start()
                    for thread in threads:
                        thread.join()
                    host = max(v[0] for v in out.values()); gpu = max(v[1] for v in out.values())
                    print(json.dumps(dict(mode=mode, K=width, chunk=chunk, threads=count,
                                          host_ms=round(1e3 * host, 3), gpu_ms=round(1e3 * gpu, 3))), flush=True)


if __name__ == '__main__':
    main()
