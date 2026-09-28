"""One chunk's device work issued eagerly against replayed from a CUDA graph.

K = 512 on four GPUs is host-bound: every chunk costs ~3.4 ms of host time
that does not parallelize across the GPUs (docs/autotune_design_20260924.md,
"A per-chunk host model"). Issuing the chunk's ~100 device ops took 0.7-0.9
ms on one thread and 2-2.5 ms each with four at once
(chunk_host_dispatch_20260927.py). A CUDA graph issues them as one launch.

Per device: the default statistics path (int8 calls to float with the
missing code, torch statistics) and the mode's tail -- min-p's reduction or
dense output's -log10 P -- on a synthetic chunk, issued (a) eagerly and (b)
as one graph replay after a copy into the graph's static input. Host
seconds are the issue wall time (no synchronization), GPU seconds come from
CUDA events; the minimum over `--repeats`. Threads run one per device at
once, like variant shards. The graph's outputs are checked against eager.

    python benchmarks/chunk_graph_dispatch_20260928.py --devices cuda:1 cuda:2 cuda:7

Found (2026-09-28): the statistics and min-p's |t| ranking capture; the
AOTInductor tail does not ("operation not permitted when stream is
capturing"). Capturing from several threads at once hung and was not
resolved; the Triton kernels (triton_scan) cut the launches instead.
"""
import argparse
import json
import threading
import time

import torch


def work(device, mode, width, chunk, n_samples):
    """The chunk's calls, its capturable core and its tail.

    The -log10 P tail is an AOTInductor graph, which refuses to run under
    stream capture ("operation not permitted when stream is capturing"); it
    is one host call in its whole-graph form, so it stays eager after the
    captured core (statistics, and min-p's |t| ranking for a complete panel).
    """
    from torchgwas.linear import _dosage_statistics
    from torchgwas.reduce import VariantReduction
    from torchgwas.tails import neg_log10_p_device
    generator = torch.Generator(device=device).manual_seed(20260928)
    calls = torch.randint(0, 3, (chunk, n_samples), device=device, dtype=torch.int8, generator=generator)
    calls[:, :7] = -9
    design = torch.randn(n_samples, width + 1, device=device, generator=generator) / n_samples ** 0.5
    phenotype_ss = torch.ones(width, device=device)
    ranking = VariantReduction('min-p')

    def core(raw):
        genotype = torch.where(raw == -9, torch.nan, raw.to(torch.float32))
        beta, t, status, variant_df = _dosage_statistics(genotype, design, phenotype_ss, width,
                                                         n_samples - 3, False, covariate_rank=1)
        if mode == 'min-p':
            return ranking.reduce(beta, t, status, variant_df, 1)
        return beta, t, status, variant_df

    def tail(outputs):
        t, variant_df = outputs[1], outputs[-1]
        return neg_log10_p_device(t, variant_df.reshape(-1, 1).expand(t.shape))
    return calls, core, tail


def measure(device, mode, width, chunk, n_samples, repeats, out, barrier):
    torch.cuda.set_device(device)
    calls, core, tail = work(device, mode, width, chunk, n_samples)
    stream = torch.cuda.Stream(device)
    static = torch.empty_like(calls)
    with torch.cuda.stream(stream):
        eager = core(calls)
        eager_logp = tail(eager)
        static.copy_(calls)
        for _ in range(2):  # warm the allocator on this stream before capture
            core(static)
    stream.synchronize()
    graph = torch.cuda.CUDAGraph()
    # thread_local: the other devices' threads keep issuing work meanwhile.
    with torch.cuda.graph(graph, stream=stream, capture_error_mode='thread_local'):
        captured = core(static)
    stream.synchronize()

    def eagerly(raw):
        tail(core(raw))

    def replay(raw):
        static.copy_(raw, non_blocking=True)
        graph.replay()
        return tail(captured)
    result = {}
    for name, issue in (('eager', eagerly), ('graph', replay)):
        host, gpu = [], []
        barrier.wait()
        with torch.cuda.stream(stream):
            for _ in range(repeats):
                start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                started = time.perf_counter()
                start.record(stream)
                issue(calls)
                stop.record(stream)
                host.append(time.perf_counter() - started)
                stop.synchronize()
                gpu.append(start.elapsed_time(stop) / 1e3)
        result[name] = (min(host), min(gpu))
    with torch.cuda.stream(stream):
        graph_logp = replay(calls)
    stream.synchronize()
    same = all(torch.equal(torch.nan_to_num(a.float()), torch.nan_to_num(b.float()))
               for a, b in zip((*eager, eager_logp), (*captured, graph_logp)))
    out[str(device)] = dict(result, same=same)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--samples', type=int, default=22_250)
    parser.add_argument('--repeats', type=int, default=20)
    args = parser.parse_args()
    from torchgwas import tails
    devices = [torch.device(name) for name in args.devices]
    for device in devices:
        tails.prepare_device_tail(device)
    for mode in ('min-p', 'dense'):
        for width in (512, 8192):
            for count in sorted({1, len(devices)}):
                out = {}
                barrier = threading.Barrier(count)
                threads = [threading.Thread(target=measure, args=(d, mode, width, 4096, args.samples,
                                                                  args.repeats, out, barrier))
                           for d in devices[:count]]
                for thread in threads:
                    thread.start()
                for thread in threads:
                    thread.join()
                row = dict(mode=mode, K=width, chunk=4096, threads=count,
                           same=all(v['same'] for v in out.values()))
                for name in ('eager', 'graph'):
                    row[f'{name}_host_ms'] = round(1e3 * max(v[name][0] for v in out.values()), 3)
                    row[f'{name}_gpu_ms'] = round(1e3 * max(v[name][1] for v in out.values()), 3)
                print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
