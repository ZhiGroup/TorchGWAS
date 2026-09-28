"""Triton's fused statistics against the Torch and native CUDA paths, one chunk at a time.

Per device, a synthetic chunk at the full-scale cohort's width (22,250
samples, 4,096 variants, 3% missing calls), for int8 calls and packed 2-bit
PGEN rows, K = 512 and 8,192, dense output and min-p's ranking:
- torch: the int8 conversion, linear._dosage_statistics, and for min-p
  VariantReduction('min-p');
- native: scan_gpu.prepare, the GEMM, scan_gpu.finish (sm80/sm90 builds only);
- triton: triton_scan.prepare, the GEMM, triton_scan.finish (or
  finish_min_p, which keeps only each variant's winner).

Checks triton against the others (beta and t, status, min-p winners), then
times each: GPU seconds by CUDA events and host issue seconds (wall time to
enqueue, no synchronization), the minimum over `--repeats`, on one thread
and on one thread per device at once.

    python benchmarks/triton_scan_20260928.py --devices cuda:5 cuda:6 cuda:7
"""
import argparse
import json
import threading
import time
import traceback

import torch


def inputs(device, n, chunk, width):
    generator = torch.Generator(device=device).manual_seed(20260928)
    calls = torch.randint(0, 3, (chunk, n), device=device, dtype=torch.int8, generator=generator)
    calls[torch.rand((chunk, n), device=device, generator=generator) < 0.03] = -9
    codes = torch.where(calls == -9, torch.full_like(calls, 3), calls).to(torch.uint8)
    width_bytes = ((n + 3) // 4 + 63) // 64 * 64
    padded = torch.zeros((chunk, width_bytes * 4), dtype=torch.uint8, device=device)
    padded[:, :n] = codes
    quads = padded.view(chunk, width_bytes, 4).to(torch.int32)
    packed = (quads[..., 0] | (quads[..., 1] << 2) | (quads[..., 2] << 4) | (quads[..., 3] << 6)).to(torch.uint8)
    design = torch.randn(n, width + 3, device=device, generator=generator) / n ** 0.5
    phenotype_ss = (design[:, :width] ** 2).sum(0)
    return calls, packed.contiguous(), design, phenotype_ss


def paths(device, n, width, calls, packed, design, phenotype_ss):
    from torchgwas import scan_gpu, triton_scan
    from torchgwas.linear import _dosage_statistics
    from torchgwas.reduce import VariantReduction
    ranking = VariantReduction('min-p')
    offset = -2.0 - 2  # rank 2 covariates (+ intercept column), less the genotype

    def torch_path(encoding, mode):
        genotype = torch.where(calls == -9, torch.nan, calls.to(torch.float32))
        beta, t, status, df = _dosage_statistics(genotype, design, phenotype_ss, width, n - 4, False,
                                                 covariate_rank=2)
        return ranking.reduce(beta, t, status, df, 1)[:4] if mode == 'min-p' else (beta, t, status)

    def fused(module, finish_min_p):
        def run(encoding, mode):
            prepared = (module.prepare(packed, encoding='pgen_2bit', n_samples=n) if encoding == 'pgen_2bit'
                        else module.prepare(calls))
            centered, ss, low, high, present = prepared
            products = centered @ design
            if mode == 'min-p' and finish_min_p is not None:
                return finish_min_p(products, ss, low, high, phenotype_ss, present, offset)
            beta, t, status = module.finish(products, ss, low, high, phenotype_ss, present, offset)
            if mode == 'min-p':
                df = present.to(torch.float32) + offset
                return ranking.reduce(beta, t, status, df, 1)[:4]
            return beta, t, status
        return run
    found = dict(torch=torch_path, triton=fused(triton_scan, triton_scan.finish_min_p))
    if scan_gpu.available(device):
        found['native'] = fused(scan_gpu, None)
    return found


def compare(results):
    """Triton against each other path: max |relative| difference of beta and t, equal status/winners."""
    out = {}
    want = results['triton']
    for name, got in results.items():
        if name == 'triton':
            continue
        diffs = []
        for a, b in zip(want[:2], got[:2]):
            scale = b.abs().clamp_min(1e-3)
            diffs.append(float(((a - b).abs() / scale).nan_to_num(0).max()))
        same_status = bool(torch.equal(want[-1] if len(want) == 3 else want[3], got[-1] if len(got) == 3 else got[3]))
        winners = float((want[2] == got[2]).float().mean()) if len(want) == 4 else None
        out[name] = dict(beta_rel=diffs[0], t_rel=diffs[1], status_equal=same_status, winner_agreement=winners)
    return out


def measure(device, n, chunk, width, repeats, out):
    try:
        torch.cuda.set_device(device)
        calls, packed, design, phenotype_ss = inputs(device, n, chunk, width)
        runs = paths(device, n, width, calls, packed, design, phenotype_ss)
        rows = {}
        for encoding in ('int8', 'pgen_2bit'):
            for mode in ('dense', 'min-p'):
                results = {name: run(encoding, mode) for name, run in runs.items()
                           if not (name == 'torch' and encoding == 'pgen_2bit')}
                torch.cuda.synchronize(device)
                row = dict(check=compare(results))
                for name, run in runs.items():
                    if name == 'torch' and encoding == 'pgen_2bit':
                        continue
                    host, gpu = [], []
                    for _ in range(repeats):
                        start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                        began = time.perf_counter()
                        start.record()
                        run(encoding, mode)
                        stop.record()
                        host.append(time.perf_counter() - began)
                        stop.synchronize()
                        gpu.append(start.elapsed_time(stop) / 1e3)
                    row[name] = dict(host_ms=round(1e3 * min(host), 3), gpu_ms=round(1e3 * min(gpu), 3))
                rows[(encoding, mode)] = row
        out[str(device)] = rows
    except Exception:  # noqa: BLE001 - reported, and the other threads carry on
        out[str(device)] = traceback.format_exc()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--samples', type=int, default=22_250)
    parser.add_argument('--chunk', type=int, default=4096)
    parser.add_argument('--repeats', type=int, default=10)
    args = parser.parse_args()
    devices = [torch.device(name) for name in args.devices]
    for width in (512, 8192):
        for count in sorted({1, len(devices)}):
            out = {}
            threads = [threading.Thread(target=measure, args=(d, args.samples, args.chunk, width, args.repeats, out))
                       for d in devices[:count]]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            failures = {k: v for k, v in out.items() if isinstance(v, str)}
            if failures:
                print(json.dumps(dict(K=width, threads=count, failures=failures)), flush=True)
                continue
            first = next(iter(out.values()))
            for key in first:
                row = dict(K=width, threads=count, encoding=key[0], mode=key[1],
                           check=first[key]['check'])
                for name in ('torch', 'native', 'triton'):
                    if name in first[key]:
                        row[name] = dict(host_ms=max(v[key][name]['host_ms'] for v in out.values()),
                                         gpu_ms=max(v[key][name]['gpu_ms'] for v in out.values()))
                print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
