"""Price variant shards against phenotype tiles from work and transfer.

When the whole phenotype panel fits on every GPU, significant pairs can be
split either way, and the GPU arithmetic is the same (2·N·M·K FLOPs divided
over the GPUs). What differs is data movement (both decode each variant
once when tiles share decode):

* variant shards: the residualized panel (N·K·4 bytes) goes from the first
  GPU to every other one before the scan (not hidden), and each GPU then
  streams only its own M/G variants;
* phenotype tiles: each GPU gets 1/G of the panel, but every GPU streams
  every variant (M·N·b bytes each), overlapped with its GEMM.

Per candidate: up-front copy + max(GEMM per GPU, stream per GPU). The hard
part to know is transfer bandwidth under contention (shared PCIe uplinks,
other users), so it is measured at startup: one pinned copy to all chosen
GPUs at once, and one peer copy from the first GPU. The GEMM rate comes from
the device's SM count and architecture (a prior; it only decides whether a
stream is hidden).
"""
from __future__ import annotations

import time


def peak_fp32_flops(device):
    """Rough dense FP32 rate: SMs x FP32 lanes x 2 x typical boost clock x 0.75."""
    import torch
    props = torch.cuda.get_device_properties(torch.device(device))
    lanes = 128 if (props.major, props.minor) in ((8, 6), (8, 9), (9, 0)) else 64
    clock = {7: 1.55e9, 8: 1.41e9, 9: 1.83e9}.get(props.major, 1.5e9)
    return 2 * lanes * props.multi_processor_count * clock * 0.75


def measure_transfer(devices, nbytes=64 << 20, repeats=3):
    """Host-to-device bytes/s per GPU with all GPUs copying at once, and peer bytes/s from devices[0]."""
    import torch
    host = torch.empty(nbytes, dtype=torch.uint8).pin_memory()
    targets = {d: torch.empty(nbytes, dtype=torch.uint8, device=d) for d in devices}
    streams = {d: torch.cuda.Stream(device=d) for d in devices}
    for d in devices:  # warm contexts and mappings
        with torch.cuda.stream(streams[d]):
            targets[d].copy_(host, non_blocking=True)
    for d in devices:
        torch.cuda.synchronize(d)
    started = time.perf_counter()
    for _ in range(repeats):
        for d in devices:
            with torch.cuda.stream(streams[d]):
                targets[d].copy_(host, non_blocking=True)
    for d in devices:
        torch.cuda.synchronize(d)
    h2d = nbytes * repeats / (time.perf_counter() - started)
    peer = {}
    root = devices[0]
    for d in devices[1:]:
        targets[d].copy_(targets[root])
        torch.cuda.synchronize(d)
        started = time.perf_counter()
        for _ in range(repeats):
            targets[d].copy_(targets[root])
        torch.cuda.synchronize(d)
        torch.cuda.synchronize(root)
        peer[d] = nbytes * repeats / (time.perf_counter() - started)
    return dict(h2d_bytes_per_second=h2d, peer_bytes_per_second=peer)


def split_costs(*, n_samples, n_traits, n_variants, genotype_bytes_per_variant, devices,
                h2d_bytes_per_second, peer_bytes_per_second, flops):
    """Seconds for variant shards and phenotype tiles over `devices` (same GPU count)."""
    count = len(devices)
    panel = float(n_samples) * n_traits * 4
    gemm = 2.0 * n_samples * n_variants * n_traits / count / min(flops.values())
    stream_all = float(n_variants) * genotype_bytes_per_variant / h2d_bytes_per_second
    copies = sum(panel / peer_bytes_per_second[d] for d in devices[1:])  # the root sends them in turn
    shards = dict(upfront=copies, gemm=gemm, stream=stream_all / count)
    tiles = dict(upfront=panel / count / h2d_bytes_per_second, gemm=gemm, stream=stream_all)
    for cost in (shards, tiles):
        cost['total'] = cost['upfront'] + max(cost['gemm'], cost['stream'])
    return dict(variant_shards=shards, phenotype_tiles=tiles)


def choose_split(costs, *, margin=0.05):
    """'traits' when tiles are clearly cheaper, else 'variants'.

    Ties go to shards: independent pipelines, each handling 1/G of the chunks,
    with no lockstep on a shared decode. Measured with device selection
    (2026-09-25), shards beat tiles at 2 and 4 GPUs on lab-2080ti
    (7.1-7.9 s vs 9.5-9.7 s at 2) and H100 (1.8-2.9 s vs 10.2-13.8 s),
    where the model priced them within 5%.
    """
    shards, tiles = costs['variant_shards']['total'], costs['phenotype_tiles']['total']
    return 'traits' if tiles < shards * (1 - margin) else 'variants'
