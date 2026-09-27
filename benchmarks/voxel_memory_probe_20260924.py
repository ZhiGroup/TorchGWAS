"""What bounds the phenotype tile at voxel scale: device residency or the host result ring?"""
import sys
sys.path.insert(0, 'src')
from torchgwas.pipeline_model import auto_trait_block, device_ring_bytes, host_pinned_bytes

GB = 1e9
cases = [
    ('frozen H100 panel', 35_365, 16_385, 27),
    ('voxel (docstring)', 33_417, 2_085_000, 27),
]
for name, n, k, c in cases:
    for label, per_variant in (('int8', float(n)), ('2-bit', float((n+3)//4))):
        for chunk in (1024, 2048):
            depth = 4
            dev_only = auto_trait_block(n_samples=n, n_traits=k, covariate_rank=c, chunk_variants=chunk, depth=depth,
                                        transfer_bytes_per_variant=per_variant, device_memory_bytes=80*2**30)
            both = auto_trait_block(n_samples=n, n_traits=k, covariate_rank=c, chunk_variants=chunk, depth=depth,
                                    transfer_bytes_per_variant=per_variant, device_memory_bytes=80*2**30,
                                    host_memory_bytes=500*2**30, trait_devices=8)
            reduced = auto_trait_block(n_samples=n, n_traits=k, covariate_rank=c, chunk_variants=chunk, depth=depth,
                                       transfer_bytes_per_variant=per_variant, device_memory_bytes=80*2**30,
                                       host_memory_bytes=500*2**30, trait_devices=8, reduction_width=1)
            w = both
            dev = device_ring_bytes(chunk_variants=chunk, depth=depth, n_samples=n, n_traits=w, covariate_rank=c,
                                    transfer_bytes_per_variant=per_variant)
            host = host_pinned_bytes(chunk_variants=chunk, depth=depth, n_traits=w, transfer_bytes_per_variant=per_variant)
            pheno = n*w*4
            print(f'{name:18s} {label:5s} chunk {chunk}: device-only fit {dev_only:>9,d}  device+host(8 GPUs) {both:>9,d}  '
                  f'reduced-output {reduced:>9,d} | at {w:,d}: device {dev/GB:6.1f} GB (phenotype {pheno/GB:5.1f}), '
                  f'host pinned {host/GB:6.1f} GB per GPU')
    print(f'{name}: whole panel as float32 = {n*k*4/GB:,.1f} GB; per GPU over 8 = {n*k*4/GB/8:,.1f} GB')
