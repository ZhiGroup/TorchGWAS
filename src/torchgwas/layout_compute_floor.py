"""Compact mandatory H2D and matrix-product work for native PGEN layouts.

Only the explicitly unpacked int8 hardcall transport is covered. The full
statistics and JAGWAS kernels do more work than the counted matrix products;
their omitted work can only raise the conditional resource floor.
"""
import math

from .analytical_plan_cache import input_identity
from .jagwas_blocks import projection_flops_per_variant
from .layout_transfer_links import transfer_link_loads


def _positive(name, value):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0):
        raise ValueError('Positive finite ' + name + ' ceiling required')
    return value


def _device_ceilings(name, values, devices):
    if not isinstance(values, dict) or set(values) != devices:
        raise ValueError('One ' + name + ' ceiling per active GPU required')
    return {device: _positive(device + ' ' + name, value)
            for device, value in values.items()}


def native_layout_compute_floor(layout, *, covariate_rank,
                                shared_h2d_bytes_per_second,
                                per_device_h2d_bytes_per_second,
                                peak_fp32_flops_per_second,
                                peak_fp64_flops_per_second=None,
                                shared_links=(),
                                transport='native_int8_hardcall'):
    """Return necessary payload/GEMM floors under declared capacity ceilings.

    `covariate_rank` is the number of Q columns; the native design also adds
    one intercept. These ceilings must be upper limits on service for the
    chosen hardware/context. Shape efficiency, kernel launch/conversion,
    setup, queues and contention beyond the supplied shared H2D ceiling are
    omitted. The result cannot select a chunk or a GPU count by itself.
    """
    if not isinstance(layout, dict) or layout.get('kind') != 'torchgwas.pgen_layout_source_floor.v1':
        raise ValueError('Typed native PGEN source layout required')
    if transport != 'native_int8_hardcall':
        raise ValueError('Only explicitly unpacked native int8 PGEN transport is covered')
    samples = layout.get('samples')
    if (type(samples) is not int or samples < 3 or type(covariate_rank) is not int or
            not 0 <= covariate_rank < samples - 2):
        raise ValueError('Bound sample count and residual covariate rank required')
    source = layout.get('input_identity')
    if not isinstance(source, dict) or 'path' not in source or input_identity(source['path']) != source:
        raise ValueError('PGEN input changed before compute-load composition')
    rows = layout.get('partitions')
    if not isinstance(rows, list) or not rows:
        raise ValueError('Nonempty fixed source partitions required')
    devices = {row.get('device') for row in rows if isinstance(row, dict)}
    if len(devices) == 0 or any(not isinstance(d, str) for d in devices):
        raise ValueError('Explicit fixed source devices required')
    shared_h2d = _positive('shared H2D', shared_h2d_bytes_per_second)
    local_h2d = _device_ceilings('H2D', per_device_h2d_bytes_per_second, devices)
    fp32 = _device_ceilings('FP32', peak_fp32_flops_per_second, devices)
    reduction = layout.get('reduction')
    if reduction not in (None, 'significant', 'jagwas'):
        raise ValueError('Unknown fixed source reduction')
    if reduction == 'jagwas':
        fp64 = _device_ceilings('FP64', peak_fp64_flops_per_second, devices)
    elif peak_fp64_flops_per_second is not None:
        raise ValueError('FP64 projection ceiling applies only to JAGWAS')
    else:
        fp64 = None
    device_work = {device: dict(h2d_bytes=0, fp32_gemm_flops=0,
                                fp64_projection_flops=0) for device in devices}
    partitions = []
    for row in rows:
        if (not isinstance(row, dict) or
                set(row) != {'id', 'device', 'variant_range', 'trait_range', 'chunks'}):
            raise ValueError('Complete fixed source partition required')
        variant, trait = row['variant_range'], row['trait_range']
        if (not isinstance(variant, list) or not isinstance(trait, list) or
                len(variant) != 2 or len(trait) != 2 or
                any(type(value) is not int for value in variant + trait) or
                not 0 <= variant[0] < variant[1] or not 0 <= trait[0] < trait[1]):
            raise ValueError('Nonempty fixed partition geometry required')
        markers, width = variant[1] - variant[0], trait[1] - trait[0]
        size, chunks = layout.get('chunk_markers'), row['chunks']
        if (type(size) is not int or size < 1 or type(chunks) is not int or
                chunks != (markers + size - 1) // size):
            raise ValueError('Fixed source chunk geometry changed')
        if reduction == 'jagwas' and trait != [0, layout.get('total_traits')]:
            raise ValueError('JAGWAS requires the full phenotype panel per GPU')
        # native_scan uploads (markers, samples) int8 and multiplies the
        # centered FP32 genotype by [phenotypes | intercept | Q].
        h2d = markers * samples
        gemm = 2 * markers * samples * (width + covariate_rank + 1)
        # The JAGWAS projection (jagwas_projection) multiplies each FP64 score
        # row by the lower-triangular inverse Cholesky factor, block by block.
        projection = markers * projection_flops_per_variant(width) if reduction == 'jagwas' else 0
        work = device_work[row['device']]
        work['h2d_bytes'] += h2d
        work['fp32_gemm_flops'] += gemm
        work['fp64_projection_flops'] += projection
        partitions.append(dict(id=row['id'], device=row['device'], markers=markers,
                               traits=width, h2d_bytes=h2d,
                               fp32_gemm_flops=gemm,
                               fp64_projection_flops=projection))
    total_h2d = sum(work['h2d_bytes'] for work in device_work.values())
    link_report = transfer_link_loads(
        {device: work['h2d_bytes'] for device, work in device_work.items()},
        shared_links, direction='h2d')
    h2d_floor = max(total_h2d / shared_h2d,
                    link_report['floor_seconds'],
                    *(work['h2d_bytes'] / local_h2d[device]
                      for device, work in device_work.items()))
    # All matrix products for a device run on its scan compute stream, so
    # their peak-rate minimum services add even when different GPUs overlap.
    compute_floor = max(
        work['fp32_gemm_flops'] / fp32[device] +
        (work['fp64_projection_flops'] / fp64[device] if fp64 is not None else 0.)
        for device, work in device_work.items())
    floor = max(h2d_floor, compute_floor)
    if not math.isfinite(floor):
        raise ValueError('Compute/transfer resource floor overflow')
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed during compute-load composition')
    return dict(kind='torchgwas.pgen_layout_compute_floor.v1',
                input_identity=dict(source), samples=samples,
                covariate_rank=covariate_rank, transport=transport,
                reduction=reduction, partitions=partitions,
                total_h2d_bytes=total_h2d, per_device_work=device_work,
                shared_h2d_link_loads=link_report,
                capacity_ceilings=dict(shared_h2d_bytes_per_second=shared_h2d,
                    per_device_h2d_bytes_per_second=local_h2d,
                    peak_fp32_flops_per_second=fp32,
                    peak_fp64_flops_per_second=fp64),
                h2d_payload_floor_seconds=h2d_floor,
                matrix_product_floor_seconds=compute_floor,
                compute_transfer_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional necessary H2D payload and matrix-product floors for fixed unpacked int8 hardcall PGEN partitions. Global, optionally declared overlapping shared links, and per-GPU H2D capacities constrain the same transfers without summing link times; each GPU has its own serial compute ceiling. Design setup/phenotype transfer, conversion, other statistics/reduction kernels, selector, output, waits and final drain are omitted. Capacity inputs must be valid upper ceilings; not elapsed-time bounds or switch authorization.')
