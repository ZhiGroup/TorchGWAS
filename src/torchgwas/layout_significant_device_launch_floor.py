"""Compact mandatory CUDA launch service for device-selected significant pairs.

This charges only count, flagged select, and conditional coordinate-scatter
launches. Other eager selector kernels, CPU dispatch and queue work remain.
"""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS, device_selection_shape


def _launch_profile(device, row):
    if (not isinstance(row, dict) or set(row) != {'torch_version', 'cuda_runtime',
            'compute_capability', 'kernel_launch_seconds', 'gpu_fraction'} or
            str(row['torch_version']).split('+')[0] != '2.5.1' or
            row['cuda_runtime'] != '12.4' or
            row['compute_capability'] not in ([8, 0], [9, 0])):
        raise ValueError('Supported installed device-selector launch policy required')
    import torch
    if (row['torch_version'] != torch.__version__ or
            row['cuda_runtime'] != torch.version.cuda or
            row['compute_capability'] !=
            list(torch.cuda.get_device_capability(device))):
        raise ValueError('Stale installed device-selector launch profile')
    launch, fraction = row['kernel_launch_seconds'], row['gpu_fraction']
    if (any(isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) for value in (launch, fraction)) or
            launch <= 0 or not 0 < fraction <= 1):
        raise ValueError('Positive independent launch price and GPU fraction required')
    return launch / fraction


def _chunk_launch_counts(markers, traits, limit):
    """Count source blocks and two-kernel count paths in four strip classes."""
    width, height, declared = device_selection_shape(markers, traits, limit)
    full_rows, last_rows = divmod(markers, height)
    full_traits, last_traits = divmod(traits, width)
    classes = ((height, width, full_rows * full_traits),
               (last_rows, width, full_traits if last_rows else 0),
               (height, last_traits, full_rows if last_traits else 0),
               (last_rows, last_traits,
                1 if last_rows and last_traits else 0))
    blocks = sum(count for _, _, count in classes)
    if blocks != declared:
        raise ValueError('Compact selection strip geometry changed')
    two_count = sum(count for rows, cols, count in classes
                    if rows * cols > 4096)
    largest = max((rows * cols for rows, cols, count in classes if count),
                  default=0)
    return blocks, two_count, largest


def native_layout_significant_device_launch_floor(source, output,
                                                   launch_profiles_by_device):
    """Price mandatory source-policy launches on each serial CUDA stream.

    The low/high endpoints vary only in how many selection blocks retain a
    pair. Count and flagged-select launches occur even for empty blocks.
    A block over 4096 cells needs a second count-reduction launch under the
    checked PyTorch 2.5.1/CUB 2.3 policy. The returned high endpoint bounds
    this launch subset, not total selector service or job completion.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            source.get('reduction') != 'significant' or
            source.get('partition_axis') != 'trait' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            output.get('significant_backend') != 'device' or
            output.get('input_identity') != source.get('input_identity')):
        raise ValueError('Matching device-selected source/output required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before device launch pricing')
    limit = output.get('device_selection_max_cells')
    if type(limit) is not int or not 1 <= limit <= DEVICE_SELECTION_MAX_CELLS:
        raise ValueError('Bound source-policy selection cell limit required')
    partitions, outputs = source.get('partitions'), output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(outputs, list) or len(outputs) != len(partitions)):
        raise ValueError('One device output row per source partition required')
    devices = {row['device'] for row in partitions}
    if (not isinstance(launch_profiles_by_device, dict) or
            set(launch_profiles_by_device) != devices):
        raise ValueError('One independent launch profile per GPU required')
    price = {device: _launch_profile(device, profile)
             for device, profile in launch_profiles_by_device.items()}
    per_device = {device: [0., 0.] for device in devices}
    rows = []
    for part, result in zip(partitions, outputs):
        markers = part['variant_range'][1] - part['variant_range'][0]
        traits = part['trait_range'][1] - part['trait_range'][0]
        size = source['chunk_markers']
        full, tail = divmod(markers, size)
        full_shape = _chunk_launch_counts(size, traits, limit) if full else (0, 0, 0)
        tail_shape = _chunk_launch_counts(tail, traits, limit) if tail else (0, 0, 0)
        blocks = full * full_shape[0] + tail_shape[0]
        extra_count = full * full_shape[1] + tail_shape[1]
        maximum_block = max(full_shape[2], tail_shape[2])
        if (result.get('id') != part['id'] or
                result.get('device') != part['device'] or
                result.get('markers') != markers or
                result.get('cells') != markers * traits or
                part['chunks'] != full + bool(tail) or
                result.get('device_selection_blocks') != blocks or
                result.get('maximum_selection_block_cells') != maximum_block):
            raise ValueError('Device launch geometry differs from bound output')
        retained = result.get('retained_rows')
        if (not isinstance(retained, list) or len(retained) != 2 or
                any(type(value) is not int or value < 0 for value in retained) or
                not retained[0] <= retained[1] <= markers * traits):
            raise ValueError('Explicit bounded retained-pair interval required')
        least_scatter = (retained[0] + maximum_block - 1) // maximum_block
        most_scatter = min(retained[1], blocks)
        # Every block launches one count and two flagged-select kernels.
        base = 3 * blocks + extra_count
        launches = [base + least_scatter, base + most_scatter]
        seconds = [value * price[part['device']] for value in launches]
        if not all(math.isfinite(value) for value in seconds):
            raise ValueError('Device selector launch service overflow')
        by_device = per_device[part['device']]
        by_device[0] += seconds[0]
        by_device[1] += seconds[1]
        rows.append(dict(id=part['id'], device=part['device'], markers=markers,
                         traits=traits, chunks=part['chunks'],
                         selection_blocks=blocks,
                         two_kernel_count_blocks=extra_count,
                         maximum_selection_block_cells=maximum_block,
                         retained_rows=list(retained),
                         coordinate_scatter_launches=[least_scatter,
                                                      most_scatter],
                         mandatory_launches=launches,
                         serial_launch_seconds=seconds))
    floor = [max(value[index] for value in per_device.values())
             for index in (0, 1)]
    if not all(math.isfinite(value) for row in per_device.values()
               for value in row):
        raise ValueError('Device selector launch service overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during device launch pricing')
    return dict(kind='torchgwas.pgen_layout_significant_device_launch_floor.v1',
                input_identity=dict(identity), reduction='significant',
                significant_backend='device', device_selection_max_cells=limit,
                output_partitions_sha256=_digest(outputs), partitions=rows,
                launch_profiles_by_device={device: dict(
                    profile, compute_capability=list(profile['compute_capability']))
                    for device, profile in launch_profiles_by_device.items()},
                per_device_serial_seconds=per_device,
                selector_launch_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Necessary per-GPU serial launch service for mandatory count/flagged-select and bounded coordinate-scatter launches under the PyTorch 2.5.1/CUB 2.3 source policy. Eager selector kernels, kernel body work, CPU dispatch, blocking transfers, queues and final drain remain outside. The high endpoint bounds this partial launch service only, not completion or a JIT switch.')