"""Compact serial count-transfer floor for device-selected significant pairs."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .reduced_output_work import significant_output_work


def _price(row):
    if (not isinstance(row, dict) or
            set(row) != {'latency_seconds', 'bytes_per_second', 'resources'} or
            not isinstance(row['resources'], (list, tuple)) or
            not all(isinstance(name, str) and name for name in row['resources']) or
            len(set(row['resources'])) != len(row['resources'])):
        raise ValueError('Explicit independent count-transfer price required')
    latency, rate = row['latency_seconds'], row['bytes_per_second']
    if (any(isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) for value in (latency, rate)) or
            latency < 0 or rate <= 0):
        raise ValueError('Finite count-transfer latency and capacity required')
    return latency + 4 / rate


def native_layout_significant_device_count_floor(source, output,
                                                  count_transfer_prices_by_device):
    """Price only unavoidable serialized CUDA-nonzero count transfers.

    The count transfer follows each block's GPU count and the host waits for
    it. Transfers on one GPU are serial, including successive phenotype tiles;
    different GPUs can overlap. All count bytes already enter the output/D2H
    floor, so this contributes latency only to the partial envelope.
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
        raise ValueError('PGEN input changed before device count pricing')
    limit = output.get('device_selection_max_cells')
    if type(limit) is not int or limit < 1:
        raise ValueError('Bound device selection cell limit required')
    partitions, outputs = source.get('partitions'), output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(outputs, list) or len(outputs) != len(partitions)):
        raise ValueError('One device output row per source partition required')
    devices = {row['device'] for row in partitions}
    if (not isinstance(count_transfer_prices_by_device, dict) or
            set(count_transfer_prices_by_device) != devices):
        raise ValueError('One independent count-transfer price per GPU required')
    seconds = {device: _price(price)
               for device, price in count_transfer_prices_by_device.items()}
    prices = {device: dict(price, resources=list(price['resources']))
              for device, price in count_transfer_prices_by_device.items()}
    by_device = {device: 0. for device in devices}
    rows = []
    for partition, result in zip(partitions, outputs):
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        traits = partition['trait_range'][1] - partition['trait_range'][0]
        geometry = significant_output_work(
            source['samples'], markers, traits, source['chunk_markers'],
            backend='device', max_selection_cells=limit, include_blocks=False)
        blocks = geometry['selection_blocks']
        if (result.get('id') != partition['id'] or
                result.get('device') != partition['device'] or
                result.get('markers') != markers or
                result.get('cells') != markers * traits or
                partition['chunks'] != geometry['source_chunks'] or
                result.get('device_selection_blocks') != blocks or
                result.get('maximum_selection_block_cells') !=
                geometry['maximum_selection_block_cells']):
            raise ValueError('Device count geometry differs from bound output')
        duration = blocks * seconds[partition['device']]
        if not math.isfinite(duration):
            raise ValueError('Device count latency overflow')
        by_device[partition['device']] += duration
        rows.append(dict(id=partition['id'], device=partition['device'],
                         markers=markers, traits=traits, chunks=partition['chunks'],
                         selection_blocks=blocks, count_d2h_bytes=4 * blocks,
                         serial_count_transfer_seconds=duration))
    floor = max(by_device.values())
    if not math.isfinite(floor):
        raise ValueError('Device count latency overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during device count pricing')
    return dict(kind='torchgwas.pgen_layout_significant_device_count_floor.v1',
                input_identity=dict(identity), reduction='significant',
                significant_backend='device', device_selection_max_cells=limit,
                output_partitions_sha256=_digest(outputs), partitions=rows,
                count_transfer_prices_by_device=prices,
                per_device_serial_seconds=by_device,
                count_barrier_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Necessary serial count-transfer latency for device-selected significant blocks. Count bytes are already charged to D2H payload and are not added again here. GPU count/select kernels, CPU dispatch, waits, queueing and final drain remain unpriced; this is not a completion ceiling or JIT switch decision.')
