"""Compact necessary host-selection work for full-panel JAGWAS shards."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest


_PRIMITIVES = ('fp32_to_fp64_view', 'finite_fp64',
               'flatnonzero_empty', 'flatnonzero_nonempty',
               'index_add', 'fp64_gather')


def _number(name, value, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError('Finite independent JAGWAS selector price required: ' + name)
    return value


def _prices(bank):
    if not isinstance(bank, dict) or set(bank) != set(_PRIMITIVES):
        raise ValueError('Complete independent JAGWAS selector prices required')
    for name, price in bank.items():
        if not isinstance(price, dict) or set(price) != {'call_cpu_seconds',
                                                         'unit_cpu_seconds'}:
            raise ValueError('Complete JAGWAS selector primitive required: ' + name)
        for field, value in price.items():
            _number(name + '.' + field, value)
    return bank


def _chunk_base(rows, bank, regimes):
    fixed = sum(bank[name]['call_cpu_seconds'] +
                rows * bank[name]['unit_cpu_seconds']
                for name in ('fp32_to_fp64_view', 'finite_fp64'))
    fixed += sum(bank[name]['call_cpu_seconds']
                 for name in ('index_add', 'fp64_gather'))
    nonzero = [bank[name]['call_cpu_seconds'] +
               rows * bank[name]['unit_cpu_seconds'] for name in regimes]
    return fixed + min(nonzero), fixed + max(nonzero)


def native_layout_jagwas_selection_floor(source, output, prices,
                                         cpu_fraction_by_device):
    """Price mandatory JAGWAS host selection without expanding source chunks.

    Each shard keeps the full phenotype panel. Retained counts are explicit
    scenarios; empty/nonempty nonzero rates are both considered when a shard
    may contain either kind of chunk. All selector work passes through the
    public single consumer. The result bounds only this conditional floor.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            source.get('reduction') != 'jagwas' or
            source.get('partition_axis') != 'variant' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            output.get('reduction') != 'jagwas' or
            output.get('input_identity') != source.get('input_identity')):
        raise ValueError('Matching bound full-panel JAGWAS source/output required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before JAGWAS selector pricing')
    bank = _prices(prices)
    partitions = source.get('partitions')
    output_rows = output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(output_rows, list) or len(output_rows) != len(partitions)):
        raise ValueError('One explicit output row per JAGWAS variant shard required')
    devices = {row['device'] for row in partitions}
    if (not isinstance(cpu_fraction_by_device, dict) or
            set(cpu_fraction_by_device) != devices):
        raise ValueError('One measured JAGWAS consumer CPU fraction per GPU required')
    fractions = {device: _number(device + ' CPU fraction', value, positive=True)
                 for device, value in cpu_fraction_by_device.items()}
    if any(value > 1 for value in fractions.values()):
        raise ValueError('JAGWAS CPU fraction exceeds one')
    cap = source['shared_capacities']
    cpu_cap = _number('shared CPU capacity', cap['cpu'], positive=True)
    dram_cap = _number('shared DRAM capacity', cap['dram'], positive=True)
    per_retained_cpu = (bank['index_add']['unit_cpu_seconds'] +
                        bank['fp64_gather']['unit_cpu_seconds'])
    rows = []
    total_cpu = [[], []]
    total_dram = [[], []]
    serial = [[], []]
    for partition, result in zip(partitions, output_rows):
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        block = source['chunk_markers']
        full, tail = divmod(markers, block)
        chunks = full + bool(tail)
        retained = result.get('retained_rows')
        if (result.get('id') != partition['id'] or
                result.get('device') != partition['device'] or
                result.get('markers') != markers or
                result.get('cells') != markers * source['total_traits'] or
                partition['trait_range'] != [0, source['total_traits']] or
                partition['chunks'] != chunks or
                not isinstance(retained, list) or len(retained) != 2 or
                any(type(value) is not int for value in retained) or
                not 0 <= retained[0] <= retained[1] <= markers):
            raise ValueError('JAGWAS selector geometry or occupancy differs from output')
        if retained[1] == 0:
            regimes = ('flatnonzero_empty',)
        elif retained[0] == markers:
            regimes = ('flatnonzero_nonempty',)
        else:
            regimes = ('flatnonzero_empty', 'flatnonzero_nonempty')
        base = [0., 0.]
        for count, width in ((full, block), (bool(tail), tail)):
            if count:
                low, high = _chunk_base(width, bank, regimes)
                base[0] += count * low
                base[1] += count * high
        cpu = [base[0] + retained[0] * per_retained_cpu,
               base[1] + retained[1] * per_retained_cpu]
        # Exact logical traffic of the five existing selection primitives:
        # FP32 view 12b, finite 9b, nonzero 1b+8h, add 16h, gather 24h.
        dram = [22 * markers + 48 * retained[0],
                22 * markers + 48 * retained[1]]
        if not all(math.isfinite(value) for value in cpu):
            raise ValueError('JAGWAS selector CPU work overflow')
        for i in (0, 1):
            total_cpu[i].append(cpu[i])
            total_dram[i].append(dram[i])
            serial[i].append(cpu[i] / fractions[partition['device']])
        rows.append(dict(id=partition['id'], device=partition['device'],
                         markers=markers, chunks=chunks, retained_rows=list(retained),
                         cpu_seconds=cpu, logical_dram_bytes=dram,
                         serial_consumer_cpu_floor_seconds=[
                             value / fractions[partition['device']] for value in cpu]))
    cpu = [math.fsum(values) for values in total_cpu]
    dram = [math.fsum(values) for values in total_dram]
    consumer = [math.fsum(values) for values in serial]
    floor = [max(cpu[i] / cpu_cap, dram[i] / dram_cap, consumer[i])
             for i in (0, 1)]
    if not all(math.isfinite(value) for value in (*consumer, *floor)):
        raise ValueError('JAGWAS selector service overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during JAGWAS selector pricing')
    return dict(kind='torchgwas.pgen_layout_jagwas_selection_floor.v1',
                input_identity=dict(identity), reduction='jagwas',
                output_partitions_sha256=_digest(output_rows),
                partitions=rows, total_cpu_seconds=cpu,
                total_logical_dram_bytes=dram,
                shared_cpu_capacity=cpu_cap, shared_dram_capacity=dram_cap,
                single_consumer_cpu_floor_seconds=consumer,
                selector_service_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional necessary full-panel JAGWAS host selector CPU, logical DRAM and one-consumer serial work from independent primitive prices and explicit retained-count scenarios. Mixed empty/nonempty chunks widen only this partial floor. GPU projection, indexed writer, queue and final drain are outside this report; no switch is authorized.')
