"""Compact conditional host-selector work for significant-pair trait tiles."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .numpy_nonzero_work import validate_host_price_protocol
from .significant_host_work import (host_significant_selection_work,
                                    host_selection_service)


def _number(name, value, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError('Finite independent host selector price required: ' + name)
    return value


def _profile(row):
    if not isinstance(row, dict):
        raise ValueError('Explicit host selector profile required')
    try:
        q = _number('cpu_fraction', row['cpu_fraction'], positive=True)
        dram = _number('shared_dram_bytes_per_second',
                       row['shared_dram_bytes_per_second'], positive=True)
    except (KeyError, TypeError):
        raise ValueError('Complete host selector CPU/DRAM profile required') from None
    if q > 1:
        raise ValueError('Host selector CPU fraction exceeds one')
    return q, dram


def _regime_price(bank, name, cells):
    price = bank[name]
    if (not isinstance(price, dict) or
            set(price) != {'call_cpu_seconds', 'unit_cpu_seconds',
                           'dram_bytes_per_unit'}):
        raise ValueError('Complete NumPy nonzero primitive required: ' + name)
    return (_number(name + '.call', price['call_cpu_seconds']) +
            cells * _number(name + '.unit', price['unit_cpu_seconds']))


def _shape_work(markers, traits, total_cells, retained, bank, profile,
                *, threshold_one):
    cells = markers * traits
    low_h = max(0, retained[0] - (total_cells - cells))
    high_h = min(cells, retained[1])
    sparse_limit = cells // 10
    regimes = []
    if low_h == 0:
        regimes.append('flatnonzero_empty')
    if sparse_limit >= 1 and low_h <= sparse_limit and high_h >= 1:
        regimes.append('flatnonzero_sparse')
    if high_h >= max(1, sparse_limit + 1):
        regimes.append('flatnonzero_dense')
    if not regimes:
        raise ValueError('No feasible host nonzero regime')
    work = host_significant_selection_work(markers, traits, 0,
        threshold_one=threshold_one, return_beta=True)
    q, dram_rate = profile
    steps = host_selection_service(work, bank, cpu_fraction=q,
        dram_bytes_per_second=dram_rate, host_serial_fraction=0.)
    empty_cpu = math.fsum(row['seconds'] * row['resources']['cpu'] for row in steps)
    empty_dram = math.fsum(row['seconds'] * row['resources']['dram'] for row in steps)
    empty_nonzero = _regime_price(bank, 'flatnonzero_empty', cells)
    bases_cpu = [empty_cpu + _regime_price(bank, name, cells) - empty_nonzero
                 for name in regimes]
    bases_dram = [empty_dram + (0 if name == 'flatnonzero_empty' else cells)
                  for name in regimes]
    return (min(bases_cpu), max(bases_cpu),
            min(bases_dram), max(bases_dram), regimes)


def native_layout_significant_host_selection_floor(source, output, prices,
                                                    profiles_by_device):
    """Bound fixed host selector primitives from whole-tile survivor ranges.

    NumPy nonzero changes algorithm at zero and 10% occupancy. The bound
    considers every regime consistent with each chunk's possible retained
    count, then adds the exact common per-retained CPU/traffic slope. It
    expands only a full and optional tail chunk shape per trait tile.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            source.get('reduction') != 'significant' or
            source.get('partition_axis') != 'trait' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            output.get('reduction') != 'significant' or
            output.get('significant_backend') != 'host' or
            type(output.get('significant_threshold_one')) is not bool or
            output.get('input_identity') != source.get('input_identity')):
        raise ValueError('Matching threshold-bound host significant source/output required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before host selector pricing')
    validate_host_price_protocol(prices)
    bank = prices.get('prices')
    if not isinstance(bank, dict):
        raise ValueError('Complete independent host selector primitive bank required')
    # The public significant API requests beta even when the indexed writer
    # saves only t. An optional t-only archive therefore does not relax scan
    # or host-selection work.
    per_h_cpu = (_number('matrix gather unit', bank['matrix_gather_flat']['unit_cpu_seconds']) * 2 +
                 _number('coordinate unit', bank['coordinate_divmod']['unit_cpu_seconds']) +
                 _number('df gather unit', bank['df_gather_row']['unit_cpu_seconds']) +
                 2 * _number('inplace add unit', bank['inplace_index_add']['unit_cpu_seconds']) +
                 _number('index cast unit', bank['index_cast']['unit_cpu_seconds']))
    per_h_dram = 128  # 8 nonzero + 32 beta/t gather + 24 divmod + 16 each df/add/cast/rebase.
    partitions = source.get('partitions')
    output_rows = output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(output_rows, list) or len(output_rows) != len(partitions)):
        raise ValueError('One output row per significant trait tile required')
    devices = {row['device'] for row in partitions}
    if not isinstance(profiles_by_device, dict) or set(profiles_by_device) != devices:
        raise ValueError('One host selector profile per active GPU required')
    profiles = {device: _profile(profile)
                for device, profile in profiles_by_device.items()}
    cpu_cap = _number('shared CPU capacity',
                      source['shared_capacities']['cpu'], positive=True)
    dram_cap = _number('shared DRAM capacity',
                       source['shared_capacities']['dram'], positive=True)
    rows = []
    cpu_totals = [[], []]
    dram_totals = [[], []]
    device_chains = {device: [[], []] for device in devices}
    for partition, result in zip(partitions, output_rows):
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        traits = partition['trait_range'][1] - partition['trait_range'][0]
        cells = markers * traits
        block = source['chunk_markers']
        full, tail = divmod(markers, block)
        chunks = full + bool(tail)
        retained = result.get('retained_rows')
        if (result.get('id') != partition['id'] or
                result.get('device') != partition['device'] or
                result.get('markers') != markers or
                result.get('cells') != cells or
                partition['chunks'] != chunks or
                not isinstance(retained, list) or len(retained) != 2 or
                any(type(value) is not int for value in retained) or
                not 0 <= retained[0] <= retained[1] <= cells):
            raise ValueError('Host selector geometry or occupancy differs from output')
        base_cpu = [0., 0.]
        base_dram = [0., 0.]
        shape_regimes = []
        profile = profiles[partition['device']]
        for count, width in ((full, block), (bool(tail), tail)):
            if count:
                low_c, high_c, low_d, high_d, regimes = _shape_work(
                    width, traits, cells, retained, bank, profile,
                    threshold_one=output['significant_threshold_one'])
                base_cpu[0] += count * low_c
                base_cpu[1] += count * high_c
                base_dram[0] += count * low_d
                base_dram[1] += count * high_d
                shape_regimes.append(dict(markers=width, chunks=count,
                                          possible_nonzero_regimes=regimes))
        cpu = [base_cpu[i] + retained[i] * per_h_cpu for i in (0, 1)]
        dram = [base_dram[i] + retained[i] * per_h_dram for i in (0, 1)]
        if not all(math.isfinite(value) and value >= 0
                   for value in (*cpu, *dram)):
            raise ValueError('Host selector service overflow')
        for i in (0, 1):
            cpu_totals[i].append(cpu[i])
            dram_totals[i].append(dram[i])
            device_chains[partition['device']][i].append(cpu[i] / profile[0])
        rows.append(dict(id=partition['id'], device=partition['device'],
                         markers=markers, traits=traits, chunks=chunks,
                         retained_rows=list(retained), cpu_seconds=cpu,
                         logical_dram_bytes=dram,
                         possible_shapes=shape_regimes))
    cpu = [math.fsum(values) for values in cpu_totals]
    dram = [math.fsum(values) for values in dram_totals]
    chains = {device: [math.fsum(values) for values in pair]
              for device, pair in device_chains.items()}
    floor = [max(cpu[i] / cpu_cap, dram[i] / dram_cap,
                 *(pair[i] for pair in chains.values())) for i in (0, 1)]
    if not all(math.isfinite(value) for value in floor):
        raise ValueError('Host selector floor overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during host selector pricing')
    return dict(kind='torchgwas.pgen_layout_significant_host_selection_floor.v1',
                input_identity=dict(identity), reduction='significant',
                significant_backend='host', return_beta=True,
                significant_threshold_one=output['significant_threshold_one'],
                output_partitions_sha256=_digest(output_rows), partitions=rows,
                total_cpu_seconds=cpu, total_logical_dram_bytes=dram,
                per_device_serial_selector_floor_seconds=chains,
                shared_cpu_capacity=cpu_cap, shared_dram_capacity=dram_cap,
                selector_service_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional host significant selector CPU, logical DRAM and per-GPU serial work using the active NumPy nonzero protocol and independent primitive prices. Bounds all chunk regimes consistent with the tile retained-count interval. Public significant scan returns beta even for t-only NPZ output. Device selection, indexed writer, queues and final drain are separate; high bounds only this partial selector floor and no JIT switch is authorized.')
