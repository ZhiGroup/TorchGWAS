"""Compact native-scan transport/output load for fixed source partitions.

This prices only mandatory payload bytes. In particular, selected-row counts
are declared scenarios rather than inferred from an early, possibly empty
significance window. File framing, selection service, CPU, GPU, queueing and
durability remain for the full finite continuation.
"""
import math

from .analytical_plan_cache import input_identity
from .binary_output_work import compact_binary_output_work
from .layout_transfer_links import transfer_link_loads
from .reduced_output_work import significant_output_work
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS


def _capacity(name, value):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0):
        raise ValueError('Positive finite ' + name + ' capacity required')
    return value


def _retained(row, limit):
    if (not isinstance(row, (tuple, list)) or len(row) != 2 or
            any(type(value) is not int for value in row) or
            not 0 <= row[0] <= row[1] <= limit):
        raise ValueError('Bounded retained-row interval required')
    return row


def native_layout_output_floor(layout, *, store_beta=True,
                               significant_backend=None, retained_ranges=None,
                               shared_d2h_bytes_per_second,
                               per_device_d2h_bytes_per_second,
                               output_bytes_per_second,
                               dense_writer_options=None,
                               jagwas_writer_fsync=None,
                               significant_writer_fsync=None,
                               significant_threshold_one=None,
                               device_selection_max_cells=DEVICE_SELECTION_MAX_CELLS,
                               shared_links=()):
    """Bound payload resource loads for a bound native PGEN layout.

    The result is an interval of *necessary payload floors*, not an interval
    for elapsed job time. Every partition is a full fixed schedule from the
    source floor; the caller still binds those partitions to the unissued
    frontier and supplies memory admission and component availability.
    """
    if not isinstance(layout, dict) or layout.get('kind') != 'torchgwas.pgen_layout_source_floor.v1':
        raise ValueError('Typed native source layout required')
    if type(store_beta) is not bool:
        raise ValueError('Boolean beta-output setting required')
    mode = layout.get('reduction')
    if mode not in (None, 'significant', 'jagwas'):
        raise ValueError('Unknown native output mode')
    if (significant_backend not in (None, 'host', 'device') or
            (mode == 'significant') != (significant_backend is not None)):
        raise ValueError('Explicit significant selection backend required only for significant output')
    if mode == 'jagwas' and store_beta is not True:
        raise ValueError('JAGWAS has no configurable beta output')
    if ((mode != 'jagwas' and jagwas_writer_fsync is not None) or
            (jagwas_writer_fsync is not None and type(jagwas_writer_fsync) is not bool)):
        raise ValueError('Explicit JAGWAS writer fsync applies only to JAGWAS')
    if ((mode != 'significant' and significant_writer_fsync is not None) or
            (significant_writer_fsync is not None and
             type(significant_writer_fsync) is not bool)):
        raise ValueError('Explicit significant writer fsync applies only to significant output')
    if ((mode != 'significant' and significant_threshold_one is not None) or
            (significant_threshold_one is not None and
             type(significant_threshold_one) is not bool)):
        raise ValueError('Explicit threshold-one selector setting applies only to significant output')
    if (type(device_selection_max_cells) is not int or
            device_selection_max_cells < 1 or
            (significant_backend != 'device' and
             device_selection_max_cells != DEVICE_SELECTION_MAX_CELLS)):
        raise ValueError('Positive device selection cell limit applies only to device-selected output')
    partitions = layout.get('partitions')
    if not isinstance(partitions, list) or not partitions:
        raise ValueError('Nonempty fixed source partitions required')
    ids = {row.get('id') for row in partitions if isinstance(row, dict)}
    if len(ids) != len(partitions) or not all(isinstance(key, str) and key for key in ids):
        raise ValueError('Unique fixed source partitions required')
    if mode is None:
        if retained_ranges is not None:
            raise ValueError('Dense output has no retained-row scenario')
        if dense_writer_options is not None and (
                not isinstance(dense_writer_options, dict) or
                set(dense_writer_options) != {'block_bytes', 'queue_depth',
                                              'borrow_chunks', 'fsync',
                                              'writeback_bytes', 'sync_file_range',
                                              'store_variant_df'}):
            raise ValueError('Complete explicit dense writer settings required')
    elif dense_writer_options is not None:
        raise ValueError('Dense writer settings do not apply to selected output')
    elif retained_ranges is not None and (not isinstance(retained_ranges, dict) or
                                           set(retained_ranges) != ids):
        raise ValueError('One retained interval per fixed partition required')
    shared_d2h = _capacity('shared D2H', shared_d2h_bytes_per_second)
    output = _capacity('shared output', output_bytes_per_second)
    devices = {row.get('device') for row in partitions}
    if (not isinstance(per_device_d2h_bytes_per_second, dict) or
            set(per_device_d2h_bytes_per_second) != devices):
        raise ValueError('One D2H capacity per active device required')
    per_device = {device: _capacity(device + ' D2H', value)
                  for device, value in per_device_d2h_bytes_per_second.items()}
    source = layout.get('input_identity')
    if not isinstance(source, dict) or 'path' not in source or input_identity(source['path']) != source:
        raise ValueError('PGEN input changed before output-load composition')
    d2h = [[], []]
    written = [[], []]
    by_device = {device: [[], []] for device in devices}
    rows = []
    writer_rows = []
    for partition in partitions:
        if (not isinstance(partition, dict) or
                set(partition) != {'id', 'device', 'variant_range', 'trait_range', 'chunks'}):
            raise ValueError('Complete fixed source partition required')
        variant, trait = partition['variant_range'], partition['trait_range']
        if (not isinstance(variant, list) or not isinstance(trait, list) or
                len(variant) != 2 or len(trait) != 2 or
                any(type(value) is not int for value in variant + trait) or
                not 0 <= variant[0] < variant[1] or not 0 <= trait[0] < trait[1]):
            raise ValueError('Nonempty fixed partition geometry required')
        markers = variant[1] - variant[0]
        cells = markers * (trait[1] - trait[0])
        chunks = partition['chunks']
        size = layout.get('chunk_markers')
        if (type(chunks) is not int or type(size) is not int or size < 1 or
                chunks != (markers + size - 1) // size):
            raise ValueError('Source partition chunk count changed')
        selected = None
        if mode is not None:
            limit = cells if mode == 'significant' else markers
            selected = ((0, limit) if retained_ranges is None else
                        _retained(retained_ranges[partition['id']], limit))
        if mode is None:
            # Native result ring also carries one status byte and FP32 df per
            # marker. Count only mandatory binary statistic payload on output;
            # df sidecars, headers and manifest add more bytes.
            transfer = [4 * (1 + int(store_beta)) * cells + 5 * markers] * 2
            payload = [4 * (1 + int(store_beta)) * cells] * 2
            if dense_writer_options is not None:
                writer = compact_binary_output_work(
                    markers, trait[1] - trait[0], size,
                    store_beta=store_beta, **dense_writer_options)
                writer_rows.append(dict(id=partition['id'], device=partition['device'],
                                        work=writer))
                payload = [writer['payload_bytes']] * 2
        elif mode == 'significant':
            payload = [(28 if store_beta else 24) * count for count in selected]
            # The public significant scan returns beta for filtering even
            # when sumstats_fields='t' omits beta from the indexed NPZ part.
            if significant_backend == 'host':
                transfer = [8 * cells + 5 * markers] * 2
            else:
                selection = significant_output_work(
                    layout['samples'], markers, trait[1] - trait[0], size,
                    backend='device', max_selection_cells=device_selection_max_cells,
                    include_blocks=False)
                # CUDA nonzero synchronizes a four-byte count for every block,
                # including blocks with zero retained pairs.
                transfer = [markers + selection['selection_count_d2h_bytes'] +
                            28 * count for count in selected]
        else:
            transfer = [17 * markers] * 2
            payload = [16 * count for count in selected]
        for i in (0, 1):
            d2h[i].append(transfer[i])
            written[i].append(payload[i])
            by_device[partition['device']][i].append(transfer[i])
        rows.append(dict(id=partition['id'], device=partition['device'],
                         markers=markers, cells=cells,
                         device_selection_blocks=(selection['selection_blocks']
                                                  if mode == 'significant' and
                                                  significant_backend == 'device' else None),
                         maximum_selection_block_cells=(
                             selection['maximum_selection_block_cells']
                             if mode == 'significant' and significant_backend == 'device'
                             else None),
                         retained_rows=None if selected is None else list(selected),
                         d2h_payload_bytes=transfer,
                         output_array_payload_bytes=payload))
    total_d2h = [sum(values) for values in d2h]
    total_output = [sum(values) for values in written]
    device_bytes = {device: [sum(values) for values in pair]
                    for device, pair in by_device.items()}
    link_loads = [transfer_link_loads(
        {device: values[i] for device, values in device_bytes.items()},
        shared_links, direction='d2h') for i in (0, 1)]
    floor = [max(total_d2h[i] / shared_d2h,
                 total_output[i] / output,
                 link_loads[i]['floor_seconds'],
                 *(values[i] / per_device[device]
                   for device, values in device_bytes.items()))
             for i in (0, 1)]
    if not all(math.isfinite(value) for value in floor):
        raise ValueError('Output payload floor overflow')
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed during output-load composition')
    return dict(kind='torchgwas.pgen_layout_output_floor.v1',
                input_identity=dict(source), reduction=mode,
                significant_backend=significant_backend,
                store_beta=store_beta, partitions=rows,
                jagwas_writer_fsync=jagwas_writer_fsync,
                significant_writer_fsync=significant_writer_fsync,
                significant_threshold_one=significant_threshold_one,
                device_selection_max_cells=(device_selection_max_cells
                                            if significant_backend == 'device' else None),
                total_d2h_payload_bytes=total_d2h,
                total_output_array_payload_bytes=total_output,
                dense_writer_work=writer_rows if dense_writer_options is not None else None,
                per_device_d2h_payload_bytes=device_bytes,
                shared_d2h_link_loads=link_loads,
                capacities=dict(shared_d2h_bytes_per_second=shared_d2h,
                                per_device_d2h_bytes_per_second=per_device,
                                output_bytes_per_second=output),
                payload_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional necessary native result D2H and output payload floors for fixed PGEN partitions. Global, optionally declared overlapping shared links, and per-GPU D2H ceilings constrain the same transfers without summing link times. Device-selected output includes the blocking four-byte nonzero count per selection block, including empty blocks. Explicit dense writer settings add exact aggregate beta/t and optional df payload, staging, minimum write and writeback request counts without expanding chunks. Selected-row intervals are declared scenarios, never inferred from early output. Dynamic selection service, metadata, writer CPU/storage service, queueing and durable commit remain unpriced. Not elapsed-time bounds or a switch decision.')
