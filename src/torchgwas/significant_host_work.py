"""Source-counted host significant selection and bounded native staging.

Complete FP32 phenotypes and native int8 PGEN only. Shared shape/setup/scan
accounting comes from the existing calculator; result occupancy never relaxes
memory admission. Timing prices are supplied by independent primitive probes.
"""
import math
import numpy as np
from .mechanistic_plan import _integer
from .trait_tiling_model import trait_tiled_shape, _tile_memory
from .decoder_work import native_reader_workspace
from .reduced_output_work import significant_execution_layout
from .numpy_nonzero_work import nonzero_regime


SELECTION_ALLOCATION_UNPRICED_TERMS = (
    'selector allocation-inclusive primitive prices require matching reference and candidate page residency; fresh-page cost cannot simply be added twice',
    'selector mask, selected-result and trait-rebase allocation routes, first-touch and final release are not independently priced',
)


def indexed_part_work(rows, *, store_beta=True):
    """Exact significant-pair NPZ bytes using the shared indexed writer ledger."""
    if type(store_beta) is not bool:
        raise ValueError('Boolean store_beta required')
    from .reduced_output_work import _indexed_array_part_work
    fields = [('variant_index', '<i8'), ('trait_index', '<i8')]
    if store_beta:
        fields.append(('beta', '<f4'))
    fields.extend([('t_stat', '<f4'), ('df', '<f4')])
    return _indexed_array_part_work(rows, fields)

def host_significant_selection_work(markers, traits, retained, *, threshold_one=False, return_beta=True):
    _integer('markers', markers)
    _integer('traits', traits)
    _integer('retained', retained, 0)
    if type(threshold_one) is not bool:
        raise ValueError('Boolean threshold_one required')
    if type(return_beta) is not bool:raise ValueError('Boolean return_beta required')
    cells = markers * traits
    if retained > cells:
        raise ValueError('Retained pairs exceed selection cells')
    from .host_significance import predicate_block_shape, host_selector, NATIVE_HOST_SELECTOR, PREDICATE_MAX_CELLS
    height, width, calls = predicate_block_shape(markers, traits)
    selector = host_selector()
    native = selector == NATIVE_HOST_SELECTOR
    if native: height, width, calls = markers, traits, 1
    # Track byte extents independently of allocation timing. NumPy may return a
    # view of its owned nonzero buffer; that storage still needs one allocation.
    # Zero-length arrays have zero data bytes, not zero dispatch/metadata cost.
    nonzero_primitive='flatnonzero_'+nonzero_regime(cells,retained)
    allocations=[dict(name='mask',bytes=cells,primitive='mask_allocate',lifetime='selector_temporary'),
        dict(name='variant_index',bytes=8*retained,primitive=nonzero_primitive,lifetime='selected_result')]
    if return_beta:
        allocations.append(dict(name='beta',bytes=4*retained,primitive='matrix_gather_flat',lifetime='selected_result'))
    allocations.extend([
        dict(name='t_stat',bytes=4*retained,primitive='matrix_gather_flat',lifetime='selected_result'),
        dict(name='trait_index',bytes=8*retained,primitive='coordinate_divmod',lifetime='selected_result'),
        dict(name='df',bytes=4*retained,primitive='df_gather_row',lifetime='selected_result'),
        dict(name='rebased_trait_index',bytes=8*retained,primitive='index_cast',lifetime='rebased_result')])
    # Variant rebasing is in place. Outer phenotype-tile rebasing copies once
    # to owned int64 storage, then adds its offset in place as well.
    return dict(markers=markers, traits=traits, cells=cells, retained=retained,
        threshold_one=threshold_one, selector=selector,
        predicate_max_cells=PREDICATE_MAX_CELLS, predicate_block_shape=[height,width],
        predicate_calls=calls, critical_lookup_rows=markers, critical_round_rows=markers,
        predicate_primitive='predicate_native' if native else 'predicate_block',
        predicate_logical_bytes=5*cells+4*markers if native else 25*cells,
        predicate_cells=cells, nonzero_cells=cells, nonzero_output_indices=retained,
        nonzero_primitive=nonzero_primitive,
        coordinate_divmod_elements=retained, matrix_gather_elements=(1+int(return_beta))*retained,
        matrix_gather_calls=1+int(return_beta), matrix_gather_primitive='matrix_gather_flat',
        df_gather_elements=retained, df_gather_primitive='df_gather_row', index_cast_elements=retained,
        inplace_index_add_elements=retained, index_add_elements=retained,
        trait_rebase_primitive='inplace_index_add',
        selected_array_bytes=(24+4*int(return_beta))*retained,
        allocations=allocations,allocation_data_bytes=sum(a['bytes'] for a in allocations),
        allocation_scope='FP32 row-df contiguous selector mask/results and outer trait rebase only; cutoff and bounded predicate scratch remain in their primitive prices. Byte extents do not imply allocator routes, fresh pages or release-thread ownership.',
        allocation_unpriced_terms=list(SELECTION_ALLOCATION_UNPRICED_TERMS),
        source='Bounded FP32 host predicate, contiguous flat payload gathers before divmod, row-only df gather and trait-coordinate rebasing')


def significant_host_memory(candidate, *, host_reserve_bytes=0, device_reserve_bytes=0,
                            device_memory_profiles=None):
    """Worst-case tensor/array staging, independently of predicted survivors.

    Reserves have the same explicit role as in the full-output calculator:
    mappings, ID/index metadata, allocator retention and library uncertainties
    are not hidden in a fitted multiplier. No claim of allocator-reserved bound.
    """
    from .decoder_work import require_regular_memory
    require_regular_memory(candidate)
    shape = trait_tiled_shape(candidate)
    if shape.get('partition_axis', 'trait') != 'trait':
        raise ValueError('Significant variant partition is not an executable API mode')
    _integer('host_reserve_bytes', host_reserve_bytes, 0)
    _integer('device_reserve_bytes', device_reserve_bytes, 0)
    layout = significant_execution_layout(shape['traits'], shape['trait_block'], shape['devices'],
        shape['reader_workers'], queue_depth=candidate['output']['queue_depth'])
    host_peaks, gpu_peaks, bins = {}, {}, {}
    selection_peaks = {}
    largest_cells = 0
    unknown = set()
    for tile in candidate['tiles']:
        device, data, profile = tile['device'], tile['data'], tile['profile']
        n, m, k, c = (data[key] for key in ['samples', 'markers', 'traits_analyzed', 'covariates'])
        b, depth = profile['chunk_markers'], profile['depth']
        rows = min(b, m)
        cells = rows * k
        largest_cells = max(largest_cells, cells)
        payload, replay = native_reader_workspace(data['encoded']['chunks'])
        host, gpu, pins, terms = _tile_memory(n, m, k, c, data.get('covariate_columns', c),
            b, depth, profile['decode_workers'], payload, None,
            (device_memory_profiles or {}).get(device), replay,return_beta=profile.get('return_beta',True))
        host_peaks[device] = max(host_peaks.get(device, 0), host)
        gpu_peaks[device] = max(gpu_peaks.get(device, 0), gpu)
        cache = bins.setdefault(device, {})
        for size, count in pins.items():
            cache[size] = max(cache.get(size, 0), count)
        # Sum distinct predicate/nonzero/gather/index-rebase arrays as if they
        # coexist, plus a whole previous iteration to cover Python references
        # surviving RHS evaluation and suspended producer frames.
        selection = 2 * (7 * cells + (56+4*int(profile.get('return_beta',True))) * cells + 64 * rows)
        selection_peaks[device] = max(selection_peaks.get(device, 0), selection)
        unknown.update(terms)
    pinned = {d:sum(size * count for size, count in counts.items()) for d, counts in bins.items()}
    payload_per_pair=max(24+4*int(tile['profile'].get('return_beta',True)) for tile in candidate['tiles'])
    queue = layout['queue_depth'] * payload_per_pair * largest_cells
    # Old consumer item may survive evaluation of next(); a newly dequeued
    # item and the old writer df copy may coexist until assignment completes.
    consumer = 2 * payload_per_pair * largest_cells + 4 * largest_cells
    # NumPy's largest member is an int64 index array. Its iterator buffer
    # is bounded independently of archive length. Do not validate ZIP file
    # offsets here: an oversized proposal must reach memory rejection before
    # output pricing is attempted.
    archive = 2 * min(8 * largest_cells, 16 << 20)
    shared_table = 8 * (shape['samples'] + 1)
    arrays = sum(host_peaks.values()) + sum(pinned.values()) + sum(selection_peaks.values())
    arrays += queue + consumer + archive + shared_table
    unknown.update(['input mmap residency, genotype/variant/NPZ-part metadata and QC require host reserve',
        'Python object/header overhead, allocator caches and fragmentation require host reserve',
        'GPU allocator/driver retention and pinned allocator fragmentation require explicit reserves'])
    return dict(host_bytes=arrays + host_reserve_bytes,
        device_bytes={d:value + device_reserve_bytes for d,value in gpu_peaks.items()},
        host_active_bytes_by_device=host_peaks, pinned_cache_bytes_by_device=pinned,
        selection_bytes_by_device=selection_peaks, queued_result_bytes=queue,
        consumer_result_bytes=consumer, archive_buffer_bytes=archive, shared_critical_table_bytes=shared_table,
        host_reserve_bytes=host_reserve_bytes, device_reserve_bytes=device_reserve_bytes,
        occupancy='all variant-trait pairs retained in every active selection block',
        unresolved_memory_terms=sorted(unknown), layout=layout)


def host_selection_service(work, prices, *, cpu_fraction, dram_bytes_per_second, host_serial_fraction):
    """Fixed dispatch plus source units; no scan-duration or density fitting."""
    if not 0 < cpu_fraction <= 1 or dram_bytes_per_second <= 0 or not 0 <= host_serial_fraction <= 1:
        raise ValueError('Explicit positive CPU/DRAM and bounded host serialization required')
    rows = []
    # DRAM terms are logical source traffic, not measured transactions. Cached
    # threshold/index accesses may cost less; candidate durations are not used.
    terms = [
        ('critical_one' if work['threshold_one'] else 'critical_lookup', 1, work['critical_lookup_rows'], None),
        ('critical_round', 1, work['critical_round_rows'], 50*work['critical_round_rows']),
        ('mask_allocate', 1, 0, 0),
        (work['predicate_primitive'], work['predicate_calls'], work['predicate_cells'], work['predicate_logical_bytes']),
        (work['nonzero_primitive'], 1, work['nonzero_cells'],
         (1 if not work['retained'] else 2)*work['nonzero_cells'] + 8*work['retained']),
        (work['matrix_gather_primitive'], work['matrix_gather_calls'], work['matrix_gather_elements'], 16*work['matrix_gather_elements']),
        ('coordinate_divmod', 1, work['coordinate_divmod_elements'], 24*work['coordinate_divmod_elements']),
        (work['df_gather_primitive'], 1, work['df_gather_elements'], 16*work['df_gather_elements']),
        ('inplace_index_add', 1, work['inplace_index_add_elements'], 16*work['inplace_index_add_elements']),
        ('index_cast', 1, work['index_cast_elements'], 16*work['index_cast_elements']),
        (work['trait_rebase_primitive'], 1, work['index_add_elements'], 16*work['index_add_elements'])]
    for name, calls, count, traffic in terms:
        if name not in prices:
            raise ValueError('Missing bounded-flat selector primitive: '+name)
        price = prices[name]
        if set(price) != {'call_cpu_seconds', 'unit_cpu_seconds', 'dram_bytes_per_unit'}:
            raise ValueError('Primitive price requires fixed CPU, unit CPU and source DRAM work')
        if any(isinstance(value, bool) or not math.isfinite(value) or value < 0 for value in price.values()):
            raise ValueError('Invalid independent primitive price')
        cpu = calls * price['call_cpu_seconds'] + count * price['unit_cpu_seconds']
        if traffic is None:
            traffic = count * (29 if work['threshold_one'] else 138)
        seconds = max(cpu / cpu_fraction, traffic / dram_bytes_per_second)
        rows.append(dict(seconds=seconds, resources=dict(cpu=cpu/seconds if seconds else 0.,
            host_serial=cpu*host_serial_fraction/seconds if seconds else 0., dram=traffic/seconds if seconds else 0.)))
    return rows
