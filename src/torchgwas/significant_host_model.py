"""Host significant candidates in the existing setup/scan/resource calculator."""
import copy
import math
from .execution_graph import ExecutionGraph
from .decoder_work import scan_chunk_count
from .mechanistic_plan import _integer
from .setup_work import setup_work
from .trait_tiling_model import trait_tiled_shape, _prepare_graph, _tile_pin_state, _cached_scan_work
from .significant_schedule import significant_trait_schedule
from .significant_host_work import (host_significant_selection_work, host_selection_service,
                                    indexed_part_work, significant_host_memory)
from .reduced_output_work import significant_execution_layout


def _positive(name, value, zero=False):
    if isinstance(value, bool) or not math.isfinite(value) or value < 0 or (not zero and value == 0):
        raise ValueError('Invalid independent price: ' + name)
    return value


def _writer_service(part, bank, profile, serial):
    if not part['rows']:
        return []
    price = bank[str(len(part['arrays']) == 5)]
    q = profile['cpu_fraction']
    fixed = _positive('archive call CPU', price['call_cpu_seconds'], True)
    bulk = _positive('archive byte CPU', price['byte_cpu_seconds'], True)
    page = _positive('pagecache copy CPU', profile['writeback_service']['pagecache_seconds_per_byte'], True)
    memcpy = _positive('NumPy copy CPU', profile['process_units']['numpy_copy_bytes'], True)
    # The writer makes one owned df copy before NPZ serialization. NPZ rates
    # come from a seekable discarding sink, excluding OS page-cache copies.
    # A seekable ZIP writer submits each local header twice: initial creation
    # and the CRC/size rewrite. Page-cache copying sees both submissions;
    # final storage extent and durable transfer count the header only once.
    submitted = part['file_bytes'] + sum(row['local_header_bytes'] for row in part['arrays'])
    cpu = fixed + bulk * part['array_payload_bytes'] + page * submitted + 4 * part['rows'] * memcpy
    traffic = 5 * part['array_payload_bytes'] + 8 * part['rows']
    traffic += 2 * (submitted - part['array_payload_bytes'])
    seconds = max(cpu / q, traffic / profile['shared_dram_bytes_per_second'])
    transfer = _positive('independent storage byte service', profile['writeback_service']['storage_seconds_per_byte'])
    return [dict(seconds=seconds, resources=dict(cpu=cpu/seconds if seconds else 0.,
                host_serial=serial*cpu/seconds if seconds else 0., dram=traffic/seconds if seconds else 0.)),
        dict(seconds=part['file_bytes']*transfer, resources=dict(output=1/transfer)),
        dict(seconds=_positive('per-file fsync', profile['fsync_seconds'], True))]


def significant_host_runtime(candidate, prices, *, occupancy, host_serial_fraction, significance_threshold=None,
                             host_serial_policy='fluid', return_graph=False,
                             max_source_chunks=10000, max_selection_blocks=100000):
    """Explicit survivor scenario, host selection and one durable indexed writer.

    occupancy is 'empty', 'dense', or one exact retained count per tile/chunk.
    These are declared scenarios, never inferred from the significance level.
    Final API metadata and startup are reported outside the priced boundary.
    """
    from .numpy_nonzero_work import validate_host_price_protocol
    validate_host_price_protocol(prices)
    shape = trait_tiled_shape(candidate)
    if shape.get('partition_axis', 'trait') != 'trait':
        raise ValueError('Significant variant partition is not an executable API mode')
    if candidate['output']['block_bytes'] is not None:
        raise ValueError('Indexed parts do not use dense block_bytes coalescing')
    if not 0 <= host_serial_fraction <= 1:
        raise ValueError('Invalid host serialization scenario')
    if (isinstance(occupancy, str) and occupancy not in ('empty', 'dense')) or not isinstance(occupancy, (str, list, tuple)):
        raise ValueError('Explicit empty/dense or per-tile survivor counts required')
    if not isinstance(occupancy, str) and len(occupancy) != len(candidate['tiles']):
        raise ValueError('One survivor-count sequence per phenotype tile required')
    layout = significant_execution_layout(shape['traits'],shape['trait_block'],shape['devices'],shape['reader_workers'],
        queue_depth=candidate['output']['queue_depth'])
    caps = dict(candidate['shared_capacities'], host_serial=1.)
    storage = candidate.get('shared_storage_bytes_per_second')
    if storage is not None:
        caps['storage'] = _positive('shared storage capacity', storage)
    links = candidate.get('shared_links', [])
    for index, link in enumerate(links):
        if not set(link['devices']) <= set(shape['devices']):
            raise ValueError('Unknown shared link device')
        for direction in ['h2d','d2h']:
            caps[f'link:{index}:{direction}'] = _positive('shared link capacity',link[direction+'_bytes_per_second'])
    tiles, reports, pins, unresolved = [], [], {}, set()
    chunk_total = sum(scan_chunk_count(tile['data'],tile['profile']) for tile in candidate['tiles'])
    if chunk_total > max_source_chunks or chunk_total > max_selection_blocks:
        raise ValueError('Significant candidate exceeds bounded graph expansion')
    for index, tile in enumerate(candidate['tiles']):
        data, profile, device = tile['data'], tile['profile'], tile['device']
        n,m,k,c = (data[key] for key in ['samples','markers','traits_analyzed','covariates'])
        work = _cached_scan_work(data, profile)
        blocks = copy.deepcopy(work['blocks'])
        if work['result_ownership'] != 'borrowed':
            raise ValueError('Significant host executor borrows its dense native ring')
        if not isinstance(occupancy, str) and len(occupancy[index]) != len(blocks):
            raise ValueError('One exact retained count per source chunk required')
        outputs = []
        payload = parts = 0
        for chunk, block in enumerate(blocks):
            rows = block['markers']
            count = 0 if occupancy == 'empty' else rows*k if occupancy == 'dense' else occupancy[index][chunk]
            _integer('retained', count, 0)
            selection = host_significant_selection_work(rows,k,count,threshold_one=significance_threshold == 1.,return_beta=profile.get('return_beta',True))
            unresolved.update(selection['allocation_unpriced_terms'])
            part = indexed_part_work(count,store_beta=candidate['output']['store_beta'])
            outputs.append([dict(cells=rows*k,retained=count,
                selection=host_selection_service(selection,prices['prices'],cpu_fraction=profile['cpu_fraction'],
                    dram_bytes_per_second=profile['shared_dram_bytes_per_second'],host_serial_fraction=host_serial_fraction),
                writer=_writer_service(part,prices['archive'],profile,host_serial_fraction))])
            payload += part['file_bytes']; parts += bool(count)
            block['host_resources']['host_serial'] = block['host_resources']['cpu'] * host_serial_fraction
            for direction in ['h2d','d2h']:
                capacity = _positive('device '+direction, profile[direction+'_bytes_per_second'])
                caps[device+':'+direction] = capacity
                resources = block.setdefault(direction+'_resources', {})
                resources[device+':'+direction] = block[direction+'_bytes']/block[direction+'_seconds'] if block[direction+'_seconds'] else 0.
                for li, link in enumerate(links):
                    if device in link['devices']:
                        resources[f'link:{li}:{direction}'] = resources[device+':'+direction]
        pages, cached = _tile_pin_state(tile,pins)
        prep, estimate = _prepare_graph(setup_work(n,k,c,reuse_observed_counts=True,
            input_contiguous=data.get('phenotype_c_contiguous'),covariate_columns=data.get('covariate_columns')),
            profile,host_serial_fraction,pages,cached)
        for resources in prep.demands.values():
            for direction in ['h2d','d2h']:
                rate = resources.pop('prep_'+direction,0.)
                resources[device+':'+direction] = rate
                for li, link in enumerate(links):
                    if device in link['devices']:
                        resources[f'link:{li}:{direction}'] = rate
        q = profile['cpu_fraction']; clean = estimate['cleanup']; seconds = clean['cpu_seconds']/q
        cleanup = [dict(seconds=seconds,resources=dict(cpu=q,host_serial=clean['serial_cpu_seconds']/seconds if seconds else 0.))]
        tiles.append(dict(device=device,backend='host',blocks=blocks,outputs=outputs,depth=work['depth'],
            decode_workers=work['workers'],cleanup=cleanup,prepare=prep))
        reports.append(dict(trait_range=tile['trait_range'],device=device,source_chunks=len(blocks),
            indexed_part_bytes=payload,parts=parts,setup=estimate,pin_fresh_pages=pages,pin_cached_calls=cached))
        unresolved.update(work['unpriced_terms']+work['allocator_unpriced_terms']+estimate['unpriced_terms'])
    queue = None
    if len(shape['devices']) > 1:
        queue_q = min(tile['profile']['cpu_fraction'] for tile in candidate['tiles'])
        queue = {kind:[dict(seconds=_positive('queue '+kind,prices['queue_cpu_seconds'][kind],True)/queue_q,
            resources=dict(cpu=queue_q,host_serial=queue_q))] for kind in ['put','get']}
    graph = significant_trait_schedule(tiles,queue_depth=layout['queue_depth'],shared_capacities=caps,
        queue_service=queue,finalize=[],return_graph=True,max_source_chunks=max_source_chunks,max_selection_blocks=max_selection_blocks)
    for resources in graph.demands.values():
        if storage is not None:
            resources['storage'] = resources.get('input',0.) + resources.get('output',0.)
    if host_serial_policy != 'fluid':
        graph = graph.with_serial_sections(host_serial_policy)
    if return_graph:
        return graph
    from .resource_balance import resource_balance
    solved = graph.solve()
    unresolved.update(['independent fixed-shape CPU primitive transfer to candidate geometry and survivor layout',
        'NPZ extent allocation and metadata fsync; modeled storage commit follows serialization',
        'kernel early writeback overlap and dirty-page CPU/DRAM pressure',
        'blocked result-queue wakeup and timeout retry CPU service',
        'API QC, shared critical-value/covariate preparation, final IDs/metadata and directory publication'])
    return dict(status='development_significant_host_candidate',estimated_tile_seconds=solved['seconds'],
        resource_balance=resource_balance(graph, solved),
        indexed_part_bytes=sum(r['indexed_part_bytes'] for r in reports),parts=sum(r['parts'] for r in reports),
        tiles=reports,layout=layout,unpriced_terms=sorted(unresolved),prediction_complete=False,
        selection_validated=False,scope='Source scan and setup plus independently priced host selection, bounded global queue and single durable part writer. No measured GWAS duration enters this estimate.')
