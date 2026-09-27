"""Indexed-reduction dependencies in the shared analytical execution graph.

All service demands must be supplied independently. This module introduces no
rates, observed GWAS durations, selectivity estimate or autotune authorization.
"""
from .execution_graph import ExecutionGraph, torch_scan_schedule
from .mechanistic_plan import _integer
from .selection_geometry import DEVICE_SELECTION_MAX_CELLS


def _steps(graph, prefix, steps, after):
    if not isinstance(steps, list):
        raise ValueError('Explicit ordered primitive steps required')
    previous = list(after)
    for index, step in enumerate(steps):
        if not isinstance(step, dict) or set(step) - {'seconds', 'resources'} or 'seconds' not in step:
            raise ValueError('Primitive steps require seconds and optional resource demands')
        name = graph.add(prefix + ':' + str(index), step['seconds'], previous, step.get('resources'))
        previous = [name]
    return graph.add(prefix + ':done', after=previous)


class _IndexedConsumer:
    def __init__(self, outputs, queue_service, cleanup, *, mode, selection_graphs=None):
        self.outputs, self.queue_service, self.cleanup = outputs, queue_service, cleanup
        self.selection_graphs = selection_graphs
        self.gpu_tail = None
        self.writes = []
        self.mode=mode
        self.queue_key=mode+':queue';self.writer_key=mode+':writer';self.fifo=mode+':fifo'
        self.consumer_selection=mode=='jagwas'

    def append(self, graph, index, after):
        previous = list(after)
        selection = None
        if self.selection_graphs is not None:
            selection = self.selection_graphs[index]
            selection_prefix = f'{self.mode}:{index}:device:'
            graph.compose(selection['graph'],selection_prefix,after=previous)
            self.gpu_tail = selection_prefix+selection['gpu_tail']
        for part, output in enumerate(self.outputs[index]):
            base = f'{self.mode}:{index}:{part}'
            selected = (selection_prefix+selection['yield_nodes'][part] if selection is not None else
                        None if self.consumer_selection and self.queue_service is not None else
                        _steps(graph, base + ':select', output['selection'], previous))
            if self.queue_service is None:
                written = _steps(graph, base + ':write', output['writer'], [selected])
                previous = [written]
            else:
                # Acquire one queue credit before put. Credit lasts until the
                # single writer starts get, not until durable part completion.
                acquire = graph.add(base + ':acquire', after=previous if selected is None else [selected])
                graph.token_actions[acquire] = {'acquire': {self.queue_key: 1}}
                put = _steps(graph, base + ':put', self.queue_service['put'], [acquire])
                get_start = graph.add(base + ':get_start', after=[put])
                graph.token_actions[get_start] = {'acquire': {self.writer_key: 1}, 'release_start': {self.queue_key: 1}}
                graph.fifo_enqueues[put] = (self.fifo, get_start)
                graph.fifo_dequeues[get_start] = self.fifo
                got = _steps(graph, base + ':get', self.queue_service['get'], [get_start])
                if self.consumer_selection:
                    got = _steps(graph, base + ':select', output['selection'], [got])
                written = _steps(graph, base + ':write', output['writer'], [got])
                graph.token_actions[written] = {'release_finish': {self.writer_key: 1}}
                previous = [put]
            if selection is not None:
                resume = selection_prefix+selection['resume_nodes'][part]
                seconds,deps = graph.nodes[resume]
                graph.nodes[resume] = (seconds,tuple(dict.fromkeys((*deps,*previous))))
                previous = [resume]
            self.writes.append(written)
        graph.add(f'consume:{index}', after=previous)

    def close(self, graph, after):
        # The producer may start another phenotype tile with old selected
        # parts still queued. Its input reads/CUDA streams must first drain.
        last = len(self.outputs) - 1
        dependencies = list(after) + [f'release:{last}', f'finish:{last}']
        if self.gpu_tail is not None: dependencies.append(self.gpu_tail)
        cleaned = _steps(graph, self.mode+':cleanup', self.cleanup, dependencies)
        graph.add(self.mode+':producer_complete', after=[cleaned])


def _indexed_reduction_schedule(tiles, *, mode, shared_prepare=None, queue_depth, shared_capacities,
                               queue_service, finalize, host_serial_policy='fluid',
                               return_graph=False, max_source_chunks=10000, max_selection_blocks=100000):
    """Shared queue/writer topology, with reduction-specific selection placement."""
    if mode not in ('significant','jagwas'):
        raise ValueError('Explicit indexed reduction mode required')
    _QUEUE,_WRITER,_FIFO=(mode+suffix for suffix in (':queue',':writer',':fifo'))
    if mode=='jagwas' and not isinstance(shared_prepare,ExecutionGraph):
        raise ValueError('Explicit shared JAGWAS preprocessing graph required')
    if not isinstance(tiles, list) or not tiles:
        raise ValueError('At least one significant tile required')
    _integer('max_source_chunks', max_source_chunks)
    _integer('max_selection_blocks', max_selection_blocks)
    if sum(len(tile['blocks']) for tile in tiles) > max_source_chunks:
        raise ValueError('Significant schedule exceeds max_source_chunks')
    if sum(len(group) for tile in tiles for group in tile['outputs']) > max_selection_blocks:
        raise ValueError('Significant schedule exceeds max_selection_blocks')
    devices = list(dict.fromkeys(tile['device'] for tile in tiles))
    if mode=='jagwas' and (len(devices)!=len(tiles) or any(tile['backend']!='host' for tile in tiles)):
        raise ValueError('JAGWAS requires one variant shard per active device and asynchronous narrow host results')
    if mode=='jagwas' and any(not isinstance(tile.get('prepare'),ExecutionGraph) for tile in tiles):
        raise ValueError('Explicit per-device factor/design preparation graph required')
    multiple = len(devices) > 1
    if multiple:
        _integer('queue_depth', queue_depth)
        if not isinstance(queue_service, dict) or set(queue_service) != {'put', 'get'}:
            raise ValueError('Independent queue put/get service required')
    elif queue_depth != 0 or queue_service is not None:
        raise ValueError('Single active device uses synchronous writing without a result queue')
    graph = ExecutionGraph()
    graph.capacities = dict(shared_capacities)
    if multiple:
        graph.token_capacities.update({_QUEUE: queue_depth, _WRITER: 1})
    prepared_shared = [] if shared_prepare is None else [graph.compose(shared_prepare,'shared_prepare:')]
    producer_done = {}
    all_writes = []
    selected_pairs = parts = selections = 0
    for index, tile in enumerate(tiles):
        if tile['backend'] not in ('host', 'device'):
            raise ValueError('Explicit significant selection backend required')
        blocks, outputs = tile['blocks'], tile['outputs']
        if not blocks or len(blocks) != len(outputs):
            raise ValueError('One output group per nonempty source chunk required')
        selection_graphs = tile.get('selection_graphs')
        if selection_graphs is not None:
            if mode!='significant' or tile['backend']!='device' or len(selection_graphs)!=len(blocks):
                raise ValueError('Source selection graphs require one graph per device significant chunk')
            for group,selection in zip(outputs,selection_graphs):
                if not isinstance(selection.get('graph'),ExecutionGraph):
                    raise ValueError('Source selection graph required')
                if selection.get('retained_per_block') != [row['retained'] for row in group]:
                    raise ValueError('Source selection occupancy differs from writer parts')
                for field in ['yield_nodes','resume_nodes']:
                    nodes = selection.get(field,[])
                    if len(nodes)!=len(group) or len(set(nodes))!=len(nodes) or any(n not in selection['graph'].nodes for n in nodes):
                        raise ValueError('Distinct source yield/resume ports required')
                if selection.get('gpu_tail') not in selection['graph'].nodes:
                    raise ValueError('Source CUDA stream tail required')
                if any(row['selection'] for row in group):
                    raise ValueError('Do not price selection twice beside its source graph')
        for block, group in zip(blocks, outputs):
            if not isinstance(group, list) or not group:
                raise ValueError('Every source chunk must yield selection blocks including empty ones')
            if block.get('consumer_seconds', 0.) or block.get('discard_seconds', 0.):
                raise ValueError('Selection steps must explicitly price consumer and owned-buffer cleanup')
            if mode=='jagwas' and len(group)!=1:
                raise ValueError('JAGWAS emits exactly one narrow result per source chunk')
            for output in group:
                if set(output) != {'cells', 'retained', 'selection', 'writer'}:
                    raise ValueError('Exact selected occupancy and primitive stages required')
                _integer('cells', output['cells'])
                _integer('retained', output['retained'], 0)
                if output['retained'] > output['cells']:
                    raise ValueError('Selected occupancy exceeds cells')
                if not output['retained'] and output['writer']:
                    raise ValueError('Empty selections do not create or fsync a part')
                if output['retained'] and not output['writer']:
                    raise ValueError('Nonempty part writer service is required')
                if not output['selection'] and selection_graphs is None:
                    raise ValueError('Selection service is required even for empty results')
                if tile['backend'] == 'device' and output['cells'] > DEVICE_SELECTION_MAX_CELLS:
                    raise ValueError('Device selection exceeds executor block bound')
                selections += 1
                selected_pairs += output['retained']
                parts += output['retained'] > 0
        consumer = _IndexedConsumer(outputs, queue_service, tile['cleanup'],mode=mode,selection_graphs=selection_graphs)
        local = torch_scan_schedule(blocks, depth=tile['depth'], decode_workers=tile['decode_workers'],
            consumer=consumer, return_graph=True, synchronous_results=tile['backend'] == 'device')
        if tile.get('prepare') is not None:
            roots = [name for name, (_, deps) in local.nodes.items() if not deps]
            prepared = local.compose(tile['prepare'], 'prepare:')
            for name in roots:
                seconds, deps = local.nodes[name]
                local.nodes[name] = (seconds, (*deps, prepared))
        if multiple:
            local.token_capacities.update({_QUEUE: queue_depth, _WRITER: 1})
        prefix = f'tile:{index}:'
        graph.compose(local, prefix, after=producer_done.get(tile['device'], prepared_shared), shared_tokens=(_QUEUE, _WRITER))
        # compose namespaces ordinary FIFOs. This executor instead has one
        # arrival-ordered result queue shared across every active tile worker.
        for name, (fifo, target) in local.fifo_enqueues.items():
            if fifo == _FIFO:
                graph.fifo_enqueues[prefix + name] = (_FIFO, prefix + target)
        for name, fifo in local.fifo_dequeues.items():
            if fifo == _FIFO:
                graph.fifo_dequeues[prefix + name] = _FIFO
        producer_done[tile['device']] = [prefix + mode+':producer_complete']
        all_writes.extend(prefix + name for name in consumer.writes)
    if multiple:
        # Each worker publishes one terminal sentinel through the same bounded
        # queue after its final tile. It creates no output part.
        for index, device in enumerate(devices):
            terminal = ExecutionGraph()
            terminal.token_capacities.update({_QUEUE: queue_depth, _WRITER: 1})
            consumer = _IndexedConsumer([[dict(selection=[], writer=[])]], queue_service, [],mode=mode)
            consumer.append(terminal, 0, [])
            prefix = f'terminal:{index}:'
            graph.compose(terminal, prefix, after=producer_done[device], shared_tokens=(_QUEUE, _WRITER))
            for name, (_, target) in terminal.fifo_enqueues.items():
                graph.fifo_enqueues[prefix + name] = (_FIFO, prefix + target)
            for name in terminal.fifo_dequeues:
                graph.fifo_dequeues[prefix + name] = _FIFO
            all_writes.extend(prefix + name for name in consumer.writes)
    _steps(graph, mode+':finalize', finalize, all_writes + [name for names in producer_done.values() for name in names])
    if host_serial_policy != 'fluid':
        graph = graph.with_serial_sections(host_serial_policy)
    if return_graph:
        return graph
    result = graph.solve()
    result.update(selected_pairs=selected_pairs, parts=parts, selections=selections, devices=devices,
        queue_depth=queue_depth, scope='Conditional primitive-service graph including shared queue, single durable writer and final drain. No fitted rates or autotune authorization.')
    return result

def significant_trait_schedule(tiles, **kwargs):
    """Selected pairs are filtered in each producer before the shared queue.

    Primitive prices and exact retained counts are mandatory. Each device may
    process multiple phenotype tiles. Host and device selection preserve their
    original native-ring and synchronous-yield behavior respectively.
    """
    return _indexed_reduction_schedule(tiles,mode='significant',**kwargs)


def jagwas_variant_schedule(shards, *, shared_prepare, **kwargs):
    """Joint per-variant results share preprocessing and one indexed writer.

    Each shard supplies factor/design preparation and narrow native scan blocks.
    The GPU joint reduction belongs in those blocks, before D2H. Host conversion,
    finite filtering and indexed writing occur after dequeue under the single
    consumer token, including empty chunks. No timing coefficient is supplied by
    this function; missing services must be rejected by the caller.
    """
    result=_indexed_reduction_schedule(shards,mode='jagwas',shared_prepare=shared_prepare,**kwargs)
    if isinstance(result,dict):
        result['retained_variants']=result.pop('selected_pairs')
        result.update(reduction='jagwas',shared_preprocessing_passes=1,
            scope='Conditional primitive-service graph: one shared preprocessing pass, per-device joint setup/scan, bounded queue and one durable indexed consumer. No timing coefficients or autotune authorization.')
    return result
