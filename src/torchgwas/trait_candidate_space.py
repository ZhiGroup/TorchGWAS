"""Bounded candidate construction for the existing detailed trait planner.

Contexts contain independently measured prices and duration-free geometry.
Only geometry is varied; worker/device context prices are never synthesized.
"""
import copy
import itertools
import math
from pathlib import Path

from .mechanistic_plan import _device, _integer
from .pgen_reader import read_header
from .linear import multigpu_variant_ranges
from .pgen_work_census import census, rechunk_census
from .trait_tiling_model import trait_tiled_shape
from .trait_tiling_plan import detailed_trait_plan


def _axis(name, values, maximum=None):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(name+' must be an explicit nonempty list')
    for value in values:
        _integer(name, value)
        if maximum is not None and value > maximum:
            raise ValueError(name+' exceeds workload size')
    if len(set(values)) != len(values):
        raise ValueError(name+' must not contain duplicates')
    return sorted(values)


def _tile_profile(template, n, m, k, c, chunk, *, work_shapes=None):
    """Copy the independent prices, retaining only this tile's compiled work.

    Geometry banks may cover many sample/chunk/trait shapes. Copying the whole
    bank into each repeated genotype pass scales with unrelated shapes. Keep
    duplicate matches so the existing coverage checks still reject ambiguity.
    Profiles and retained kernel records remain independent mutable copies.
    """
    sizes={min(chunk,m)}
    if m%chunk:sizes.add(m%chunk)
    if work_shapes is not None:
        sizes=set(work_shapes)
    view=dict(template)
    if 'kernel_geometry' in view:
        view['kernel_geometry']=[row for row in view['kernel_geometry']
            if row.get('N')==n and row.get('B') in sizes and row.get('K')==k
            and row.get('C')==c and row.get('validate_range') is True]
    if 'joint_kernel_geometry' in view:
        view['joint_kernel_geometry']=[row for row in view['joint_kernel_geometry']
            if row.get('N')==n and row.get('B') in sizes and row.get('K')==k
            and row.get('compute_dtype')=='float32']
    return copy.deepcopy(view)


def prepare_trait_candidates(workload, contexts, *, chunks, trait_blocks, output,
                             max_candidates=1000, max_candidate_tiles=10000,
                             max_census_chunks=1000000, partition_axes=('trait',), reduction=None,
                             census_chunk_size=None, memory_only=False,
                             _header_receiver=None, _prepared_header=None,
                             _index_receiver=None, _compact_memory=False,
                             _compact_cache_dir=None, _compact_cache_receiver=None,
                             _defer_compact_bases=False):
    """Expand finite axes, sharing exact counts across aligned chunk grids.

    Each named context supplies devices, per-device profiles and shared
    capacities, with optional shared links/storage. Profile worker allocation,
    depth and prices remain those of the measured context. Only chunk_markers
    changes. Full-file encoded censuses are shared read-only across tiles;
    profiles and output settings are independent copies.

    This prepares proposals; missing geometry is reported explicitly and
    prevents bounded_trait_plan from selecting a partial supported subset.
    File/header reads occur only after validating all combinatorial budgets.
    """
    required = {'genotype', 'samples', 'markers', 'traits', 'covariates',
                'matching_sample_order', 'complete_phenotypes'}
    if set(workload)-{'phenotype_c_contiguous','covariate_columns'} != required:
        raise ValueError('Explicit genotype, dimensions and complete/matching input declarations required')
    if 'phenotype_c_contiguous' in workload and type(workload['phenotype_c_contiguous']) is not bool:
        raise ValueError('phenotype_c_contiguous must be boolean when supplied')
    n, m, k, c = [workload[key] for key in ['samples', 'markers', 'traits', 'covariates']]
    for name, value, minimum in [('samples', n, 32), ('markers', m, 1), ('traits', k, 1), ('covariates', c, 0)]:
        _integer(name, value, minimum)
    columns=workload.get('covariate_columns',c)
    _integer('covariate_columns',columns,c)
    if columns>=n-2 or workload['matching_sample_order'] is not True or workload['complete_phenotypes'] is not True:
        raise ValueError('Detailed trait context requires positive df and complete matching samples/phenotypes')
    path = workload['genotype']
    if not isinstance(path, str) or not path:
        raise ValueError('Explicit PGEN path required')
    if (not isinstance(partition_axes,(list,tuple)) or not partition_axes
        or any(axis not in ('trait','variant') for axis in partition_axes)
        or len(set(partition_axes))!=len(partition_axes)):
        raise ValueError('Explicit unique trait/variant partition axes required')
    chunks = _axis('chunks', chunks)
    if type(memory_only) is not bool or (memory_only and census_chunk_size is None):
        raise ValueError('Memory-only construction requires an explicit source grid')
    if type(_compact_memory) is not bool or (_compact_memory and not memory_only):
        raise ValueError('Compact source vectors require memory-only construction')
    if type(_defer_compact_bases) is not bool or (_defer_compact_bases and not _compact_memory):
        raise ValueError('Deferred base indexing requires compact memory admission')
    if _header_receiver is not None and not callable(_header_receiver):
        raise ValueError('Prepared header receiver must be callable')
    if _index_receiver is not None and (not memory_only or not callable(_index_receiver)):
        raise ValueError('Validated index receiver requires memory-only construction')
    census_sizes=chunks
    if census_chunk_size is not None:
        _integer('census_chunk_size',census_chunk_size)
        if any(chunk%census_chunk_size for chunk in chunks):
            raise ValueError('Census grid must divide every candidate chunk size')
        census_sizes=sorted(set(chunks)|{census_chunk_size})
    widths = _axis('trait_blocks', trait_blocks, k)
    if reduction not in (None, 'jagwas'):
        raise ValueError('Unsupported candidate-space reduction')
    if reduction == 'jagwas':
        if tuple(partition_axes) != ('variant',) or widths != [k]:
            raise ValueError('JAGWAS requires variant-only partitioning and the full phenotype panel')
        if k > n-c-1:
            raise ValueError('JAGWAS trait count exceeds residual phenotype rank')
    for name, value in [('max_candidates', max_candidates), ('max_candidate_tiles', max_candidate_tiles),
                        ('max_census_chunks', max_census_chunks)]:
        _integer(name, value)
    if not isinstance(contexts, (list, tuple)) or not contexts:
        raise ValueError('Nonempty measured context list required')
    total = len(chunks)*len(contexts)*((len(widths) if 'trait' in partition_axes else 0)+('variant' in partition_axes))
    if total > max_candidates:
        raise ValueError('Search exceeds max_candidates before census')
    trait_tile_count=sum(math.ceil(k/width) for width in widths)*len(chunks)*len(contexts) if 'trait' in partition_axes else 0
    if trait_tile_count > max_candidate_tiles:
        raise ValueError('Search exceeds max_candidate_tiles before census')
    if sum(math.ceil(m/chunk) for chunk in census_sizes) > max_census_chunks:
        raise ValueError('Search exceeds max_census_chunks before census')
    names = set()
    for context in contexts:
        if set(context)-{'name', 'devices', 'profiles', 'shared_capacities',
                         'shared_transfer_capacities', 'shared_links',
                         'shared_storage_bytes_per_second'}:
            raise ValueError('Unknown measured context fields')
        transfer=context.get('shared_transfer_capacities')
        if transfer is not None and (not isinstance(transfer,dict) or
                set(transfer)!={'h2d','d2h'} or any(
                    isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or value<=0
                    for value in transfer.values())):
            raise ValueError('Positive shared H2D/D2H context capacities required')
        name = context['name']
        if not isinstance(name, str) or not name or name in names:
            raise ValueError('Unique nonempty context names required')
        names.add(name)
        devices = context['devices']
        if not isinstance(devices, (list, tuple)) or not devices or len(set(devices)) != len(devices):
            raise ValueError('Unique explicit context devices required')
        for device in devices:
            _device(device)
        if set(context['profiles']) != set(devices):
            raise ValueError('Exactly one independent profile per active context device required')
        if reduction == 'jagwas':
            for profile in context['profiles'].values():
                if (profile.get('reduction') != 'jagwas' or profile.get('result_ownership') != 'owned'
                        or profile.get('borrow_results', False)
                        or profile.get('compute_dtype', 'float32') != 'float32'):
                    raise ValueError('JAGWAS requires independently priced native FP32 owned-result profiles')
    range_requests=set()
    if 'variant' in partition_axes:
        if trait_tile_count+len(chunks)*sum(len(c['devices']) for c in contexts)>max_candidate_tiles:
            raise ValueError('Search exceeds max_candidate_tiles before census')
        for context,chunk in itertools.product(contexts,chunks):
            spans=multigpu_variant_ranges(m,chunk,len(context['devices']))
            if len(spans)==len(context['devices']):
                range_requests.update((chunk,lo,hi) for lo,hi in spans if (lo,hi)!=(0,m))
        if sum(math.ceil(m/chunk) for chunk in census_sizes)+sum(math.ceil((hi-lo)/chunk) for chunk,lo,hi in range_requests)>max_census_chunks:
            raise ValueError('Search exceeds max_census_chunks before census')
    if set(output) != {'block_bytes', 'queue_depth', 'store_beta', 'fsync'}:
        raise ValueError('Explicit complete writer settings required')
    if output['block_bytes'] is not None:
        _integer('block_bytes', output['block_bytes'])
    _integer('queue_depth', output['queue_depth'])
    if type(output['store_beta']) is not bool or output['fsync'] is not True:
        raise ValueError('Boolean fields and durable fsync output required')
    if reduction == 'jagwas' and (output['block_bytes'] is not None or output['store_beta'] is not False):
        raise ValueError('JAGWAS requires indexed output without beta or dense block coalescing')
    path = str(Path(path).resolve(strict=True))
    from .analytical_plan_cache import input_identity as full_input_identity
    full_identity=full_input_identity(path)
    def file_identity():
        stat = Path(path).stat()
        return dict(path=path, bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                    device=stat.st_dev, inode=stat.st_ino)
    input_identity = file_identity()
    if _prepared_header is None:
        header = read_header(path)
    else:
        from .pgen_reader import PgenHeader
        if (not isinstance(_prepared_header,tuple) or len(_prepared_header)!=2
                or _prepared_header[0]!=full_identity
                or not isinstance(_prepared_header[1],PgenHeader)):
            raise ValueError('Prepared PGEN header differs from candidate input')
        header = _prepared_header[1]
    if (header.sample_ct, header.variant_ct) != (n, m):
        raise ValueError('Declared workload differs from PGEN header')
    encoded={};file_passes=0;validated_bases=[];compact_source=None
    if memory_only:
        from .pgen_memory_layout import memory_layout,rechunk_memory_layout,compact_rechunk_memory_layout
    admission_cache=None
    if _compact_memory and _compact_cache_dir is not None:
        from .pgen_admission_cache import PgenAdmissionCache
        admission_cache=PgenAdmissionCache(_compact_cache_dir,full_identity,header,census_chunk_size)
    if _compact_memory:
        cached=None if admission_cache is None else admission_cache.load()
        if cached is None:
            compact_source=memory_layout(path,census_chunk_size,header=header,compact=True,
                _validated_bases_receiver=(validated_bases.append if _index_receiver is not None
                    and not _defer_compact_bases else None),
                _defer_bases=_defer_compact_bases)
        else:
            compact_source,bases=cached
            validated_bases.append(bases)
        encoded={chunk:compact_rechunk_memory_layout(compact_source,chunk) for chunk in census_sizes}
        ranged={(chunk,lo,hi):compact_rechunk_memory_layout(compact_source,chunk,(lo,hi))
            for chunk,lo,hi in sorted(range_requests)}
    else:
        regroup=rechunk_memory_layout if memory_only else rechunk_census
        for chunk in census_sizes:
            finer=next((size for size in reversed(encoded) if chunk%size==0),None)
            if finer is None:
                if memory_only:encoded[chunk]=memory_layout(path,chunk,header=header,
                    _validated_bases_receiver=validated_bases.append if _index_receiver is not None else None)
                else:encoded[chunk]=census(path,chunk,include_chunks=True);file_passes+=1
            else:encoded[chunk]=regroup(encoded[finer],chunk)
        ranged = {(chunk,lo,hi):regroup(encoded[chunk],chunk,variant_range=(lo,hi))
                  for chunk,lo,hi in sorted(range_requests)}
    if file_identity() != input_identity:
        raise ValueError('PGEN changed during candidate preparation')
    candidates, assignments, excluded, required_geometry, missing_geometry = [], [], [], [], []
    geometry_seen = set()
    for context, chunk, width in itertools.product(contexts, chunks, widths if 'trait' in partition_axes else []):
        devices = list(context['devices'])
        count = math.ceil(k/width)
        identity = dict(context=context['name'], chunk_size=chunk, trait_block=width)
        if len(devices) > count:
            excluded.append(dict(**identity, reason='idle_devices'))
            continue
        tiles = []
        for index, start in enumerate(range(0, k, width)):
            size, device = min(width, k-start), devices[index % len(devices)]
            profile = _tile_profile(context['profiles'][device],n,m,size,c,chunk)
            profile['return_beta']=output['store_beta'] if reduction is None else True
            profile['chunk_markers'] = chunk
            data = dict(samples=n, markers=m, traits_analyzed=size, covariates=c,
                        matching_sample_order=True, encoded=encoded[chunk])
            if columns!=c:data['covariate_columns']=columns
            if 'phenotype_c_contiguous' in workload:
                data['phenotype_c_contiguous']=workload['phenotype_c_contiguous'] and size==k
            tiles.append(dict(trait_range=[start, start+size], device=device, data=data, profile=profile))
            blocks = {min(chunk, m)}
            if m % chunk:
                blocks.add(m % chunk)
            for b in sorted(blocks):
                key = (context['name'], device, n, b, size, c)
                if key in geometry_seen:
                    continue
                geometry_seen.add(key)
                request = dict(context=context['name'], device=device, shape=[n, b, size, c], validate_range=True)
                required_geometry.append(request)
                matches = [row for row in profile.get('kernel_geometry', [])
                    if [row.get(field) for field in ['N', 'B', 'K', 'C']] == [n, b, size, c]
                    and row.get('validate_range') is True]
                if len(matches) != 1:
                    missing_geometry.append(dict(**request, reason='missing' if not matches else 'ambiguous'))
        candidate = dict(tiles=tiles, trait_block=width, devices=devices,
                         shared_capacities=copy.deepcopy(context['shared_capacities']), output=copy.deepcopy(output))
        for field in ['shared_links', 'shared_storage_bytes_per_second']:
            if field in context:
                candidate[field] = copy.deepcopy(context[field])
        trait_tiled_shape(candidate)
        assignments.append(dict(candidate_index=len(candidates), **identity))
        candidates.append(candidate)
    if 'variant' in partition_axes:
        for context,chunk in itertools.product(contexts,chunks):
            devices=list(context['devices']);spans=multigpu_variant_ranges(m,chunk,len(devices))
            proposal=dict(context=context['name'],chunk_size=chunk,trait_block=None,partition_axis='variant')
            if len(spans)!=len(devices):
                excluded.append(dict(**proposal,reason='idle_devices'));continue
            tiles=[]
            for device,(lo,hi) in zip(devices,spans):
                profile=_tile_profile(context['profiles'][device],n,hi-lo,k,c,chunk)
                profile['return_beta']=output['store_beta'] if reduction is None else True
                profile['chunk_markers']=chunk
                data=dict(samples=n,markers=hi-lo,traits_analyzed=k,covariates=c,matching_sample_order=True,
                    encoded=encoded[chunk] if (lo,hi)==(0,m) else ranged[chunk,lo,hi])
                if columns!=c:data['covariate_columns']=columns
                if 'phenotype_c_contiguous' in workload:data['phenotype_c_contiguous']=workload['phenotype_c_contiguous']
                if reduction == 'jagwas':data['phenotype_complete']=True
                tiles.append(dict(trait_range=[0,k],variant_range=[lo,hi],device=device,data=data,profile=profile))
                sizes={min(chunk,hi-lo)}
                if (hi-lo)%chunk:sizes.add((hi-lo)%chunk)
                for b in sorted(sizes):
                    key=(context['name'],device,n,b,k,c)
                    if key in geometry_seen:continue
                    geometry_seen.add(key)
                    request=dict(context=context['name'],device=device,shape=[n,b,k,c],validate_range=True)
                    required_geometry.append(request)
                    matches=[row for row in profile.get('kernel_geometry',[]) if
                        [row.get(field) for field in ['N','B','K','C']]==[n,b,k,c] and row.get('validate_range') is True]
                    if len(matches)!=1:missing_geometry.append(dict(**request,reason='missing' if not matches else 'ambiguous'))
                    if reduction == 'jagwas':
                        joint_request=dict(context=context['name'],device=device,shape=[n,b,k],
                                           reduction='jagwas',compute_dtype='float32')
                        required_geometry.append(joint_request)
                        joint_matches=[row for row in profile.get('joint_kernel_geometry',[]) if
                            [row.get(field) for field in ['N','B','K']]==[n,b,k]
                            and row.get('compute_dtype')=='float32']
                        if len(joint_matches)!=1:
                            missing_geometry.append(dict(**joint_request,
                                reason='missing' if not joint_matches else 'ambiguous'))
            candidate=dict(tiles=tiles,trait_block=k,devices=devices,partition_axis='variant',
                shared_capacities=copy.deepcopy(context['shared_capacities']),output=copy.deepcopy(output))
            for field in ['shared_links','shared_storage_bytes_per_second']:
                if field in context:candidate[field]=copy.deepcopy(context[field])
            trait_tiled_shape(candidate,reduction=reduction)
            assignments.append(dict(candidate_index=len(candidates),**proposal));candidates.append(candidate)
    if file_identity()!=input_identity or full_input_identity(path)!=full_identity:
        raise ValueError('PGEN changed during candidate preparation')
    if _header_receiver is not None or _index_receiver is not None:
        for array in (header.vblock_offsets,header.vrtypes,header.record_offsets,header.record_lengths):
            array.flags.writeable=False
    if _header_receiver is not None:
        _header_receiver((full_identity,header))
    compact_bases=None
    if _index_receiver is not None:
        if len(validated_bases)>1 or (not validated_bases and not _defer_compact_bases):
            raise ValueError('Memory admission did not certify one source index')
        if validated_bases:
            # PGEN variant indices are uint32 on disk; retaining int64 locators
            # would double the persistent JIT metadata for no added range.
            compact_bases=validated_bases[0].astype('uint32',copy=False)
            compact_bases.flags.writeable=False
            _index_receiver((full_identity,header,compact_bases))
    if admission_cache is not None:
        admission_cache.stage(compact_source,compact_bases)
        if _compact_cache_receiver is not None:_compact_cache_receiver(admission_cache)
    if not candidates:
        raise ValueError('No active candidate in the requested context/trait space')
    result=dict(candidates=candidates, assignments=assignments, excluded=excluded,input_file_identity=input_identity,
        required_geometry=required_geometry, missing_geometry=missing_geometry,
        raw_candidates=total, census_passes=file_passes, census_views=len(encoded)+len(ranged), census_chunks=sum(x.get('logical_chunks',len(x['chunks'])) for x in [*encoded.values(),*ranged.values()]),
        scope='Explicit bounded axes over independently priced device/worker contexts. Each chunk-size census is '
              'computed once; every tile retains a full genotype pass and exact tails. No association timing or '
              'coarse-score fallback. Complete-input, source and live calibration context checks remain required before execution.')
    if 'variant' in partition_axes:
        result['scope']='Bounded trait/variant partitions over independently priced device contexts. Trait tiles reread the complete genotype input; variant shards census and process disjoint ranges with the complete phenotype panel on each device. No association timing or coarse-score fallback. Input/source/live-context checks remain required.'
    if reduction == 'jagwas':
        result['reduction']='jagwas'
        result['scope']='Bounded variant partitions with a full phenotype panel and independent joint factor on every device. Exact PGEN censuses and separate statistics/projection geometry; source, live context and memory admission remain required.'
    if census_chunk_size is not None:
        result['source_census']=compact_source if _compact_memory else encoded[census_chunk_size]
    if memory_only:
        result['memory_only']=True
        result['scope']='Header-only reader extents for memory admission; genotype record payloads and decoder operation counts were not read. This candidate cannot be runtime-ranked until populated with exact work evidence.'
    return result


def bounded_trait_plan(workload, contexts, *, bounds, joint, output):
    """Construct proposals then use the existing detailed minimax planner."""
    space = prepare_trait_candidates(workload, contexts, output=output, **bounds)
    # The detailed planner checks reader/memory budgets before geometry.
    # An infeasible whole-panel candidate must not demand an impossible capture.
    # Missing geometry on any feasible candidate remains an error below.
    result = detailed_trait_plan(space['candidates'], **joint)
    unsupported = [row for row in result['rejected'] if row['reason'] == 'unsupported_model_context']
    if unsupported:
        raise ValueError('Requested feasible geometry lacks model coverage: '+str(unsupported))
    result['search_space'] = {key: value for key, value in space.items() if key != 'candidates'}
    return result
