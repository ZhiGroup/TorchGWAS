"""Independent NumPy allocator CPU prices; numeric-copy service is separate."""
import math
import statistics


def number(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError('Finite nonnegative '+name+' required')
    return value


def allocator_service_prices(probe,*,workers,numpy_version,libc,cpu_affinity):
    if (probe.get('context_verified') is not True or probe.get('numpy_version')!=numpy_version
        or probe.get('libc')!=list(libc) or probe.get('affinity')!=list(cpu_affinity)
        or probe.get('numpy_madvise_hugepage') is not False):
        raise ValueError('Allocator runtime, affinity or advice context mismatch')
    extent,page=probe['bytes'],probe['page_bytes']
    if extent!=1<<26 or page!=4096:raise ValueError('Expected fixed 64 MiB / 4 KiB allocator probe')
    contexts=[row for row in probe['results'] if row['worker_count']==workers]
    if len(contexts)!=1 or sorted(w['worker'] for w in contexts[0]['workers'])!=list(range(workers)):
        raise ValueError('Complete unique allocator worker context required')
    groups={};partitions={};coverage=[];mappings=set()
    for worker in contexts[0]['workers']:
        if worker.get('context_verified') is not True:raise ValueError('Unverified allocator worker')
        for name,intervals in [('native',1),('python',0)]:
            control=worker['controls'][name]
            cpu=number(control['cpu_seconds'],'control CPU');detached=number(control['detached_cpu_seconds'],'control detached CPU')
            if (control['balance_errors']!=0 or control['detached_intervals']!=intervals or not cpu
                or detached>cpu or (intervals and detached/cpu<=.9) or (not intervals and detached)):
                raise ValueError('Allocator GIL control failed')
        meters={}
        for row in worker['meter_rows']:
            pair=meters.setdefault(row['repeat'],{})
            recording=row['recording']
            if (not isinstance(recording,bool) or recording in pair or row['intervals']!=20000
                or row['detached_intervals']!=(20000 if recording else 0) or row['balance_errors']!=0):
                raise ValueError('Allocator meter control failed')
            pair[recording]=number(row['cpu_seconds'],'meter CPU')
        if len(meters)<3 or any(set(pair)!={False,True} for pair in meters.values()):
            raise ValueError('Paired allocator meter controls required')
        seen=set()
        for row in worker['rows']:
            repeat,iteration,recording,size=(row[key] for key in ['repeat','iteration','recording','bytes'])
            if any(type(v) is not int or v<0 for v in (repeat,iteration)) or type(recording) is not bool or size not in (0,extent) or row['worker']!=worker['worker']:
                raise ValueError('Invalid allocator sample coordinates')
            key=(repeat,iteration,recording,size)
            if key in seen:raise ValueError('Duplicate allocator sample')
            seen.add(key);empty=number(row['empty_cpu_seconds'],'empty timer CPU')
            for phase in ['allocate','copy','release']:
                sample=row[phase];cpu=number(sample['cpu_seconds'],phase+' CPU')
                detached=number(sample['detached_cpu_seconds'],phase+' detached CPU')
                intervals=sample['detached_intervals']
                if type(intervals) is not int or intervals<0 or sample['balance_errors']!=0 or detached>cpu:
                    raise ValueError('Unbalanced allocator GIL sample')
                if not recording and (intervals or detached):raise ValueError('Disabled allocator meter reports intervals')
                if phase!='copy' and (intervals or detached):raise ValueError('Allocator service is not wholly GIL-held')
                if recording and size and phase=='copy' and (intervals<1 or detached/cpu<=.9):
                    raise ValueError('Numeric copy positive control failed')
                if phase!='copy':groups.setdefault((repeat,recording,size,phase),[]).append(max(0.,cpu-empty))
            if workers==1:
                mapping=row['mapping']
                if size:
                    if (mapping['allocate_count']!=1 or mapping['release_count']!=-1
                        or mapping['allocate_bytes']!=extent+page or mapping['release_bytes']!=-extent-page):
                        raise ValueError('Allocator mmap/unmap evidence differs')
                    mappings.add(mapping['allocate_bytes'])
                elif any(mapping.values()):raise ValueError('Zero-element control unexpectedly maps pages')
        repeats={key[0] for key in seen};iterations={key[1] for key in seen}
        if len(repeats)<3 or not iterations or seen!={(r,i,recording,size) for r in repeats for i in iterations for recording in [False,True] for size in [0,extent]}:
            raise ValueError('Incomplete paired allocator sample coverage')
        coverage.append(seen)
    if any(keys!=coverage[0] for keys in coverage):raise ValueError('Allocator worker coverage differs')
    # Arithmetic means conserve expensive observations; no median of individual
    # calls is multiplied by source counts. Zero-element CPU removes already
    # priced tiny finish allocation/destruction bookkeeping once per array.
    for phase in ['allocate','release']:
        records=[];deltas=[];paired=[]
        for repeat in sorted({key[0] for key in groups}):
            pair={record:statistics.fmean(groups[repeat,record,extent,phase])-statistics.fmean(groups[repeat,record,0,phase]) for record in [False,True]}
            # A difference of two noisy positive measurements can be negative.
            # Retain it in the repeat estimator; neither discard nor clamp it.
            # Only the aggregate used as a physical service must be nonnegative.
            paired.append(dict(recording_enabled=pair[True],recording_disabled=pair[False]))
            records.append(pair[True])
            deltas.append(pair[True]/pair[False]-1 if pair[False]>0 else None)
        estimate=statistics.median(records)
        if estimate<0:raise ValueError('Negative aggregate additional allocator service')
        partitions[phase]=dict(additional_cpu_seconds=estimate,repeat_mean_seconds=records,
            paired_relative_deltas=deltas,paired_repeat_mean_seconds=paired)
    return dict(worker_count=workers,page_bytes=page,probe_bytes=extent,mapping_bytes=extent+page,
        allocate_cpu_seconds_per_mmap=partitions['allocate']['additional_cpu_seconds'],
        release_cpu_seconds_per_mapped_page=partitions['release']['additional_cpu_seconds']/((extent+page)//page),
        observations=partitions,serial_fraction=1.,
        scope='Additional CPU beyond zero-element NumPy calls; median of complete repeat arithmetic means. Allocation per mapping and release per mapped page are explicit fixed-extent scaling assumptions, not a latency bound. Numeric copy is not charged here.')


def owned_allocator_service(work,prices,geometry):
    """Attach measured mmap work to source array extents and lifetime roles."""
    page=prices['page_bytes'];rows=geometry['array_bytes'];arrays={};unpriced=set()
    for name,size in work['array_bytes'].items():
        if str(size) not in rows:raise ValueError('Missing untimed allocator geometry for array bytes '+str(size))
        route=rows[str(size)]
        allocate=release=0.
        if route['route']=='mmap':
            mapping=route['mapped_bytes']
            if type(mapping) is not int or mapping<size or mapping%page:raise ValueError('Invalid mapped extent')
            allocate=number(prices['allocate_cpu_seconds_per_mmap'],'allocation price')
            release=(mapping//page)*number(prices['release_cpu_seconds_per_mapped_page'],'release page price')
        elif route['route'] in ('arena','variable'):
            unpriced.add('arena reuse, allocation/free service beyond the tiny finish baseline')
        else:raise ValueError('Unknown allocator route')
        arrays[name]=dict(bytes=size,route=route['route'],allocate_cpu_seconds=allocate,release_cpu_seconds=release,
            release_owner='consumer' if name in work.get('consumer_arrays',('beta','t')) else 'finish_worker')
    return dict(arrays=arrays,allocate_cpu_seconds=sum(row['allocate_cpu_seconds'] for row in arrays.values()),
        worker_release_cpu_seconds=sum(row['release_cpu_seconds'] for row in arrays.values() if row['release_owner']=='finish_worker'),
        consumer_release_cpu_seconds=sum(row['release_cpu_seconds'] for row in arrays.values() if row['release_owner']=='consumer'),
        unpriced_terms=sorted(unpriced|{'allocator route can change with subsequent history',
            'fixed-64MiB allocator scaling, munmap call overhead and CPU placement beyond measured context'}))


def validate_allocator_geometry(geometry,*,numpy_version,libc,cpu_affinity,allocator_environment):
    if (geometry.get('numpy_version')!=numpy_version or geometry.get('libc')!=list(libc)
        or geometry.get('affinity')!=list(cpu_affinity) or geometry.get('numpy_madvise_hugepage') is not False
        or geometry.get('allocator_handler')!='default_allocator'
        or geometry.get('allocator_environment')!=allocator_environment):
        raise ValueError('Allocator geometry context mismatch')
    if not geometry['array_bytes']:raise ValueError('Allocator geometry is empty')
    for size,row in geometry['array_bytes'].items():
        size=int(size);observations=row['observations'];routes=set();mappings=set()
        if size<1 or len(observations)!=4 or {obs['repeat'] for obs in observations}!={0,1,2,3}:
            raise ValueError('Incomplete untimed allocator geometry')
        for obs in observations:
            if obs['route']=='mmap':
                amount=obs['allocate_bytes']
                if (type(amount) is not int or amount<size or amount%4096 or obs['allocate_count']!=1
                    or obs['release_count']!=-1 or obs['release_bytes']!=-amount):
                    raise ValueError('Invalid allocator geometry mmap counters')
                mappings.add(amount)
            elif obs['route']=='arena':
                if any(obs[key] for key in ['allocate_bytes','allocate_count','release_bytes','release_count']):
                    raise ValueError('Invalid allocator geometry arena counters')
            else:raise ValueError('Unknown allocator geometry observation')
            routes.add(obs['route'])
        expected=next(iter(routes)) if len(routes)==1 else 'variable'
        if row['route']!=expected or row['mapped_bytes']!=(next(iter(mappings)) if routes=={'mmap'} and len(mappings)==1 else None):
            raise ValueError('Allocator geometry summary does not match observations')
    return geometry
