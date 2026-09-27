"""Source host-allocation lifetimes and independent fixed-extent services.

No association timings. Geometry is a duration-free allocator-route census;
prices are one fixed 64-MiB copy/download experiment in a named worker context.
"""
import math
import statistics


def _number(value,name,minimum=0.):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<minimum:
        raise ValueError('Invalid '+name)
    return value


def setup_host_allocations(work):
    """Phenotype-sized CPU arrays; covariate/basis and Python metadata excluded."""
    arrays=[];n,k=work['samples'],work['traits']
    residual=[(i,p) for i,p in enumerate(work['phases']) if p['phase']=='residual_block']
    def add(name,kind,size,phase,touches,release):
        if size:
            arrays.append(dict(name=name,kind=kind,bytes=size,allocate_phase=phase,
                touches=touches,release_phase=release))
    for i,p in enumerate(work['phases']):
        amount=p['input_contiguous_copy_bytes']
        add('contiguous:'+str(i),'numpy',amount,i,[dict(phase=i,bytes=amount)],i)
        if p['phase']=='residual_block':
            amount=4*n*p['traits']
            add('download:'+str(i),'torch',amount,i,[dict(phase=i,bytes=amount)],
                'after_scan' if len(residual)==1 else i)
    if len(residual)>1:
        add('assembled_result','numpy',4*n*k,0,
            [dict(phase=i,bytes=4*n*p['traits']) for i,p in residual],'after_scan')
    for a in arrays:
        if sum(t['bytes'] for t in a['touches'])!=a['bytes']:
            raise ValueError('Host allocation touch bytes not conserved')
    return arrays


def pageable_host_prices(probe,*,devices,device,numpy_version,torch_version,python_version,libc,cpu_affinity,source_sha256):
    """Median paired fixed-extent observations; all raw repeats remain visible."""
    if (probe.get('bytes')!=64<<20 or probe.get('page_bytes')!=4096
        or probe['numpy_version']!=numpy_version or probe['torch_version']!=torch_version
        or probe['python_version']!=python_version or probe['libc']!=list(libc)
        or probe['affinity']!=list(cpu_affinity) or probe['torch_threads']!=4
        or probe['numpy_madvise_hugepage'] is not False or probe['source_sha256']!=source_sha256):
        raise ValueError('Pageable host primitive source/runtime context differs')
    contexts=[c for c in probe['contexts'] if c['devices']==list(devices)]
    if len(contexts)!=1 or device not in devices:raise ValueError('Missing exact pageable worker context')
    workers=contexts[0]['workers']
    if len(set(devices))!=len(devices) or sorted(w['device'] for w in workers)!=sorted(devices):raise ValueError('Incomplete pageable workers')
    worker=next(w for w in workers if w['device']==device)
    prices={};size=probe['bytes'];pages=size//probe['page_bytes'];warnings=[]
    for operation,kind in [('numpy_copy','numpy'),('pageable_d2h','torch')]:
        rows=[r for r in worker['rows'] if r['operation']==operation]
        if len(rows)!=7 or {r['repeat'] for r in rows}!=set(range(7)):
            raise ValueError('Seven complete pageable repetitions required')
        repeats=[]
        for row in rows:
            if row['bytes']!=size or row['mapping']['heap'] or row['mapping']['mapping_bytes']<size:
                raise ValueError('Pageable extent/mapping mismatch')
            for phase in ['first','warm','release']+(['allocation'] if kind=='numpy' else []):
                for key in ['seconds','cpu_seconds']:_number(row[phase][key],phase+' '+key)
            for phase in ['first','warm']:
                for key in ['minor_faults','major_faults']:
                    if type(row[phase][key]) is not int or row[phase][key]<0:raise ValueError('Invalid page counter')
                if row[phase]['major_faults']:raise ValueError('Major faults in independent host primitive')
            delta=row['first']['cpu_seconds']-row['warm']['cpu_seconds']
            faults=row['first']['minor_faults']-row['warm']['minor_faults']
            if not math.isclose(delta,row['signed_first_less_warm_cpu'],abs_tol=1e-15) or faults!=row['signed_first_less_warm_faults']:
                raise ValueError('Pageable signed difference corrupted')
            if row['warm']['minor_faults']>4 or abs(row['first']['minor_faults']-pages)>64:
                warnings.append(dict(operation=operation,repeat=row['repeat'],reason='disturbed first/reused-page fault control'))
            repeats.append(dict(repeat=row['repeat'],warm_cpu_seconds_per_byte=row['warm']['cpu_seconds']/size,
                first_touch_cpu_seconds_per_page=delta/pages,
                release_cpu_seconds_per_page=row['release']['cpu_seconds']/pages,
                allocation_cpu_seconds=row.get('allocation',{}).get('cpu_seconds',0.),
                first_minor_faults=row['first']['minor_faults'],warm_minor_faults=row['warm']['minor_faults']))
        values={key:statistics.median(r[key] for r in repeats) for key in
            ['warm_cpu_seconds_per_byte','first_touch_cpu_seconds_per_page','release_cpu_seconds_per_page','allocation_cpu_seconds']}
        if min(values.values())<0:raise ValueError('Negative aggregate host service')
        prices[kind]=dict(values,repeat_prices=repeats)
    return dict(page_bytes=probe['page_bytes'],probe_bytes=size,devices=list(devices),device=device,
        prices=prices,control_warnings=warnings,
        context={k:probe[k] for k in ['numpy_version','torch_version','python_version','libc','affinity','torch_threads','numpy_madvise_hugepage','allocator_environment']},
        source_sha256=source_sha256,
        scope='Fixed-64-MiB per-byte/per-page scaling of complete repeat medians, including disturbed controls. Not a universal allocator or latency bound. Fresh-page increments are distinct from warm copy/transfer service.')


def validate_pageable_geometry(geometry,prices):
    context=prices['context']
    for key in ['numpy_version','torch_version','python_version','libc','affinity','torch_threads','numpy_madvise_hugepage','allocator_environment']:
        if geometry.get(key)!=context[key]:raise ValueError('Pageable allocator geometry context differs: '+key)
    if geometry.get('durations_recorded') is not False or geometry.get('page_bytes')!=prices['page_bytes']:
        raise ValueError('Untimed matching-page geometry required')
    for kind,bank in geometry['arrays'].items():
        if kind not in ('numpy','torch'):raise ValueError('Unknown pageable allocator')
        for size,row in bank.items():
            amount=int(size);obs=row['observations']
            if str(amount)!=size or amount<=0 or len(obs)!=4 or {x['repeat'] for x in obs}!=set(range(4)):
                raise ValueError('Incomplete pageable allocator geometry')
            routes=set();extents=set()
            for x in obs:
                if x['route']=='mmap':
                    if (x['allocate_count']!=1 or x['release_count']!=-1 or x['allocate_bytes']<amount
                        or x['allocate_bytes']%prices['page_bytes'] or x['release_bytes']!=-x['allocate_bytes']):
                        raise ValueError('Invalid pageable mmap counters')
                    extents.add(x['allocate_bytes'])
                elif x['route']=='arena':
                    if any(x[k] for k in ['allocate_count','allocate_bytes','release_count','release_bytes']):
                        raise ValueError('Invalid pageable arena counters')
                else:raise ValueError('Unknown pageable observation route')
                routes.add(x['route'])
            route=next(iter(routes)) if len(routes)==1 else 'variable'
            mapped=next(iter(extents)) if routes=={'mmap'} and len(extents)==1 else None
            if route=='mmap' and mapped is None:raise ValueError('Unstable pageable mapped extent')
            if row['route']!=route or row['mapped_bytes']!=mapped:raise ValueError('Pageable geometry summary differs')
    return geometry


def setup_host_service(work,inputs):
    """Place allocation/touch/destruction at source lifetime boundaries.

    Explicit arena freshness is a caller scenario, never inferred from GWAS.
    mmap bytes use observed geometry; future allocator history remains uncertain.
    """
    prices=inputs['prices'];geometry=validate_pageable_geometry(inputs['geometry'],prices)
    arena=inputs['arena_fresh_fraction']
    if isinstance(arena,bool) or not isinstance(arena,(int,float)) or not 0<=arena<=1:
        raise ValueError('Explicit arena_fresh_fraction in [0,1] required')
    page=prices['page_bytes'];rows=[dict(touch_cpu_seconds=0.,allocation_cpu_seconds=0.,
        release_cpu_seconds=0.,serial_cpu_seconds=0.,zero_bytes=0.) for _ in work['phases']]
    cleanup=dict(cpu_seconds=0.,serial_cpu_seconds=0.);records=[];unpriced=set()
    for a in setup_host_allocations(work):
        kind,size=a['kind'],a['bytes']
        if str(size) not in geometry['arrays'].get(kind,{}):raise ValueError('Missing pageable geometry: '+kind+' '+str(size))
        route=geometry['arrays'][kind][str(size)];p=prices['prices'][kind]
        if route['route']=='mmap':fresh=1.;release_pages=route['mapped_bytes']/page
        else:
            fresh=float(arena);release_pages=0.
            unpriced.add('arena allocation/release and future allocator route history remain conditional')
        # Only bytes beyond the fixed N32/T1 CPU download baseline are new.
        baseline=128 if kind=='torch' else 0
        touch_pages=max(0,(size-baseline)/page)*fresh
        allocate=p['allocation_cpu_seconds'] if route['route']=='mmap' and kind=='numpy' else 0.
        rows[a['allocate_phase']]['allocation_cpu_seconds']+=allocate
        rows[a['allocate_phase']]['serial_cpu_seconds']+=allocate
        for t in a['touches']:
            fraction=t['bytes']/size;row=rows[t['phase']]
            row['touch_cpu_seconds']+=touch_pages*fraction*p['first_touch_cpu_seconds_per_page']
            if kind=='numpy' and t['bytes']//4<=500:
                row['serial_cpu_seconds']+=touch_pages*fraction*p['first_touch_cpu_seconds_per_page']
            row['zero_bytes']+=touch_pages*fraction*page
        release=release_pages*p['release_cpu_seconds_per_page']
        serial=release if kind=='numpy' else 0.
        if a['release_phase']=='after_scan':
            cleanup['cpu_seconds']+=release;cleanup['serial_cpu_seconds']+=serial
        else:
            row=rows[a['release_phase']];row['release_cpu_seconds']+=release;row['serial_cpu_seconds']+=serial
        records.append(dict(a,route=route['route'],fresh_fraction=fresh,touch_pages=touch_pages,
            allocation_cpu_seconds=allocate,release_cpu_seconds=release))
    return dict(phases=rows,cleanup=cleanup,allocations=records,
        warm_copy_cpu_seconds_per_byte=prices['prices']['numpy']['warm_cpu_seconds_per_byte'],
        warm_d2h_cpu_seconds_per_byte=prices['prices']['torch']['warm_cpu_seconds_per_byte'],
        unpriced_terms=sorted(unpriced|{'fixed-64MiB page/copy/release scaling and allocator-history changes',
            'multi-block assembled-result first touches distributed by written bytes',
            'pageable H2D staging, covariate/basis host allocations and NUMA migration faults'}))
