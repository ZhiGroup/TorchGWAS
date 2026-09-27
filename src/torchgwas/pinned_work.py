"""Pinned allocations for the unreduced native int8 pipeline.

PyTorch 2.5 CachingHostAllocatorImpl::allocate uses PowerOf2Ceil(size)
for a fresh allocation. Driver backing, cached blocks and fragmentation may
consume additional memory; this is a request ledger, not an OS memory bound.
"""
def fresh_pin_prices(probe, *, torch_version, cpu_affinity, device, torch_threads=4):
    """Price fixed generic fresh 64 MiB allocations from complete retained runs.

    No regression or scan input. Within-run arithmetic means preserve total
    service; the median across fresh processes defines the supplied scenario.
    First allocation and other sizes are retained as diagnostics, not fitted.
    """
    import math
    import statistics
    runs=probe['runs']
    if len(runs)<3 or sorted(r['repeat'] for r in runs)!=list(range(len(runs))):
        raise ValueError('Complete unique fresh pinned repeats required')
    prices=[]
    for run in runs:
        if (run['torch_version']!=torch_version or run['affinity']!=list(cpu_affinity)
            or run['device']!=str(device) or run.get('torch_threads')!=torch_threads
            or run.get('all_allocations_retained') is not True):
            raise ValueError('Fresh pinned primitive context differs')
        rows=run['rows']
        if ([r['bytes'] for r in rows]!=[4096]*4+[1<<20]*4+[64<<20]*4
            or [r['index'] for r in rows]!=list(range(12))
            or len({r['pointer'] for r in rows})!=12):
            raise ValueError('Fresh pinned coverage or retained pointer audit failed')
        for r in rows:
            if r['pages']!=r['bytes']//4096 or any(isinstance(r[k],bool) or not math.isfinite(r[k]) or r[k]<=0 for k in ['cpu_seconds','wall_seconds']):
                raise ValueError('Invalid fresh pinned observation')
        bulk=rows[-4:];pages=sum(r['pages'] for r in bulk)
        cpu=sum(r['cpu_seconds'] for r in bulk)/pages
        signed=sum(r['wall_seconds']-r['cpu_seconds'] for r in bulk)/pages
        prices.append(dict(repeat=run['repeat'],cpu_seconds_per_page=cpu,
                           signed_non_cpu_seconds_per_page=signed))
    non_cpu=statistics.median(r['signed_non_cpu_seconds_per_page'] for r in prices)
    return dict(pin_cpu_seconds_per_page=statistics.median(r['cpu_seconds_per_page'] for r in prices),
        pin_driver_seconds_per_page=max(0.,non_cpu),signed_non_cpu_seconds_per_page=non_cpu,repeat_prices=prices,
        scope='Fresh fixed 64 MiB retained allocations; mean per registered 4 KiB page within each process, '
              'median across processes. No scan timings or fitted intercept. First-call/size transfer and concurrent registration remain unresolved.')


def pinned_scan_work(samples, chunk, traits, depth, *, transfer_bytes_per_variant=None, reduction=None, return_beta=True):
    for name,value in [('samples',samples),('chunk',chunk),('traits',traits),('depth',depth)]:
        if isinstance(value,bool) or not isinstance(value,int) or value<1:
            raise ValueError(name+' must be a positive integer')
    if reduction not in (None,'jagwas','device_significant'):
        raise ValueError('Unsupported pinned reduction layout')
    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced pinned layout')
    row_bytes=samples if transfer_bytes_per_variant is None else transfer_bytes_per_variant
    if isinstance(row_bytes,bool) or not isinstance(row_bytes,int) or row_bytes<1:
        raise ValueError('transfer_bytes_per_variant must be a positive integer')
    sizes={'genotype':row_bytes*chunk,'beta':4*chunk*traits,
           'tstat':4*chunk*traits,'flags':chunk,'df':4*chunk}
    if reduction=='jagwas':
        sizes.update(beta=4*chunk,tstat=4*chunk,trait_index=4*chunk)
    elif reduction=='device_significant':
        sizes={'genotype':row_bytes*chunk}
    if not return_beta:sizes.pop('beta')
    allocations=[]
    for name,size in sizes.items():
        rounded=1<<(size-1).bit_length()
        allocations.append(dict(name=name,count=depth,requested_bytes=size,
            allocator_bytes=rounded,pages=(rounded+4095)//4096))
    return dict(allocations=allocations,allocation_count=len(sizes)*depth,
        requested_bytes=sum(r['count']*r['requested_bytes'] for r in allocations),
        allocator_bytes=sum(r['count']*r['allocator_bytes'] for r in allocations),
        allocation_pages=sum(r['count']*r['pages'] for r in allocations),
        scope='Fresh PyTorch power-of-two host allocation requests; '+('unreduced' if reduction is None else reduction)+' host-native scan with explicit transfer row bytes. Driver backing, cache history and allocator service are not bounded.')
