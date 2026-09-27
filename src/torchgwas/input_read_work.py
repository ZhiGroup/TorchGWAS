"""CPU/copy demand inside buffered input reads, separately from storage service."""
from .first_principles import positive


def buffered_read_service(payload_bytes, calls, prices, *, storage_bytes_per_second,
                          dram_bytes_per_second, cpu_fraction):
    """Fluid overlap scenario for complete preadv requests into resident pages.

    Byte service and fixed-call prices must come from independent generic reads.
    CPU work is charged during the read node; it is not added again to decode.
    Logical kernel-to-user copying is two memory accesses per payload byte.
    """
    amount=positive('payload_bytes',payload_bytes,True)
    if isinstance(calls,bool) or not isinstance(calls,int) or calls<0:
        raise ValueError('Read call count must be a nonnegative integer')
    if amount and not calls:raise ValueError('Nonempty read needs a syscall')
    if isinstance(cpu_fraction,bool) or not 0<cpu_fraction<=1:
        raise ValueError('Positive CPU fraction at most one required')
    storage=positive('storage_bytes_per_second',storage_bytes_per_second)
    dram=positive('dram_bytes_per_second',dram_bytes_per_second)
    byte_cpu=positive('cpu_seconds_per_byte',prices['cpu_seconds_per_byte'],True)
    fixed_cpu=positive('cpu_seconds_per_call',prices['cpu_seconds_per_call'],True)
    cpu=amount*byte_cpu+calls*fixed_cpu
    seconds=max(amount/storage,2*amount/dram,cpu/cpu_fraction)
    resources={name:value/seconds if seconds else 0. for name,value in [('input',amount),('dram',2*amount),('cpu',cpu)]}
    return dict(seconds=seconds,cpu_seconds=cpu,resources=resources,logical_copy_bytes=2*amount,calls=calls,
                scope='Supplied complete-read call-count scenario with fluid IO/CPU/copy overlap. Destination allocation/faults, short-read retries and detailed page-cache fill traffic remain unpriced.')


def direct_read_capacity(probe, *, workers, cpu_affinity, filesystem):
    """Load independently measured aggregate storage service, not buffered CPU rate.

    O_DIRECT observations must transfer the declared bytes in process I/O
    accounting while leaving file-page residency zero. This is a measured
    service in the recorded context, not a physical maximum or universal bound.
    Buffered read CPU/copy work is priced separately by buffered_read_service.
    """
    import math
    import statistics
    if isinstance(workers,bool) or not isinstance(workers,int) or workers not in (1,4):
        raise ValueError('Unmeasured direct-read worker context')
    if probe.get('affinity')!=list(cpu_affinity):
        raise ValueError('Direct-read CPU affinity mismatch')
    inputs=probe['input']
    if inputs.get('filesystem')!=filesystem:
        raise ValueError('Direct-read filesystem mismatch')
    block=inputs.get('block_bytes');count=inputs.get('blocks');amount=inputs.get('bytes')
    if (block,count,amount)!=(16*1024*1024,4,64*1024*1024):
        raise ValueError('Unverified fixed direct-read extent')
    rows=[row for row in probe['rows'] if row.get('mode')=='direct' and row.get('workers')==workers]
    if len(rows)<3 or len({row['repeat'] for row in rows})!=len(rows):
        raise ValueError('Insufficient unique direct-read repeats')
    rates=[]
    for row in rows:
        if (row.get('cache')!='bypass' or row.get('bytes')!=amount or row.get('resident_at_launch')!=0
                or row.get('resident_after')!=0 or row.get('process_io_delta',{}).get('read_bytes')!=amount):
            raise ValueError('Direct-read bypass or transferred-byte evidence failed')
        seconds=positive('direct read elapsed',row['seconds'])
        calls=row['calls']
        if len(calls)!=count or sorted(call['part'] for call in calls)!=list(range(count)):
            raise ValueError('Direct-read request coverage mismatch')
        for call in calls:
            if (call['bytes']!=block or call['offset']!=call['part']*block
                    or call['worker']!=call['part']%workers):
                raise ValueError('Direct-read request geometry mismatch')
            start=positive('read start',call['start'],True)
            end=positive('read end',call['end'])
            if end<=start or not math.isclose(call['wall_seconds'],end-start,rel_tol=1e-8,abs_tol=1e-10):
                raise ValueError('Direct-read interval mismatch')
        if not math.isclose(seconds,max(call['end'] for call in calls),rel_tol=1e-12):
            raise ValueError('Direct-read timing boundary mismatch')
        rates.append(amount/seconds)
    return dict(bytes_per_second=statistics.median(rates),workers=workers,
        repeat_bytes_per_second=rates,filesystem=filesystem,
        service_kind='independent_direct_read_aggregate',
        unpriced_terms=['Transfer of generic direct-read service to buffered files, fragmentation and device/controller-cache state'],
        scope='Fixed16MiB requests,64MiB total, page-aligned O_DIRECT destinations and verified process read bytes. Aggregate wall service at the requested reader context. Buffered CPU work is separate; no GWAS input, elapsed residual or fitted multiplier.')