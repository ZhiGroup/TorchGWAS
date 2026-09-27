"""The blocking uint8 status copy and NumPy QC before device selection."""
import math
from .first_principles import positive
from .mechanistic_plan import _integer


def native_status_service(rows, prices, *, cpu_fraction, dram_bytes_per_second):
    """Independent API partitions and fixed/unit CPU prices, never scan time.

    The successful scan has no malformed status. Missing and invariant rows
    still execute both comparisons and sums. Copy API CPU excludes runtime
    waiting and returned-array destruction; those have separate lifetimes.
    """
    _integer('status rows', rows)
    q=positive('CPU fraction',cpu_fraction)
    if q>1: raise ValueError('CPU fraction must not exceed one')
    dram=positive('DRAM capacity',dram_bytes_per_second)
    if not isinstance(prices,dict) or prices.get('dtype')!='uint8':
        raise ValueError('Independent uint8 status-copy and QC services required')
    copy=prices.get('copy',{})
    if set(copy)!={'before_cpu_seconds','after_cpu_seconds'}:
        raise ValueError('Status copy CPU must be partitioned around the barrier')
    before=positive('status copy before',copy['before_cpu_seconds'],True)
    after=positive('status copy after',copy['after_cpu_seconds'],True)
    transfer=prices.get('transfer',{})
    if set(transfer)!={'latency_seconds','bytes_per_second'}:
        raise ValueError('Explicit pageable status transfer latency and capacity required')
    latency=positive('status latency',transfer['latency_seconds'],True)
    bandwidth=positive('status bandwidth',transfer['bytes_per_second'])
    serial=prices.get('host_serial_fraction');wait=prices.get('wait_cpu_fraction')
    for name,value in [('host serialization',serial),('status wait CPU',wait)]:
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('Explicit bounded '+name+' scenario required')
    qc=prices.get('qc',{})
    if set(qc)!={'malformed_empty','count_status'}:
        raise ValueError('Independent empty malformed check and status count prices required')
    phases=[dict(seconds=after/q,host_serial_fraction=serial)]
    cpu=before+after
    # uint8 input read, bool comparison write, then predicate read. No
    # malformed coordinates are allocated on the successful path.
    for name,calls in [('malformed_empty',1),('count_status',2)]:
        price=qc[name]
        if set(price)!={'call_cpu_seconds','row_cpu_seconds'}:
            raise ValueError('QC prices require independent fixed and per-row CPU work')
        fixed=positive(name+' call',price['call_cpu_seconds'],True)
        bulk=positive(name+' row',price['row_cpu_seconds'],True)
        work=calls*(fixed+rows*bulk);traffic=3*calls*rows
        seconds=max(work/q,traffic/dram);cpu+=work
        phases.append(dict(seconds=seconds,host_serial_fraction=serial,
            resources=dict(cpu=work/seconds if seconds else 0.,dram=traffic/seconds if seconds else 0.)))
    return dict(status_bytes=rows,status_submit_seconds=before/q,
        status_submit_resources=dict(cpu=q,host_serial=q*serial),
        d2h_seconds=latency+rows/bandwidth,finish_seconds=sum(p['seconds'] for p in phases),
        finish_operations=phases,result_wait_resources=dict(cpu=q*wait),cpu_seconds=cpu,
        qc_work=dict(malformed_checks=1,status_counts=2,predicate_rows=3*rows,logical_dram_bytes=9*rows),
        unpriced_terms=['status scalar/counter bookkeeping and retained NumPy status-array release',
            'pageable status transfer staging/DRAM traffic beyond supplied independent transfer service'],
        scope='Successful source status path before device selection; no dense result copy or finish future.')
