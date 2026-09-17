"""Input-derived work counts; no elapsed-time fitting or machine constants."""
from dataclasses import dataclass
import math
from .pgen_reader import ld_safe_start

@dataclass(frozen=True)
class CpuOperatingConditions:
    """CPU capacity available to this job in its affinity/cgroup domain.

    This is an explicit scenario or an independent measurement, in cores,
    not the host load average (which also counts I/O-blocked tasks).
    """
    available_cores: float
    source: str

    def capacity(self, runnable_workers):
        if not math.isfinite(self.available_cores) or self.available_cores<=0:
            raise ValueError('available_cores must be finite and positive')
        if not self.source:
            raise ValueError('operating-condition source is required')
        return min(float(runnable_workers),self.available_cores)


def native_pgen_work(vrtypes, record_lengths, samples, ranges, difference_entries=None):
    """Price the records actually visited by NativePgenReader._decode_into.

    Every requested range independently rebuilds its LD base, replaying the
    prefix from ld_safe_start. This is the current native implementation,
    not a claim that all PGEN readers replay the same records. Exact sparse
    patch counts require supplied per-record difference-list lengths.
    Headers and metadata parsing are outside these genotype-record counts.
    """
    types=vrtypes;lengths=record_lengths
    if samples<=0 or len(types)!=len(lengths) or any(v<0 or int(v)!=v for v in lengths):
        raise ValueError('invalid sample count or record lengths')
    if any(v not in (0,1,2,3,4,6,7) for v in types):
        raise ValueError('work counting currently supports biallelic hardcall records only')
    if len(types) and types[0] in (2,3):raise ValueError('LD record lacks a preceding base')
    if difference_entries is not None:
        difference_entries=[int(v) for v in difference_entries]
        if len(difference_entries)!=len(types) or any(v<0 or v>samples for v in difference_entries):
            raise ValueError('invalid per-record difference counts')
    counts={v:0 for v in (0,1,2,3,4,6,7)}
    requested=decoded=read_bytes=replay_bytes=prefix_records=read_calls=0
    patches=0 if difference_entries is not None else None
    base_ranges=[]
    for start,end in ranges:
        if not 0<=start<=end<=len(types):raise ValueError('invalid decode range')
        if start==end:continue
        safe=ld_safe_start(types,start)
        if types[safe] in (2,3):raise ValueError('LD record lacks a preceding base')
        read_calls+=1;requested+=end-start;decoded+=end-safe;prefix_records+=start-safe
        read_bytes+=int(sum(lengths[safe:end]));replay_bytes+=int(sum(lengths[safe:start]))
        for i in range(safe,end):
            counts[int(types[i])]+=1
            if patches is not None and types[i]!=0:patches+=difference_entries[i]
        base_ranges.append((safe,start,end))
    packed_width=(samples+3)//4
    return dict(samples=samples,requested_records=requested,decoded_records=decoded,
                replay_records=prefix_records,record_read_bytes=read_bytes,
                replay_read_bytes=replay_bytes,minimum_record_read_calls=read_calls,
                decoded_by_form=counts,packed_materialization_bytes=decoded*packed_width,
                difference_patches=patches,decode_ranges=base_ranges,
                scope='NativePgenReader per-range LD-prefix replay; no timing coefficient')


def form_service_cpu_seconds(work, seconds_per_record_by_form):
    """Sum count * independently supplied form service; no default rates.

    A form average is conditional on sample size and its difference-list
    distribution. Use primitive patch/byte costs for transfer across those
    distributions; a form label alone is not sufficient provenance.
    """
    total=0.0
    for form,count in work['decoded_by_form'].items():
        if not count:continue
        if form not in seconds_per_record_by_form:raise ValueError(f'missing rate for record form {form}')
        rate=seconds_per_record_by_form[form]
        if not math.isfinite(rate) or rate<0:raise ValueError('invalid service rate')
        total+=count*rate
    return total
