"""Resource accounting for independently measured source operations.

CPU service means summed on-CPU seconds, including the cache/memory stalls of
its measured working-set regime. Scheduling fraction converts service to wall;
load average is never used. Memory traffic supplies an additional capacity lower
bound, not a second additive charge for stalls already present in CPU service.
A materially different memory regime requires fresh component rates; a generic
bandwidth number cannot uniquely decompose compute and memory stall time.
"""
from dataclasses import asdict,dataclass
import math


def positive(name,value):
    if not math.isfinite(value) or value<=0:raise ValueError(f'{name} must be finite and positive')

@dataclass(frozen=True)
class ResourceState:
    cpu_fraction_serial: float
    cpu_fraction_workers: float
    cpu_capacity_cores: float
    memory_bytes_per_second: float
    storage_read_bytes_per_second: float
    storage_write_bytes_per_second: float
    input_cache_hit_fraction: float
    gpu_service_fraction: float=1.
    cpu_fraction_barrier_worker: float=1.
    def __post_init__(self):
        for name,value in asdict(self).items():
            if name=='input_cache_hit_fraction':
                if not math.isfinite(value) or not 0<=value<=1:raise ValueError('cache hit fraction must be in [0,1]')
            else:positive(name,value)
        for key in ('cpu_fraction_serial','cpu_fraction_workers','gpu_service_fraction','cpu_fraction_barrier_worker'):
            if getattr(self,key)>1:raise ValueError(f'{key} must not exceed 1')
    def cpu_wall(self,service_seconds,threads=1):
        if service_seconds<0 or not math.isfinite(service_seconds):raise ValueError('invalid CPU service')
        positive('threads',threads)
        fraction=self.cpu_fraction_serial if threads==1 else self.cpu_fraction_workers
        return service_seconds/min(threads*fraction,self.cpu_capacity_cores)
    def input_io_wall(self,bytes_):
        if bytes_<0:raise ValueError('negative bytes')
        return bytes_*(1-self.input_cache_hit_fraction)/self.storage_read_bytes_per_second

@dataclass(frozen=True)
class StageDemand:
    name: str
    cpu_seconds: float=0.
    threads: int=1
    memory_bytes: float=0.
    gpu_seconds: float=0.
    read_bytes: float=0.
    write_bytes: float=0.
    serial_wall_seconds: float=0.
    def __post_init__(self):
        if self.threads<1:raise ValueError('positive thread count required')
        for key,value in asdict(self).items():
            if key not in ('name','threads') and (not math.isfinite(value) or value<0):raise ValueError(f'invalid {key}')
    def seconds(self,state):
        # Competing limits within an overlapped stage, plus truly serial work.
        return self.serial_wall_seconds+max(state.cpu_wall(self.cpu_seconds,self.threads),self.memory_bytes/state.memory_bytes_per_second,self.gpu_seconds/state.gpu_service_fraction,state.input_io_wall(self.read_bytes),self.write_bytes/state.storage_write_bytes_per_second)


def shared_resource_floor(stages,state):
    """Whole-work capacity bounds; do not pretend stages have separate DRAM."""
    return max(sum(s.cpu_seconds for s in stages)/state.cpu_capacity_cores,
               sum(s.memory_bytes for s in stages)/state.memory_bytes_per_second,
               sum(s.gpu_seconds for s in stages)/state.gpu_service_fraction,
               state.input_io_wall(sum(s.read_bytes for s in stages)),
               sum(s.write_bytes for s in stages)/state.storage_write_bytes_per_second)


def interpolate_unit_service(samples,request_n,unit_count):
    """Interpolate CPU/GPU seconds per source work unit, not GWAS runtime.

    Samples are (N, independently measured component seconds). Work counts
    must come from that component's source. Interpolation is only within the
    measured range, in log working-set size. This empirical primitive-rate
    interpolation must itself pass withheld-size component checks.
    """
    samples=sorted(samples)
    if not samples or request_n<samples[0][0] or request_n>samples[-1][0]:raise ValueError('N outside independently profiled shape domain')
    positive('N',request_n)
    for n,seconds in samples:
        positive('component seconds',seconds);positive('work units',unit_count(n))
        if n==request_n:return seconds
    for (lo,a),(hi,b) in zip(samples,samples[1:]):
        if lo<request_n<hi:
            weight=math.log(request_n/lo)/math.log(hi/lo)
            return ((1-weight)*a/unit_count(lo)+weight*b/unit_count(hi))*unit_count(request_n)
    raise ValueError('unbracketed N')


def bounded_gpu_pipeline(variants,chunk_variants,workers,depth,fill_seconds,
                         controller_seconds,host_seconds,copy_seconds,
                         gpu_seconds,result_seconds,reader_init_seconds=0.):
    """Deterministic service-time schedule for the source's bounded rings.

    Ordered producer emission; one CPU fill per worker; host slots released
    after H2D; device slots gated by prior compute; result futures limit host
    lookahead. All supplied times are per full chunk except reader init.
    This predicts the stated constant-service scenario, not stochastic tails.
    """
    import heapq
    if variants<1 or chunk_variants<1 or workers<1 or depth<2:raise ValueError('invalid pipeline geometry')
    services=(fill_seconds,controller_seconds,host_seconds,copy_seconds,gpu_seconds,result_seconds,reader_init_seconds)
    if any(not math.isfinite(v) or v<0 for v in services):raise ValueError('invalid service time')
    workers=min(workers,depth);chunks=math.ceil(variants/chunk_variants)
    worker_free=[(0.,i) for i in range(workers)];heapq.heapify(worker_free)
    initialized=set();slot_free=[0.]*depth;fills=[];ends=[];results=[]
    controller=host=copy=compute=result=0.;submitted=0
    for j in range(chunks):
        while submitted<min(chunks,j+workers):
            index=submitted;fraction=min(chunk_variants,variants-index*chunk_variants)/chunk_variants
            worker_ready,worker=heapq.heappop(worker_free)
            controller=max(controller,slot_free[index%depth])+controller_seconds
            start=max(controller,worker_ready)
            end=start+fill_seconds*fraction+(reader_init_seconds if worker not in initialized else 0.)
            initialized.add(worker);fills.append(end);heapq.heappush(worker_free,(end,worker));submitted+=1
        controller=max(controller,fills[j])
        fraction=min(chunk_variants,variants-j*chunk_variants)/chunk_variants
        host=max(host,controller,results[j-depth] if j>=depth else 0.)+host_seconds
        copy=max(copy,host,ends[j-depth] if j>=depth else 0.)+copy_seconds*fraction
        slot_free[j%depth]=copy
        compute=max(compute,copy)+gpu_seconds*fraction;ends.append(compute)
        result=max(result,compute)+result_seconds*fraction;results.append(result)
    return dict(seconds=result,result_ready_seconds=results,chunks=chunks,
                scope='source ring constraints with deterministic independently measured component services')


def resource_state_with_cpu_probes(state, probes, worker_threads=4, affinity_threads=None):
    """Apply independent scheduling probes; never infer load from GWAS runtime.

    A full-affinity probe is required to update shared CPU capacity. Without it,
    retain the separately supplied capacity, rather than inventing one from a
    smaller team. Mean worker share is the explicit migrating-worker scenario.
    This does not capture frequency, SMT, cache contention or future load changes.
    """
    from dataclasses import replace
    by_threads = {}
    for probe in probes:
        threads = probe['threads']
        if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
            raise ValueError('CPU probe thread counts must be positive integers')
        if threads in by_threads:
            raise ValueError('Duplicate CPU probe thread count')
        fraction = probe['scheduled_fraction']
        positive('CPU scheduling fraction', fraction)
        by_threads[threads] = min(1., fraction)
    if 1 not in by_threads or worker_threads not in by_threads:
        raise ValueError('Serial and worker-team CPU probes are required')
    updates = dict(cpu_fraction_serial=by_threads[1],
                   cpu_fraction_workers=by_threads[worker_threads],
                   cpu_fraction_barrier_worker=by_threads[worker_threads])
    if affinity_threads is not None:
        if affinity_threads not in by_threads:
            raise ValueError('Full-affinity CPU probe required for shared capacity')
        updates['cpu_capacity_cores'] = affinity_threads * by_threads[affinity_threads]
    return replace(state, **updates)
