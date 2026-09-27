"""Source work accounting and resource service, without workload timing fits.

This module deliberately distinguishes necessary resource service from an
end-to-end runtime prediction. A roofline floor is not a predicted elapsed
runtime, and two floors cannot establish a crossover. Unknown work remains
visible. Units and provenance are part of the interface.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
from typing import Mapping


def positive(name, value, allow_zero=False):
    value = float(value)
    if not math.isfinite(value) or value < 0 or (not allow_zero and value == 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return value


@dataclass(frozen=True)
class Cohort:
    samples: int
    markers: int
    covariates: int = 8
    traits: int = 1
    missing_rate: float = .001
    stored_genotype_bytes: int = 0
    metadata_bytes: int = 0
    output_bytes: int | None = None

    def __post_init__(self):
        for key in ('samples', 'markers', 'traits'):
            value = getattr(self, key)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f'{key} must be a positive integer')
        for key in ('covariates', 'stored_genotype_bytes', 'metadata_bytes'):
            value = getattr(self, key)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f'{key} must be a nonnegative integer')
        positive('missing_rate', self.missing_rate, True)
        if self.missing_rate > 1:
            raise ValueError('missing_rate must be at most one')
        if self.output_bytes is not None:
            positive('output_bytes', self.output_bytes, True)


@dataclass(frozen=True)
class Term:
    name: str
    amount: float
    unit: str
    source: str
    interpretation: str

    def __post_init__(self):
        positive(self.name, self.amount, True)
        if not self.source or not self.interpretation:
            raise ValueError('Every work term needs its derivation and interpretation')


@dataclass(frozen=True)
class Availability:
    """Available capacity is resource-specific, never inferred from loadavg.

    CPU fractions describe scheduling only. They do not describe CPU frequency,
    SMT competition or cache/bandwidth losses. Those alter primitive capacity.
    Bandwidths below must already represent bandwidth available to this job;
    multiplying them by CPU scheduling fraction would double-count contention.
    """
    serial_cpu_fraction: float
    worker_cpu_fraction: float
    cpu_capacity_cores: float
    memory_bytes_per_second: float
    storage_read_bytes_per_second: float
    storage_write_bytes_per_second: float
    gpu_fraction: float
    pcie_h2d_bytes_per_second: float
    pcie_d2h_bytes_per_second: float
    input_ram_hit_fraction: float = 0.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            positive(name, value, True)
        for name in ('serial_cpu_fraction', 'worker_cpu_fraction',
                     'gpu_fraction', 'input_ram_hit_fraction'):
            if getattr(self, name) > 1:
                raise ValueError(f'{name} must lie in [0,1]')

    def cpu_cores(self, workers):
        if not isinstance(workers, int) or workers < 1:
            raise ValueError('workers must be a positive integer')
        share = self.serial_cpu_fraction if workers == 1 else self.worker_cpu_fraction
        return min(workers * share, self.cpu_capacity_cores)


def service(work, available_rate):
    """Conservation of work: zero capacity cannot perform positive work."""
    positive('work', work, True)
    positive('available_rate', available_rate, True)
    return 0.0 if work == 0 else (work / available_rate if available_rate else math.inf)


@dataclass(frozen=True)
class CapacityEpoch:
    duration_seconds: float
    rate_per_second: float

    def __post_init__(self):
        positive('duration_seconds', self.duration_seconds)
        positive('rate_per_second', self.rate_per_second, True)


def time_varying_service(work, epochs, start=0.0):
    """Finish time from integral_start^finish capacity(t) dt >= work.

    Exhausting a finite resource forecast returns None, never an extrapolated
    stationary load. A zero-capacity epoch is a pause, not a divide-by-zero or
    silent replacement by an idle-server value.
    """
    positive('work', work, True)
    positive('start', start, True)
    if work == 0:
        return start
    now = 0.0
    left = float(work)
    for epoch in epochs:
        end = now + epoch.duration_seconds
        begin = max(start, now)
        if end > begin:
            capacity = (end - begin) * epoch.rate_per_second
            if epoch.rate_per_second and left <= capacity:
                return begin + left / epoch.rate_per_second
            left -= capacity
        now = end
    return None


def pgen_structure(samples, forms: Mapping[int, int], difflist_entries=None,
                   varint_bytes=None, ld_reference_replays=0):
    """Data-dependent decode work, not a per-form elapsed-time lookup.

    Counts are for the basic biallelic hardcall record. The low record-type bits
    alone do not price multiallelic/dosage/phase payloads. Entries and varint
    lengths must come from a census, not a constant inferred from compressed
    file size. LD references need not be the immediately preceding variant.
    """
    if not isinstance(samples, int) or samples <= 0:
        raise ValueError('samples must be positive integer')
    if not isinstance(ld_reference_replays, int) or ld_reference_replays < 0:
        raise ValueError('ld_reference_replays must be nonnegative integer')
    if any(k not in {0,1,2,3,4,6,7} for k in forms):
        raise ValueError('Unsupported PGEN record form or additional payload bits')
    if any(not isinstance(v, int) or v < 0 for v in forms.values()):
        raise ValueError('Record counts must be nonnegative integers')
    m = sum(forms.values())
    packed = (samples + 3) // 4
    diff_records = m - forms.get(0, 0)
    unknown = []
    if diff_records and difflist_entries is None:
        unknown.append('exact difflist entry count')
    if diff_records and varint_bytes is None:
        unknown.append('exact variable-integer byte count')
    for key, val in [('difflist_entries', difflist_entries), ('varint_bytes', varint_bytes)]:
        if val is not None: positive(key, val, True)
    return dict(record_counts=dict(forms), packed_bytes_per_variant=packed,
        plain_copy_bytes=packed * forms.get(0,0),
        ld_base_copy_bytes=packed * (forms.get(2,0)+forms.get(3,0)),
        ld_invert_calls=samples * forms.get(3,0),
        one_bit_calls=samples * forms.get(1,0),
        constant_background_calls=samples * sum(forms.get(f,0) for f in [4,6,7]),
        difflist_entries=difflist_entries, varint_bytes=varint_bytes,
        ld_reference_replays=ld_reference_replays,
        unknown=unknown+(['record census and lengths of replayed LD bases'] if ld_reference_replays else []))


def plink_missing_branches(n, m, missing_rate, restart_segments=1):
    """Expected branches under explicitly iid missingness, not a timing fit.

    A complete variant after a missing variant (or at each worker/block start)
    opens the fast dense path. Consecutive complete variants can reuse the
    precomputed covariate products. Exact observed branch counts should replace
    these expectations when genotype missingness clusters across variants.
    """
    if not isinstance(restart_segments, int) or not 1 <= restart_segments <= m:
        raise ValueError('restart_segments must be in [1,M]')
    if not 0 <= missing_rate <= 1 or n <= 0 or m <= 0:
        raise ValueError('Invalid dimensions or missingness')
    p = math.exp(n*math.log1p(-missing_rate)) if missing_rate < 1 else 0.0
    gram = m * (1-p)
    opener = restart_segments*p + (m-restart_segments)*p*(1-p)
    sparse = (m-restart_segments)*p*p
    return dict(gram=gram, opener=opener, sparse=sparse,
                complete_probability=p, assumption='iid across subjects and markers')


def work_ledger(method, cohort: Cohort, *, branch_counts=None, sparse_carriers=None):
    n,m,c,k = cohort.samples,cohort.markers,cohort.covariates,cohort.traits
    terms=[]
    def add(name, amount, unit, source, interpretation='source-level count'):
        terms.append(Term(name,amount,unit,source,interpretation))
    common_unknown=['process/import/dynamic-loader and library initialization service',
        'metadata parse, sample matching and covariate preprocessing service',
        'PGEN record census, replay dependencies and decoder instruction service',
        'allocator, synchronization, controller and teardown service',
        'physical cache-line traffic and instruction scheduling of compiled kernels']
    add('genotype_input',cohort.stored_genotype_bytes,'bytes','file length; one hardcall pass')
    add('metadata_input',cohort.metadata_bytes,'bytes','actual input file lengths')
    if method=='fastGWA':
        if k!=1: raise ValueError('This ledger covers one GCTA process/phenotype; no new K model')
        p=c+1
        src='GCTA1.95.3 FastFAM.cpp:1363,2433,2436-2437'
        add('projection_flops',4*n*p*m,'flop',src,'leading multiply/add count of two GEMVs')
        add('association_dot_flops',4*n*m,'flop',src,'leading multiply/add count of two dot products')
        add('projection_matrix_reads',16*n*p*m,'logical_bytes',src,'X and H are each read once per marker; not necessarily DRAM traffic')
        add('projection_vector_traffic',(56*n+32*p)*m,'logical_bytes',src,
            'Eigen 3.4 alias-safe temporary: copy y to temporary and back (32N), H*y input (8N), temporary GEMV read/update (16N), and Hy initialize/read/update/read (32P); tiled implementation can add L1 references')
        add('projection_temporary_allocations',2*m,'calls','Eigen 3.4 ProductEvaluators.h:184-215 and nested product evaluation')
        add('projection_temporary_allocated_bytes',8*(n+p)*m,'bytes_allocated','N-vector assignment temporary plus P-vector nested H*y')
        add('expand_count_dot_traffic',m*(2*((n+3)//4)+40*n),'logical_bytes',
            'Geno.cpp:977-983; PgenReader.cpp:413,418,558; FastFAM.cpp:2436-2437')
        add('genotype_vector_allocations',m,'calls','Geno.cpp:977')
        add('genotype_vector_allocated_bytes',8*n*m,'bytes_allocated','Geno.cpp:977')
        add('chisquare_tail',m,'calls','FastFAM.cpp:2448')
        add('numeric_output_fields',5*m,'fields','FastFAM.cpp:2796-2801')
        add('string_concatenations',8*m,'calls','Marker.cpp:1353-1358 and FastFAM output')
        topology=dict(analysis_workers=1,reader_workers=1,block_markers=1024,reader_slots=3,
            build='official 1.95.3 Linux AppImage a82f1b737ab2ce9087202ff405233e7105689b8b5c62f3a882efdfe5c920f29b',
            note='Build-specific serial analysis loop, MKL_NUM_THREADS=1; no assumed OpenMP speedup')
        common_unknown += ['numeric formatting and chi-square library instruction service']
    elif method=='PLINK2':
        if k!=1: raise ValueError('This ledger deliberately covers the validated K=1 branch only')
        p=c+2
        if branch_counts is None:
            branch_counts=plink_missing_branches(n,m,cohort.missing_rate)
        gram,opener,sparse=(float(branch_counts[t]) for t in ['gram','opener','sparse'])
        if any(not math.isfinite(v) or v<0 for v in [gram,opener,sparse]) or not math.isclose(gram+opener+sparse,m,abs_tol=1e-7):
            raise ValueError('PLINK branch counts must be nonnegative and sum to M')
        # N_observed differs conditionally: all missing calls belong to Gram
        # variants. An exact census can pass gram_observed_sample_sum.
        observed_sum=branch_counts.get('gram_observed_sample_sum',n*(gram-m*cohort.missing_rate))
        if not 0<=observed_sum<=n*gram+1e-7: raise ValueError('Invalid observed sample sum')
        src='plink-ng ca0f464: plink2_glm_linear.cc Gram branch; plink2_matrix.cc dsyrk'
        add('gram_flops',observed_sum*p*(p+1),'flop',src,'One triangular SYRK, not a full 2NP^2 product')
        add('xty_flops',2*observed_sum*p,'flop','plink2_matrix.h LinearRegressionInv')
        add('opener_crossproducts',2*n*(c+2)*opener,'flop','plink2_glm_linear.cc fast dense path')
        add('gather_elements',(c+1)*(observed_sum+n*opener),'elements','BitIter1 covariate/phenotype gather')
        add('small_matrix_inversions',gram,'calls','LinearRegressionInv')
        add('rank1_inverse_updates',opener+sparse,'calls','InvertRank1Symm')
        add('vif_checks',gram,'calls','Gram-path VIF check')
        if sparse_carriers is None and sparse>0:
            common_unknown.append('number of minor-allele carriers on sparse-eligible records')
        elif sparse_carriers is not None:
            positive('sparse_carriers',sparse_carriers,True)
            add('sparse_crossproducts',2*sparse_carriers*(c+1),'flop','sparse covariate/phenotype carrier loop')
        topology=dict(block_markers=65536,main='read next block; join; launch compute; format previous block',
            branch_counts=branch_counts,compute_workers='explicit requested/build-dependent count')
        common_unknown += ['small-matrix factorization/update, VIF, p-value and dtoa instruction service',
            'worker/block restart census when using expected missingness branches']
    elif method=='torchGWAS':
        src='native_scan.py dosage_cuda_iterator and linear.py _dosage_statistics'
        add('gemm_flops',2*n*m*(k+c+1),'flop',src,'Full FP32 design has phenotype, intercept and C covariates')
        add('h2d',n*m,'bytes',src,'int8 dosage transfer; one GPU, no multi-GPU speedup')
        add('d2h',m*(8*k+5),'bytes',src,'FP32 beta and t per cell plus uint8 status and FP32 df per marker')
        add('binary_output',8*m*k,'bytes','api.py _write_linear_binary_streaming','beta/t binary matrix payload; headers and sidecars additional')
        from .binary_output_work import binary_output_work
        binary=binary_output_work(m,k,2048)
        add('binary_staging_zero_bytes',binary['zero_initialization_bytes'],'logical_bytes','sumstats.py _BlockStream.__init__','bytearray initialization of queue_depth+1 staging blocks for each array')
        add('binary_staging_copy_traffic',binary['staging_logical_memory_bytes'],'logical_bytes','sumstats.py _BlockStream.append','source reads and staging writes; no borrowing for sub-block chunks')
        add('binary_payload_queued_at_close',binary['payload_queued_only_at_close_bytes'],'bytes','sumstats.py _BlockStream.close','these bytes cannot overlap earlier scan execution; partial buffers are emitted at close')
        add('binary_fsync_calls',binary['fsync_calls'],'calls','api.py default sumstats_fsync=True; sumstats.py close','native beta/t payload commit, before non-fsynced sidecars')
        # Eager Torch materializes zeros/full_like and intermediate arrays.
        # These are logical accesses, NOT HBM-byte claims or elapsed times.
        # CUDA bool sum materializes int64: 1N read + 8N write + 8N reduce read.
        # Verified by a duration-free kernel census on all three GPU families.
        # Row reductions count their input reads; small outputs listed separately.
        traffic={'int8_missing_compare':2,'int8_to_fp32':5,'conversion_where':9,
            'isnan':5,'observed_invert':2,'present_reduce':17,
            'zero_missing_fill_and_where':17,'minimum_fill_where_reduce':21,
            'maximum_fill_where_reduce':21,'mean_reduce':4,
            'center_subtract_fill_where':25,'square_and_reduce':12}
        for label,b in traffic.items():add(label,n*m*b,'logical_bytes',src,'Eager, non-fused FP32 path; reuse/caching not priced as HBM')
        add('gemm_unique_input',4*n*m+4*n*(k+c+1),'logical_bytes',src,'Unique matrix input volume, not tile-reload traffic')
        add('result_host_copies',2*m*(8*k+5),'logical_bytes','native_scan.py finish','Read and write for numpy.copy')
        topology=dict(gpus=1,reader_workers=4,chunk_markers=2048,host_slots=4,device_slots=4,
            backend='int8 CPU decode; eager Torch FP32 statistics; binary beta/t output',
            binary_writer=dict(block_bytes=binary['block_bytes'],queue_depth=binary['queue_depth'],fsync_calls=binary['fsync_calls'],payload_queued_at_close=binary['payload_queued_only_at_close_bytes']))
        common_unknown += ['elementwise/reduction instruction counts and small result-array traffic',
            'actual GEMM tile selection, finite-wave utilization and kernel launch service',
            'bounded producer/copy/compute/result queue synchronization service',
            'binary headers, marker/trait sidecars and durability boundary']
    else:
        raise ValueError('method must be torchGWAS, PLINK2 or fastGWA')
    if method!='torchGWAS':
        if cohort.output_bytes is None: common_unknown.append('native text output bytes from data value/identifier statistics')
        else: add('native_output',cohort.output_bytes,'bytes','explicit output value/identifier statistics')
    return dict(method=method,cohort=asdict(cohort),terms=[asdict(t) for t in terms],
                topology=topology,unpriced_mechanisms=common_unknown,
                prediction_seconds=None,status='source work ledger; incomplete elapsed-time model')


def necessary_resource_service(ledger, availability: Availability, *, cpu_flops_per_core_second,
                               gpu_flops_per_second, cpu_workers=1):
    """Partial resource floor using declared capacity ceilings, not a prediction.

    Logical memory accesses are intentionally not divided by DRAM bandwidth:
    doing so would silently assume a cache-miss rate. Omitted work can only
    increase elapsed time. Native buffered output bytes are not priced as
    durable storage writes without a matching fsync timing boundary.
    """
    positive('cpu_flops_per_core_second',cpu_flops_per_core_second)
    positive('gpu_flops_per_second',gpu_flops_per_second)
    terms=ledger['terms'];method=ledger['method']
    total=lambda name:sum(t['amount'] for t in terms if t['name']==name)
    flops=sum(t['amount'] for t in terms if t['unit']=='flop')
    rates={}
    rates['cpu_compute']=service(0 if method=='torchGWAS' else flops,
        cpu_flops_per_core_second*availability.cpu_cores(cpu_workers))
    rates['gpu_compute']=service(flops if method=='torchGWAS' else 0,
        gpu_flops_per_second*availability.gpu_fraction)
    read_bytes=total('genotype_input')+total('metadata_input')
    rates['input_storage']=service(read_bytes*(1-availability.input_ram_hit_fraction),availability.storage_read_bytes_per_second)
    rates['h2d']=service(total('h2d'),availability.pcie_h2d_bytes_per_second)
    rates['d2h']=service(total('d2h'),availability.pcie_d2h_bytes_per_second)
    return dict(partial_resource_floor_seconds=max(rates.values()),resource_seconds=rates,
        prediction_seconds=None,unpriced_mechanisms=ledger['unpriced_mechanisms'],
        interpretation='Necessary service under declared resource ceilings; not a runtime estimate or validated crossover')


def compare_runtime_intervals(torch_interval, competitor_interval):
    """Only full elapsed-time intervals can establish a guaranteed sign."""
    if torch_interval is None or competitor_interval is None:
        return 'undetermined: full runtime intervals unavailable'
    for interval in (torch_interval,competitor_interval):
        if len(interval)!=2 or any(not math.isfinite(x) or x<0 for x in interval) or interval[0]>interval[1]:
            raise ValueError('Need ordered finite nonnegative full-runtime intervals')
    if torch_interval[1]<competitor_interval[0]: return 'torchGWAS faster throughout supplied intervals'
    if torch_interval[0]>competitor_interval[1]: return 'torchGWAS slower throughout supplied intervals'
    return 'undetermined: intervals overlap'


@dataclass(frozen=True)
class AffineServiceEnvelope:
    """Conditional bounds from a closed source/primitive service model.

    At fixed N and resource state, startup and steady per-marker costs may be
    represented as A + B*M only after pipeline fill/drain and finite blocks have
    been bounded. This class does not fit A or B to observed association times.
    """
    startup_seconds: tuple[float, float]
    seconds_per_marker: tuple[float, float]
    provenance: str
    unpriced_mechanisms: tuple[str, ...] = ()

    def __post_init__(self):
        if not self.provenance:
            raise ValueError('A source/primitive derivation is required')
        for interval in (self.startup_seconds, self.seconds_per_marker):
            if len(interval)!=2 or any(not math.isfinite(x) or x<0 for x in interval) or interval[0]>interval[1]:
                raise ValueError('Expected ordered finite nonnegative conditional bounds')


def conservative_marker_boundary(torch_model, competitor_model):
    """Solve upper(torch) < lower(competitor), without a fitted hyperbola.

    A finite boundary is conditional on every service envelope, including load.
    If compute or I/O availability can vanish, no finite universal bound follows.
    A piecewise resource/cache/block model must solve each valid piece separately.
    """
    missing=tuple(torch_model.unpriced_mechanisms)+tuple(competitor_model.unpriced_mechanisms)
    if missing:
        return dict(status='incomplete',marker_boundary=None,unpriced_mechanisms=missing)
    a=torch_model.startup_seconds[1]-competitor_model.startup_seconds[0]
    advantage=competitor_model.seconds_per_marker[0]-torch_model.seconds_per_marker[1]
    if advantage<=0:
        return dict(status='no guaranteed eventual torchGWAS advantage from these envelopes',
                    marker_boundary=None,steady_advantage_seconds_per_marker=advantage)
    return dict(status='conditional on supplied complete service envelopes',
                marker_boundary=max(0.0,a/advantage),
                inequality='M strictly greater than marker_boundary',
                steady_advantage_seconds_per_marker=advantage,
                derivation='max(0, (A_t_upper-A_c_lower)/(B_c_lower-B_t_upper))')
