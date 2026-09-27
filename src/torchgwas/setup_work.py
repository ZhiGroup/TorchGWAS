"""Source-derived FP32 complete-phenotype preprocessing and design work.

No association timings. Phase services are priced from fixed tiny generic
operations plus explicit additional arithmetic and transfer work. Logical byte
counts are not claims of physical HBM traffic or allocator-reserved peaks.
"""
import hashlib
from pathlib import Path
from .first_principles import positive


def setup_reference_shape(rank):
    if type(rank) is not int or rank<0:raise ValueError('Invalid covariate rank')
    return [max(32,rank+3),1,rank]


def setup_primitive_bank(profile,rank):
    """Exact independent tiny reference; never interpolate across ranks."""
    bank=profile.get('setup_primitives_by_rank',{}).get(str(rank))
    if bank is None:
        bank=profile.get('setup_primitives',{})
    expected=setup_reference_shape(rank)
    phases={'residual_common','residual_block','design_common','design_block'}
    if set(bank)!=phases or any(row.get('reference_shape')!=expected for row in bank.values()):
        raise ValueError('Missing independent setup primitives for covariate rank '+str(rank))
    for row in bank.values():
        positive('cpu_seconds',row['cpu_seconds'],True)
        positive('non_cpu_seconds',row['non_cpu_seconds'],True)
    return bank

def phase_work(phase,n,t,c=8):
    nc=n*c;nt=n*t
    if phase=='residual_common':
        return dict(h2d_bytes=4*nc,d2h_bytes=0,gemm_flops=0,gemm_bytes=0,
                    vector_ops=0,vector_bytes=0,host_copy_bytes=0)
    if phase=='residual_block':
        return dict(h2d_bytes=4*nt,d2h_bytes=4*nt,gemm_flops=4*n*c*t,
                    gemm_bytes=8*(nc+nt+c*t) if c else 0,vector_ops=(9 if c else 8)*nt,
                    vector_bytes=(36 if c else 24)*nt+38*t,host_copy_bytes=0)
    if phase=='design_common':
        return dict(h2d_bytes=4*nc,d2h_bytes=0,gemm_flops=0,gemm_bytes=0,
                    vector_ops=n,vector_bytes=4*n+(16 if c else 8)*n*(c+1),host_copy_bytes=0)
    if phase=='design_block':
        return dict(h2d_bytes=4*nt,d2h_bytes=0,gemm_flops=0,gemm_bytes=0,
                    vector_ops=(2*n-1)*t,vector_bytes=20*nt+4*t,host_copy_bytes=0)
    raise ValueError('Unknown setup phase: '+phase)


def setup_memory(samples, traits, covariates, *, sm_count=None, max_threads_per_sm=None):
    """Live FP32 setup tensor bounds; separate residual and design lifetimes."""
    from .preprocess import _device_trait_block
    from .native_scan import design_column_block
    for name,value,minimum in [('samples',samples,1),('traits',traits,1),('covariates',covariates,0)]:
        if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
            raise ValueError('Invalid '+name)
    n,k,c=samples,traits,covariates
    residual_width=min(k,_device_trait_block(n,4,None))
    design_width=design_column_block(n,k)
    residual=4*n*c+12*n*residual_width+4*c*residual_width+16*residual_width
    design=4*n*(k+c+1)+4*k+4*n+4*n*(c+1)+4*n*c+8*n*design_width
    workspace=None;residual_workspace=None
    if (sm_count is None)!=(max_threads_per_sm is None):
        raise ValueError('Both GPU reduction properties are required')
    if sm_count is not None:
        from .reduction_memory import column_sum_workspace,column_std_workspace
        # A final partial block can choose a different reduction geometry.
        widths={design_width,k%design_width} - {0}
        workspace=max(column_sum_workspace(n,w,sm_count=sm_count,max_threads_per_sm=max_threads_per_sm)['workspace_bytes']+column_sum_workspace(n,w,sm_count=sm_count,max_threads_per_sm=max_threads_per_sm)['semaphore_bytes'] for w in widths)
        design+=workspace
        residual_widths={residual_width,k%residual_width} - {0}
        residual_workspace=max(column_std_workspace(n,w,sm_count=sm_count,max_threads_per_sm=max_threads_per_sm)['workspace_bytes']+column_std_workspace(n,w,sm_count=sm_count,max_threads_per_sm=max_threads_per_sm)['semaphore_bytes'] for w in residual_widths)
        # The sequential mean workspace is smaller than Welford's. Add the
        # maximum to the conservative live tensor sum, not both capacities.
        residual+=residual_workspace
    return dict(design_reduction_workspace_bytes=workspace,residual_reduction_workspace_bytes=residual_workspace,
        unresolved_memory_terms=[] if workspace is not None else ['CUDA design and residual mean/std reduction workspaces require GPU properties'],
        residual_device_live_bytes_upper=residual,
        design_device_live_bytes_upper=design,
        device_live_bytes_upper=max(residual,design))

def setup_work(samples,traits,covariates=8, *, reuse_observed_counts=False, input_contiguous=None, covariate_columns=None):
    from .preprocess import _device_trait_block
    from .native_scan import design_column_block
    if any(type(v) is not int or v<1 for v in (samples,traits)) or type(covariates) is not int or covariates<0 or samples<=covariates+2:
        raise ValueError('Setup census requires positive N,K and valid covariate rank/df')
    columns=covariates if covariate_columns is None else covariate_columns
    if type(columns) is not int or not covariates<=columns<samples-2:
        raise ValueError('Invalid covariate column count')
    if type(reuse_observed_counts) is not bool:
        raise ValueError('reuse_observed_counts must be boolean')
    if input_contiguous is not None and type(input_contiguous) is not bool:
        raise ValueError('input_contiguous must be boolean or None')
    n,k,c=samples,traits,covariates
    residual_width=min(k,_device_trait_block(n,4,None))
    design_width=design_column_block(n,k)
    phases=[dict(phase='residual_common',traits=k,**phase_work('residual_common',n,k,c))]
    for start in range(0,k,residual_width):
        t=min(k-start,residual_width);phases.append(dict(phase='residual_block',traits=t,**phase_work('residual_block',n,t,c)))
    phases.append(dict(phase='design_common',traits=k,**phase_work('design_common',n,k,c)))
    for start in range(0,k,design_width):
        t=min(k-start,design_width);phases.append(dict(phase='design_block',traits=t,**phase_work('design_block',n,t,c)))
    # Copies depend on the actual block layout, not simply on element count.
    # The residual download is returned directly for a single block. An input
    # layout not supplied by the caller is a conservative contiguous-copy bound.
    for p in phases:
        t=p['traits'];input_copy=assembly=0
        if p['phase']=='residual_block':
            input_copy=4*n*t if input_contiguous is not True or (n>1 and t<k) else 0
            assembly=4*n*t if k>residual_width else 0
        elif p['phase']=='design_block':
            input_copy=4*n*t if n>1 and t<k else 0
        p.update(input_contiguous_copy_bytes=input_copy,result_assembly_copy_bytes=assembly,
                 host_copy_bytes=input_copy+assembly)
    paths=[Path(__file__).with_name(name) for name in ['preprocess.py','native_scan.py']]
    result=dict(samples=n,traits=k,covariates=c,residual_block_traits=residual_width,
        design_block_traits=design_width,phases=phases,
        observed_counts_reused=reuse_observed_counts,input_contiguous=input_contiguous,
        host_copy_bytes=sum(p['host_copy_bytes'] for p in phases),
        residual_result_assembly_bytes=sum(p['result_assembly_copy_bytes'] for p in phases),
        observation_cpu_work=dict(phenotype_isnan_cells=0 if reuse_observed_counts else n*k,
            boolean_not_cells=0 if reuse_observed_counts else n*k,
            column_count_reduction_cells=0 if reuse_observed_counts else n*k,
            count_vector_cells=k),
        h2d_bytes=sum(p['h2d_bytes'] for p in phases),d2h_bytes=4*n*k,
        residual_result_host_bytes=4*n*k,
        **setup_memory(n,k,c),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        scope='Actual residualization and design-upload blocking; complete FP32 phenotypes, C8. '
              'Logical transfer/arithmetic work and a conservative live-tensor setup bound. '
              'Excludes allocator reservations, scan temporaries, input arrays and writer memory.')
    if columns!=c:result['covariate_columns']=columns
    if c!=8:result['scope']=result['scope'].replace('C8','explicit covariate rank '+str(c))
    return result


def setup_service(work,profile):
    bank=setup_primitive_bank(profile,work['covariates']);q=positive('cpu_fraction',profile['cpu_fraction'])
    gpu=profile['gpu_resources'];g=positive('gpu_fraction',gpu['gpu_fraction'])
    bw=positive('hbm_bytes_per_second',gpu['hbm_bytes_per_second'])*g
    flops=positive('fp32_flops_per_second',gpu['fp32_flops_per_second'])*g
    h2d=positive('h2d_bytes_per_second',profile['h2d_bytes_per_second'])
    d2h=positive('d2h_bytes_per_second',profile['d2h_bytes_per_second'])
    copy=positive('numpy_copy_bytes',profile['process_units']['numpy_copy_bytes'],True)
    host_service=None
    if profile.get('pageable_host_service') is not None:
        from .pageable_host_service import setup_host_service
        host_service=setup_host_service(work,profile['pageable_host_service'])
        copy=host_service['warm_copy_cpu_seconds_per_byte']
    rows=[]
    for index,p in enumerate(work['phases']):
        name=p['phase']
        if name not in bank:raise ValueError('Missing independent setup primitive: '+name)
        r=bank[name]
        if work['samples']<r['reference_shape'][0]:raise ValueError('Workload smaller than independent setup reference')
        ref=phase_work(name,*r['reference_shape'])
        extra={key:max(0,p[key]-value) for key,value in ref.items()}
        fixed=positive('cpu_seconds',r['cpu_seconds'],True)/q+positive('non_cpu_seconds',r['non_cpu_seconds'],True)
        bulk=max(extra['gemm_flops']/flops,extra['gemm_bytes']/bw)
        bulk+=max(extra['vector_ops']/flops,extra['vector_bytes']/bw)
        reduction=None
        if name=='design_block' and profile.get('reduction_gpu_properties') is not None:
            from .reduction_memory import column_sum_workspace,column_std_workspace
            reduction=column_sum_workspace(work['samples'],p['traits'],**profile['reduction_gpu_properties'])
            bulk+=reduction['partial_sum_read_write_bytes']/bw
        if name=='residual_block' and profile.get('reduction_gpu_properties') is not None:
            from .reduction_memory import column_sum_workspace,column_std_workspace
            reduction={label:fn(work['samples'],p['traits'],**profile['reduction_gpu_properties'])
                for label,fn in [('mean',column_sum_workspace),('std',column_std_workspace)]}
            bulk+=sum(r['partial_sum_read_write_bytes'] for r in reduction.values())/bw
        host_extra=dict(touch_cpu_seconds=0.,allocation_cpu_seconds=0.,release_cpu_seconds=0.,serial_cpu_seconds=0.,zero_bytes=0.)
        warm_d2h=0.
        if host_service is not None:
            host_extra=host_service['phases'][index]
            warm_d2h=extra['d2h_bytes']*host_service['warm_d2h_cpu_seconds_per_byte']/q
        # The warm download CPU price includes transfer waiting. CPU and link
        # service are simultaneous demands, not two consecutive transfers.
        transfer=extra['h2d_bytes']/h2d+max(extra['d2h_bytes']/d2h,warm_d2h)
        pages=sum(host_extra[key] for key in ['touch_cpu_seconds','allocation_cpu_seconds','release_cpu_seconds'])/q
        host=extra['host_copy_bytes']*copy/q
        rows.append(dict(phase=name,traits=p['traits'],fixed_seconds=fixed,
                         fixed_cpu_seconds=r['cpu_seconds']/q,fixed_non_cpu_seconds=r['non_cpu_seconds'],bulk_seconds=bulk,
                         transfer_seconds=transfer,host_copy_seconds=host,pageable_d2h_cpu_seconds=warm_d2h,
                         host_page_seconds=pages,host_page_serial_seconds=host_extra['serial_cpu_seconds']/q,
                         host_dram_bytes=2*p['host_copy_bytes']+host_extra['zero_bytes']+(p['h2d_bytes']+p['d2h_bytes'] if host_service is not None else 0),
                         reduction=reduction,seconds=fixed+bulk+transfer+host+pages))
    return dict(seconds=sum(r['seconds'] for r in rows),phases=rows,observation_cpu_work=work['observation_cpu_work'],
        pageable_host=host_service,
        cleanup=dict(cpu_seconds=0.,serial_cpu_seconds=0.) if host_service is None else host_service['cleanup'],
        unpriced_terms=(['prevalidated count-vector validation service'] if work['observed_counts_reused'] else ['CPU phenotype observed-count scans require independent prices'])+['setup GEMM shape efficiency and allocator-reservation service','CUDA reduction synchronization and occupancy service']+(['first-touch host pages, CPU allocation/release and pageable CUDA transfer service'] if host_service is None else host_service['unpriced_terms'])+(['input layout unknown; setup copy count uses a conservative bound'] if work['input_contiguous'] is None else [])+([] if profile.get('reduction_gpu_properties') is not None else ['global setup-reduction partial traffic requires GPU properties']),
        scope='Independent warm tiny phase costs plus additional source work. '
              'Peak-throughput/logical-HBM approximation, not calibrated shape timings. '
              'Cold library initialization must be added separately.')
