"""Conservative eager tensor-storage accounting, without measured multipliers."""
from functools import lru_cache
from .tensor_work import eager_statistics_work

@lru_cache(maxsize=128)
def eager_scan_memory(samples, chunk, traits, covariates, depth, *, reduction=None, return_beta=True):
    if any(isinstance(v,bool) or not isinstance(v,int) or v<1 for v in (samples,chunk,traits,depth)):
        raise ValueError('Positive integer scan dimensions required')
    if isinstance(covariates,bool) or not isinstance(covariates,int) or covariates<0:
        raise ValueError('Nonnegative integer covariate rank required')
    if reduction not in (None,'jagwas','device_significant'):
        raise ValueError('Unsupported eager reduction memory layout')
    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced memory layout')
    work=eager_statistics_work(samples,chunk,traits,covariates,validate_range=True)
    initial={s['storage'] for s in work['initial_storages']}
    temporaries={}
    for step in work['steps']:
        for s in step['inputs']+step['outputs']:
            if s['device']=='meta' and s['storage'] not in initial:
                temporaries[s['storage']]=max(temporaries.get(s['storage'],0),s['storage_bytes'])
    # Sum every distinct intermediate storage as if it coexisted. Views alias
    # storage and add nothing. This intentionally avoids speculative early frees.
    rounded=lambda size:((size+511)//512)*512
    temporary_bytes=sum(rounded(size) for size in temporaries.values())
    outputs={s['storage']:s['storage_bytes'] for s in work['result_storages']}
    reduction_work=None;factor=0;previous_result=0
    if not return_beta:
        # Beta is still computed. Its previous local survives the next RHS,
        # but no beta storage remains pending on the asynchronous D2H stream.
        beta_storage=work['result_storages'][0]
        previous_result=rounded(beta_storage['storage_bytes'])
        outputs.pop(beta_storage['storage'])
    if reduction=='jagwas':
        from .reduction_tensor_work import jagwas_tensor_work
        reduction_work=jagwas_tensor_work(samples,chunk,traits,phase='reduce')
        # The dense statistics are still live during reduction. Add every
        # reduction temporary, but retain only narrow results across D2H slots.
        temporary_bytes+=reduction_work['distinct_temporary_bytes']
        outputs={s['storage']:s['storage_bytes'] for s in reduction_work['result_storages']}
        factor=rounded(8*traits*traits)
        # staged_values in the prior loop iteration may survive the next RHS.
        previous_result=sum(rounded(size) for size in outputs.values())
    if reduction=='device_significant':
        # Dense output locals survive evaluation of the next statistics RHS,
        # but there are no outstanding dense D2H slots in this synchronous path.
        previous_result=sum(rounded(size) for size in outputs.values())
    retained=0 if reduction=='device_significant' else (depth-1)*sum(rounded(size) for size in outputs.values())
    resident=4*samples*(traits+covariates+1)+4*traits
    # native_scan retains intercept/covariates while the generator is alive.
    critical=rounded(4*(samples+1)) if reduction=='device_significant' else 0
    resident+=4*samples*(covariates+2)+factor+critical
    staging=depth*samples*chunk
    # The caller's old genotype reference can survive evaluation of its next
    # conversion before assignment replaces it.
    previous=4*samples*chunk
    live=eager_live_storage(samples,chunk,traits,covariates)
    joint_temporary=0 if reduction_work is None else reduction_work['distinct_temporary_bytes']
    result=dict(tensor_storage_budget=resident+staging+temporary_bytes+retained+previous+previous_result,
        lifetime_candidate_bytes=resident+staging+live['temporary_live_bytes']+joint_temporary+retained+previous+previous_result,
        lifetime_trace=live,
        resident_bytes=resident,staging_bytes=staging,
        distinct_temporary_bytes=temporary_bytes,retained_output_bytes=retained,
        previous_genotype_bytes=previous,distinct_temporary_storages=len(temporaries),
        source_sha256=work['source_sha256'],
        scope='Conservative all-intermediate storage sum for native int8 eager FP32 scan. Includes hidden bool-to-int64 promotion and asynchronous output retention. Not an allocator-reserved or library-workspace guarantee.')
    if reduction=='device_significant':
        result.update(reduction=reduction,critical_device_bytes=critical,
            previous_dense_result_bytes=previous_result,selection_memory_included=False)
    if not return_beta:
        result.update(return_beta=False,previous_beta_bytes=previous_result)
    if reduction_work is not None:
        result.update(reduction='jagwas',reduction_work=reduction_work,
            persistent_factor_bytes=factor,previous_reduced_result_bytes=previous_result)
        result['source_sha256']=dict(result['source_sha256'],**reduction_work['source_sha256'])
        result['distinct_temporary_storages']+=reduction_work['distinct_temporary_storages']
    return result

def eager_live_storage(samples, chunk, traits, covariates):
    """Observe Python tensor lifetimes on meta; never time association work."""
    import weakref
    import torch
    from torch.utils._python_dispatch import TorchDispatchMode
    from torch.utils._pytree import tree_flatten
    from .linear import _dosage_statistics
    storage_objects={};sizes={};counts={};references={};excluded=set()
    peak=0;peak_op=None
    def register(t):
        if t.device.type!='meta' or id(t) in references:return
        storage=t.untyped_storage();key=storage._cdata;identity=id(t)
        storage_objects[key]=storage  # prevent storage-handle identity reuse
        sizes[key]=((storage.nbytes()+511)//512)*512
        counts[key]=counts.get(key,0)+1
        def gone(ref):
            counts[key]-=1
            references.pop(identity,None)
        references[identity]=weakref.ref(t,gone)
    def update(name,extra=0):
        nonlocal peak,peak_op
        value=sum(sizes[key] for key,count in counts.items() if count and key not in excluded)+extra
        if value>peak:peak=value;peak_op=name
    class Observe(TorchDispatchMode):
        def __torch_dispatch__(self,func,types,args=(),kwargs=None):
            kwargs=kwargs or {}
            inputs=[v for v in tree_flatten((args,kwargs))[0] if isinstance(v,torch.Tensor)]
            for t in inputs:register(t)
            output=func(*args,**kwargs)
            outputs=[v for v in tree_flatten(output)[0] if isinstance(v,torch.Tensor)]
            for t in outputs:register(t)
            extra=0
            if str(func)=='aten.sum.dim_IntList' and inputs[0].dtype==torch.bool:
                extra=((inputs[0].numel()*8+511)//512)*512
            update(str(func),extra)
            return output
    x=torch.empty((chunk,samples),device='meta',dtype=torch.int8)
    design=torch.empty((samples,traits+covariates+1),device='meta')
    ss=torch.empty(traits,device='meta')
    for t in (x,design,ss):
        register(t);excluded.add(t.untyped_storage()._cdata)
    with Observe():
        genotype=torch.where(x==-9,torch.nan,x.to(torch.float32))
        result=_dosage_statistics(genotype,design,ss,traits,samples-covariates-2,
            True,covariate_rank=covariates)
    # Keep outputs and the caller's converted input live through the trace.
    assert len(result)==4
    return dict(temporary_live_bytes=peak,peak_operation=peak_op,
        scope='Meta tensor Python-reference lifetimes with alias deduplication and hidden bool-sum promotion. Excludes asynchronous deferred frees and library workspace.')


def eager_memory_plan(samples, chunk, traits, covariates, depth, device_profile,
                      *, preprocessing_traits=None, reduction=None, return_beta=True):
    """Compose FP32 native tensor and library requests under explicit resources."""
    from .setup_work import setup_memory
    from .cublas_memory import cublas_workspace
    required={'torch_version','sm_count','max_threads_per_sm','compute_capability',
              'cublas_workspace_config','cublas_handle_stream_pairs'}
    missing=required-set(device_profile)
    if missing:raise ValueError('Missing device memory profile: '+', '.join(sorted(missing)))
    if device_profile['torch_version'].split('+')[0]!='2.5.1':
        raise ValueError('Library memory rules require PyTorch 2.5.1')
    scan=eager_scan_memory(samples,chunk,traits,covariates,depth,reduction=reduction,return_beta=return_beta)
    setup=setup_memory(samples,traits,covariates,sm_count=device_profile['sm_count'],
        max_threads_per_sm=device_profile['max_threads_per_sm'])
    full=setup_memory(samples,traits if preprocessing_traits is None else preprocessing_traits,covariates,
        sm_count=device_profile['sm_count'],max_threads_per_sm=device_profile['max_threads_per_sm'])
    blas=cublas_workspace(device_profile['compute_capability'],device_profile['cublas_workspace_config'],
        handle_stream_pairs=device_profile['cublas_handle_stream_pairs'])
    library=blas['total_bytes']
    factor_setup=None;factor_workspace=None
    setup_peak=max(setup['design_device_live_bytes_upper'],full['residual_device_live_bytes_upper'])+library
    if reduction=='device_significant':
        setup_peak+=scan['critical_device_bytes']
    if reduction=='jagwas':
        if preprocessing_traits is not None and preprocessing_traits!=traits:
            raise ValueError('JAGWAS requires the whole phenotype panel')
        from .reduction_tensor_work import jagwas_tensor_work
        factor_setup=jagwas_tensor_work(samples,chunk,traits,phase='prepare')
        joint=4*samples*traits+factor_setup['distinct_temporary_bytes']
        from .reduction_tensor_work import jagwas_cutoff_method
        method=factor_setup['method'] if 'method' in factor_setup else jagwas_cutoff_method()
        # Each cutoff has its own census: xpotrf for the rounding cutoff's
        # Cholesky, xsyevd and xgeqrf for the default eigen factor.
        if method=='rounding' and 'jagwas_factor_workspace_census' in device_profile:
            from .cusolver_memory import jagwas_factor_workspace
            factor_workspace=jagwas_factor_workspace(
                device_profile['jagwas_factor_workspace_census'],traits,device_profile)
            joint+=factor_workspace['device_rounded_bytes']
        elif method=='eigen' and 'jagwas_eigen_workspace_census' in device_profile:
            from .cusolver_memory import jagwas_eigen_factor_workspace
            factor_workspace=jagwas_eigen_factor_workspace(
                device_profile['jagwas_eigen_workspace_census'],traits,device_profile)
            joint+=factor_workspace['device_rounded_bytes']
        setup_peak=max(full['residual_device_live_bytes_upper'],joint,
            setup['design_device_live_bytes_upper']+scan['persistent_factor_bytes'])+library
    result=dict(status='incomplete_memory_candidate',prediction_complete=False,
        device_bytes=max(setup_peak,scan['tensor_storage_budget']+library),
        lifetime_candidate_bytes=max(setup_peak,scan['lifetime_candidate_bytes']+library),
        setup_bytes=setup_peak,scan=scan,setup=setup,cublas=blas,
        unresolved_memory_terms=['other CUDA reduction temporaries and library workspaces',
            'allocator reservations, cache history and driver allocations'],
        scope='Composed source memory accounting; conservative storage sum remains the selection budget. Tight lifetime candidate is diagnostic. Host memory is accounted separately.')
    if reduction=='device_significant':
        result['unresolved_memory_terms'].append('Device selection intermediates, CUB scratch and owned host payloads require the selection ledger')
    if factor_setup is not None:
        result.update(reduction='jagwas',factor_setup=factor_setup)
        if factor_workspace is None and factor_setup.get('method')=='eigen':
            result['unresolved_memory_terms'].append('JAGWAS FP64 eigh (syevd) and QR (geqrf) workspace and per-row reduction workspace')
        elif factor_workspace is None:
            result['unresolved_memory_terms'].append('JAGWAS FP64 Cholesky/triangular-solve workspace and per-row reduction workspace')
        else:
            result['factor_workspace']=factor_workspace
            result['factor_host_workspace_bytes']=factor_workspace['host_requested_bytes']
            result['unresolved_memory_terms'].append(
                'JAGWAS factor info tensors and per-row reduction workspace' if factor_workspace.get('method')=='eigen'
                else 'JAGWAS factor info/error tensors, triangular-solve workspace and per-row reduction workspace')
    return result
