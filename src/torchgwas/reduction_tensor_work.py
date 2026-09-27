"""Duration-free JAGWAS tensor work traced from the executable reduction.

Meta execution allocates no phenotype, correlation or genotype data. Counts
are logical tensor requests, not physical traffic or allocator/RSS bounds.
"""
import hashlib
from pathlib import Path
from .mechanistic_plan import _integer


def jagwas_cutoff_method(method=None):
    """'eigen' (the default eigen truncation) or 'rounding' (rcond=0: the rounding cutoff over traits).

    None takes the configured default (TORCHGWAS_JAGWAS_RCOND). The two issue
    different factor operations and a differently oriented projection, so
    every ledger below is per method. A trait-dropping threshold
    (min_residual) adds one reduction to the rounding path; it is unpriced.
    """
    if method is None:
        from .jagwas_blocks import cutoff_settings
        method='eigen' if cutoff_settings()[0] is not None else 'rounding'
    if method not in ('eigen','rounding'):
        raise ValueError('Explicit eigen/rounding JAGWAS cutoff method required')
    return method


def jagwas_tensor_work(samples, markers, traits, *, phase, compute_dtype="float32", method=None):
    """Source tensor operations for factor preparation or one scan reduction.

    Initial inputs and every distinct intermediate are reported separately.
    The all-intermediate sum is conservative about Python reference lifetimes;
    cuBLAS/cuSOLVER workspace and asynchronous retention are explicitly absent.
    No timing coefficient or association measurement enters this function.
    The eigen factor keeps k <= K rows, known only once R's spectrum is read:
    the trace takes k = K.
    """
    from .structural_tensor_cache import cached_tensor_trace
    from .jagwas_projection import JagwasReduction, JagwasRankSelection, gram_rows, triangular_blocks
    method=jagwas_cutoff_method(method)
    return cached_tensor_trace('jagwas.v3',dict(samples=samples,markers=markers,traits=traits,
        phase=phase,compute_dtype=compute_dtype,method=method),
        lambda:_trace_jagwas_tensor_work(samples,markers,traits,phase=phase,compute_dtype=compute_dtype,method=method),
        implementation=(_trace_jagwas_tensor_work,JagwasReduction.__init__,JagwasReduction.prepare,
            JagwasReduction._prepare,JagwasReduction._eigen_factor,JagwasReduction._set_upper_factor,
            JagwasReduction._pivoted_selection,JagwasReduction._subset_factor,JagwasReduction._set_factor,
            JagwasReduction.reduce,JagwasRankSelection.__init__,JagwasRankSelection.trace_limit,
            gram_rows,triangular_blocks))


def _trace_jagwas_tensor_work(samples,markers,traits,*,phase,compute_dtype,method):
    for name,value in [('samples',samples),('markers',markers),('traits',traits)]:
        _integer(name,value)
    if phase not in ('prepare','reduce'):
        raise ValueError('Explicit prepare/reduce JAGWAS phase required')
    if compute_dtype not in ('float32','float64'):
        raise ValueError('Explicit float32/float64 compute dtype required')
    import torch
    from torch.utils._python_dispatch import TorchDispatchMode
    from torch.overrides import TorchFunctionMode
    from torch.utils._pytree import tree_flatten
    from .jagwas_projection import JagwasReduction
    from .jagwas_blocks import DEFAULT_RCOND
    reduction=JagwasReduction(rcond=DEFAULT_RCOND if method=='eigen' else 0)
    dtype=getattr(torch,compute_dtype)
    storage_ids={}; keep=[]; steps=[];host_calls=[];active_call=[]
    def describe(tensor):
        storage=tensor.untyped_storage(); key=storage._cdata
        if key not in storage_ids:storage_ids[key]=len(storage_ids)
        return dict(storage=storage_ids[key],shape=list(tensor.shape),stride=list(tensor.stride()),
            dtype=str(tensor.dtype),device=str(tensor.device),bytes=tensor.numel()*tensor.element_size(),
            storage_bytes=storage.nbytes(),offset_bytes=tensor.storage_offset()*tensor.element_size())
    class HostRecord(TorchFunctionMode):
        def __torch_function__(self,func,types,args=(),kwargs=None):
            call_id=len(host_calls);kwargs=kwargs or {}
            inputs=[value for value in tree_flatten((args,kwargs))[0] if isinstance(value,torch.Tensor)]
            host_calls.append(dict(id=call_id,name=getattr(func,'__name__',str(func)),
                input_dtypes=[str(value.dtype) for value in inputs],input_shapes=[list(value.shape) for value in inputs],
                kwargs={key:str(value) for key,value in kwargs.items() if not isinstance(value,torch.Tensor)},step_indices=[]))
            active_call.append(call_id)
            try:return func(*args,**kwargs)
            finally:active_call.pop()
    class Record(TorchDispatchMode):
        def __torch_dispatch__(self,func,types,args=(),kwargs=None):
            kwargs=kwargs or {}
            inputs=[value for value in tree_flatten((args,kwargs))[0] if isinstance(value,torch.Tensor)]
            before=[describe(value) for value in inputs]
            result=func(*args,**kwargs)
            outputs=[value for value in tree_flatten(result)[0] if isinstance(value,torch.Tensor)]
            after=[describe(value) for value in outputs]
            keep.extend(inputs);keep.extend(outputs)
            name=str(func)
            alias=bool(after) and all(value['storage'] in {row['storage'] for row in before} for value in after) and not func._schema.is_mutable
            allocate_only=name.startswith(('aten.empty','aten.new_empty'))
            fill=name.startswith(('aten.full_like','aten.zeros_like','aten.ones_like','aten.eye','aten.empty_like'))
            if name.endswith('.out') and isinstance(kwargs.get('out'),torch.Tensor):
                # An out= product writes its output view; it does not read it.
                before=[describe(value) for value in tree_flatten((args,{k:v for k,v in kwargs.items() if k!='out'}))[0]
                        if isinstance(value,torch.Tensor)]
            read={} if fill else {value['storage']:value for value in before}
            write={value['storage']:value for value in after}
            reads=0 if alias or allocate_only else sum(value['bytes'] for value in read.values())
            writes=0 if alias or allocate_only else sum(value['bytes'] for value in write.values())
            row=dict(op=name,inputs=before,outputs=after,alias_only=alias,allocation_only=allocate_only,
                shape_only_inputs=fill,phase=phase,host_call_id=active_call[-1] if active_call else None,
                read_bytes=reads,write_bytes=writes,logical_bytes=reads+writes)
            if name in ('aten.mm.default','aten.mm.out'):
                row['matmul_flops']=2*before[0]['shape'][0]*before[0]['shape'][1]*before[1]['shape'][1]
                row['matmul_dtype']=before[0]['dtype']
            elif name=='aten.addmm_.default':
                # Accumulating product (the blocked FP64 Gram): self += mat1 @ mat2.
                mat1,mat2=[describe(value) for value in args[1:3]]
                row['matmul_flops']=2*mat1['shape'][0]*mat1['shape'][1]*mat2['shape'][1]
                row['matmul_dtype']=mat1['dtype']
            steps.append(row)
            return result
    if phase=='prepare':
        phenotype=torch.empty((samples,traits),device='meta',dtype=dtype)
        inputs=[phenotype]
    else:
        beta=torch.empty((markers,traits),device='meta',dtype=dtype)
        stat=torch.empty_like(beta)
        status=torch.empty(markers,device='meta',dtype=torch.uint8)
        df=torch.empty(markers,device='meta',dtype=dtype)  # the scan's per-variant df has its compute dtype
        factor=torch.empty((traits,traits),device='meta',dtype=torch.float64)
        if method=='eigen':reduction._set_upper_factor(factor,traits)
        else:reduction._set_factor(factor)
        inputs=[beta,stat,status,df,reduction._inverse_cholesky]
    initial=[describe(value) for value in inputs]
    with HostRecord(),Record():
        if phase=='prepare':
            reduction.prepare(phenotype,device='meta')
            result=[reduction._inverse_cholesky]
        else:
            result=list(reduction.reduce(beta,stat,status,df,1))
    for index,step in enumerate(steps):
        if step['host_call_id'] is not None:host_calls[step['host_call_id']]['step_indices'].append(index)
    unpriced_host_calls=[call for call in host_calls if not call['step_indices']]
    host_calls=[call for call in host_calls if call['step_indices']]
    initial_ids={value['storage'] for value in initial}
    temporaries={}
    for step in steps:
        for value in step['inputs']+step['outputs']:
            if value['storage'] not in initial_ids:
                temporaries[value['storage']]=max(temporaries.get(value['storage'],0),value['storage_bytes'])
    rounded=lambda size:((size+511)//512)*512
    sources=['jagwas_projection.py','jagwas_blocks.py']
    return dict(phase=phase,compute_dtype=compute_dtype,method=method,samples=samples,markers=markers,traits=traits,steps=steps,
        host_calls=host_calls,unpriced_host_calls=unpriced_host_calls,
        initial_storages=initial,result_storages=[describe(value) for value in result],
        distinct_temporary_bytes=sum(rounded(value) for value in temporaries.values()),
        distinct_temporary_storages=len(temporaries),logical_bytes=sum(step['logical_bytes'] for step in steps),
        tensor_dispatches=len(steps),
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources},
        torch_version=torch.__version__,
        unpriced_terms=['Actual compiled kernel geometry and CPU launch service',
            'FP64 GEMM, factorization (eigh and QR, or Cholesky and triangular solve) service and library workspace',
            'CUDA allocator retention, asynchronous pending outputs, host staging and driver allocations'],
        scope='Source JAGWAS meta tensor trace with alias deduplication and rounded all-intermediate storage requests. It does not authorize timing prediction or memory admission without surrounding scan/setup and library terms.')


def jagwas_factor_memory_floor(samples, traits, *, compute_dtype='float32', method=None):
    """Necessary explicit tensor capacity of factor preparation.

    Eigen (default), at k = K kept rows: the phenotype and correlation stay
    live throughout (JagwasReduction._prepare holds both). eigh adds the
    eigenvectors and values; scaling them adds Lambda^-1/2 U' before the
    eigenvectors are released; QR then holds that and R. Rounding (rcond=0):
    the phenotype, correlation, Cholesky factor, identity and inverse factor
    coexist during the triangular solve (the full-rank path; the collinear
    fallback factors a smaller kept block). Either way the FP64 Gram holds the
    phenotype, the FP64 correlation and one FP64 sample block before that; the
    largest stage is the floor. This excludes workspace, allocator retention
    and the scan, so passing is not admission.
    """
    _integer('samples',samples);_integer('traits',traits)
    if compute_dtype not in ('float32','float64'):
        raise ValueError('Explicit float32/float64 compute dtype required')
    method=jagwas_cutoff_method(method)
    from .jagwas_blocks import gram_rows
    itemsize=4 if compute_dtype=='float32' else 8
    square=8*traits*traits
    if method=='eigen':
        arrays=dict(phenotype=itemsize*samples*traits,correlation=square,eigenvectors=square,
            eigenvalues=8*traits,scaled_eigenvectors=square,factor=square)
        base=arrays['phenotype']+arrays['correlation']
        factor_bytes=max(base+arrays['eigenvectors']+arrays['eigenvalues']+arrays['scaled_eigenvectors'],
            base+arrays['scaled_eigenvectors']+arrays['factor'])
        persistent=arrays['factor']
    else:
        arrays=dict(phenotype=itemsize*samples*traits,correlation=square,
            cholesky=square,identity=square,inverse_cholesky=square)
        factor_bytes=sum(arrays.values());persistent=arrays['inverse_cholesky']
    # An FP64 panel's blocks are views; an FP32 panel's are FP64 copies.
    arrays['gram_block']=(0 if itemsize==8 else 8*gram_rows(samples,traits)*traits)
    gram_bytes=arrays['phenotype']+arrays['correlation']+arrays['gram_block']
    sources=['jagwas_projection.py','jagwas_blocks.py']
    return dict(samples=samples,traits=traits,compute_dtype=compute_dtype,method=method,
        arrays=arrays,explicit_live_bytes=max(factor_bytes,gram_bytes),
        persistent_factor_bytes=persistent,
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources},
        scope='Necessary explicit factor-preparation tensor capacity only. Passing does not establish scan, workspace, allocator or host-memory feasibility.')


def require_jagwas_factor_capacity(samples, traits, devices, *, compute_dtype='float32', method=None,
                                   group_sizes=None):
    """Refuse a definitely oversized joint factor on any active CUDA device.

    group_sizes (JagwasGroups): groups factor one at a time, so the floor is
    the whole panel, the largest group's live factor matrices, and every other
    group's retained FP64 factor (at most k x k).
    """
    import torch
    if group_sizes:
        sizes=[int(size) for size in group_sizes]
        largest=max(sizes)
        work=dict(jagwas_factor_memory_floor(samples,largest,compute_dtype=compute_dtype,method=method))
        panel_largest=4*int(samples)*largest
        work['explicit_live_bytes']=(4*int(samples)*int(traits)+work['explicit_live_bytes']-panel_largest
                                     +8*(sum(size*size for size in sizes)-largest*largest))
        work['group_sizes']=sizes
    else:
        work=jagwas_factor_memory_floor(samples,traits,compute_dtype=compute_dtype,method=method)
    for name in devices:
        device=torch.device(name)
        if device.type!='cuda':continue
        free,_=torch.cuda.mem_get_info(device)
        # Cached but unallocated blocks can satisfy future requests. Fragmented
        # blocks are not guaranteed reusable; this is only an impossibility gate.
        available=int(free)+int(torch.cuda.memory_reserved(device))-int(torch.cuda.memory_allocated(device))
        if work['explicit_live_bytes']>available:
            raise ValueError(
                f"reduce='jagwas' requires at least {work['explicit_live_bytes']} bytes "
                f"for the full phenotype and live factor matrices on {device}, "
                f"against {available} bytes available. Workspace and scan buffers "
                "need additional memory. The joint statistic cannot be trait-tiled; "
                "reduce='significant' supports trait tiling.")
    return work


def jagwas_host_primitive_name(call):
    """Typed fixed-operation keys for independent joint-reduction CPU probes."""
    name=call['name'];dtypes=call['input_dtypes'];kwargs=call['kwargs']
    width='fp32' if dtypes[:1]==['torch.float32'] else 'fp64'
    if name=='double':return 'joint_fp32_to_fp64' if dtypes==['torch.float32'] else 'joint_fp64_noop'
    if name=='isfinite':return 'joint_isfinite_'+width
    if name=='all':return 'joint_all_bool_rows'
    if name=='__invert__':return 'joint_invert_bool'
    if name=='nan_to_num':return 'joint_nan_to_num_'+width
    # Score form z = t / sqrt(1 + t^2 / df), in the scan precision.
    if name=='square':return 'joint_square_out_'+width
    if name=='div_':return 'joint_div_rows_inplace_'+width
    if name=='add_':return 'joint_add_scalar_inplace_'+width
    if name=='rsqrt_':return 'joint_rsqrt_inplace_'+width
    if name=='mul_':return 'joint_mul_inplace_'+width
    if name=='__get__':return 'joint_transpose_view'
    if name=='mm':return 'joint_gemm_out_fp64'
    if name=='empty':return 'joint_empty_fp64'
    if name=='__getitem__':return 'joint_slice_view'
    if name=='square_':return 'joint_square_inplace_fp64'
    if name=='sum':return 'joint_sum_fp64_columns'
    if name in ('ne','__ne__'):return 'joint_status_compare'
    if name=='__or__':return 'joint_or_bool'
    if name=='unsqueeze':
        return {'torch.float64':'joint_unsqueeze_view_fp64','torch.float32':'joint_unsqueeze_view_fp32'}.get(
            dtypes[0] if dtypes else None,'joint_unsqueeze_view')
    if name=='masked_fill':return 'joint_masked_fill_fp64'
    if name=='full_like':return 'joint_fill_fp64' if kwargs.get('dtype')=='torch.float64' else 'joint_fill_fp32'
    if name=='to':return 'joint_fp64_to_fp32' if kwargs.get('dtype')=='torch.float32' or dtypes==['torch.float64'] else 'joint_cast_unknown'
    if name=='zeros_like':return 'joint_fill_int32'
    raise ValueError('Unpriced joint host API '+name)


def jagwas_projection_service(samples, markers, traits, resources, kernels, *,
                              compute_dtype='float32', initial_cache=None, host_primitives):
    """Use the shared tensor/cache model with explicit FP64 arithmetic rates.

    Kernel geometry must describe the exact projection shape and math mode.
    Fixed host API prices are independent tiny operations. No reduction/scan
    elapsed timings or shape-indexed service table enter this calculation.
    """
    from .tensor_service import tensor_stage_service
    from .jagwas_projection import projection_gemm_dimensions
    work=jagwas_tensor_work(samples,markers,traits,phase='reduce',compute_dtype=compute_dtype)
    dimensions=projection_gemm_dimensions(traits,markers,upper=work['method']=='eigen')
    result=tensor_stage_service(work,resources,kernels,initial_cache=initial_cache,
        host_primitives=host_primitives,gemm_dimensions=dimensions[0] if len(dimensions)==1 else dimensions,
        gemm_dtype='float64',host_primitive_resolver=jagwas_host_primitive_name)
    if result.get('status')=='zero_available_capacity':return result
    result['unpriced_terms']+=['Joint row-reduction workspace/internal traffic and operation-specific instruction latency',
        'Source no-op/property host calls and allocator reservation service']
    result.update(reduction='jagwas',prediction_complete=False,
        compute_dtype=compute_dtype,unpriced_host_calls=work['unpriced_host_calls'])
    return result


def append_jagwas_projection(scan, projection, cpu_fraction):
    """Compose eager host submission and same-stream GPU work before D2H.

    The two component cache traces have independent storage namespaces. The
    projection therefore starts from its explicitly empty cache, not from an
    invented reuse fraction. This boundary approximation remains unpriced.
    """
    from .first_principles import positive
    q=positive('CPU fraction',cpu_fraction)
    if q>1:raise ValueError('CPU fraction must be no greater than one')
    if projection.get('reduction')!='jagwas':raise ValueError('Explicit JAGWAS projection required')
    host_offset=scan['host_dispatch_cpu_seconds']/q
    id_offset=1+max((call['call_id'] for call in scan.get('host_calls',[])),default=-1)
    calls=list(scan.get('host_calls',[]))+[
        dict(call,call_id=call['call_id']+id_offset,submit_finish=call['submit_finish']+host_offset)
        for call in projection['host_calls']]
    operations=[dict(op) for op in scan['operations']]+[
        dict(op,phase='joint_projection',host_submit_finish=op['host_submit_finish']+host_offset)
        for op in projection['operations']]
    end=0.
    for op in operations:
        end=max(end,op['host_submit_finish'])+op['kernel_service_seconds']
        op['gpu_finish']=end
    total_host=scan['host_dispatch_cpu_seconds']+projection['host_dispatch_cpu_seconds']
    return dict(scan,operations=operations,host_calls=calls,
        host_dispatch_cpu_seconds=total_host,estimated_span_seconds=max(end,total_host/q),
        kernel_service_seconds=scan['kernel_service_seconds']+projection['kernel_service_seconds'],
        kernel_count=scan.get('kernel_count',len(scan['operations']))+projection['kernel_count'],
        modeled_hbm_bytes=scan.get('modeled_hbm_bytes',0)+projection['modeled_hbm_bytes'],
        logical_bytes=scan.get('logical_bytes',0)+projection['logical_bytes'],
        source_sha256=dict(scan.get('source_sha256',{}),**projection['source_sha256']),
        reduction_service_seconds=projection['kernel_service_seconds'],
        reduction_component=projection,reduction='jagwas',prediction_complete=False,
        unpriced_terms=list(scan.get('unpriced_terms',[]))+list(projection['unpriced_terms'])+[
            'Statistics-to-joint projection cache state: the projection uses an explicitly empty modeled cache'])


def jagwas_factor_arithmetic_work(samples, traits, *, compute_dtype='float32', method=None):
    """Mathematical operation ledger for the source's factor algorithm, per cutoff method."""
    method=jagwas_cutoff_method(method)
    if method=='eigen':
        return _eigen_factor_arithmetic_work(samples,traits,compute_dtype=compute_dtype)
    return _rounding_factor_arithmetic_work(samples,traits,compute_dtype=compute_dtype)


def _gram_expected(samples,traits,fp32):
    from .jagwas_blocks import gram_rows
    blocks=-(-samples//gram_rows(samples,traits))
    expected=[('correlation','aten.zeros.default')]
    for _ in range(blocks):
        expected+=([('cast','aten._to_copy.default')] if fp32 else [])+[('correlation','aten.addmm_.default')]
    return blocks,expected+[('correlation','aten.div_.Tensor')]


def _eigen_factor_arithmetic_work(samples,traits,*,compute_dtype):
    """Eigen truncation at k = K kept rows: FP64 Gram, eigh, one spectrum read, scaling and QR.

    Conventional counts, not cuSOLVER's issued instructions:
    - eigh (syevd with vectors): tridiagonal reduction (2/3)K^3 and
      back-transformation of the K eigenvectors K^3 multiply-adds; the
      divide-and-conquer on the tridiagonal matrix depends on deflation and is
      unpriced;
    - the spectrum read: K FP64 values to the host, which picks k;
    - scaling: K square roots and k K divisions for Lambda_k^-1/2 U_k';
    - QR of that k x K matrix (Householder, R only): k^2 (K - k/3)
      multiply-adds, one square root and K - j divisions per reflector j.
    k is data; pricing at k = K bounds the scaling and QR from above.
    """
    work=jagwas_tensor_work(samples,1,traits,phase='prepare',compute_dtype=compute_dtype,method='eigen')
    active=[step for step in work['steps'] if not step['alias_only'] and not step['allocation_only']]
    fp32=compute_dtype=='float32'
    n,k=samples,traits;itemsize=4 if fp32 else 8
    blocks,expected=_gram_expected(n,k,fp32)
    expected+=[('eigh','aten._linalg_eigh.default'),('scale','aten.sqrt.default'),('scale','aten.div.Tensor'),
               ('qr','aten.linalg_qr.default')]
    if [step['op'] for step in active]!=[op for _,op in expected]:
        raise ValueError('Factor source operations changed; reconcile arithmetic ledger')
    phase_bytes={}
    for (name,_),step in zip(expected,active):
        phase_bytes[name]=phase_bytes.get(name,0)+step['logical_bytes']
    gram=[step for step in active if step['op'] in ('aten.mm.default','aten.addmm_.default')]
    if sum(step['matmul_flops'] for step in gram)!=2*n*k*k:
        raise ValueError('Factor correlation dimensions differ from dense Gram matrix')
    square=k*k
    correlation=dict(dtype='float64',useful_flops=2*n*square,normalization_divisions=square,sample_blocks=blocks)
    eigh=dict(dtype='float64',multiply_add_flops=(2*k**3)//3+k**3)
    spectrum=dict(d2h_bytes=8*k)
    scale=dict(dtype='float64',square_roots=k,divisions=square)
    qr=dict(dtype='float64',multiply_add_flops=k**3-(k**3)//3,square_roots=k,divisions=k*(k+1)//2)
    return dict(samples=n,traits=k,compute_dtype=compute_dtype,method='eigen',
        h2d_bytes=itemsize*n*k,persistent_factor_bytes=8*square,
        correlation=correlation,eigh=eigh,spectrum=spectrum,scale=scale,qr=qr,
        fp64_multiply_add_flops=correlation['useful_flops']+eigh['multiply_add_flops']+qr['multiply_add_flops'],
        fp64_divisions=scale['divisions']+qr['divisions'],fp64_square_roots=scale['square_roots']+qr['square_roots'],
        tensor_operations=[dict(op=step['op'],logical_bytes=step['logical_bytes']) for step in active],
        phase_logical_bytes=phase_bytes,
        source_sha256=work['source_sha256'],prediction_complete=False,
        unpriced_terms=['Closed-source eigh (tridiagonalization, divide and conquer, back-transformation) and QR issued work',
            'Divide-and-conquer work, which depends on deflation in the spectrum',
            'Factorization dependency latency, library workspace traffic and synchronization',
            'Host upload, pageable staging, dispatch, allocation and factor release service',
            'The kept count k <= K: scaling, QR and the factor are priced at k = K'],
        scope='Source-checked useful arithmetic per independently prepared device factor; conventional algorithm counts, not physical traffic or execution time.')


def _rounding_factor_arithmetic_work(samples,traits,*,compute_dtype):
    """Mathematical operation ledger for the rounding cutoff's full-rank factor algorithm (rcond=0).

    Counts describe conventional dot-product Cholesky and a dense triangular
    solve with K right-hand sides. They are useful arithmetic, not cuSOLVER's
    issued instructions, dependency/occupancy model, or a runtime prediction.
    In particular, the identity RHS does not turn the dense TRSM API into a
    sparse inverse algorithm. No association or factor timings enter here.
    The rank-revealing factor forms R with an FP64 Gram in sample blocks
    (each FP32 block cast to FP64 first) and adds the rank check:
    ||L^-1||_F (tr R^-1) and one read of two FP64 scalars. Its collinear-panel
    fallback (host dpstrf, the prefix traces, then the kept block) is
    data-dependent and unpriced.
    """
    work=jagwas_tensor_work(samples,1,traits,phase='prepare',compute_dtype=compute_dtype,method='rounding')
    active=[step for step in work['steps'] if not step['alias_only'] and not step['allocation_only']]
    fp32=compute_dtype=='float32'
    n,k=samples,traits;itemsize=4 if fp32 else 8
    blocks,expected=_gram_expected(n,k,fp32)
    expected+=[('cholesky','aten.linalg_cholesky_ex.default'),
               ('identity','aten.eye.default'),('solve','aten.linalg_solve_triangular.default')]
    # info cast to FP64, ||L^-1||_F, one stacked read.
    expected+=[('rank_check',op) for op in ['aten._to_copy.default','aten.linalg_vector_norm.default','aten.stack.default']]
    if [step['op'] for step in active]!=[op for _,op in expected]:
        raise ValueError('Factor source operations changed; reconcile arithmetic ledger')
    phase_bytes={}
    for (name,_),step in zip(expected,active):
        phase_bytes[name]=phase_bytes.get(name,0)+step['logical_bytes']
    gram=[step for step in active if step['op'] in ('aten.mm.default','aten.addmm_.default')]
    if sum(step['matmul_flops'] for step in gram)!=2*n*k*k:
        raise ValueError('Factor correlation dimensions differ from dense Gram matrix')
    square=k*k
    correlation=dict(dtype='float64',useful_flops=2*n*square,normalization_divisions=square,sample_blocks=blocks)
    cholesky=dict(multiply_add_flops=(k**3-k)//3,divisions=k*(k-1)//2,square_roots=k)
    solve=dict(right_hand_sides=k,multiply_add_flops=square*(k-1),divisions=square)
    check=dict(dtype='float64',multiply_add_flops=2*square,square_roots=1,d2h_bytes=16)
    gram_fp64=correlation['useful_flops']
    return dict(samples=n,traits=k,compute_dtype=compute_dtype,method='rounding',
        h2d_bytes=itemsize*n*k,persistent_factor_bytes=8*square,
        correlation=correlation,
        cholesky=dict(dtype='float64',**cholesky),triangular_solve=dict(dtype='float64',**solve),rank_check=check,
        fp64_multiply_add_flops=gram_fp64+cholesky['multiply_add_flops']+solve['multiply_add_flops']+check['multiply_add_flops'],
        fp64_divisions=cholesky['divisions']+solve['divisions'],fp64_square_roots=k+check['square_roots'],
        tensor_operations=[dict(op=step['op'],logical_bytes=step['logical_bytes']) for step in active],
        phase_logical_bytes=phase_bytes,
        source_sha256=work['source_sha256'],prediction_complete=False,
        unpriced_terms=['Closed-source blocked factorization/solve issued work and scalar versus Tensor Core split',
            'Factorization dependency latency, library workspace traffic and synchronization',
            'Host upload, pageable staging, dispatch, allocation and factor release service']
            +['Collinear-panel fallback: host FP64 copy of R, LAPACK dpstrf, prefix traces and the kept-block factor'],
        scope='Source-checked useful arithmetic per independently prepared device factor; conventional algorithm counts, not physical traffic or execution time.')
