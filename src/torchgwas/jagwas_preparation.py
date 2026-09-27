"""JAGWAS preparation from source work and independent tiny phase services.

The factor model is an explicit blocking-phase approximation. It uses one
fixed reference, resource capacities and mathematical work, never a timing grid.
Closed vendor-kernel details and loaded-context transfer remain unresolved.
"""
import copy
import math
from .execution_graph import ExecutionGraph
from .jagwas_candidate import jagwas_candidate_shape
from .pinned_work import pinned_scan_work
from .reduction_tensor_work import jagwas_cutoff_method,jagwas_factor_arithmetic_work
from .setup_work import setup_work,setup_service

FACTOR_REFERENCE=[65,32]
# Eigen truncation (the default) and the rounding cutoff over traits (rcond=0).
FACTOR_PHASES_BY_METHOD=dict(
    eigen=('upload','correlation','cast','eigh','spectrum','scale','qr'),
    rounding=('upload','correlation','cast','cholesky','identity','solve','rank_check'))
FACTOR_PHASES=FACTOR_PHASES_BY_METHOD['eigen']


def factor_phases(method=None):
    return FACTOR_PHASES_BY_METHOD[jagwas_cutoff_method(method)]


def _number(name,value,*,zero=False,fraction=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (not zero and value==0) or (fraction and value>1):
        raise ValueError('Invalid independent '+name)
    return float(value)


def factor_phase_work(samples,traits,method=None):
    """FP64 Gram in sample blocks (the cast phase is each FP32 block's FP64
    copy), then per method: eigen (default) eigh, the K-value spectrum read,
    scaling and QR; rounding (rcond=0) Cholesky, identity, the inverse-factor
    solve and the rank check."""
    work=jagwas_factor_arithmetic_work(samples,traits,method=method)
    ops=work['phase_logical_bytes'];phases=FACTOR_PHASES_BY_METHOD[work['method']]
    n,k=samples,traits
    rows={name:dict(h2d_bytes=0,d2h_bytes=0,fp32_ops=0,fp64_flops=0,fp64_special_ops=0,logical_gpu_bytes=0) for name in phases}
    rows['upload']['h2d_bytes']=4*n*k
    if work['correlation']['dtype']=='float64':
        rows['correlation'].update(fp64_flops=2*n*k*k,fp64_special_ops=k*k,logical_gpu_bytes=ops['correlation'])
    else:
        rows['correlation'].update(fp32_ops=2*n*k*k+k*k,logical_gpu_bytes=ops['correlation'])
    rows['cast'].update(logical_gpu_bytes=ops.get('cast',0))
    if work['method']=='eigen':
        rows['eigh'].update(fp64_flops=work['eigh']['multiply_add_flops'],logical_gpu_bytes=ops['eigh'])
        rows['spectrum']['d2h_bytes']=work['spectrum']['d2h_bytes']
        rows['scale'].update(fp64_special_ops=work['scale']['square_roots']+work['scale']['divisions'],
            logical_gpu_bytes=ops['scale'])
        rows['qr'].update(fp64_flops=work['qr']['multiply_add_flops'],
            fp64_special_ops=work['qr']['square_roots']+work['qr']['divisions'],logical_gpu_bytes=ops['qr'])
        return dict(samples=n,traits=k,method='eigen',phases=[dict(phase=name,**rows[name]) for name in phases],
            arithmetic=work)
    # Cholesky info is read with ||L^-1||_F in the rank check.
    rows['cholesky'].update(fp64_flops=work['cholesky']['multiply_add_flops'],
        fp64_special_ops=work['cholesky']['divisions']+work['cholesky']['square_roots'],
        logical_gpu_bytes=ops['cholesky'])
    rows['identity']['logical_gpu_bytes']=ops['identity']
    rows['solve'].update(fp64_flops=work['triangular_solve']['multiply_add_flops'],
        fp64_special_ops=work['triangular_solve']['divisions'],logical_gpu_bytes=ops['solve'])
    rows['rank_check'].update(fp64_flops=work['rank_check']['multiply_add_flops'],
        fp64_special_ops=work['rank_check']['square_roots'],logical_gpu_bytes=ops['rank_check'],
        d2h_bytes=work['rank_check']['d2h_bytes'])
    return dict(samples=n,traits=k,method='rounding',phases=[dict(phase=name,**rows[name]) for name in phases],
        arithmetic=work)


def factor_service(samples,traits,profile,*,library_arithmetic):
    """Conditional scalar/Tensor Core scenarios for opaque library arithmetic.

    The chosen capacity prices multiply/add work only. Division/sqrt remain
    scalar-equivalent work with unresolved instruction latency. Neither choice
    is asserted to bound elapsed time. CPU wait is a supplied scenario.
    """
    if library_arithmetic not in ('scalar','tensor'):
        raise ValueError('Explicit factor library arithmetic scenario required')
    if samples<FACTOR_REFERENCE[0] or traits<FACTOR_REFERENCE[1]:
        raise ValueError('Factor shape smaller than independent fixed reference')
    bank=profile['jagwas_factor_primitives']
    if bank.get('reference_shape')!=FACTOR_REFERENCE or bank.get('compute_dtype')!='float32' or bank.get('boundary')!='call_then_device_synchronize':
        raise ValueError('Independent fixed-reference factor phase boundary required')
    method=jagwas_cutoff_method()
    if bank.get('method','rounding')!=method or set(bank.get('phases',{}))!=set(FACTOR_PHASES_BY_METHOD[method]):
        raise ValueError('Complete independent factor phase bank for the '+method+' factor required')
    q=_number('cpu_fraction',profile['cpu_fraction'],fraction=True)
    wait=_number('factor wait fraction',profile['factor_wait_cpu_fraction'],zero=True,fraction=True)
    gpu=profile['gpu_resources'];g=_number('gpu_fraction',gpu['gpu_fraction'],fraction=True)
    hbm=_number('HBM capacity',gpu['hbm_bytes_per_second'])*g
    fp32=_number('FP32 capacity',gpu['fp32_flops_per_second'])*g
    scalar=_number('FP64 scalar capacity',gpu['fp64_flops_per_second'])*g
    math_rate=(scalar if library_arithmetic=='scalar' else _number('FP64 tensor capacity',gpu['fp64_tensor_flops_per_second'])*g)
    links={direction:_number(direction+' capacity',profile[direction+'_bytes_per_second']) for direction in ('h2d','d2h')}
    work=factor_phase_work(samples,traits,method);reference=factor_phase_work(*FACTOR_REFERENCE,method=method);rows=[]
    context=profile['jagwas_factor_context']
    if not isinstance(context,dict) or not context or bank.get('context')!=context:
        raise ValueError('Factor primitive execution context differs')
    if bank.get('source_sha256')!=work['arithmetic']['source_sha256']:
        raise ValueError('Factor primitive executor source differs')
    for phase,ref in zip(work['phases'],reference['phases']):
        name=phase['phase'];fixed=bank['phases'][name]
        cpu=_number(name+' fixed CPU',fixed['cpu_seconds'],zero=True)
        noncpu=_number(name+' fixed non-CPU',fixed['non_cpu_seconds'],zero=True)
        if not cpu+noncpu:raise ValueError('Positive factor work requires nonzero reference service')
        extra={key:max(0,phase[key]-ref[key]) for key in phase if key!='phase'}
        arithmetic=extra['fp32_ops']/fp32+extra['fp64_flops']/math_rate+extra['fp64_special_ops']/scalar
        bulk=max(arithmetic,extra['logical_gpu_bytes']/hbm)
        transfer=sum(extra[d+'_bytes']/links[d] for d in links)
        # A blocking reference includes its own wait; only additional device
        # service receives additional CPU wait demand. It overlaps that work.
        extra_wait_cpu=(bulk+transfer)*q*wait
        seconds=cpu/q+noncpu+bulk+transfer
        rows.append(dict(phase=name,seconds=seconds,cpu_seconds=cpu+extra_wait_cpu,
            fixed_cpu_seconds=cpu,extra_wait_cpu_seconds=extra_wait_cpu,
            h2d_bytes=phase['h2d_bytes'],d2h_bytes=phase['d2h_bytes'],
            host_dram_bytes=phase['h2d_bytes']+phase['d2h_bytes'],
            extra_gpu_seconds=bulk,extra_transfer_seconds=transfer,work=phase))
    return dict(phases=rows,seconds=sum(row['seconds'] for row in rows),work=work,
        library_arithmetic=library_arithmetic,prediction_complete=False,
        unpriced_terms=list(work['arithmetic']['unpriced_terms'])+[
            'Factor blocking-phase reference adds synchronization between normally asynchronous APIs',
            'Reference-to-large factor internal launch count and host/library dispatch scaling',
            'Division and sqrt instruction latency beyond scalar-equivalent operation counts',
            'Single-worker fixed primitive transfer to concurrently prepared factors',
            'Factor pageable upload staging, first touch and allocator release beyond reference',
            'One-time serialized 2x2 CUDA linalg initialization before concurrent factors'],
        scope='One fixed N65/K32 reference plus additional source arithmetic, logical bytes and independent resource capacities. Explicit library arithmetic/wait scenario; not a timing fit or bound.')


def build_jagwas_preparation(candidate,*,host_serial_fraction,library_arithmetic,shared_cpu_steps,finalize):
    """Build shared residualization and independent factor/design/pinned graphs.

    Callers must supply common CPU setup/count/basis and final commit service;
    neither silently defaults to zero. Repeated device work comes from source.
    """
    from .indexed_schedule import _steps
    shape=jagwas_candidate_shape(candidate)
    serial=_number('host serialization fraction',host_serial_fraction,zero=True,fraction=True)
    if not isinstance(shared_cpu_steps,list) or not shared_cpu_steps or not isinstance(finalize,list) or not finalize:
        raise ValueError('Explicit common CPU preparation and final commit services required')
    if not isinstance(library_arithmetic,dict) or set(library_arithmetic)!=set(shape['devices']):
        raise ValueError('Factor arithmetic scenario required for every device')
    shared=ExecutionGraph();last=_steps(shared,'common_cpu',shared_cpu_steps,[])
    devices={};cleanup={};unknown=set();reports={}
    n,k,c=shape['samples'],shape['traits'],shape['covariates']
    columns=shape['covariate_columns']
    def setup_subset(tile,prefix):
        ledger=setup_work(n,k,c,reuse_observed_counts=prefix=='design',
            input_contiguous=tile['data'].get('phenotype_c_contiguous'),covariate_columns=columns)
        ledger['phases']=[phase for phase in ledger['phases'] if phase['phase'].startswith(prefix)]
        return ledger,setup_service(ledger,tile['profile'])
    def append_setup(graph,last,device,ledger,estimate,profile):
        q=profile['cpu_fraction']
        for i,(phase,cost) in enumerate(zip(ledger['phases'],estimate['phases'])):
            seconds=cost['seconds'];cpu=q*sum(cost[key] for key in ['fixed_cpu_seconds','host_copy_seconds','host_page_seconds','pageable_d2h_cpu_seconds'])
            held=q*((cost['fixed_cpu_seconds']+(cost['host_copy_seconds'] if n*phase['traits']<=500 else 0))*serial+cost['host_page_serial_seconds'])
            resources={}
            for key,amount in [('cpu',cpu),('host_serial',held),('dram',cost['host_dram_bytes'])]+[(device+':'+d,phase[d+'_bytes']) for d in ['h2d','d2h']]:
                if amount:
                    if not seconds:raise ValueError('Positive setup work has zero service')
                    resources[key]=amount/seconds
            last=graph.add(phase['phase']+':'+str(i),seconds,[last],resources)
        unknown.update(estimate['unpriced_terms'])
        return last
    first=candidate['tiles'][0];ledger,estimate=setup_subset(first,'residual')
    last=append_setup(shared,last,first['device'],ledger,estimate,first['profile'])
    shared_cleanup=estimate['cleanup'];reports['shared']=estimate
    for tile in candidate['tiles']:
        device,profile=tile['device'],tile['profile'];q=profile['cpu_fraction']
        graph=ExecutionGraph();last=graph.add('begin')
        factor=factor_service(n,k,profile,library_arithmetic=library_arithmetic[device])
        for row in factor['phases']:
            seconds=row['seconds'];amounts=dict(cpu=row['cpu_seconds'],
                host_serial=serial*row['fixed_cpu_seconds'],dram=row['host_dram_bytes'])
            amounts.update({device+':'+d:row[d+'_bytes'] for d in ['h2d','d2h']})
            if not seconds and any(amounts.values()):raise ValueError('Positive factor work has zero service')
            last=graph.add('factor:'+row['phase'],seconds,[last],{key:value/seconds for key,value in amounts.items() if value})
        ledger,estimate=setup_subset(tile,'design')
        last=append_setup(graph,last,device,ledger,estimate,profile)
        pins=pinned_scan_work(n,shape['chunk_size'],k,shape['depth'],reduction='jagwas')
        pin_cpu=pins['allocation_pages']*_number('pin CPU/page',profile['pin_cpu_seconds_per_page'],zero=True)
        pin_driver=pins['allocation_pages']*_number('pin driver/page',profile['pin_driver_seconds_per_page'],zero=True)
        last=graph.add('pin_cpu',pin_cpu/q,[last],dict(cpu=q,host_serial=q*serial))
        graph.add('pin_driver',pin_driver,[last]);devices[device]=graph
        clean=estimate['cleanup'];seconds=clean['cpu_seconds']/q
        cleanup[device]=[dict(seconds=seconds,resources=dict(cpu=q,host_serial=clean['serial_cpu_seconds']/seconds if seconds else 0.))]
        reports[device]=dict(factor=factor,design=estimate,pins=pins)
        unknown.update(factor['unpriced_terms'])
    final=copy.deepcopy(finalize)
    if shared_cleanup['cpu_seconds']:
        q=first['profile']['cpu_fraction'];seconds=shared_cleanup['cpu_seconds']/q
        final.insert(0,dict(seconds=seconds,resources=dict(cpu=q,host_serial=shared_cleanup['serial_cpu_seconds']/seconds)))
    unknown.update(['Fresh pinned-allocation scenario; process cache history may reuse completed buffers',
        'Shared CPU count/basis/cast services and final commit are caller-supplied independent components',
        'Factor-to-design GPU stream overlap and host object release placement'])
    contract=dict(dimensions=[n,k,c,columns,shape['chunk_size'],shape['depth']],shared_graph=shared,
        device_graphs=devices,cleanup=cleanup,finalize=final,unpriced_terms=sorted(unknown))
    return dict(preparation=contract,components=reports,prediction_complete=False)


def attach_factor_calibration(profile,bank,capacity,*,device,factor_wait_cpu_fraction):
    """Attach complete fixed-reference observations without claiming transfer.

    Context, repetitions, work, rates and instruction families must agree.
    All signed reference partitions are retained in the observation identity;
    a negative aggregate is refused rather than clipped or discarded.
    """
    import hashlib
    import json
    import statistics
    _number('factor wait fraction',factor_wait_cpu_fraction,zero=True,fraction=True)
    from pathlib import Path
    expected_source=jagwas_factor_arithmetic_work(*FACTOR_REFERENCE)['source_sha256']
    if bank.get('source_sha256')!=expected_source:raise ValueError('Factor primitive executor source differs')
    method=jagwas_cutoff_method();names=FACTOR_PHASES_BY_METHOD[method]
    if bank.get('method','rounding')!=method:raise ValueError('Factor primitive bank is for another cutoff method')
    context=bank['context']
    for key in ['device','device_name','torch_version','cuda_version','affinity','torch_threads']:
        if context.get(key)!=capacity.get(key):raise ValueError('Factor/capacity context differs: '+key)
    if context['device']!=device or context.get('allow_tf32') is not False:
        raise ValueError('Factor calibration device or precision context differs')
    if context.get('boundary')!='single_worker_call_then_device_synchronize':
        raise ValueError('Factor worker measurement boundary differs')
    if bank.get('reference_shape')!=FACTOR_REFERENCE or bank.get('compute_dtype')!='float32' or bank.get('boundary')!='call_then_device_synchronize':
        raise ValueError('Fixed factor reference contract differs')
    records=bank['records'];expected={(name,r) for name in names for r in range(5)}
    if len(records)!=len(expected) or {(r['phase'],r['repeat']) for r in records}!=expected:
        raise ValueError('Five complete unique factor repeats required')
    for row in records:
        if row['calls']!=32:raise ValueError('Fixed factor call count differs')
        cpu=_number('factor observed CPU',row['cpu_seconds'],zero=True)
        wall=_number('factor observed wall',row['wall_seconds'])
        signed=row['non_cpu_seconds']
        if isinstance(signed,bool) or not isinstance(signed,(int,float)) or not math.isfinite(signed) or not math.isclose(cpu+signed,wall,rel_tol=1e-12,abs_tol=1e-15):
            raise ValueError('Factor signed partition does not conserve elapsed time')
    phases={name:{key:statistics.median(r[key] for r in records if r['phase']==name)
        for key in ['cpu_seconds','non_cpu_seconds']} for name in names}
    if phases!=bank.get('phases') or any(v<0 for phase in phases.values() for v in phase.values()):
        raise ValueError('Factor aggregate service differs or is negative')
    rows=capacity['rows'];expected={(kind,r) for kind in ['scalar','tensor'] for r in range(5)}
    if len(rows)!=10 or {(r['arithmetic'],r['repeat']) for r in rows}!=expected:
        raise ValueError('Five complete paired FP64 capacity repeats required')
    for kind,instruction in [('scalar',1),('tensor',8)]:
        selected=capacity['selections'][kind]
        if any(type(selected[key]) is not int or selected[key]<0 for key in ['instruction_flags','preference_mask','workspace_bytes']):
            raise ValueError('Invalid FP64 instruction/workspace metadata')
        if selected['instruction_flags']&0xff!=instruction or selected['preference_mask']&0xff!=instruction or not 0<=selected['workspace_bytes']<=64<<20:
            raise ValueError('FP64 instruction family or workspace contract differs')
        kernels=capacity['kernels'][kind]
        names=[row['name'] for row in kernels]
        tensor=any('tensorop_d884gemm' in name or 'tensorop_d1688gemm' in name for name in names)
        scalar=any('dgemm' in name and 'tensorop' not in name for name in names)
        if (kind=='scalar' and (tensor or not scalar)) or (kind=='tensor' and not tensor):
            raise ValueError('FP64 kernel census differs from instruction family')
        if _number('CPU-reference absolute error',capacity['max_abs_cpu_reference_error'][kind],zero=True)>1e-9:
            raise ValueError('Fixed FP64 primitive CPU reference check differs')
    for row in rows:
        if (row['dimension'],row['loops'],row['useful_flops'])!=(4096,16,2*4096**3):
            raise ValueError('Independent fixed FP64 work differs')
        seconds=_number('observed GPU seconds',row['gpu_seconds'])
        rate=_number('observed FP64 rate',row['flops_per_second'])
        if not math.isclose(rate,row['useful_flops']/seconds,rel_tol=1e-12):
            raise ValueError('FP64 rate does not conserve work')
    rates={field:statistics.median(r['flops_per_second'] for r in rows if r['arithmetic']==kind)
        for field,kind in [('fp64_flops_per_second','scalar'),('fp64_tensor_flops_per_second','tensor')]}
    if rates!=capacity.get('resources'):raise ValueError('FP64 summary differs from complete repeats')
    result=copy.deepcopy(profile)
    result['gpu_resources'].update(rates)
    result['jagwas_factor_context']=copy.deepcopy(context)
    keys=['reference_shape','compute_dtype','boundary','context','source_sha256']
    result['jagwas_factor_primitives']={key:copy.deepcopy(bank[key]) for key in keys}
    result['jagwas_factor_primitives']['method']=method
    result['jagwas_factor_primitives']['phases']=phases
    result['factor_wait_cpu_fraction']=factor_wait_cpu_fraction
    result['factor_observation_sha256']={name:hashlib.sha256(json.dumps(record,sort_keys=True,allow_nan=False).encode()).hexdigest()
        for name,record in [('fixed_phases',bank),('fp64_capacity',capacity)]}
    result['factor_calibration_transfer_qualified']=False
    return result
