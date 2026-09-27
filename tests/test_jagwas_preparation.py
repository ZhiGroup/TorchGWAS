"""Preparation ordering and resource conservation, using synthetic primitives."""
import copy
import pytest
from torchgwas.jagwas_preparation import factor_service,factor_phase_work,build_jagwas_preparation,FACTOR_PHASES,FACTOR_REFERENCE
from torchgwas.jagwas_candidate import jagwas_candidate_runtime
from torchgwas.reduction_tensor_work import jagwas_factor_arithmetic_work
from test_jagwas_actual_candidate import input_path,actual_candidate,writer_prices


def profile(p):
    p=copy.deepcopy(p)
    context=dict(device_name='synthetic resource-control context')
    p.update(jagwas_factor_context=context,factor_wait_cpu_fraction=1.,
        pin_cpu_seconds_per_page=1e-7,pin_driver_seconds_per_page=2e-8)
    p['jagwas_factor_primitives']=dict(reference_shape=FACTOR_REFERENCE,compute_dtype='float32',
        boundary='call_then_device_synchronize',context=context,
        source_sha256=jagwas_factor_arithmetic_work(65,32)['source_sha256'],method='eigen',
        phases={name:dict(cpu_seconds=1e-6,non_cpu_seconds=2e-6) for name in FACTOR_PHASES})
    p['setup_primitives']={phase:dict(reference_shape=[32,1,2],cpu_seconds=1e-6,non_cpu_seconds=2e-6)
        for phase in ['residual_common','residual_block','design_common','design_block']}
    return p


def choice(path,count=2):
    c=actual_candidate(path,256,count)
    for tile in c['tiles']:tile['profile']=profile(tile['profile'])
    return c


def build(c,**kwargs):
    options=dict(host_serial_fraction=.5,library_arithmetic={d:'scalar' for d in c['devices']},
        shared_cpu_steps=[dict(seconds=.001,resources=dict(cpu=1.,host_serial=.5))],
        finalize=[dict(seconds=.001,resources=dict(cpu=1.,host_serial=.5))])
    options.update(kwargs)
    return build_jagwas_preparation(c,**options)


def test_reference_is_exact_and_extra_work_is_not_a_shape_timing_table(input_path):
    p=choice(input_path,1)['tiles'][0]['profile'];before=copy.deepcopy(p)
    tiny=factor_service(65,32,p,library_arithmetic='scalar')
    assert tiny['seconds']==pytest.approx(3e-6*len(FACTOR_PHASES))
    assert all(r['extra_gpu_seconds']==r['extra_transfer_seconds']==0 for r in tiny['phases'])
    scalar=factor_service(2049,512,p,library_arithmetic='scalar')
    tensor=factor_service(2049,512,p,library_arithmetic='tensor')
    assert scalar['seconds']>tensor['seconds']>tiny['seconds']
    assert sum(r['h2d_bytes'] for r in scalar['phases'])==4*2049*512
    # The eigen factor's one read of R's K-value FP64 spectrum.
    assert sum(r['d2h_bytes'] for r in scalar['phases'])==8*512
    assert [r['phase'] for r in scalar['phases']]==list(FACTOR_PHASES)
    assert p==before and not scalar['prediction_complete']
    p['factor_wait_cpu_fraction']=0.
    sleeping=factor_service(2049,512,p,library_arithmetic='scalar')
    assert sleeping['seconds']==scalar['seconds']
    assert sum(r['cpu_seconds'] for r in sleeping['phases'])<sum(r['cpu_seconds'] for r in scalar['phases'])


@pytest.mark.parametrize('count',[1,2])
def test_source_preparation_has_one_residual_pass_and_per_device_factor(input_path,count):
    c=choice(input_path,count);before=copy.deepcopy(c);built=build(c);prep=built['preparation']
    assert c==before
    assert len([n for n in prep['shared_graph'].nodes if n.startswith('residual_block:')])==1
    assert all(not n.startswith('design') for n in prep['shared_graph'].nodes)
    for device,graph in prep['device_graphs'].items():
        assert not any(n.startswith('residual') for n in graph.nodes)
        assert len([n for n in graph.nodes if n.startswith('factor:')])==len(FACTOR_PHASES)
        assert graph.nodes['design_common:0'][1]==('factor:'+FACTOR_PHASES[-1],)
        assert graph.nodes['pin_cpu'][1]==('design_block:1',)
        pins=built['components'][device]['pins'];assert pins['allocation_count']==12
    graph=jagwas_candidate_runtime(c,writer_prices(),preparation=prep,occupancy='dense',host_serial_fraction=.5,return_graph=True)
    solved=graph.solve()
    for i in range(count):
        assert solved['start'][f'tile:{i}:prepare:factor:upload']>=solved['end']['shared_prepare:complete']
        assert solved['start'][f'tile:{i}:submit_decode:0']>=solved['end'][f'tile:{i}:prepare:complete']
    for resource,capacity in graph.capacities.items():
        amount=sum(graph.nodes[n][0]*d.get(resource,0.) for n,d in graph.demands.items())
        assert solved['seconds']+1e-10>=amount/capacity
    uploads=sum(graph.nodes[n][0]*sum(d.get(device+':h2d',0.) for device in c['devices'])
                for n,d in graph.demands.items() if 'prepare:' in n)
    assert uploads==pytest.approx(4*2049*(512+2)+count*4*2049*(2*512+2))


@pytest.mark.parametrize('fault',['reference','boundary','context','source','missing_phase','negative','zero','fp64','wait','library','small'])
def test_invalid_factor_pricing_refused(input_path,fault):
    p=choice(input_path,1)['tiles'][0]['profile'];bank=p['jagwas_factor_primitives'];n,k=2049,512;scenario='scalar'
    if fault=='reference':bank['reference_shape']=[64,32]
    if fault=='boundary':bank['boundary']='scan_time_fit'
    if fault=='context':bank['context']={}
    if fault=='source':bank['source_sha256']={}
    if fault=='missing_phase':bank['phases'].pop('qr')
    if fault=='negative':bank['phases']['qr']['cpu_seconds']=-1
    if fault=='zero':bank['phases']['qr']=dict(cpu_seconds=0.,non_cpu_seconds=0.)
    if fault=='fp64':p['gpu_resources']['fp64_flops_per_second']=0.
    if fault=='wait':p['factor_wait_cpu_fraction']=1.1
    if fault=='library':scenario='automatic'
    if fault=='small':k=7
    with pytest.raises(ValueError):factor_service(n,k,p,library_arithmetic=scenario)


def test_required_shared_and_final_services_cannot_be_omitted(input_path):
    c=choice(input_path)
    for kwargs in [dict(shared_cpu_steps=[]),dict(finalize=[]),dict(library_arithmetic={'cuda:0':'scalar'})]:
        with pytest.raises(ValueError):build(c,**kwargs)
