"""Device candidates share native setup, exact status/selection and one writer."""
import copy
from unittest.mock import patch
import pytest
from test_trait_tiling_model import candidate,input_path
from test_mechanistic_shapes import component
from test_native_status_work import status_prices
from test_significant_host_model import bank
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.device_significance_service import device_selection_host_primitives
from torchgwas.reduced_output_work import significant_output_work
from torchgwas.significant_device_model import significant_device_runtime


def fixture(path,occupancy='dense',**kwargs):
    c=candidate(path,block_bytes=None,**kwargs);services=[]
    for tile in c['tiles']:
        data,p=tile['data'],tile['profile'];data['phenotype_complete']=True
        p.update(reduction='device_significant',result_ownership='owned',device_status_service=status_prices())
        p.pop('result_finish_service');n,m,k=data['samples'],data['markers'],data['traits_analyzed']
        rows=[]
        for start in range(0,m,p['chunk_markers']):
            b=min(m-start,p['chunk_markers'])
            ledger=significant_output_work(n,b,k,b,backend='device')
            counts=[0 if occupancy=='empty' else 1 if occupancy=='single' else x['cells'] for x in ledger['blocks']]
            work=device_significant_tensor_work(n,b,k,counts)
            host={name:({'before_cpu_seconds':1e-6,'after_cpu_seconds':1e-6}
                if name.startswith(('copy_','nonzero_')) else {'cpu_seconds':1e-6})
                for name in device_selection_host_primitives(work)['counts']}
            operations={str(i):dict(op=s['op'],seconds=1e-6+s['gpu_logical_bytes']/1e12)
                for i,s in enumerate(work['steps']) if not(s['alias_only'] or s['allocation_only'] or s['host_copy'] or s['op']=='aten.nonzero.default')}
            nonzero=[dict(count=dict(seconds=1e-6),flagged_select=dict(seconds=2e-6),
                **({'coordinate_scatter':dict(seconds=1e-6)} if count else {})) for count in counts]
            rows.append(dict(retained_per_block=counts,host_prices=host,operation_services=operations,nonzero_services=nonzero,
                transfer_prices={name:dict(latency_seconds=1e-6,bytes_per_second=1e8,resources=['d2h']) for name in ['count','payload']},
                yield_cpu_seconds=1e-6,wait_cpu_fraction=1.))
        services.append(rows)
    prices=bank()
    prices['emit']={name:dict(call_cpu_seconds=1e-6,row_cpu_seconds=1e-9) for name in ['index_cast','inplace_index_add']}
    prices['critical']={d:dict(round_call_cpu_seconds=1e-6,round_row_cpu_seconds=1e-9,
        copy=dict(before_cpu_seconds=1e-6,after_cpu_seconds=1e-6),
        transfer=dict(latency_seconds=1e-6,bytes_per_second=1e8),wait_cpu_fraction=1.) for d in c['devices']}
    return c,prices,services


def run(c,prices,services,**kwargs):
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return significant_device_runtime(c,prices,selection_services=services,host_serial_fraction=.5,**kwargs)


@pytest.mark.parametrize('occupancy',['empty','single','dense'])
@pytest.mark.parametrize('count',[1,2])
def test_complete_device_candidate_conserves_payloads_and_input_passes(input_path,occupancy,count):
    c,prices,services=fixture(input_path,occupancy,count=count)
    original=copy.deepcopy((c,prices,services))
    report=run(c,prices,services)
    selected=0 if occupancy=='empty' else 9 if occupancy=='single' else 50
    assert sum(r['selected_payload_d2h_bytes'] for r in report['tiles'])==28*selected
    assert sum(r['status_d2h_bytes'] for r in report['tiles'])==30
    assert sum(r['nonzero_count_d2h_bytes'] for r in report['tiles'])==36
    assert report['parts']==(0 if occupancy=='empty' else 9)
    assert (c,prices,services)==original
    assert not report['prediction_complete'] and not report['selection_validated']
    assert report['resource_balance']['resource_lower_bound_seconds']<=report['estimated_tile_seconds']


def test_shared_transfer_work_includes_status_count_payload_and_setup(input_path):
    c,prices,services=fixture(input_path)
    c['shared_links']=[dict(devices=c['devices'],h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8)]
    c['shared_storage_bytes_per_second']=1e8
    report=run(c,prices,services)
    assert report['resource_balance']['resources']['link:0:d2h']['work']==pytest.approx(4*32*5+30+36+20*50)
    assert report['resource_balance']['resources']['storage']['work']>report['indexed_part_bytes']


def test_critical_upload_and_emit_rebase_follow_source_order(input_path):
    c,prices,services=fixture(input_path)
    graph=run(c,prices,services,return_graph=True);solved=graph.solve()
    assert solved==graph._solve_shared_python()
    for tile in range(3):
        base=f'tile:{tile}:'
        assert solved['start'][base+'prepare:critical:round']>=solved['end'][base+'prepare:phase:1']
        assert solved['start'][base+'prepare:phase:2']>=solved['end'][base+'prepare:critical:ready']
        for chunk in range(3):
            part=base+f'significant:{chunk}:'
            assert solved['start'][part+'0:put:0']>=solved['end'][part+'device:emit:0:inplace_index_add']
    assert solved['start']['tile:2:prepare:critical:round']>=solved['end']['tile:0:significant:producer_complete']
    intervals=sorted((solved['start'][n],solved['end'][n]) for n in graph.nodes if n.endswith(':write:1'))
    assert len(intervals)==9
    assert all(a[1]<=b[0] for a,b in zip(intervals,intervals[1:]))


@pytest.mark.parametrize('policy',['held-first','held-last'])
def test_declared_host_serial_policies_preserve_full_candidate(policy,input_path):
    c,prices,services=fixture(input_path,'single')
    report=run(c,prices,services,host_serial_policy=policy)
    assert report['parts']==9 and report['estimated_tile_seconds']>0


def test_beta_storage_changes_archive_bytes_but_not_selected_transport(input_path):
    c,prices,services=fixture(input_path)
    before=run(c,prices,services);c['output']['store_beta']=False
    after=run(c,prices,services)
    assert before['indexed_part_bytes']>after['indexed_part_bytes']
    assert [r['selected_payload_d2h_bytes'] for r in before['tiles']]==[r['selected_payload_d2h_bytes'] for r in after['tiles']]


def test_missing_extents_and_silent_service_overrides_are_refused(input_path):
    c,prices,services=fixture(input_path)
    with pytest.raises(ValueError,match='per tile'):run(c,prices,services[:1])
    with pytest.raises(ValueError,match='bounded graph'):run(c,prices,services,max_source_chunks=1)
    with pytest.raises(ValueError,match='bounded selection'):run(c,prices,services,max_selection_blocks=1)
    services[0][0]['cpu_fraction']=.5
    with pytest.raises(ValueError,match='independent primitive'):run(c,prices,services)


def test_emit_rejects_old_allocating_add_price(input_path):
    c,prices,services=fixture(input_path)
    prices['emit']['index_add']=prices['emit'].pop('inplace_index_add')
    with pytest.raises(ValueError,match='inplace_index_add'):run(c,prices,services)
