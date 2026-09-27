"""Source branch boundaries, separate prices and rejection of stale protocols."""
from copy import deepcopy
import numpy as np
import pytest
from torchgwas.numpy_nonzero_work import nonzero_protocol, nonzero_regime, validate_host_price_protocol
from torchgwas.significant_host_work import host_significant_selection_work, host_selection_service
from test_significant_host_model import bank
from test_window_model import configuration, evaluate
from test_trait_tiling_model import input_path


@pytest.mark.parametrize('cells,retained,regime',[(0,0,'empty'),(99,0,'empty'),
    (99,9,'sparse'),(99,10,'dense'),(100,10,'sparse'),(100,11,'dense'),
    (1<<20,104857,'sparse'),(1<<20,104858,'dense')])
def test_source_branch_boundary(cells,retained,regime):
    assert nonzero_regime(cells,retained)==regime


@pytest.mark.parametrize('cells,retained',[(-1,0),(10,-1),(10,11),(True,0),(10,1.5)])
def test_invalid_counts_rejected(cells,retained):
    with pytest.raises(ValueError,match='bounded retained'): nonzero_regime(cells,retained)


def test_mask_occupancy_selects_distinct_measured_services():
    prices=bank()['prices']
    for name,cost in [('empty',.001),('sparse',.01),('dense',.1)]:
        prices['flatnonzero_'+name]=dict(call_cpu_seconds=0.,unit_cpu_seconds=cost,dram_bytes_per_unit=2)
    for count,regime,cost in [(0,'empty',.001),(10,'sparse',.01),(11,'dense',.1)]:
        work=host_significant_selection_work(10,10,count)
        assert work['nonzero_primitive']=='flatnonzero_'+regime
        result=host_selection_service(work,prices,cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
        assert result[4]['seconds']==100*cost
    legacy=deepcopy(prices);legacy['flatnonzero_nonempty']=legacy.pop('flatnonzero_sparse')
    with pytest.raises(ValueError,match='flatnonzero_sparse'):
        host_selection_service(host_significant_selection_work(10,10,1),legacy,
            cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)


def test_protocol_is_version_bound(monkeypatch):
    prices=bank(); validate_host_price_protocol(prices)
    prices.pop('nonzero_protocol')
    with pytest.raises(ValueError,match='nonzero protocol'): validate_host_price_protocol(prices)
    prices=bank();prices['nonzero_protocol']['numpy_version']='different'
    with pytest.raises(ValueError,match='nonzero protocol'): validate_host_price_protocol(prices)
    monkeypatch.setattr(np,'__version__','99.0')
    with pytest.raises(ValueError,match='Unqualified NumPy'): nonzero_protocol()


def test_prepared_window_rejects_unsplit_legacy_price_bank(input_path):
    windows,options=configuration(input_path,reduction='significant')
    options['prices'].pop('nonzero_protocol')
    with pytest.raises(ValueError,match='nonzero protocol'): evaluate(windows,options)


def test_prepared_significant_window_discloses_unpriced_allocation_transfer(input_path):
    windows,options=configuration(input_path,reduction='significant')
    report=evaluate(windows,options)
    assert any('page residency' in term for term in report['unpriced_terms'])
    assert not report['prediction_complete'] and not report['selection_validated']


@pytest.mark.parametrize('return_beta',[False,True])
def test_flat_and_row_gathers_have_distinct_prices_and_source_traffic(return_beta):
    prices=bank()['prices'];work=host_significant_selection_work(13,19,37,return_beta=return_beta)
    prices['matrix_gather_flat']=dict(call_cpu_seconds=.002,unit_cpu_seconds=.003,dram_bytes_per_unit=16)
    prices['df_gather_row']=dict(call_cpu_seconds=.005,unit_cpu_seconds=.007,dram_bytes_per_unit=16)
    steps=host_selection_service(work,prices,cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
    payloads=1+int(return_beta)
    assert steps[5]['seconds']==payloads*(.002+37*.003)
    assert steps[7]['seconds']==.005+37*.007
    assert steps[5]['seconds']*steps[5]['resources']['dram']==pytest.approx(16*37*payloads)
    assert steps[7]['seconds']*steps[7]['resources']['dram']==pytest.approx(16*37)
    legacy=deepcopy(prices);legacy['matrix_gather']=legacy.pop('matrix_gather_flat')
    with pytest.raises(ValueError,match='matrix_gather_flat'):
        host_selection_service(work,legacy,cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)


@pytest.mark.parametrize('old_selector',['bounded_flat_v1','native_row_flat_v1'])
def test_two_dimensional_gather_prices_cannot_be_relabelled(old_selector):
    prices=bank();prices['host_selector']=old_selector
    with pytest.raises(ValueError,match='host selector'):validate_host_price_protocol(prices)
