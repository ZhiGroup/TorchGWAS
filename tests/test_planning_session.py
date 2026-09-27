"""Incremental work reuse must not reuse changed evidence or consume startup."""
import copy
from concurrent.futures import ThreadPoolExecutor
import threading
from unittest.mock import patch

import pytest

from torchgwas import planning_session as session
from torchgwas.planning_session import PlanningWorkCache,IncrementalPlanningBudget,cached_planning_work
from torchgwas.mechanistic_torch import torch_scan_work
from test_mechanistic_shapes import fixture,component


def compute():return {'rows':[1,2,3]}


def get(key=None,builder=compute):
    return cached_planning_work('test',{} if key is None else key,builder,implementation=(compute,))


def test_scoped_cache_survives_steps_but_not_changed_bindings_or_implementation():
    cache=PlanningWorkCache()
    with patch(__name__+'.compute',wraps=compute) as build:
        with cache.activate():
            first=get({'rate':1.},build);first['rows'].append(99)
        with cache.activate():
            assert get({'rate':1.},build)=={'rows':[1,2,3]}
            get({'rate':2.},build)
            assert build.call_count==2
        get({'rate':1.},build)
        assert build.call_count==3
    report=cache.snapshot();assert report['hits']==1 and report['entries']==2
    cache.close();assert cache.snapshot()['estimated_retained_bytes']==0
    with pytest.raises(ValueError),cache.activate():pass


def test_limits_failure_nested_scopes_and_non_json_inputs():
    outer=PlanningWorkCache(max_entries=2);inner=PlanningWorkCache()
    with outer.activate():
        for i in range(3):get({'i':i})
        assert outer.snapshot()['evictions']==1
        with pytest.raises(RuntimeError),inner.activate():
            get();raise RuntimeError('cancelled step')
        get({'i':2});assert outer.snapshot()['hits']==1
        get({'tuple':(1,2)});assert outer.snapshot()['bypasses']==1
    tiny=PlanningWorkCache(max_bytes=1024)
    with tiny.activate():
        get(builder=lambda:{'large':'x'*10000})
        assert tiny.snapshot()['entries']==0 and tiny.snapshot()['bypasses']==1
    bounded=PlanningWorkCache(max_bytes=4096)
    with bounded.activate():
        for i in range(3):get({'i':i},lambda:{'large':'x'*1800})
    assert bounded.snapshot()['evictions']>=1
    assert bounded.snapshot()['estimated_retained_bytes']<=4096


def test_failed_builder_is_not_cached_and_closed_scope_cannot_recompute():
    cache=PlanningWorkCache()
    with cache.activate():
        with pytest.raises(RuntimeError):get(builder=lambda:(_ for _ in ()).throw(RuntimeError()))
        assert cache.snapshot()['entries']==0
        cache.close()
        with pytest.raises(ValueError):get({'tuple':(1,2)})


@pytest.mark.parametrize('damage',['rate','geometry','host','shape','serial','dtype'])
def test_cached_tensor_components_match_uncached_with_changed_inputs(damage):
    import torch
    data,profile=fixture();before=copy.deepcopy((data,profile));cache=PlanningWorkCache()
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component) as build:
        with cache.activate():
            first=torch_scan_work(data,profile);assert build.call_count==2
            second=torch_scan_work(data,profile);assert second==first and build.call_count==2
            second['components'][8]['operations'][0]['kernel_service_seconds']=123.
            assert torch_scan_work(data,profile)==first
            changed=copy.deepcopy(profile)
            if damage=='rate':changed['gpu_resources']['hbm_bytes_per_second']/=2
            elif damage=='geometry':changed['kernel_geometry'][0]['kernels']=[{'identity':'changed'}]
            elif damage=='host':changed['host_primitives']={'changed':1.}
            elif damage=='shape':
                changed['chunk_markers']=2
                changed['kernel_geometry']=[r for r in changed['kernel_geometry'] if r['B']==2]
            elif damage=='serial':changed['cpu_fraction']=.5
            previous=torch.get_default_dtype()
            try:
                if damage=='dtype':torch.set_default_dtype(torch.float64)
                predicted=torch_scan_work(data,changed)
                # B=2 was already priced, but its per-chunk reader/control work
                # is rebuilt for the different schedule and agrees uncached.
                if damage!='shape':assert build.call_count>2
                with PlanningWorkCache().activate():assert torch_scan_work(data,changed)==predicted
            finally:torch.set_default_dtype(previous)
    assert (data,profile)==before


def test_implementation_replacement_misses_even_with_same_prices():
    data,profile=fixture();cache=PlanningWorkCache()
    with cache.activate():
        with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component) as a:
            torch_scan_work(data,profile);assert a.call_count==2
        with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component) as b:
            torch_scan_work(data,profile);assert b.call_count==2


class Clock:
    def __init__(self):self.wall=0.;self.cpu=0.
    def advance(self,wall,cpu):self.wall+=wall;self.cpu+=cpu;return 'calculated'


@pytest.fixture
def clock():
    value=Clock()
    with patch.object(session.time,'perf_counter',lambda:value.wall),patch.object(session.time,'thread_time',lambda:value.cpu):
        yield value


def step(budget,callback=lambda:None,**kwargs):
    values=dict(remaining_seconds=10.,expected_cpu_seconds=.01,expected_wall_seconds=.1)
    values.update(kwargs)
    return budget.run_step(callback,**values)


def test_no_planning_before_useful_output_or_after_forecast_payback_window(clock):
    budget=IncrementalPlanningBudget(max_cpu_seconds=1.)
    callback=lambda:pytest.fail('planning must not start')
    assert step(budget,callback)['reason']=='no_useful_output_yet'
    assert budget.start_after_first_output()==0.
    assert step(budget,callback,remaining_seconds=.1)['reason']=='insufficient_remaining_horizon'
    assert step(budget,callback,expected_gain_seconds=.1)['reason']=='forecast_gain_does_not_repay_step'
    assert step(budget,callback,expected_gain_seconds=-1.)['reason']=='forecast_gain_does_not_repay_step'
    assert budget.snapshot()['steps']==[]


def test_switching_must_fit_the_window_and_gain_cannot_exceed_remaining_work(clock):
    budget=IncrementalPlanningBudget(max_cpu_seconds=1.,max_window_seconds=1.)
    budget.start_after_first_output()
    assert step(budget,switching_seconds=.95)['reason']=='window_budget'
    with pytest.raises(ValueError):step(budget,remaining_seconds=1.,expected_gain_seconds=2.)


def test_actual_step_costs_are_charged_and_start_does_not_reset_deadline(clock):
    budget=IncrementalPlanningBudget(max_cpu_seconds=.1,max_window_seconds=2.,max_steps=2)
    budget.start_after_first_output()
    result=step(budget,lambda:clock.advance(.2,.02))
    assert result['evaluated'] and result['usable_for_decision'] and result['value']=='calculated'
    assert result['cpu_seconds']==.02 and result['wall_seconds']==.2
    clock.advance(1.75,0.)
    assert budget.start_after_first_output()==0.
    assert step(budget)['reason']=='window_budget'


def test_cumulative_planning_and_deferred_costs_gate_the_next_step(clock):
    budget=IncrementalPlanningBudget(max_cpu_seconds=1.);budget.start_after_first_output()
    first=step(budget,lambda:clock.advance(.2,.01))
    assert first['cumulative_planning_wall_seconds']==pytest.approx(.2)
    callback=lambda:pytest.fail('must check all costs before calculation')
    assert step(budget,callback,expected_gain_seconds=.31,publication_seconds=.02)['reason']=='forecast_gain_does_not_repay_step'
    assert step(budget,callback,remaining_seconds=.15,publication_seconds=.06)['reason']=='insufficient_remaining_horizon'
    second=step(budget,lambda:clock.advance(.3,.01),expected_gain_seconds=.5,publication_seconds=.02)
    assert not second['usable_for_decision'] and second['total_tuning_cost_seconds']==pytest.approx(.52)


@pytest.mark.parametrize('options',[dict(publication_seconds=-1.),dict(reserve_seconds=float('nan')),
    dict(publication_seconds=True),dict(publication_seconds=1e308,reserve_seconds=1e308)])
def test_invalid_deferred_tuning_costs_are_rejected(options):
    with pytest.raises(ValueError):step(IncrementalPlanningBudget(),**options)


@pytest.mark.parametrize('overrun',['cpu','window','horizon','gain','close'])
def test_cooperative_overrun_never_authorizes_a_late_switch(clock,overrun):
    budget=IncrementalPlanningBudget(max_cpu_seconds=.1,max_window_seconds=1.)
    budget.start_after_first_output()
    def work():
        if overrun=='close':budget.finish()
        return clock.advance(2. if overrun=='window' else .5,.2 if overrun=='cpu' else .01)
    options=dict(remaining_seconds=.3) if overrun=='horizon' else dict(expected_gain_seconds=.3) if overrun=='gain' else {}
    result=step(budget,work,**options)
    assert result['evaluated'] and not result['usable_for_decision']
    assert budget.snapshot()['steps'][0]['wall_seconds']>0


def test_exception_counts_and_finite_step_limit(clock):
    budget=IncrementalPlanningBudget(max_steps=1);budget.start_after_first_output()
    def bad():clock.advance(.2,.02);raise RuntimeError('failed calculation')
    with pytest.raises(RuntimeError):step(budget,bad)
    report=budget.snapshot()
    assert not report['step_in_flight'] and report['steps'][0]['error']=='RuntimeError'
    assert step(budget)['reason']=='step_budget'


def test_only_one_planning_step_may_be_in_flight():
    entered=threading.Event();release=threading.Event()
    budget=IncrementalPlanningBudget(max_cpu_seconds=1.,max_window_seconds=10.)
    budget.start_after_first_output()
    def work():entered.set();assert release.wait(5.);return 1
    with ThreadPoolExecutor(max_workers=1) as pool:
        future=pool.submit(step,budget,work)
        assert entered.wait(5.)
        try:assert step(budget)['reason']=='step_in_flight'
        finally:release.set()
        assert future.result()['evaluated']


@pytest.mark.parametrize('options',[{'max_steps':True},{'max_steps':0},{'max_cpu_seconds':float('nan')},
    {'max_cpu_seconds':0.},{'max_window_seconds':-1.}])
def test_invalid_budget_rejected(options):
    with pytest.raises(ValueError):IncrementalPlanningBudget(**options)
