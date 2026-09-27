"""Planner-cost reuse remains conservative, bounded and immutable."""
import time

import pytest

from torchgwas.productive_planning_cost import ProductivePlanningCostHistory


def options(tmp_path, **changes):
    return dict(cache_dir=str(tmp_path), max_age_seconds=3600.,
                publication_seconds=.05, **changes)


def forecasts():
    return dict(expected_cpu_seconds=.01, expected_wall_seconds=.02,
                publication_seconds=.03, reserve_seconds=.04)


def test_prior_productive_cost_only_raises_gate_and_does_not_renew_age(tmp_path):
    dependencies=dict(source='fixed', workload='jagwas')
    first=ProductivePlanningCostHistory(options(tmp_path),dependencies,max_steps=2)
    observed=time.time()-3
    first.observe(dict(evaluated=True,cpu_seconds=.3,wall_seconds=1.2),
                  observed_unix_seconds=observed)
    first.observe(dict(evaluated=False),observed_unix_seconds=time.time())
    stored=first.finish(successful=True)
    assert stored['publication']['status']=='stored'
    second=ProductivePlanningCostHistory(options(tmp_path),dependencies,max_steps=2)
    costs=second.costs(forecasts())
    assert second.lookup_report['hit']
    assert costs['expected_cpu_seconds']==.3
    assert costs['expected_wall_seconds']==1.2
    assert costs['publication_seconds']==pytest.approx(.08)
    assert costs['reserve_seconds']==.04
    assert second.prior['observed_unix_seconds']==observed
    assert second.prior['age_seconds']>=3
    assert second.finish(successful=True)['publication']['status']=='not_stored'
    third=ProductivePlanningCostHistory(options(tmp_path),dependencies,max_steps=2)
    assert third.lookup()['hit']
    assert third.prior['observed_unix_seconds']==observed
    assert third.prior['record_sha256']==second.prior['record_sha256']


def test_stale_and_changed_dependencies_cannot_supply_prior_cost(tmp_path):
    config=options(tmp_path)
    first=ProductivePlanningCostHistory(config,dict(source='old'),max_steps=1)
    first.observe(dict(evaluated=True,cpu_seconds=.3,wall_seconds=1.),
                  observed_unix_seconds=time.time()-10)
    first.finish(successful=True)
    changed=ProductivePlanningCostHistory(config,dict(source='new'),max_steps=1)
    assert not changed.lookup()['hit']
    assert changed.costs(forecasts())['expected_wall_seconds']==.02
    stale=ProductivePlanningCostHistory(dict(config,max_age_seconds=1.),dict(source='old'),max_steps=1)
    assert not stale.lookup()['hit']
    assert stale.lookup_report['reason']=='expired'
    assert stale.costs(forecasts())['expected_wall_seconds']==.02


def test_bounded_history_validation_and_preparation_charge(tmp_path):
    with pytest.raises(ValueError,match='step count'):
        ProductivePlanningCostHistory(options(tmp_path),dict(source='x'),max_steps=33)
    history=ProductivePlanningCostHistory(options(tmp_path),dict(source='x'),max_steps=1)
    with pytest.raises(ValueError,match='Positive'):
        history.observe(dict(evaluated=True,cpu_seconds=.1,wall_seconds=0.),
                        observed_unix_seconds=time.time())
    history.charge_preparation(.25)
    assert history.snapshot()['preparation_wall_seconds']==.25
    assert history.finish(successful=False)['publication']['status']=='not_stored'
