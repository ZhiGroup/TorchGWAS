import math
import pytest
from torchgwas.primitive_prices import cpu_service_samples, mean_cpu_service_prices


def test_additive_cpu_keeps_expensive_tail_and_conserves_work():
    row = cpu_service_samples([1., 1., 10.], [.1, .2, .3])
    assert row['cpu_seconds_per_call'] == pytest.approx(3.8)
    assert row['median_cpu_seconds_per_call'] == pytest.approx(.8)
    assert row['sample_count'] == 3
    assert row['raw_cpu_sum_seconds'] == 12.
    assert 3 * (row['cpu_seconds_per_call'] + row['timer_overhead_seconds']) == pytest.approx(12.)


def test_below_timer_floor_has_zero_service():
    assert cpu_service_samples([.1], [.2])['cpu_seconds_per_call'] == 0.


@pytest.mark.parametrize('values,timers', [([], [0.]), ([0.], []), ([math.nan], [0.]),
                        ([math.inf], [0.]), ([-1.], [0.]), ([True], [0.]), ([0.], ['1'])])
def test_invalid_samples_rejected(values, timers):
    with pytest.raises(ValueError):
        cpu_service_samples(values, timers)


def rows():
    return [dict(cpu_service_samples([1., 1., 10.], [.2]), primitive=name, device=1, repeat=repeat)
            for name in ['copy', 'submit'] for repeat in range(3)]


def test_validated_repeat_means_keep_call_tail():
    assert mean_cpu_service_prices(rows(), device=1) == {'copy': 3.8, 'submit': 3.8}


@pytest.mark.parametrize('field,value', [('cpu_seconds_per_call', .8), ('raw_cpu_seconds_per_call', 1.),
                         ('raw_cpu_sum_seconds', 3.), ('sample_count', True), ('timer_overhead_seconds', float('nan'))])
def test_loader_rejects_inconsistent_mean_or_median_substitution(field, value):
    observations = rows()
    observations[0][field] = value
    with pytest.raises(ValueError):
        mean_cpu_service_prices(observations, device=1)


def test_missing_or_duplicate_context_rejected():
    observations = rows()
    for bad, device in [(observations, 2), (observations[:-1], 1), (observations + observations[:1], 1)]:
        with pytest.raises(ValueError):
            mean_cpu_service_prices(bad, device=device)
