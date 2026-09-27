import copy
import pytest
from torchgwas.pinned_work import fresh_pin_prices


def probe():
    runs=[]
    for repeat in range(5):
        rows=[]
        for index,size in enumerate([4096]*4+[1<<20]*4+[64<<20]*4):
            pages=size//4096;cpu=(repeat+1)*pages*1e-6
            rows.append(dict(index=index,bytes=size,pages=pages,pointer=100+index,
                             cpu_seconds=cpu,wall_seconds=cpu+pages*1e-7))
        runs.append(dict(repeat=repeat,rows=rows,torch_version='test',affinity=[1,2],device='1',
                         torch_threads=4,all_allocations_retained=True))
    return dict(runs=runs)


def price(value):
    return fresh_pin_prices(value,torch_version='test',cpu_affinity=[1,2],device=1)


def test_complete_process_mean_prices_preserve_signed_controls():
    value=probe();before=copy.deepcopy(value)
    result=price(value)
    assert value==before
    assert result['pin_cpu_seconds_per_page']==pytest.approx(3e-6)
    assert result['pin_driver_seconds_per_page']==pytest.approx(1e-7)
    for run in value['runs']:
        for row in run['rows']:row['wall_seconds']=row['cpu_seconds']-row['pages']*1e-7
    result=price(value)
    assert result['pin_driver_seconds_per_page']==0
    assert result['signed_non_cpu_seconds_per_page']==pytest.approx(-1e-7)


@pytest.mark.parametrize('mutate',[
    lambda p:p['runs'].pop(1),
    lambda p:p['runs'][0].update(device='2'),
    lambda p:p['runs'][0].update(all_allocations_retained=False),
    lambda p:p['runs'][0]['rows'][0].update(pointer=101),
    lambda p:p['runs'][0]['rows'][0].update(pages=2),
    lambda p:p['runs'][0]['rows'][0].update(cpu_seconds=float('nan')),
])
def test_context_and_raw_observation_guards(mutate):
    value=probe()
    mutate(value)
    with pytest.raises(ValueError):price(value)
