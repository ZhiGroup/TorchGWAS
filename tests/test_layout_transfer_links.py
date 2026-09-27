"""Nested transfer bottlenecks use one shared resource load per link."""
import pytest

from torchgwas.layout_transfer_links import transfer_link_loads


def test_nested_links_take_maximum_concurrent_constraint():
    work = {'cuda:0': 100, 'cuda:1': 200, 'cuda:2': 300}
    report = transfer_link_loads(work, [
        dict(devices=['cuda:0', 'cuda:1', 'cuda:2'],
             h2d_bytes_per_second=300., d2h_bytes_per_second=600.),
        dict(devices=['cuda:1', 'cuda:2'],
             h2d_bytes_per_second=100., d2h_bytes_per_second=200.)],
        direction='h2d')
    assert [row['bytes'] for row in report['links']] == [600, 500]
    assert report['floor_seconds'] == 5.


@pytest.mark.parametrize('links', [
    [dict(devices=['cuda:1'], h2d_bytes_per_second=1., d2h_bytes_per_second=1.)],
    [dict(devices=['cuda:0', 'cuda:0'], h2d_bytes_per_second=1., d2h_bytes_per_second=1.)],
    [dict(devices=['cuda:0'], h2d_bytes_per_second=0., d2h_bytes_per_second=1.)],
    [dict(devices=['cuda:0'], h2d_bytes_per_second=1.)],
])
def test_link_binding_rejects_unknown_or_invalid_resources(links):
    with pytest.raises(ValueError):
        transfer_link_loads({'cuda:0': 100}, links, direction='d2h')
