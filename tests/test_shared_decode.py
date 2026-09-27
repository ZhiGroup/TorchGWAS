"""One decode pass must reach every tile scan intact, and never stall the others."""
from contextlib import contextmanager
import threading

import numpy as np
import pytest
import torch

from torchgwas.shared_decode import SharedDecodeHub


class FakeSource:
    allows_direct_native_fill = True
    native_dtype = np.int8
    decode_workers = 2

    def __init__(self, values, fail_at=None):
        self.values, self.fail_at = values, fail_at
        self.shape = values.shape
        self.fills = []
        self.lock = threading.Lock()

    @contextmanager
    def native_reader_session(self):
        def fill(start, end, out):
            if self.fail_at is not None and start >= self.fail_at:
                raise RuntimeError('decode failed')
            with self.lock:
                self.fills.append((start, end))
            np.copyto(out, self.values[:, start:end].T)
        yield fill


@pytest.fixture(autouse=True)
def cuda():
    if not torch.cuda.is_available():
        pytest.skip('pinned memory needs CUDA')


def consume(subscriber, seen, stop_after=None):
    try:
        for count, item in enumerate(subscriber):
            index, data, start, end = item[:4]
            if stop_after is not None and count == stop_after:
                break
            if len(item) > 5:  # GPU fan-out: wait for the root copy, then read the device tensor
                item[5].synchronize()
                seen.append((start, end, data.cpu().numpy()))
            else:
                seen.append((start, end, data.numpy().copy()))
            subscriber.release(index)
    finally:
        subscriber.close()


def values(m=101, n=16):
    return np.random.default_rng(3).integers(0, 3, (n, m)).astype(np.int8)


def test_every_subscriber_receives_each_chunk_from_one_decode():
    source = FakeSource(values())
    hub = SharedDecodeHub(source, 8, 3, 2, subscribers=3)
    seen = [[] for _ in range(3)]
    threads = [threading.Thread(target=consume, args=(hub.subscriber(i), seen[i])) for i in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert not any(t.is_alive() for t in threads)
    assert len(source.fills) == 13 == hub.chunks  # decoded once, not three times
    for rows in seen:
        assert [(s, e) for s, e, _ in rows] == [(s, min(s+8, 101)) for s in range(0, 101, 8)]
        for s, e, data in rows:
            np.testing.assert_array_equal(data, source.values[:, s:e].T)
    with pytest.raises(ValueError):
        hub.subscriber(0)


def test_gpu_fanout_copies_each_chunk_once_to_the_root_gpu():
    # One PCIe copy into a ring on cuda:0; every subscriber reads that ring.
    source = FakeSource(values(m=200))
    hub = SharedDecodeHub(source, 8, 3, 2, subscribers=2, fanout_device='cuda:0')
    seen = [[], []]
    threads = [threading.Thread(target=consume, args=(hub.subscriber(i), seen[i])) for i in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert not any(t.is_alive() for t in threads)
    assert hub.audit()['transfer'] == 'pcie_to_root_then_gpu_peer' and len(source.fills) == 25
    for rows in seen:
        assert len(rows) == 25
        for s, e, data in rows:
            np.testing.assert_array_equal(data, source.values[:, s:e].T)


def test_a_scan_that_stops_early_does_not_stall_the_others():
    source = FakeSource(values(m=400))
    hub = SharedDecodeHub(source, 8, 2, 2, subscribers=2)
    early, full = [], []
    threads = [threading.Thread(target=consume, args=(hub.subscriber(0), early, 3)),
               threading.Thread(target=consume, args=(hub.subscriber(1), full))]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert not any(t.is_alive() for t in threads)
    assert len(early) == 3 and len(full) == 50


def test_decode_errors_reach_every_scan_and_untaken_shares_are_released():
    # Subscribers must consume concurrently (as tile scans do); one that is
    # live but idle holds buffer shares and paces the decoder.
    hub = SharedDecodeHub(FakeSource(values(), fail_at=40), 8, 3, 2, subscribers=2)
    errors = [None, None]

    def run(i):
        try:
            consume(hub.subscriber(i), [])
        except RuntimeError as error:
            errors[i] = error
    threads = [threading.Thread(target=run, args=(i,)) for i in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert not any(t.is_alive() for t in threads)
    assert all(e is not None and 'decode failed' in str(e) for e in errors)
    # A tile that will never start holds buffer shares until abandoned; after
    # that the remaining tile runs to completion and the decoder shuts down.
    idle = SharedDecodeHub(FakeSource(values()), 8, 3, 2, subscribers=2)
    started = idle.subscriber(0)
    idle.abandon_untaken()
    rows = []
    consume(started, rows)
    assert len(rows) == 13 and idle._loader_closed


def test_pcie_root_port_reads_the_sysfs_tree(tmp_path):
    from torchgwas.shared_decode import pcie_root_port
    devices = tmp_path/'sys'/'devices'
    bus = tmp_path/'sys'/'bus'/'pci'/'devices'
    bus.mkdir(parents=True)
    # Two GPUs behind one switch under root port 0000:00:01.1 (lab-a100 GPUs 0 and 1).
    for gpu, port in [('0000:07:00.0', '0000:04:00.0'), ('0000:0a:00.0', '0000:04:10.0')]:
        leaf = devices/'pci0000:00'/'0000:00:01.1'/'0000:01:00.0'/port/gpu
        leaf.mkdir(parents=True)
        (bus/gpu).symlink_to(leaf)
    other = devices/'pci0000:40'/'0000:40:01.1'/'0000:47:00.0'
    other.mkdir(parents=True)
    (bus/'0000:47:00.0').symlink_to(other)
    root = [pcie_root_port(b, sysfs=str(bus)) for b in ('00000000:07:00.0', '00000000:0A:00.0', '00000000:47:00.0')]
    assert root == ['pci0000:00/0000:00:01.1', 'pci0000:00/0000:00:01.1', 'pci0000:40/0000:40:01.1']
    assert pcie_root_port('00000000:99:00.0', sysfs=str(bus)) is None  # unknown device
    assert pcie_root_port('not-a-bus-id', sysfs=str(bus)) is None
