"""Trait-sharded reduced scans stop every device once one of them fails."""
import threading
import time

import numpy as np
import pytest

from torchgwas.api import _trait_blocked_reduced_chunks


class _KeepLatest:
    def merge(self, previous, incoming, offset):
        return incoming


def test_a_failed_device_stops_the_others():
    chunks_per_block = 100
    seen = {'b': 0}
    started, failed = threading.Event(), threading.Event()
    closed = []

    def scan_once(offset, width, device):
        def iterate():
            try:
                if device == 'a':
                    # Fail once b is mid-block, so b has a scan to abandon.
                    started.wait(5)
                    failed.set()
                    raise RuntimeError('device a failed')
                started.set()
                failed.wait(5)
                for start in range(chunks_per_block):
                    seen['b'] += 1
                    time.sleep(0.001)
                    yield (start, start + 1, np.zeros((1, 1)), np.zeros((1, 1)), None, np.zeros((1, 1), np.int64))
            finally:
                closed.append(device)
        return iterate()

    with pytest.raises(RuntimeError, match='device a failed'):
        list(_trait_blocked_reduced_chunks(scan_once, _KeepLatest(), 4, 1, devices=['a', 'b']))
    # Device b owns two blocks of 100 chunks; it stops within a few chunks.
    assert seen['b'] < 10
    assert closed.count('b') == 1
