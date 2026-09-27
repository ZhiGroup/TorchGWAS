"""NumPy's hugepage advice is off during a run and restored after."""
import threading

import numpy as np
import pytest

from torchgwas import host_pages
from torchgwas.host_pages import numpy_hugepage_advice, without_numpy_hugepages

pytestmark = pytest.mark.skipif(numpy_hugepage_advice() is None, reason='NumPy without the hugepage switch')


@pytest.fixture(autouse=True)
def advice_on(monkeypatch):
    monkeypatch.delenv('TORCHGWAS_NUMPY_HUGEPAGE', raising=False)
    setter = host_pages._setter()
    before = setter(True)
    yield
    setter(before)


def test_off_inside_and_restored_after():
    seen = []

    @without_numpy_hugepages
    def run():
        seen.append(numpy_hugepage_advice())
        np.empty(1 << 21)

    run()
    assert seen == [False] and numpy_hugepage_advice() is True


def test_restored_after_a_failure_and_reentrant():
    @without_numpy_hugepages
    def inner():
        assert numpy_hugepage_advice() is False
        raise RuntimeError('scan failed')

    @without_numpy_hugepages
    def outer():
        with pytest.raises(RuntimeError):
            inner()
        return numpy_hugepage_advice()

    assert outer() is False
    assert numpy_hugepage_advice() is True


def test_concurrent_runs_restore_only_when_the_last_ends():
    first_inside, release = threading.Event(), threading.Event()

    @without_numpy_hugepages
    def long_run():
        first_inside.set()
        release.wait(5)

    @without_numpy_hugepages
    def short_run():
        return numpy_hugepage_advice()

    worker = threading.Thread(target=long_run)
    worker.start()
    first_inside.wait(5)
    assert short_run() is False
    assert numpy_hugepage_advice() is False
    release.set()
    worker.join()
    assert numpy_hugepage_advice() is True


def test_environment_keeps_numpy_setting(monkeypatch):
    monkeypatch.setenv('TORCHGWAS_NUMPY_HUGEPAGE', '1')
    assert without_numpy_hugepages(numpy_hugepage_advice)() is True
