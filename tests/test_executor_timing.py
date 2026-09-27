"""A fake event clock checks real eager and worker-lazy factor boundaries."""
import threading
from types import SimpleNamespace

import pytest

from torchgwas.executor_timing import (
    SETUP_SCAN_WRITE_BOUNDARY, SETUP_SCAN_WRITE_METRIC, streaming_timing,
)


def test_same_endpoint_preserves_writer_interval_and_names_setup():
    row = streaming_timing(10., 17., 22.)
    assert row['scan_setup_seconds'] == 7.
    assert row['scan_and_write_seconds'] == 5.
    assert row[SETUP_SCAN_WRITE_METRIC] == 12.
    assert row['setup_scan_and_write_boundary'] == SETUP_SCAN_WRITE_BOUNDARY


@pytest.mark.parametrize('bounds', [(17., 10., 22.), (10., 22., 17.)])
def test_reversed_boundaries_are_not_reported_as_valid_timings(bounds):
    with pytest.raises(ValueError, match='ordered'):
        streaming_timing(*bounds)


@pytest.mark.parametrize('device_count', [1, 2, 3])
def test_joint_executor_clock_includes_every_factor_and_writer_publication(
        tmp_path, monkeypatch, device_count):
    from test_jagwas_variant_devices import fixture
    from torchgwas import api, sumstats_indexed
    from torchgwas.jagwas_projection import JagwasReduction

    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    path, _, y, c = fixture(tmp_path, 'pgen', missing=False)
    clock = [0.]
    lock = threading.Lock()
    def advance(seconds):
        with lock:
            clock[0] += seconds
    monkeypatch.setattr(api, 'time', SimpleNamespace(perf_counter=lambda: clock[0]))

    original_qc = api.prepare_inputs_for_prep
    def qc(*args, **kwargs):
        result = original_qc(*args, **kwargs)
        advance(11.)  # Input QC must be outside the new executor interval.
        return result
    monkeypatch.setattr(api, 'prepare_inputs_for_prep', qc)

    original_prepare = JagwasReduction.prepare
    def prepare(self, *args, **kwargs):
        result = original_prepare(self, *args, **kwargs)
        advance(7.)
        return result
    monkeypatch.setattr(JagwasReduction, 'prepare', prepare)

    original_writer = sumstats_indexed.write_indexed_sumstats
    def writer(*args, **kwargs):
        result = original_writer(*args, **kwargs)
        assert (tmp_path / 'out' / 'sumstats' / 'manifest.json').exists()
        advance(5.)
        return result
    monkeypatch.setattr(sumstats_indexed, 'write_indexed_sumstats', writer)

    original_phases = api._phase_breakdown
    def phases(*args, **kwargs):
        advance(23.)  # Later run bookkeeping must not move the endpoint.
        return original_phases(*args, **kwargs)
    monkeypatch.setattr(api, '_phase_breakdown', phases)

    result = api.run_linear_gwas(path, y, c,
        variant_devices=[f'cpu:{i}' for i in range(device_count)],
        output_dir=tmp_path / 'out', genotype_format='pgen',
        compute_dtype='float32', reduce='jagwas', chunk_size=4,
        reader_workers=3, sumstats_queue_depth=1, sumstats_fsync=True)
    row = result.run_metadata['sumstats_write']
    assert row[SETUP_SCAN_WRITE_METRIC] == 7. * device_count + 5.
    assert row['scan_setup_seconds'] == (7. if device_count == 1 else 0.)
    assert row['scan_and_write_seconds'] == (5. if device_count == 1 else 7. * device_count + 5.)
    assert row['scan_setup_seconds'] + row['scan_and_write_seconds'] == row[SETUP_SCAN_WRITE_METRIC]
    assert row['setup_scan_and_write_boundary'] == SETUP_SCAN_WRITE_BOUNDARY
    assert clock[0] == 11. + 7. * device_count + 5. + 23.
