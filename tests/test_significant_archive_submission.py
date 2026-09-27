"""Significant archives charge actual submitted bytes separately from extent."""
import io
import numpy as np
import pytest
from torchgwas.significant_host_work import indexed_part_work
from torchgwas.significant_host_model import _writer_service


class CountingArchive(io.BytesIO):
    def __init__(self):
        super().__init__()
        self.submitted = 0
    def write(self, data):
        self.submitted += len(data)
        return super().write(data)


@pytest.mark.parametrize('rows', [1, 37, 4096])
@pytest.mark.parametrize('store_beta', [False, True])
def test_npz_rewrites_charge_pagecache_and_dram_once_per_submission(rows, store_beta):
    part = indexed_part_work(rows, store_beta=store_beta)
    values = {r['field']: np.arange(rows, dtype=r['dtype']) for r in part['arrays']}
    with CountingArchive() as stream:
        np.savez(stream, **values)
        assert len(stream.getvalue()) == part['file_bytes']
        submitted = stream.submitted
        assert submitted == part['file_bytes'] + sum(r['local_header_bytes'] for r in part['arrays'])
        stream.seek(0)
        with np.load(stream, allow_pickle=False) as reopened:
            for name, value in values.items():
                np.testing.assert_array_equal(reopened[name], value)
    profile = dict(cpu_fraction=.5, shared_dram_bytes_per_second=1e9,
        writeback_service=dict(pagecache_seconds_per_byte=2e-9, storage_seconds_per_byte=3e-9),
        process_units=dict(numpy_copy_bytes=4e-9), fsync_seconds=.01)
    bank = {str(store_beta): dict(call_cpu_seconds=1e-5, byte_cpu_seconds=5e-9)}
    copy, transfer, commit = _writer_service(part, bank, profile, .3)
    expected_cpu = 1e-5 + 5e-9*part['array_payload_bytes'] + 2e-9*submitted + 4e-9*4*rows
    expected_dram = 5*part['array_payload_bytes'] + 8*rows + 2*(submitted-part['array_payload_bytes'])
    assert copy['seconds'] * copy['resources']['cpu'] == pytest.approx(expected_cpu)
    assert copy['seconds'] * copy['resources']['dram'] == pytest.approx(expected_dram)
    assert transfer['seconds'] * transfer['resources']['output'] == pytest.approx(part['file_bytes'])
    assert commit['seconds'] == .01


def test_no_archive_submission_for_empty_selection():
    assert _writer_service(indexed_part_work(0), {}, {}, 0.) == []
