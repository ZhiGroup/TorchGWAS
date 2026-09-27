"""Autotuned chunk switching on BED, BGEN (CPU and GPU decode) and disk-backed inputs.

Each test runs the same input once with a fixed chunk and once with autotune
(trials forced on), requires the tuner to have switched sizes and committed,
and compares every statistic with the fixed run.
"""
import json
import sqlite3
import struct
import zlib

import numpy as np
import pytest
import torch

N, M, K = 96, 6144, 12
# start_chunk=16: these tests exercise probing and switching; the default
# start is the largest size (api.py), which probes nothing above itself.
OPTIONS = dict(chunk_sizes=[16, 32, 64], warmup_fraction=.02, trial_fraction=.2, min_job_seconds=0,
               start_chunk=16)
RING = dict(prefetch_chunks=2, reader_workers=4)


@pytest.fixture(autouse=True)
def native(monkeypatch):
    if not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    for key, value in dict(TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_SCAN_PROFILE='0').items():
        monkeypatch.setenv(key, value)


def dosages(seed=11):
    rng = np.random.default_rng(seed)
    return rng.integers(0, 3, (N, M)).astype(np.float32)


def phenotypes(tmp_path, calls):
    rng = np.random.default_rng(12)
    y = rng.normal(size=(N, K)).astype(np.float32)
    y[:, :2] += .5*calls[:, :2]
    np.save(tmp_path/'y.npy', y)
    np.save(tmp_path/'c.npy', rng.normal(size=(N, 2)).astype(np.float32))


def run(tmp_path, name, genotype, **kwargs):
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats import open_binary_sumstats
    out = tmp_path/name
    run_linear_gwas(genotype=genotype, phenotype=tmp_path/'y.npy', covariates=tmp_path/'c.npy',
                    compute_dtype='float32', device='cuda:0', output_dir=out, **kwargs)
    beta, t, _logp, _ = open_binary_sumstats(out/'sumstats')
    return json.loads((out/'run.json').read_text()), np.asarray(t, dtype=np.float64)


def check(tmp_path, genotype, path, **kwargs):
    _, fixed = run(tmp_path, 'fixed', genotype, chunk_size=32, **kwargs)
    meta, tuned = run(tmp_path, 'tuned', genotype, autotune=True, autotune_options=OPTIONS, **RING, **kwargs)
    chunk = meta['autotune']['chunk']
    assert meta['autotune']['control_path'] == path
    assert meta['autotune']['tuner_disabled_reason'] is None, meta['autotune']['tuner_disabled_reason']
    assert chunk['state'] == 'committed' and chunk['choice'] in (16, 32, 64)
    # Every size from the start up really ran (switching works on this path);
    # smaller sizes are never probed (model_autotune).
    assert {s['size'] for s in chunk['segments']} == {size for size in chunk['sizes'] if size >= chunk['initial']}
    assert len({s['size'] for s in chunk['segments']}) >= 2
    finite = np.isfinite(fixed) & np.isfinite(tuned)
    assert np.array_equal(np.isfinite(fixed), np.isfinite(tuned))
    np.testing.assert_allclose(tuned[finite], fixed[finite], rtol=2e-5, atol=2e-5)
    return meta


def test_bed_packed_path(tmp_path, monkeypatch):
    from test_statistics import _write_bed
    from torchgwas import scan_gpu
    if not scan_gpu.available(torch.device('cuda:0')):
        pytest.skip('native BED kernel not built for this GPU (chunk switching needs it)')
    calls = dosages()
    phenotypes(tmp_path, calls)
    path = _write_bed(tmp_path/'input', calls)
    check(tmp_path, str(path), 'packed')


def write_bgen(path, calls):
    """Layout 2, zlib, 8-bit probabilities; dosage d -> (p_AA, p_AB) exactly."""
    n, m = calls.shape
    samples = [f's{i}' for i in range(n)]
    sampleblock = (struct.pack('<II', 8+sum(2+len(x) for x in samples), n)
                   + b''.join(struct.pack('<H', len(x))+x.encode() for x in samples))
    data = bytearray(struct.pack('<IIII4sI', 20+len(sampleblock), 20, m, n, b'bgen', 1 | (2 << 2) | (1 << 31))
                     + sampleblock)
    code = {0: (255, 0), 1: (0, 255), 2: (0, 0)}
    rows = []
    for i in range(m):
        probs = b''.join(bytes(code[int(d)]) for d in calls[:, i])
        raw = struct.pack('<IHBB', n, 2, 2, 2)+bytes([2]*n)+bytes([0, 8])+probs
        compressed = zlib.compress(raw)

        def text(value, width='H'):
            return struct.pack('<'+width, len(value))+value.encode()
        record = (text(f'id{i}')+text(f'rs{i}')+text('1')+struct.pack('<IH', 100+i, 2)
                  + text('A', 'I')+text('G', 'I')+struct.pack('<II', len(compressed)+4, len(raw))+compressed)
        rows.append(('1', 100+i, f'rs{i}', 2, 'A', 'G', len(data), len(record)))
        data += record
    path.write_bytes(bytes(data))
    with sqlite3.connect(str(path)+'.bgi') as db:
        db.execute('CREATE TABLE Variant (chromosome TEXT,position INTEGER,rsid TEXT,number_of_alleles INTEGER,'
                   'allele1 TEXT,allele2 TEXT,file_start_position INTEGER,size_in_bytes INTEGER)')
        db.executemany('INSERT INTO Variant VALUES (?,?,?,?,?,?,?,?)', rows)


@pytest.mark.parametrize('backend', ['cpu', 'gpu'])
def test_bgen_cpu_and_gpu_decode(tmp_path, backend):
    calls = dosages()
    phenotypes(tmp_path, calls)
    path = tmp_path/'input.bgen'
    write_bgen(path, calls)
    if backend == 'gpu':
        from torchgwas.bgen import BgenGenotype
        probe = BgenGenotype(str(path), decode_backend='auto')
        if probe.resolve_decode_backend(torch.device('cuda:0')) != 'gpu':
            pytest.skip(f'GPU BGEN decoder unavailable: {probe.backend_reason}')
    meta = check(tmp_path, str(path), 'device' if backend == 'gpu' else 'range', bgen_decode_backend=backend)
    assert meta.get('decode_backend_used', backend) == backend


def test_disk_backed_generic_path(tmp_path):
    from torchgwas.io import DiskBackedGenotype
    calls = dosages()
    phenotypes(tmp_path, calls)
    np.save(tmp_path/'geno.npy', calls)
    genotype = DiskBackedGenotype(tmp_path/'geno.npy', np.array([f's{i}' for i in range(N)]),
                                  np.array([f'v{i}' for i in range(M)]))
    check(tmp_path, genotype, 'range')
