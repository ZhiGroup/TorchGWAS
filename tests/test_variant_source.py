"""Streaming stores record their input and variant indices, not the variant IDs."""
import json

import numpy as np
import pytest

from torchgwas.api import run_linear_gwas
from torchgwas.variant_source import store_variant_ids, store_variants, variant_digest


def write_input(tmp_path, n=64, m=23, k=3):
    from test_pgen_native_reader import write_pgen
    rng = np.random.default_rng(1207)
    calls = rng.integers(0, 3, size=(n, m)).astype(np.uint8)
    path = tmp_path / 'input.pgen'
    write_pgen(path, calls.T)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'1\t{100 + i}\trs{9000 + i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    np.save(tmp_path / 'y.npy', rng.normal(size=(n, k)).astype(np.float32))
    return path, [f'rs{9000 + i}' for i in range(m)]


@pytest.fixture(autouse=True)
def native(monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.setenv('TORCHGWAS_PGEN_PACKED', '0')


OPTIONS = dict(genotype_format='pgen', pgen_mode='hardcall', compute_dtype='float32', device='cpu',
               chunk_size=4, reader_workers=2, prefetch_chunks=2)


def test_dense_store_records_its_input_and_resolves_the_ids(tmp_path):
    path, ids = write_input(tmp_path)
    run_linear_gwas(path, tmp_path / 'y.npy', output_dir=tmp_path / 'dense', variant_range=(2, 19), **OPTIONS)
    store = tmp_path / 'dense' / 'sumstats'
    manifest = json.loads((store / 'manifest.json').read_text())
    source = manifest['variant_source']
    assert not (store / 'variant_ids.txt').exists()
    assert source['format'] == 'pgen' and source['variant_offset'] == 2 and source['n_variants'] == 17
    assert source['genotype'] == str(path.resolve())
    np.testing.assert_array_equal(store_variant_ids(store), ids[2:19])
    _, metadata = store_variants(store)
    assert isinstance(metadata, dict)


def test_indexed_store_records_its_input_and_resolves_the_ids(tmp_path):
    path, ids = write_input(tmp_path)
    run_linear_gwas(path, tmp_path / 'y.npy', output_dir=tmp_path / 'pairs', reduce='significant',
                    significance_threshold=1.0, **OPTIONS)
    store = tmp_path / 'pairs' / 'sumstats'
    manifest = json.loads((store / 'manifest.json').read_text())
    assert 'variant_ids' not in manifest and not (store / 'variant_ids.npy').exists()
    assert manifest['variant_source']['n_variants'] == len(ids)
    np.testing.assert_array_equal(store_variant_ids(store), ids)


def test_a_changed_variant_list_is_refused(tmp_path):
    path, ids = write_input(tmp_path)
    run_linear_gwas(path, tmp_path / 'y.npy', output_dir=tmp_path / 'dense', **OPTIONS)
    pvar = path.with_suffix('.pvar')
    pvar.write_text(pvar.read_text().replace('rs9005', 'rs1'))
    with pytest.raises(ValueError, match='does not match'):
        store_variant_ids(tmp_path / 'dense' / 'sumstats')


def test_ids_can_still_be_embedded(tmp_path):
    path, ids = write_input(tmp_path)
    run_linear_gwas(path, tmp_path / 'y.npy', output_dir=tmp_path / 'dense', sumstats_variant_ids=True, **OPTIONS)
    run_linear_gwas(path, tmp_path / 'y.npy', output_dir=tmp_path / 'pairs', reduce='significant',
                    significance_threshold=1.0, sumstats_variant_ids=True, **OPTIONS)
    assert (tmp_path / 'dense' / 'sumstats' / 'variant_ids.txt').read_text().split() == ids
    assert np.load(tmp_path / 'pairs' / 'sumstats' / 'variant_ids.npy').tolist() == ids
    # Embedded IDs are read without opening the genotype.
    path.with_suffix('.pvar').unlink()
    np.testing.assert_array_equal(store_variant_ids(tmp_path / 'dense' / 'sumstats'), ids)


def test_digest_is_order_sensitive_and_cheap_at_scale():
    ids = np.asarray([f'rs{i}' for i in range(100_000)])
    assert variant_digest(ids) == variant_digest(ids.copy())
    assert variant_digest(ids) != variant_digest(ids[::-1])
    assert variant_digest(ids) != variant_digest(ids[:-1])
