from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import special

from torchgwas.sumstats import (
    BinarySumstatsWriter,
    open_binary_sumstats,
    read_manifest,
)
from torchgwas.tails import upper_tail_log10_from_t


def _write_store(tmp_path, n_variants=57, n_traits=5, block_bytes=1024, seed=0):
    rng = np.random.default_rng(seed)
    beta = rng.normal(size=(n_variants, n_traits)).astype(np.float32)
    tstat = rng.normal(scale=3.0, size=(n_variants, n_traits)).astype(np.float32)
    beta[3] = np.nan
    tstat[3] = np.nan
    directory = tmp_path / "sumstats"
    writer = BinarySumstatsWriter(
        directory=directory,
        n_variants=n_variants,
        trait_names=[f"trait_{i}" for i in range(n_traits)],
        n_samples=1000,
        df=990,
        block_bytes=block_bytes,
        queue_depth=2,
    )
    step = 7
    for start in range(0, n_variants, step):
        end = min(start + step, n_variants)
        writer.write_chunk(start, end, beta[start:end], tstat[start:end])
    summary = writer.close()
    return directory, beta, tstat, summary


def test_roundtrip_is_bit_exact(tmp_path):
    directory, beta, tstat, summary = _write_store(tmp_path)
    stored_beta, stored_t, stored_logp, manifest = open_binary_sumstats(directory)
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)
    np.testing.assert_allclose(
        np.asarray(stored_logp),
        -np.log10(2.0 * special.stdtr(990, -np.abs(tstat))),
        rtol=1e-6, atol=1e-8, equal_nan=True)
    assert manifest["shape"] == list(beta.shape)
    assert manifest["df"] == 990
    assert summary["payload_bytes"] == beta.nbytes + tstat.nbytes + stored_logp.nbytes


def test_binary_logp_remains_finite_above_one_thousand(tmp_path):
    tstat = np.asarray([[100.0]], dtype=np.float32)
    logp = upper_tail_log10_from_t(tstat, 20000.0)
    writer = BinarySumstatsWriter(
        directory=tmp_path / "extreme",
        n_variants=1,
        trait_names=["trait"],
        n_samples=20002,
        df=20000,
        fsync=False,
    )
    writer.write_chunk(0, 1, np.ones_like(tstat), tstat, logp)
    writer.close()
    _beta, _tstat, stored_logp, _manifest = open_binary_sumstats(
        tmp_path / "extreme")
    assert np.isfinite(stored_logp[0, 0])
    assert stored_logp[0, 0] > 1000.0
    assert stored_logp[0, 0] == pytest.approx(float(logp[0, 0]), rel=1e-7)


def test_manifest_preserves_trait_specific_df(tmp_path):
    n_variants, n_traits = 4, 3
    trait_df = [87, 89, 86]
    writer = BinarySumstatsWriter(
        directory=tmp_path / "sumstats",
        n_variants=n_variants,
        trait_names=[f"trait_{index}" for index in range(n_traits)],
        n_samples=100,
        df=trait_df,
        store_beta=False,
        fsync=False,
    )
    tstat = np.arange(n_variants * n_traits, dtype=np.float32).reshape(
        n_variants, n_traits
    )
    writer.write_chunk(0, n_variants, None, tstat)
    writer.close()

    manifest = read_manifest(tmp_path / "sumstats")
    assert manifest["df"] == trait_df
    assert "p_value" not in manifest
    assert "neg_log10_p" in manifest["significance"]


def test_payload_is_twelve_bytes_per_cell(tmp_path):
    directory, beta, _tstat, summary = _write_store(tmp_path)
    assert summary["payload_bytes"] == 12 * beta.size
    assert summary["cells"] == beta.size


@pytest.mark.parametrize("block_bytes", [64, 1000, 1 << 20])
def test_block_size_does_not_change_bytes(tmp_path, block_bytes):
    directory, beta, tstat, _ = _write_store(
        tmp_path / str(block_bytes), block_bytes=block_bytes
    )
    stored_beta, stored_t, _stored_logp, _ = open_binary_sumstats(directory)
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)


def _write_geometry(tmp_path, n_variants, n_traits, chunk, block_bytes, borrow, seed=11):
    """Drive the writer with a chosen chunk/block geometry and return the data."""
    rng = np.random.default_rng(seed)
    beta = rng.normal(size=(n_variants, n_traits)).astype(np.float32)
    tstat = rng.normal(scale=3.0, size=(n_variants, n_traits)).astype(np.float32)
    beta[0] = np.nan
    tstat[0] = np.nan
    if n_variants > 3:
        beta[n_variants - 2] = np.nan
        tstat[n_variants - 2] = np.nan
    writer = BinarySumstatsWriter(
        directory=tmp_path,
        n_variants=n_variants,
        trait_names=[f"trait_{i}" for i in range(n_traits)],
        n_samples=5000,
        df=4970,
        block_bytes=block_bytes,
        queue_depth=2,
        borrow_chunks=borrow,
    )
    for start in range(0, n_variants, chunk):
        end = min(start + chunk, n_variants)
        # Owned per-chunk copies, matching the scan contract.
        writer.write_chunk(start, end, beta[start:end].copy(), tstat[start:end].copy())
    return beta, tstat, writer.close()


def test_pass_through_fires_when_a_chunk_fills_a_block(tmp_path):
    """Chunks at or above the block size must skip the staging memcpy.

    Regression guard for a geometry where block_bytes exceeds the per-chunk
    payload: the writer then silently coalesces every chunk through staging and
    pays a full memcpy of the result stream. Block counts are the only visible
    difference, so assert on them.
    """
    n_variants, n_traits, chunk = 4000, 64, 500
    chunk_payload = chunk * n_traits * 4  # 128000 bytes
    block_bytes = 4096
    assert chunk_payload >= block_bytes

    beta, tstat, summary = _write_geometry(
        tmp_path / "borrow", n_variants, n_traits, chunk, block_bytes, borrow=True
    )
    chunks = n_variants // chunk
    # One queued item per chunk per array, with no staging coalescing.
    assert summary["blocks"] == 3 * chunks
    assert summary["payload_bytes"] == 3 * tstat.nbytes

    stored_beta, stored_t, _stored_logp, _ = open_binary_sumstats(tmp_path / "borrow")
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)


def test_borrow_disabled_coalesces_through_staging(tmp_path):
    """The same geometry with borrowing off must stage, and still be exact."""
    n_variants, n_traits, chunk = 4000, 64, 500
    block_bytes = 4096
    beta, tstat, summary = _write_geometry(
        tmp_path / "staged", n_variants, n_traits, chunk, block_bytes, borrow=False
    )
    per_array = beta.nbytes
    expected_blocks = 3 * -(-per_array // block_bytes)
    assert summary["blocks"] == expected_blocks
    assert summary["blocks"] > 2 * (n_variants // chunk)
    assert summary["payload_bytes"] == 3 * tstat.nbytes

    stored_beta, stored_t, _stored_logp, _ = open_binary_sumstats(tmp_path / "staged")
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)


def test_oversized_block_silently_disables_pass_through(tmp_path):
    """Document the misconfiguration: block larger than a chunk always stages.

    This is the shape the full-scale run hit — chunk 10000 x 128 traits is
    4.88 MiB against a 32 MiB block, so every chunk took the copy.
    """
    n_variants, n_traits, chunk = 4000, 64, 500
    chunk_payload = chunk * n_traits * 4
    block_bytes = 4 * chunk_payload
    _beta, _tstat, summary = _write_geometry(
        tmp_path / "oversized", n_variants, n_traits, chunk, block_bytes, borrow=True
    )
    # Borrowing is requested but cannot engage: four chunks per block.
    assert summary["blocks"] == 3 * (n_variants // chunk) // 4


@pytest.mark.parametrize(
    "n_variants,n_traits,chunk,block_bytes",
    [
        (4000, 64, 500, 4096),        # pass-through, exact multiple
        (4000, 64, 500, 1 << 20),     # staging, block far above chunk
        (4001, 64, 500, 4096),        # ragged tail below the block
        (997, 7, 333, 8192),          # non-aligned chunk and trait width
        (64, 3, 1, 64),               # one variant per chunk
        (1500, 128, 700, 2048),       # tail chunk smaller than a full chunk
    ],
)
@pytest.mark.parametrize("borrow", [True, False])
def test_geometries_are_bit_exact(tmp_path, n_variants, n_traits, chunk, block_bytes, borrow):
    directory = tmp_path / f"{n_variants}_{chunk}_{block_bytes}_{borrow}"
    beta, tstat, summary = _write_geometry(
        directory, n_variants, n_traits, chunk, block_bytes, borrow
    )
    stored_beta, stored_t, stored_logp, manifest = open_binary_sumstats(directory)
    assert manifest["shape"] == [n_variants, n_traits]
    assert summary["payload_bytes"] == beta.nbytes + tstat.nbytes + stored_logp.nbytes
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)
    assert np.isnan(np.asarray(stored_t)[0]).all()


def test_borrowed_payload_is_not_read_after_the_caller_mutates_it(tmp_path):
    """Borrowing must copy when the caller keeps the array.

    With borrow_chunks=False the writer owns a staged copy, so mutating the
    source array after the call cannot corrupt the store.
    """
    n_variants, n_traits, chunk = 2000, 64, 500
    rng = np.random.default_rng(5)
    beta = rng.normal(size=(n_variants, n_traits)).astype(np.float32)
    tstat = rng.normal(size=(n_variants, n_traits)).astype(np.float32)
    directory = tmp_path / "reused"
    writer = BinarySumstatsWriter(
        directory=directory,
        n_variants=n_variants,
        trait_names=[f"t{i}" for i in range(n_traits)],
        n_samples=100,
        df=95,
        block_bytes=4096,
        queue_depth=2,
        borrow_chunks=False,
    )
    scratch_beta = np.empty((chunk, n_traits), dtype=np.float32)
    scratch_t = np.empty((chunk, n_traits), dtype=np.float32)
    for start in range(0, n_variants, chunk):
        end = start + chunk
        scratch_beta[:] = beta[start:end]
        scratch_t[:] = tstat[start:end]
        writer.write_chunk(start, end, scratch_beta, scratch_t)
        scratch_beta.fill(np.nan)  # a reusing iterator would do this next
        scratch_t.fill(np.nan)
    writer.close()
    stored_beta, stored_t, _stored_logp, _ = open_binary_sumstats(directory)
    np.testing.assert_array_equal(np.asarray(stored_beta), beta)
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)


def test_t_only_store_omits_beta_and_keeps_logp(tmp_path):
    n_variants, n_traits, chunk = 600, 16, 100
    rng = np.random.default_rng(19)
    beta = rng.normal(size=(n_variants, n_traits)).astype(np.float32)
    tstat = rng.normal(scale=2.0, size=(n_variants, n_traits)).astype(np.float32)
    directory = tmp_path / "screen"
    writer = BinarySumstatsWriter(
        directory=directory,
        n_variants=n_variants,
        trait_names=[f"trait_{i}" for i in range(n_traits)],
        n_samples=900,
        df=880,
        block_bytes=4096,
        store_beta=False,
    )
    for start in range(0, n_variants, chunk):
        end = start + chunk
        writer.write_chunk(start, end, beta[start:end].copy(), tstat[start:end].copy())
    summary = writer.close()

    assert summary["payload_bytes"] == 8 * n_variants * n_traits
    assert summary["payload_bytes"] == 2 * tstat.nbytes
    assert not (directory / "beta.f32").exists()

    stored_beta, stored_t, stored_logp, manifest = open_binary_sumstats(directory)
    assert stored_beta is None
    assert "beta" not in manifest["arrays"]
    np.testing.assert_array_equal(np.asarray(stored_t), tstat)
    assert stored_logp.shape == tstat.shape




def test_t_only_writer_with_no_chunks_creates_no_beta(tmp_path):
    """The zero-chunk close path must respect the configured field set."""
    directory = tmp_path / "empty"
    writer = BinarySumstatsWriter(
        directory=directory,
        n_variants=0,
        trait_names=["a", "b"],
        n_samples=10,
        df=5,
        store_beta=False,
    )
    summary = writer.close()
    assert summary["payload_bytes"] == 0
    assert not (directory / "beta.f32").exists()
    assert (directory / "tstat.f32").exists()
    stored_beta, _stored_t, _stored_logp, manifest = open_binary_sumstats(directory)
    assert stored_beta is None
    assert "beta" not in manifest["arrays"]


def test_t_only_is_rejected_when_output_is_disabled(tmp_path):
    from torchgwas.api import run_linear_gwas

    genotype, phenotype, covariates = _toy_linear_inputs()
    with pytest.raises(ValueError, match="applies only to"):
        run_linear_gwas(
            genotype=genotype,
            phenotype=phenotype,
            covariates=covariates,
            device="cpu",
            output_dir=tmp_path / "bad",
            sumstats_format="none",
            sumstats_fields="t",
        )


def test_out_of_order_chunk_is_rejected(tmp_path):
    writer = BinarySumstatsWriter(
        directory=tmp_path / "s",
        n_variants=10,
        trait_names=["a"],
        n_samples=100,
        df=95,
        block_bytes=256,
    )
    try:
        writer.write_chunk(0, 4, np.zeros((4, 1), np.float32), np.zeros((4, 1), np.float32))
        with pytest.raises(ValueError, match="in order"):
            writer.write_chunk(6, 10, np.zeros((4, 1), np.float32), np.zeros((4, 1), np.float32))
    finally:
        writer.abort()


def test_incomplete_store_is_rejected(tmp_path):
    writer = BinarySumstatsWriter(
        directory=tmp_path / "s",
        n_variants=10,
        trait_names=["a"],
        n_samples=100,
        df=95,
        block_bytes=256,
    )
    writer.write_chunk(0, 4, np.zeros((4, 1), np.float32), np.zeros((4, 1), np.float32))
    with pytest.raises(ValueError, match="expected 10"):
        writer.close()


def test_wrong_trait_width_is_rejected(tmp_path):
    writer = BinarySumstatsWriter(
        directory=tmp_path / "s",
        n_variants=4,
        trait_names=["a", "b"],
        n_samples=100,
        df=95,
        block_bytes=256,
    )
    try:
        with pytest.raises(ValueError, match="expected chunk shape"):
            writer.write_chunk(0, 4, np.zeros((4, 3), np.float32), np.zeros((4, 3), np.float32))
    finally:
        writer.abort()






def _toy_linear_inputs(seed=3, n_samples=200, n_markers=97, n_traits=4):
    rng = np.random.default_rng(seed)
    genotype = rng.integers(0, 3, size=(n_samples, n_markers)).astype(np.float32)
    genotype[:, 5] = 1.0  # invariant: excluded by the scan contract
    covariates = rng.normal(size=(n_samples, 3))
    phenotype = rng.normal(size=(n_samples, n_traits))
    return genotype, phenotype, covariates


def test_binary_store_matches_in_memory_statistics(tmp_path):
    from torchgwas.api import run_linear_gwas

    genotype, phenotype, covariates = _toy_linear_inputs()
    marker_ids = [f"rs{i}" for i in range(genotype.shape[1])]
    trait_columns = [f"trait_{i}" for i in range(phenotype.shape[1])]
    shared = dict(
        genotype=genotype,
        phenotype=phenotype,
        covariates=covariates,
        marker_ids=marker_ids,
        trait_columns=trait_columns,
        chunk_size=11,
        device="cpu",
        compute_dtype="float64",
    )
    binary_dir = tmp_path / "binary"
    reference = run_linear_gwas(**shared)
    run_linear_gwas(output_dir=binary_dir, sumstats_format="binary", **shared)

    rows = reference.table
    beta, tstat, logp, manifest = open_binary_sumstats(binary_dir / "sumstats")
    assert manifest["traits"] == trait_columns

    # The in-memory path drops the invariant column, so the store covers the
    # retained markers and variant_ids.txt names exactly those.
    stored_ids = (binary_dir / "sumstats" / "variant_ids.txt").read_text().split()
    assert "rs5" not in stored_ids
    assert stored_ids == [name for name in marker_ids if name != "rs5"]
    assert manifest["shape"] == [len(stored_ids), phenotype.shape[1]]

    marker_position = {name: index for index, name in enumerate(stored_ids)}
    trait_position = {name: index for index, name in enumerate(trait_columns)}
    assert rows, "in-memory reference produced no rows"
    assert len(rows) == len(stored_ids) * len(trait_columns)
    assert all("p_value" not in row for row in rows)
    for row in rows:
        variant = marker_position[row["marker_id"]]
        trait = trait_position[row["trait"]]
        assert float(beta[variant, trait]) == pytest.approx(float(row["beta"]), rel=1e-6)
        assert float(tstat[variant, trait]) == pytest.approx(float(row["t_stat"]), rel=1e-6)
        assert float(logp[variant, trait]) == pytest.approx(float(row["-log10_p"]), rel=1e-6)


def test_dropped_columns_do_not_shift_marker_labels(tmp_path):
    """Regression: an invariant column must not relabel every later marker."""
    from torchgwas.api import run_linear_gwas

    genotype, phenotype, covariates = _toy_linear_inputs()
    marker_ids = [f"rs{i}" for i in range(genotype.shape[1])]
    result = run_linear_gwas(
        genotype=genotype,
        phenotype=phenotype,
        covariates=covariates,
        marker_ids=marker_ids,
        chunk_size=11,
        device="cpu",
        compute_dtype="float64",
    )
    reported = {row["marker_id"] for row in result.table}
    assert "rs5" not in reported
    assert reported == set(marker_ids) - {"rs5"}
    # The last marker must still be named rs96, not rs95.
    assert marker_ids[-1] in reported


def test_sumstats_format_none_writes_no_table(tmp_path):
    from torchgwas.api import run_linear_gwas

    genotype, phenotype, covariates = _toy_linear_inputs()
    out = tmp_path / "none"
    run_linear_gwas(
        genotype=genotype,
        phenotype=phenotype,
        covariates=covariates,
        chunk_size=11,
        device="cpu",
        compute_dtype="float64",
        output_dir=out,
        sumstats_format="none",
    )
    assert not (out / "results.tsv.gz").exists()
    assert not (out / "sumstats").exists()
    run_metadata = json.loads((out / "run.json").read_text())
    assert run_metadata["sumstats_format"] == "none"


def test_threshold_selector_uses_indexed_binary(tmp_path):
    from torchgwas.api import run_linear_gwas
    from binary_helpers import binary_rows
    genotype, phenotype, covariates = _toy_linear_inputs()
    run_linear_gwas(
        genotype,
        phenotype,
        covariates=covariates,
        device="cpu",
        output_dir=tmp_path,
        p_value_threshold=1.0,
    )
    assert len(binary_rows(tmp_path)) == (genotype.shape[1] - 1) * phenotype.shape[1]
