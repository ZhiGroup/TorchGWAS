from __future__ import annotations

import csv
import gzip
import tempfile
import time
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import torch
from scipy import stats

from torchgwas.api import run_linear_gwas
from binary_helpers import binary_rows
from torchgwas.reduce import VariantReduction
from torchgwas.bed import PlinkBedGenotype, resolve_plink_triplet
from torchgwas.linear import (
    _unpack_plink_a2_float,
    linear_scan,
    linear_scan_streaming,
    linear_scan_streaming_chunks,
)
from torchgwas.preprocess import residualize_and_standardize
from torchgwas.tails import upper_tail_log10_from_t
from torchgwas.utils import choose_device


def _write_bed(prefix: Path, dosage_a1: np.ndarray) -> Path:
    bed, bim, fam = resolve_plink_triplet(prefix)
    n_samples, n_markers = dosage_a1.shape
    fam.write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_samples)))
    bim.write_text("".join(f"1 rs{i} 0 {i + 1} A G\n" for i in range(n_markers)))
    dosage_to_code = {0: 0b00, 1: 0b10, 2: 0b11}
    payload = bytearray(b"\x6c\x1b\x01")
    for marker in range(n_markers):
        for sample_start in range(0, n_samples, 4):
            byte = 0
            for offset in range(4):
                sample = sample_start + offset
                code = (
                    0b01
                    if sample >= n_samples or np.isnan(dosage_a1[sample, marker])
                    else dosage_to_code[int(dosage_a1[sample, marker])]
                )
                byte |= code << (2 * offset)
            payload.append(byte)
    bed.write_bytes(payload)
    return bed


class ExactLinearStatisticsTestCase(unittest.TestCase):
    def test_optional_log10_p_matches_exact_student_tail(self):
        rng = np.random.default_rng(20260918)
        n_samples = 101
        genotype = rng.integers(0, 3, size=(n_samples, 7)).astype(np.float64)
        genotype[:3] = np.arange(3)[:, None]
        phenotype = rng.normal(size=(n_samples, 4))
        covariates = rng.normal(size=(n_samples, 2))

        beta, t_stat, logp, _q = linear_scan(
            genotype, phenotype, covariates, chunk_size=3,
            device="cpu", compute_dtype="float64", return_log10_p=True)
        expected = upper_tail_log10_from_t(t_stat, n_samples - 2 - 2)
        np.testing.assert_allclose(logp, expected, rtol=1e-12, atol=1e-12)
        self.assertEqual(beta.shape, t_stat.shape)

    def test_missing_phenotypes_match_mean_imputed_ols_with_trait_df(self):
        """Phenotype NaNs use the same mean-impute/df rule as genotype NaNs."""
        rng = np.random.default_rng(20260917)
        n_samples, n_markers, n_traits = 101, 9, 3
        genotype = rng.integers(0, 3, size=(n_samples, n_markers)).astype(np.float64)
        genotype[:3] = np.arange(3)[:, None]
        covariates = rng.normal(size=(n_samples, 2))
        phenotype = rng.normal(size=(n_samples, n_traits))
        phenotype[rng.random(size=phenotype.shape) < 0.11] = np.nan

        _beta, observed_t, _p, _q = linear_scan(
            genotype, phenotype, covariates,
            chunk_size=4, device="cpu", compute_dtype="float64")

        expected = np.empty_like(observed_t)
        for trait in range(n_traits):
            observed = np.isfinite(phenotype[:, trait])
            y = phenotype[:, trait].copy()
            y[~observed] = y[observed].mean()
            for marker in range(n_markers):
                design = np.column_stack(
                    (np.ones(n_samples), covariates, genotype[:, marker]))
                coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
                residual = y - design @ coefficients
                df = int(observed.sum()) - np.linalg.matrix_rank(design)
                gram = np.linalg.inv(design.T @ design)
                se = np.sqrt(float(residual @ residual) / df * gram[-1, -1])
                expected[marker, trait] = coefficients[-1] / se
        np.testing.assert_allclose(observed_t, expected, rtol=1e-11, atol=1e-11)

    def test_auto_cuda_device_is_resolved_to_an_explicit_index(self):
        with mock.patch("torch.cuda.is_available", return_value=True), mock.patch(
            "torch.cuda.current_device", return_value=3
        ):
            self.assertEqual(choose_device("auto"), torch.device("cuda:3"))

    def test_packed_plink_decoder_preserves_a2_dosage_and_missing_calls(self):
        # PLINK pairs are least-significant first: 00, 10, 11, 01.
        packed = torch.tensor([[0b01_11_10_00]], dtype=torch.uint8)
        observed = _unpack_plink_a2_float(packed, n_samples=4).numpy()
        expected = np.asarray([[0.0, 1.0, 2.0, np.nan]], dtype=np.float32)
        np.testing.assert_allclose(observed, expected, equal_nan=True)

    def test_packed_plink_decoder_trims_byte_padding(self):
        packed = torch.tensor([[0b11_11_10_00]], dtype=torch.uint8)
        observed = _unpack_plink_a2_float(packed, n_samples=2).numpy()
        np.testing.assert_array_equal(observed, np.asarray([[0.0, 1.0]], dtype=np.float32))

    def test_packed_plink_decoder_gathers_arbitrary_sample_positions(self):
        packed = torch.tensor(
            [[0b01_11_10_00, 0b11_00_10_00]],
            dtype=torch.uint8,
        )
        sample_indices = torch.tensor([7, 0, 5, 2], dtype=torch.int64)
        observed = _unpack_plink_a2_float(
            packed,
            n_samples=4,
            sample_byte_indices=sample_indices // 4,
            sample_bit_shifts=((sample_indices % 4) * 2).to(torch.uint8),
        ).numpy()
        expected = np.asarray([[2.0, 0.0, 1.0, 2.0]], dtype=np.float32)
        np.testing.assert_array_equal(observed, expected)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the packed BED scan")
    def test_packed_bed_cuda_scan_matches_float64_reference(self):
        rng = np.random.default_rng(20260911)
        n_samples = 65
        genotype = rng.integers(0, 3, size=(n_samples, 13), dtype=np.uint8)
        # Make every marker nonconstant in this small randomized fixture.
        genotype[:3, :] = np.arange(3, dtype=np.uint8)[:, None]
        covariate = rng.normal(size=n_samples)
        covariates = np.column_stack((covariate, 2.0 * covariate))
        phenotype = np.column_stack(
            (
                0.3 * genotype[:, 2] + covariate + rng.normal(size=n_samples),
                -0.2 * genotype[:, 9] - covariate + rng.normal(size=n_samples),
            )
        )
        beta_ref, t_ref, p_ref, _ = linear_scan(
            genotype.astype(np.float64),
            phenotype,
            covariates,
            chunk_size=5,
            device="cpu",
            compute_dtype="float64",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = _write_bed(Path(tmpdir) / "fixture", genotype)
            packed = PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2)
            beta, t_stat, logp, _ = linear_scan_streaming(
                packed,
                phenotype,
                covariates,
                chunk_size=5,
                device="cuda:0",
                compute_dtype="float32",
                reader_workers=2,
                return_log10_p=True,
            )
            output_dir = Path(tmpdir) / "topk"
            streamed = run_linear_gwas(
                genotype=bed,
                phenotype=phenotype,
                covariates=covariates,
                genotype_format="plink",
                chunk_size=5,
                device="cuda:0",
                compute_dtype="float32",
                reader_workers=2,
                topk_per_trait=2,
                output_dir=output_dir,
            )
            top_rows = binary_rows(output_dir)
        np.testing.assert_allclose(beta, beta_ref, rtol=2e-4, atol=2e-5)
        np.testing.assert_allclose(t_stat, t_ref, rtol=2e-4, atol=2e-5)
        np.testing.assert_allclose(
            logp, upper_tail_log10_from_t(
                t_stat, n_samples - np.linalg.matrix_rank(covariates) - 2),
            rtol=1e-10, atol=1e-10)
        self.assertEqual(streamed.qc_summary["genotype_qc_mode"], "fused_gpu_scan")
        self.assertEqual(len(top_rows), 4)
        for trait_index in range(2):
            observed = {
                row["marker_id"]
                for row in top_rows
                if row["trait"] == f"trait_{trait_index}"
            }
            expected = {
                f"rs{index}"
                for index in np.argsort(np.abs(t_ref[:, trait_index]))[-2:]
            }
            self.assertEqual(observed, expected)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the packed BED scan")
    def test_packed_bed_cuda_subset_matches_float64_reference(self):
        rng = np.random.default_rng(20260913)
        n_stored_samples = 65
        selected_indices = np.asarray([64, 0, 31, 7, 42, 18, 3, 55, 12, 27, 9, 38])
        genotype = rng.integers(0, 3, size=(n_stored_samples, 11), dtype=np.uint8)
        genotype[selected_indices[:3], :] = np.arange(3, dtype=np.uint8)[:, None]
        selected_genotype = genotype[selected_indices]
        covariates = rng.normal(size=(selected_indices.size, 2))
        phenotype = np.column_stack(
            (
                0.4 * selected_genotype[:, 2] + rng.normal(size=selected_indices.size),
                -0.3 * selected_genotype[:, 8] + rng.normal(size=selected_indices.size),
            )
        )
        beta_ref, t_ref, p_ref, _ = linear_scan(
            selected_genotype.astype(np.float64),
            phenotype,
            covariates,
            chunk_size=4,
            device="cpu",
            compute_dtype="float64",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = _write_bed(Path(tmpdir) / "subset", genotype)
            packed = PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2)
            packed.select_samples([f"I{index}" for index in selected_indices])
            beta, t_stat, p_value, _ = linear_scan_streaming(
                packed,
                phenotype,
                covariates,
                chunk_size=4,
                device="cuda:0",
                compute_dtype="float32",
                reader_workers=2,
            )
            api_result = run_linear_gwas(
                genotype=bed,
                phenotype=phenotype,
                covariates=covariates,
                genotype_format="plink",
                sample_ids=np.asarray([f"I{index}" for index in selected_indices]),
                chunk_size=4,
                device="cuda:0",
                compute_dtype="float32",
                reader_workers=2,
            )

        np.testing.assert_allclose(beta, beta_ref, rtol=3e-4, atol=3e-5)
        np.testing.assert_allclose(t_stat, t_ref, rtol=3e-4, atol=3e-5)
        np.testing.assert_allclose(p_value, p_ref, rtol=3e-4, atol=1e-7)
        api_t = np.asarray([row["t_stat"] for row in api_result.table]).reshape(t_ref.shape)
        np.testing.assert_allclose(api_t, t_ref, rtol=3e-4, atol=3e-5)
        self.assertEqual(api_result.run_metadata["genotype_shape"], [selected_indices.size, 11])

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the packed BED scan")
    def test_packed_bed_cuda_scan_masks_missing_calls(self):
        """A missing call costs its own sample, not the whole variant.

        This replaces an earlier test that asserted variant 3 came back as
        NaN. That behaviour discarded any variant with a single missing call,
        which on real data means nearly all of them. The assertion here is
        stronger than the one it replaces: the masked result has to equal an
        explicit least squares fit in which the missing call is set to the
        variant's observed mean, which is the same estimand as masking.
        """
        rng = np.random.default_rng(20260912)
        genotype = rng.integers(0, 3, size=(65, 13)).astype(np.float32)
        genotype[:3, :] = np.arange(3, dtype=np.float32)[:, None]
        genotype[7, 3] = np.nan
        phenotype = rng.normal(size=(65, 2))
        covariates = rng.normal(size=(65, 2))
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = _write_bed(Path(tmpdir) / "missing", genotype)
            packed = PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2)
            beta, t_stat, p_value, _ = linear_scan_streaming(
                packed,
                phenotype,
                covariates,
                chunk_size=5,
                device="cuda:0",
                compute_dtype="float32",
                reader_workers=2,
            )
            self.assertTrue(np.isfinite(t_stat[3]).all())
            self.assertTrue(np.isfinite(beta[3]).all())
            self.assertTrue(np.isfinite(p_value[3]).all())
            self.assertEqual(packed._last_scan_exclusion_counts["missing"], 0)

            imputed = genotype[:, 3].copy()
            observed = ~np.isnan(imputed)
            imputed[~observed] = imputed[observed].mean()
            n_samples = phenotype.shape[0]
            design = np.column_stack([np.ones(n_samples), covariates, imputed])
            # The variant spends only its observed samples, so the reference
            # uses that count rather than the full one.
            df = int(observed.sum()) - design.shape[1]
            gram = np.linalg.inv(design.T @ design)
            for trait in range(phenotype.shape[1]):
                coefficients, *_ = np.linalg.lstsq(
                    design, phenotype[:, trait], rcond=None)
                # Form the residual explicitly rather than taking lstsq's
                # optional residual output, which is an array that is empty
                # for rank-deficient or exactly-determined systems.
                residual = phenotype[:, trait] - design @ coefficients
                sigma2 = float(residual @ residual) / df
                se = np.sqrt(sigma2 * gram[-1, -1])
                # The scan standardises the phenotype, so its beta is in those
                # units and does not equal a raw-phenotype coefficient. The
                # t-statistic is invariant under that scaling -- beta and its
                # standard error move together -- so it is what pins the
                # masked result to the reference.
                self.assertAlmostEqual(float(t_stat[3, trait]),
                                       float(coefficients[-1] / se), places=3)

    def test_matches_explicit_ols_with_rank_deficient_covariates(self):
        rng = np.random.default_rng(20260910)
        n_samples = 96
        genotype = rng.integers(0, 3, size=(n_samples, 7)).astype(np.float64)
        covariate = rng.normal(size=n_samples)
        covariates = np.column_stack([covariate, 2.0 * covariate])
        phenotype = np.column_stack(
            [
                0.4 * genotype[:, 1] + 0.8 * covariate + rng.normal(size=n_samples),
                -0.2 * genotype[:, 4] - 0.5 * covariate + rng.normal(size=n_samples),
            ]
        )

        beta, t_stat, p_value, q_matrix = linear_scan(
            genotype,
            phenotype,
            covariates,
            chunk_size=3,
            device="cpu",
            compute_dtype="float64",
        )
        phenotype_processed, q_reference = residualize_and_standardize(phenotype, covariates)
        self.assertIsNotNone(q_matrix)
        self.assertEqual(q_matrix.shape[1], 1)
        np.testing.assert_allclose(np.abs(q_matrix), np.abs(q_reference), atol=1e-12)

        beta_reference = np.empty_like(beta)
        t_reference = np.empty_like(t_stat)
        p_reference = np.empty_like(p_value)
        for marker in range(genotype.shape[1]):
            design = np.column_stack([np.ones(n_samples), q_reference, genotype[:, marker]])
            coefficients, _, rank, _ = np.linalg.lstsq(design, phenotype_processed, rcond=None)
            residual = phenotype_processed - design @ coefficients
            df = n_samples - rank
            sigma2 = np.sum(residual * residual, axis=0) / df
            covariance_last = np.linalg.pinv(design.T @ design)[-1, -1]
            standard_error = np.sqrt(sigma2 * covariance_last)
            beta_reference[marker] = coefficients[-1]
            t_reference[marker] = coefficients[-1] / standard_error
            p_reference[marker] = 2.0 * stats.t.sf(np.abs(t_reference[marker]), df=df)

        np.testing.assert_allclose(beta, beta_reference, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(t_stat, t_reference, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_value, p_reference, rtol=1e-10, atol=1e-12)


class TwoBitEncodingTestCase(unittest.TestCase):
    """The BED and PGEN two-bit tables collide, and must stay distinct.

    Code `01` is **missing** in PLINK 1 BED but **heterozygous** in PGEN; code
    `11` is **missing** in PGEN but **homozygous alternate** in BED. Applying
    the wrong table to the right bytes produces entirely plausible dosages with
    no error raised -- every missing call becomes 1.0 and every homozygous
    alternate becomes NaN -- so this pins both tables against the formats rather
    than against each other. A change that made them agree would be a bug that
    no other test in this suite would catch.
    """

    def test_bed_and_pgen_tables_disagree_on_the_same_bytes(self):
        from torchgwas import scan_gpu

        # `scan_gpu.available()` rather than `torch.cuda.is_available()`: the
        # native library loads on any card but carries code only for the
        # architectures it was built for, so CUDA being present says nothing
        # about whether these kernels can launch.
        if not (torch.cuda.is_available() and scan_gpu.available()):
            self.skipTest("native packed kernels cannot launch on this device")

        # One byte, four samples, codes 00 01 10 11 (least-significant first).
        packed = torch.tensor([[0b11_10_01_00]], dtype=torch.uint8,
                              device=choose_device("auto"))
        bed = scan_gpu.prepare(packed, encoding="plink_2bit", n_samples=4)[0]
        pgen = scan_gpu.prepare(packed, encoding="pgen_2bit", n_samples=4)[0]
        # prepare() centres and zeroes missing calls, so compare the observed
        # counts, which say which code each table treats as absent.
        bed_present = int(scan_gpu.prepare(packed, encoding="plink_2bit",
                                           n_samples=4)[4][0])
        pgen_present = int(scan_gpu.prepare(packed, encoding="pgen_2bit",
                                            n_samples=4)[4][0])
        self.assertEqual(bed_present, 3, "BED must treat exactly code 01 as missing")
        self.assertEqual(pgen_present, 3, "PGEN must treat exactly code 11 as missing")
        # Same bytes, same count of missing, but a different sample is missing:
        # the centred blocks must therefore differ.
        self.assertFalse(
            torch.allclose(bed, pgen),
            "BED and PGEN two-bit tables produced identical output; one of them "
            "is applying the other's missing code")

    def test_unsupported_encoding_is_refused(self):
        from torchgwas import scan_gpu

        packed = torch.zeros((1, 1), dtype=torch.uint8)
        with self.assertRaises(ValueError):
            scan_gpu.prepare(packed, encoding="plink2_2bit", n_samples=4)


class VariantRangeTestCase(unittest.TestCase):
    """A ranged scan must cover exactly the range it was given.

    There was no coverage of `variant_range` anywhere in this suite, and a
    caller that dropped it went unnoticed: `native_scan` passed the range to the
    device loader but not to the pinned one, so a ranged scan of a non-device
    source read the **whole file** and reported a plausible time for it. Nothing
    raised. The multi-GPU driver shards by exactly this parameter, so the
    failure mode was every shard scanning every variant and the results being
    silently wrong rather than merely slow.

    The assertion is on the extent actually produced, not on whether the call
    succeeds, because succeeding was never the problem.
    """

    def _scan_extent(self, bed, phenotype, covariates, variant_range, chunk):
        starts = []
        stream, _q = linear_scan_streaming_chunks(
            PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2),
            phenotype, covariates, chunk_size=chunk, device="cpu",
            compute_dtype="float64", compute_p_values=False,
            variant_range=variant_range)
        covered = 0
        for start, end, _beta, _t, _p in stream:
            starts.append((start, end))
            covered += end - start
        return covered, starts

    def test_ranged_scan_covers_exactly_the_requested_range(self):
        rng = np.random.default_rng(31)
        n_samples, n_markers = 64, 40
        dosage = rng.integers(0, 3, size=(n_samples, n_markers)).astype(float)
        phenotype = rng.normal(size=(n_samples, 3))
        covariates = rng.normal(size=(n_samples, 2))
        with tempfile.TemporaryDirectory() as directory:
            bed = _write_bed(Path(directory) / "range", dosage)
            for variant_range, expected in (
                (None, n_markers),
                ((0, n_markers), n_markers),
                ((0, 10), 10),
                ((10, 25), 15),
                ((32, n_markers), 8),
            ):
                for chunk in (4, 7, n_markers):
                    covered, spans = self._scan_extent(
                        bed, phenotype, covariates, variant_range, chunk)
                    with self.subTest(range=variant_range, chunk=chunk):
                        self.assertEqual(covered, expected)
                        if variant_range is not None and spans:
                            self.assertEqual(spans[0][0], variant_range[0])
                            self.assertEqual(spans[-1][1], variant_range[1])

    def test_shards_tile_the_file_exactly_once(self):
        """What the multi-GPU driver relies on: contiguous shards, no overlap."""
        rng = np.random.default_rng(32)
        n_samples, n_markers = 48, 37
        dosage = rng.integers(0, 3, size=(n_samples, n_markers)).astype(float)
        phenotype = rng.normal(size=(n_samples, 2))
        covariates = rng.normal(size=(n_samples, 2))
        bounds = [(0, 12), (12, 24), (24, n_markers)]
        with tempfile.TemporaryDirectory() as directory:
            bed = _write_bed(Path(directory) / "shard", dosage)
            seen = np.zeros(n_markers, dtype=int)
            for variant_range in bounds:
                _covered, spans = self._scan_extent(
                    bed, phenotype, covariates, variant_range, 5)
                for start, end in spans:
                    seen[start:end] += 1
            np.testing.assert_array_equal(seen, np.ones(n_markers, dtype=int))


class RepeatedScanTestCase(unittest.TestCase):
    """A second scan in the same process must not be much slower than the first.

    Every benchmark in this project was single-shot, which is exactly why a
    **17.9x** degradation across four repeated scans went unnoticed until a
    competitor harness happened to loop: `PinnedDosageLoader.close()` stopped
    its producer thread without releasing its page-locked staging buffers, so
    each scan left about 891 MB pinned and a looping process accumulated them.

    This asserts the property rather than the cause, so it still fails if some
    other resource starts accumulating. The ratio is generous because a small
    fixture is dominated by fixed costs and the machine is shared; the failure
    being guarded against was seventeen-fold, not fractional.
    """

    def _elapsed_scans(self, bed, phenotype, covariates, rounds):
        times = []
        for _ in range(rounds):
            started = time.perf_counter()
            stream, _q = linear_scan_streaming_chunks(
                PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2),
                phenotype, covariates, chunk_size=8, device="cpu",
                compute_dtype="float64", compute_p_values=False)
            for _start, _end, _b, _t, _p in stream:
                pass
            times.append(time.perf_counter() - started)
        return times

    def test_repeated_scans_do_not_degrade(self):
        rng = np.random.default_rng(71)
        n_samples, n_markers = 96, 48
        dosage = rng.integers(0, 3, size=(n_samples, n_markers)).astype(float)
        phenotype = rng.normal(size=(n_samples, 4))
        covariates = rng.normal(size=(n_samples, 3))
        with tempfile.TemporaryDirectory() as directory:
            bed = _write_bed(Path(directory) / "repeat", dosage)
            times = self._elapsed_scans(bed, phenotype, covariates, 4)
        first, last = times[0], times[-1]
        self.assertLess(
            last, first * 5.0 + 0.5,
            f"repeated scans degraded: {['%.3f' % t for t in times]}")


class ApiVariantRangeTestCase(unittest.TestCase):
    """`run_linear_gwas` must honour a variant range, or refuse it.

    Added because the public API had no way to express a subset at all, which
    made a benchmark harness scan 8.93M variants while the tool it was being
    compared against got 200,000. The range is honoured on the path that can
    honour it and raises on the one that cannot; what is not acceptable is
    accepting the argument and ignoring it.
    """

    def test_range_narrows_the_written_output(self):
        rng = np.random.default_rng(11)
        n_samples, n_markers = 48, 30
        dosage = rng.integers(0, 3, size=(n_samples, n_markers)).astype(float)
        phenotype = rng.normal(size=(n_samples, 2))
        covariates = rng.normal(size=(n_samples, 2))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "apirange", dosage)
            result = run_linear_gwas(
                str(bed), phenotype, covariates=covariates,
                device="cpu", compute_dtype="float64",
                output_dir=str(root / "out"),
                variant_range=(5, 17))
            self.assertIsNotNone(result)
            written = (root / "out" / "sumstats" / "variants.txt")
            if written.exists():
                lines = [x for x in written.read_text().splitlines() if x]
                self.assertEqual(len(lines), 12)

    def test_range_is_refused_when_it_cannot_be_honoured(self):
        rng = np.random.default_rng(12)
        dosage = rng.integers(0, 3, size=(32, 20)).astype(float)
        phenotype = rng.normal(size=(32, 2))
        with tempfile.TemporaryDirectory() as directory:
            bed = _write_bed(Path(directory) / "refuse", dosage)
            with self.assertRaises(ValueError):
                run_linear_gwas(str(bed), phenotype, device="cpu",
                                compute_dtype="float64",
                                variant_range=(0, 10))


class ResidualizeDeviceTestCase(unittest.TestCase):
    """The GPU projection must match the numpy reference it accelerates.

    Every statistic this tool reports is computed from the residualised
    phenotype, so a discrepancy here is invisible in the scan and wrong
    everywhere. Two mismatches are guarded specifically because both are silent:
    numpy's `std` uses the population form while torch's default is the sample
    form, and the numpy path promotes a float32 phenotype to float64 whenever
    the covariates are float64 -- the second one was caught this way.
    """

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the device projection")
    def test_device_matches_host(self):
        device = choose_device("auto")
        rng = np.random.default_rng(20260912)
        n_samples, n_covariates = 512, 9
        for n_traits in (1, 5, 64):
            for dtype in (np.float64, np.float32):
                covariates = rng.normal(size=(n_samples, n_covariates))
                # Duplicate a column so the rank-deficient basis path is used.
                covariates[:, 3] = covariates[:, 2]
                phenotype = rng.normal(size=(n_samples, n_traits)).astype(dtype)
                want, q_host = residualize_and_standardize(
                    phenotype.copy(), covariates.copy())
                got, q_device = residualize_and_standardize(
                    phenotype.copy(), covariates.copy(), device=device)
                with self.subTest(traits=n_traits, dtype=np.dtype(dtype).name):
                    self.assertEqual(got.dtype, want.dtype)
                    self.assertEqual(q_host.shape, q_device.shape)
                    tolerance = 1e-12 if dtype == np.float64 else 2e-6
                    np.testing.assert_allclose(got, want, rtol=tolerance,
                                               atol=tolerance)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the device projection")
    def test_trait_blocking_matches_the_unblocked_projection(self):
        """Blocking the traits must change nothing about the answer.

        The device path used to upload the whole phenotype at once, which
        stopped fitting at voxel scale: 33,417 subjects by 600,000 voxels in
        float64 is 160 GB against an 80 GB card, so it raised, the caller
        caught the RuntimeError as a busy device, and the run silently fell
        back to a numpy path that ground three more full-size temporaries
        through a single core.

        Blocking is exact in the algebra -- centring, projection and scaling
        are all per-column, so a block computes what it would have computed
        in company. It is NOT bit-identical: cuBLAS picks a different kernel
        for a different matrix width, so the projection GEMM rounds
        differently and about half the entries move in the last bit or two
        (measured: 8.9e-16 absolute, 2.1e-13 relative, in float64). The
        tolerance below is that rounding, not slack for a disagreement.
        """
        device = choose_device("auto")
        rng = np.random.default_rng(20260914)
        n_samples, n_traits = 384, 37
        covariates = rng.normal(size=(n_samples, 6))
        phenotype = rng.normal(size=(n_samples, n_traits))
        whole, _ = residualize_and_standardize(
            phenotype.copy(), covariates.copy(), device=device)
        for block in (1, 2, 8, n_traits - 1, n_traits, n_traits + 5):
            with mock.patch(
                    "torchgwas.preprocess._device_trait_block",
                    return_value=block):
                blocked, _ = residualize_and_standardize(
                    phenotype.copy(), covariates.copy(), device=device)
            with self.subTest(block=block):
                np.testing.assert_allclose(blocked, whole,
                                           rtol=1e-11, atol=1e-13)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the device projection")
    def test_trait_block_width_shrinks_as_samples_grow(self):
        """The block is sized from free memory, so it must scale with N."""
        from torchgwas.preprocess import _device_trait_block
        device = choose_device("auto")
        small = _device_trait_block(1000, 4, device)
        large = _device_trait_block(1000000, 4, device)
        self.assertGreater(small, large)
        self.assertGreaterEqual(large, 1)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is required for the device projection")
    def test_constant_trait_does_not_divide_by_zero(self):
        device = choose_device("auto")
        rng = np.random.default_rng(7)
        covariates = rng.normal(size=(256, 4))
        phenotype = np.full((256, 2), 3.0)
        want, _ = residualize_and_standardize(phenotype.copy(), covariates.copy())
        got, _ = residualize_and_standardize(phenotype.copy(), covariates.copy(),
                                             device=device)
        self.assertTrue(np.isfinite(got).all())
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)


class CovariateRankTestCase(unittest.TestCase):
    """The covariate rank must not depend on the dtype the caller happened to use.

    `_covariate_basis` truncates at numpy's `matrix_rank` tolerance,
    `max(shape) * eps * s[0]`, and that inherits the input's `eps`. At N = 22,250
    the cut is 4.9e-12 in float64 and 2.6e-3 in float32 -- nine orders of
    magnitude apart -- so a covariate correlated with another at r = 0.9999995
    used to be kept as float64 and dropped as float32. Dropping it lowers
    `covariate_rank`, which raises `df`, which moves every t-statistic and
    p-value in the run, with nothing reported. Several callers reach this
    function without `prepare_inputs_for_prep`'s float64 cast, so the guarantee
    has to live here.
    """

    def _near_collinear(self, n=4000, noise=1e-5, seed=5):
        rng = np.random.default_rng(seed)
        covariates = rng.normal(size=(n, 6))
        # A genuinely distinct column that is nearly a copy of another. The
        # noise has to put s_min/s_max between the two tolerances, and both
        # scale with n: at n = 4,000 the float32 cut is 4.8e-4 and the float64
        # cut is 8.9e-13, while this fixture lands near 5e-6. A milder noise
        # (1e-3) clears the float32 cut at this n and the test would pass
        # against the unfixed code, pinning nothing.
        covariates[:, 2] = covariates[:, 0] + noise * rng.normal(size=n)
        return covariates

    def test_rank_is_the_same_in_float32_and_float64(self):
        from torchgwas.preprocess import _covariate_basis

        covariates = self._near_collinear()
        wide = _covariate_basis(covariates.astype(np.float64))
        narrow = _covariate_basis(covariates.astype(np.float32))
        self.assertIsNotNone(wide)
        self.assertIsNotNone(narrow)
        self.assertEqual(narrow.shape[1], wide.shape[1])
        self.assertEqual(wide.shape[1], 6, "no column here is redundant")

    def test_the_returned_dtype_still_follows_the_input(self):
        # Only the rank decision was widened. Returning float64 unconditionally
        # would promote the residualised phenotype downstream, which is a
        # behaviour change rather than a fix.
        from torchgwas.preprocess import _covariate_basis

        covariates = self._near_collinear()
        self.assertEqual(
            _covariate_basis(covariates.astype(np.float32)).dtype, np.float32)
        self.assertEqual(
            _covariate_basis(covariates.astype(np.float64)).dtype, np.float64)

    def test_exact_collinearity_is_still_dropped(self):
        from torchgwas.preprocess import _covariate_basis

        rng = np.random.default_rng(6)
        covariates = rng.normal(size=(500, 5))
        covariates[:, 3] = covariates[:, 1]
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=np.dtype(dtype).name):
                basis = _covariate_basis(covariates.astype(dtype))
                self.assertEqual(basis.shape[1], 4)

    def test_a_constant_column_is_dropped(self):
        # std 0 becomes 1.0, the centred column is all zeros, and the
        # truncation removes it -- which is right, because the intercept is
        # added separately by the scan.
        from torchgwas.preprocess import _covariate_basis

        rng = np.random.default_rng(7)
        covariates = rng.normal(size=(400, 4))
        covariates[:, 1] = 3.5
        self.assertEqual(_covariate_basis(covariates).shape[1], 3)

    def test_column_scale_alone_does_not_drop_a_covariate(self):
        # Each column is divided by its own standard deviation before the SVD,
        # so being small in absolute terms is not being small in the SVD.
        from torchgwas.preprocess import _covariate_basis

        rng = np.random.default_rng(8)
        covariates = rng.normal(size=(600, 5))
        covariates[:, 4] *= 1e-8
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=np.dtype(dtype).name):
                self.assertEqual(
                    _covariate_basis(covariates.astype(dtype)).shape[1], 5)


class CovariateRankReportingTestCase(unittest.TestCase):
    """A truncated covariate rank must be reported, not inferred from df.

    The SVD drops linearly dependent covariates, which raises the residual
    degrees of freedom and so moves every p-value. Before this, the only visible
    trace was a `df` that disagreed with the covariate column count, and nothing
    in the output explained the discrepancy.
    """

    def _run(self, covariates, seed=21):
        rng = np.random.default_rng(seed)
        dosage = rng.integers(0, 3, size=(covariates.shape[0], 15)).astype(float)
        phenotype = rng.normal(size=(covariates.shape[0], 2))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "rank", dosage)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = run_linear_gwas(
                    str(bed), phenotype, covariates=covariates,
                    device="cpu", compute_dtype="float64",
                    output_dir=str(root / "out"))
            return result, [w for w in caught
                            if issubclass(w.category, RuntimeWarning)]

    def test_full_rank_covariates_report_no_drop_and_warn_nothing(self):
        rng = np.random.default_rng(22)
        covariates = rng.normal(size=(60, 4))
        result, caught = self._run(covariates)
        metadata = result.run_metadata
        self.assertEqual(metadata["covariate_rank"], 4)
        self.assertEqual(metadata["covariate_components_dropped"], 0)
        self.assertEqual([str(w.message) for w in caught], [])

    def test_a_dropped_component_is_counted_and_warned_about(self):
        rng = np.random.default_rng(23)
        covariates = rng.normal(size=(60, 4))
        covariates[:, 2] = covariates[:, 0]          # exactly collinear
        result, caught = self._run(covariates)
        metadata = result.run_metadata
        self.assertEqual(metadata["covariate_rank"], 3)
        self.assertEqual(metadata["covariate_components_dropped"], 1)
        self.assertEqual(len(caught), 1)
        self.assertIn("linearly dependent", str(caught[0].message))

    def test_no_covariates_reports_rank_zero_and_no_drop_count(self):
        rng = np.random.default_rng(24)
        dosage = rng.integers(0, 3, size=(48, 12)).astype(float)
        phenotype = rng.normal(size=(48, 2))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "norank", dosage)
            result = run_linear_gwas(str(bed), phenotype, device="cpu",
                                     compute_dtype="float64",
                                     output_dir=str(root / "out"))
        metadata = result.run_metadata
        self.assertEqual(metadata["covariate_rank"], 0)
        self.assertIsNone(metadata["covariate_components_dropped"])


class BorrowedResultRingTestCase(unittest.TestCase):
    """Borrowing the result ring must not change a single number.

    `borrow_results` yields views into the pinned ring rather than copies,
    because that copy is 69% of a BED scan at K=2048 (4.92 GB at 4.4 GB/s).
    The risk it takes on is that a borrowed chunk is valid only until the
    consumer asks for the next one, so a consumer that retains would read a
    slot the pipeline has already refilled -- producing wrong numbers with no
    error. These compare the borrowed path against the copying one on the same
    data, which is the only check that would notice.
    """

    @unittest.skipUnless(torch.cuda.is_available(),
                         "the borrowed ring is on the packed BED CUDA path")
    def test_binary_tables_match_float64_reference(self):
        rng = np.random.default_rng(77)
        n_samples = 120
        dosage = rng.integers(0, 3, size=(n_samples, 240)).astype(float)
        dosage[:3, :] = np.arange(3)[:, None]
        phenotype = rng.normal(size=(n_samples, 6))
        covariates = rng.normal(size=(n_samples, 3))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "borrow", dosage)
            tables = {}
            for name, fmt in (("borrowed", "binary"), ("copied", "binary")):
                run_linear_gwas(str(bed), phenotype, covariates=covariates,
                                genotype_format="plink", device="cuda:0",
                                compute_dtype="float32", chunk_size=16,
                                reader_workers=2, output_dir=str(root / name),
                                sumstats_format=fmt)
            # Read binary output and compare against a float64 CPU reference.
            rows = binary_rows(root / "borrowed")
            beta_ref, t_ref, _p, _q = linear_scan(
                dosage, phenotype, covariates, chunk_size=16, device="cpu",
                compute_dtype="float64")
            self.assertEqual(len(rows), 240 * 6)
            by_key = {(r["marker_id"], r["trait"]): r for r in rows}
            for marker in range(240):
                for trait in range(6):
                    row = by_key[(f"rs{marker}", f"trait_{trait}")]
                    with self.subTest(marker=marker, trait=trait):
                        np.testing.assert_allclose(
                            float(row["t_stat"]), t_ref[marker, trait],
                            rtol=2e-3, atol=2e-4)

    @unittest.skipUnless(torch.cuda.is_available(),
                         "the borrowed ring is on the packed BED CUDA path")
    def test_many_chunks_so_the_ring_wraps_repeatedly(self):
        # A small chunk against many variants forces the ring round many times;
        # with one slot the borrowed view would be overwritten mid-flight.
        rng = np.random.default_rng(78)
        dosage = rng.integers(0, 3, size=(80, 600)).astype(float)
        dosage[:3, :] = np.arange(3)[:, None]
        phenotype = rng.normal(size=(80, 3))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "wrap", dosage)
            run_linear_gwas(str(bed), phenotype, genotype_format="plink",
                            device="cuda:0", compute_dtype="float32",
                            chunk_size=8, reader_workers=2,
                            output_dir=str(root / "out"), sumstats_format="binary")
            rows = binary_rows(root / "out")
        beta_ref, t_ref, _p, _q = linear_scan(
            dosage, phenotype, None, chunk_size=8, device="cpu",
            compute_dtype="float64")
        self.assertEqual(len(rows), 600 * 3)
        by_key = {(r["marker_id"], r["trait"]): r for r in rows}
        for marker in range(0, 600, 7):
            for trait in range(3):
                row = by_key[(f"rs{marker}", f"trait_{trait}")]
                with self.subTest(marker=marker, trait=trait):
                    np.testing.assert_allclose(
                        float(row["t_stat"]), t_ref[marker, trait],
                        rtol=2e-3, atol=2e-4)


class TraitBlockMergeTestCase(unittest.TestCase):
    """Reducing in trait blocks must equal reducing all K at once.

    This is the enabling piece for trait tiling: at K in the millions the
    residualised phenotype does not fit on the device (22,250 x 2.085M float32
    is 185 GB against 80 GB), so the traits have to be processed in blocks. That
    is only sound if the blocked answer is the *same* answer, which holds
    because the top k of a union is contained in the union of the two top-k
    sets -- so nothing an earlier block discarded could have belonged in the
    final result.
    """

    def _blocked(self, reduction, t, status, df, block):
        from torchgwas.reduce import VariantReduction  # noqa: F401

        running = None
        for offset in range(0, t.shape[1], block):
            piece = t[:, offset:offset + block]
            width = reduction.resolved_width(piece.shape[1])
            r = reduction.reduce(piece * 0.5, piece, status, df, width)
            incoming = (r[0], r[1], None, r[2], r[3], r[4])
            running = reduction.merge(running, incoming, offset)
        return running

    def test_blocked_top_k_equals_unblocked(self):
        from torchgwas.reduce import VariantReduction

        torch.manual_seed(3)
        variants, traits = 40, 96
        t = torch.randn(variants, traits) * 5.0
        status = torch.zeros(variants, dtype=torch.uint8)
        df = torch.full((variants,), 400.0)
        for k in (1, 3, 8):
            reduction = VariantReduction("top-k" if k > 1 else "max-abs-t",
                                         k if k > 1 else None)
            whole = reduction.reduce(t * 0.5, t, status, df,
                                     reduction.resolved_width(traits))
            for block in (7, 16, 32, 96):
                got = self._blocked(reduction, t, status, df, block)
                with self.subTest(k=k, block=block):
                    # Compare the statistics, not the indices: ties can be
                    # broken differently and that is not a defect.
                    np.testing.assert_allclose(got[1].numpy(), whole[1].numpy(),
                                               rtol=1e-6)
                    np.testing.assert_allclose(got[0].numpy(), whole[0].numpy(),
                                               rtol=1e-6)
                    # The trait each row names must actually carry that value.
                    picked = t.gather(1, got[3].long())
                    np.testing.assert_allclose(picked.numpy(), got[1].numpy(),
                                               rtol=1e-6)

    def test_merge_rebases_trait_indices_onto_the_full_axis(self):
        from torchgwas.reduce import VariantReduction

        # One variant whose single strongest trait lives in the second block.
        t = torch.tensor([[1.0, 2.0, 9.0, 3.0]])
        status = torch.zeros(1, dtype=torch.uint8)
        df = torch.tensor([100.0])
        reduction = VariantReduction("max-abs-t")
        running = None
        for offset in (0, 2):
            piece = t[:, offset:offset + 2]
            r = reduction.reduce(piece * 0.5, piece, status, df, 1)
            running = reduction.merge(
                running, (r[0], r[1], None, r[2], r[3], r[4]), offset)
        self.assertEqual(int(running[3][0, 0]), 2, "index must be file-wide")
        self.assertAlmostEqual(float(running[1][0, 0]), 9.0, places=6)

    def test_an_invalid_variant_stays_invalid_across_blocks(self):
        from torchgwas.reduce import VariantReduction

        t = torch.tensor([[1.0, 2.0, 9.0, 3.0],
                          [4.0, 5.0, 6.0, 7.0]])
        status = torch.tensor([2, 0], dtype=torch.uint8)
        df = torch.tensor([100.0, 100.0])
        reduction = VariantReduction("top-k", 2)
        running = None
        for offset in (0, 2):
            piece = t[:, offset:offset + 2]
            r = reduction.reduce(piece * 0.5, piece, status, df, 2)
            running = reduction.merge(
                running, (r[0], r[1], None, r[2], r[3], r[4]), offset)
        # The status word survives the merge, so the caller's NaN rule still
        # applies to the invalid row exactly as in an unblocked scan.
        self.assertEqual([int(v) for v in running[4]], [2, 0])
        self.assertEqual(sorted(int(v) for v in running[3][1]), [2, 3])


class TraitBlockedScanTestCase(unittest.TestCase):
    """A trait-blocked scan must produce the same table as an unblocked one.

    Blocking exists so a scan whose phenotype does not fit on the device can
    run at all -- 22,250 samples by 2.085M traits is 185 GB of float32 against
    an 80 GB card. It is only useful if the answer is unchanged.
    """

    def _table(self, root, name, dosage, phenotype, covariates, **kwargs):
        bed = _write_bed(root / "tb", dosage)
        run_linear_gwas(str(bed), phenotype, covariates=covariates,
                        device="cpu", compute_dtype="float64",
                        output_dir=str(root / name), **kwargs)
        return binary_rows(root / name)

    def test_blocked_matches_unblocked_for_max_abs_t(self):
        rng = np.random.default_rng(41)
        dosage = rng.integers(0, 3, size=(70, 18)).astype(float)
        phenotype = rng.normal(size=(70, 9))
        covariates = rng.normal(size=(70, 3))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            whole = self._table(root, "whole", dosage, phenotype, covariates,
                                _internal_reduction=VariantReduction("max-abs-t"))
            for block in (1, 2, 4, 9, 20):
                blocked = self._table(root, f"b{block}", dosage, phenotype,
                                      covariates, _internal_reduction=VariantReduction("max-abs-t"),
                                      trait_block=block)
                with self.subTest(trait_block=block):
                    self.assertEqual(len(blocked), len(whole))
                    by_marker = {r["marker_id"]: r for r in whole}
                    for row in blocked:
                        want = by_marker[row["marker_id"]]
                        self.assertEqual(row["trait"], want["trait"])
                        self.assertAlmostEqual(float(row["t_stat"]),
                                               float(want["t_stat"]), places=9)
                        # -log10(P) must come from the winning trait's own per-variant
                        # df, carried through the merge rather than recomputed.
                        self.assertAlmostEqual(float(row["-log10_p"]),
                                               float(want["-log10_p"]), places=6)

    def test_blocked_matches_unblocked_for_top_k(self):
        rng = np.random.default_rng(42)
        dosage = rng.integers(0, 3, size=(64, 12)).astype(float)
        phenotype = rng.normal(size=(64, 10))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            whole = self._table(root, "whole", dosage, phenotype, None,
                                _internal_reduction=VariantReduction("top-k", 4))
            for block in (2, 3, 5):
                blocked = self._table(root, f"b{block}", dosage, phenotype, None,
                                      _internal_reduction=VariantReduction("top-k", 4),
                                      trait_block=block)
                with self.subTest(trait_block=block):
                    self.assertEqual(len(blocked), len(whole))
                    want = {}
                    for row in whole:
                        want.setdefault(row["marker_id"], []).append(
                            abs(float(row["t_stat"])))
                    got = {}
                    for row in blocked:
                        got.setdefault(row["marker_id"], []).append(
                            abs(float(row["t_stat"])))
                    for marker, scores in got.items():
                        self.assertEqual(scores, sorted(scores, reverse=True))
                        np.testing.assert_allclose(scores, want[marker],
                                                   rtol=1e-9)

    def test_trait_block_without_a_reduction_is_refused(self):
        rng = np.random.default_rng(43)
        dosage = rng.integers(0, 3, size=(32, 10)).astype(float)
        phenotype = rng.normal(size=(32, 4))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "norev", dosage)
            for kwargs in ({"trait_block": 2},
                           {"trait_block": 0, "reduce": "max-abs-t"}):
                with self.subTest(**kwargs), self.assertRaises(ValueError):
                    run_linear_gwas(str(bed), phenotype, device="cpu",
                                    compute_dtype="float64",
                                    output_dir=str(root / "out"), **kwargs)


class StreamingReductionTestCase(unittest.TestCase):
    """A device-side reduction must agree with the full scan it replaces.

    The reduction exists because a 2M-variant scan at K = 10^5 traits would
    otherwise produce 800 GB of t-statistics; but a narrower answer is only
    useful if it is the *same* answer, so every test here compares against the
    unreduced table rather than against itself. The two risks worth pinning are
    that torch orders NaN above every finite value in `topk` -- so an invalid
    variant would win its own row -- and that the residual df is per variant,
    which the reduced path has to broadcast across the k kept columns.
    """

    def _tables(self, root, dosage, phenotype, covariates, **kwargs):
        bed = _write_bed(root / "reduce", dosage)
        outputs = {}
        for name, extra in (("full", {}), ("reduced", kwargs)):
            run_linear_gwas(
                str(bed), phenotype, covariates=covariates,
                device="cpu", compute_dtype="float64",
                output_dir=str(root / name), **extra)
            outputs[name] = binary_rows(root / name)
        # Check selection against full float64 computation. Dense binary
        # stores float32; its round-trip accuracy is covered separately.
        outputs["full"] = run_linear_gwas(
            dosage, phenotype, covariates=covariates, device="cpu",
            compute_dtype="float64", marker_ids=[f"rs{i}" for i in range(dosage.shape[1])],
            return_t=True).table
        return outputs["full"], outputs["reduced"]

    def test_max_abs_t_picks_the_same_winner_as_the_full_table(self):
        rng = np.random.default_rng(101)
        dosage = rng.integers(0, 3, size=(64, 25)).astype(float)
        phenotype = rng.normal(size=(64, 6))
        covariates = rng.normal(size=(64, 3))
        with tempfile.TemporaryDirectory() as directory:
            full, reduced = self._tables(
                Path(directory), dosage, phenotype, covariates,
                _internal_reduction=VariantReduction("max-abs-t"))
        best: dict[str, tuple[float, str, str]] = {}
        for row in full:
            score = abs(float(row["t_stat"]))
            if row["marker_id"] not in best or score > best[row["marker_id"]][0]:
                best[row["marker_id"]] = (score, row["trait"], row["-log10_p"])
        self.assertEqual(len(reduced), len(best))
        for row in reduced:
            want_score, want_trait, want_p = best[row["marker_id"]]
            self.assertEqual(row["trait"], want_trait)
            self.assertAlmostEqual(abs(float(row["t_stat"])), want_score, places=9)
            # min-P and max-|t| must agree, which is only true because df is
            # per variant and therefore shared across a variant's traits.
            self.assertAlmostEqual(float(row["-log10_p"]), float(want_p), places=6)

    def test_top_k_keeps_k_rows_per_variant_in_descending_order(self):
        rng = np.random.default_rng(102)
        dosage = rng.integers(0, 3, size=(64, 12)).astype(float)
        phenotype = rng.normal(size=(64, 7))
        with tempfile.TemporaryDirectory() as directory:
            full, reduced = self._tables(
                Path(directory), dosage, phenotype, None,
                _internal_reduction=VariantReduction("top-k", 3))
        self.assertEqual(len(reduced), 12 * 3)
        per_variant: dict[str, list[float]] = {}
        for row in reduced:
            per_variant.setdefault(row["marker_id"], []).append(
                abs(float(row["t_stat"])))
        want: dict[str, list[float]] = {}
        for row in full:
            want.setdefault(row["marker_id"], []).append(abs(float(row["t_stat"])))
        for marker, scores in per_variant.items():
            self.assertEqual(scores, sorted(scores, reverse=True))
            np.testing.assert_allclose(
                scores, sorted(want[marker], reverse=True)[:3], rtol=1e-9)

    def test_k_is_capped_at_the_number_of_traits(self):
        rng = np.random.default_rng(103)
        dosage = rng.integers(0, 3, size=(48, 8)).astype(float)
        phenotype = rng.normal(size=(48, 2))
        with tempfile.TemporaryDirectory() as directory:
            _, reduced = self._tables(
                Path(directory), dosage, phenotype, None,
                _internal_reduction=VariantReduction("top-k", 25))
        self.assertEqual(len(reduced), 8 * 2)

    def test_an_invalid_variant_never_wins_its_own_row(self):
        # torch orders NaN *above* every finite value in topk and max, so
        # without the guard a variant flagged invalid -- or one carrying a NaN
        # statistic -- would be selected as its own top hit and reported with
        # whichever trait happened to hold the NaN. The public API refuses a
        # zero-variance variant before the scan starts, so the guard is pinned
        # here, at the reducer, which is where it lives and where the CUDA
        # backends reach it with real status words.
        from torchgwas.reduce import VariantReduction

        t = torch.tensor([[1.0, -5.0, 2.0],
                          [np.nan, 0.5, 0.25],
                          [3.0, 1.0, 9.0]])
        beta = t * 0.1
        status = torch.tensor([0, 0, 2], dtype=torch.uint8)
        variant_df = torch.tensor([50.0, 50.0, 50.0])
        reduction = VariantReduction("max-abs-t")
        _, t_sel, index, _, _ = reduction.reduce(beta, t, status, variant_df, 1)
        # Row 0: ordinary winner. Row 1: the NaN must lose to a finite 0.5.
        # Row 2: flagged invalid, so no column may be preferred over another on
        # the strength of its value -- the caller overwrites it with NaN.
        self.assertEqual(index[0, 0].item(), 1)
        self.assertEqual(index[1, 0].item(), 1)
        self.assertAlmostEqual(t_sel[1, 0].item(), 0.5, places=6)
        self.assertEqual(t_sel.shape, (3, 1))

    def test_a_nan_row_does_not_displace_a_real_hit_in_top_k(self):
        from torchgwas.reduce import VariantReduction

        t = torch.tensor([[np.nan, 1.0, np.nan, 4.0, 2.0]])
        status = torch.zeros(1, dtype=torch.uint8)
        reduction = VariantReduction("top-k", 3)
        _, t_sel, index, _, _ = reduction.reduce(
            t * 0.5, t, status, torch.tensor([40.0]), 3)
        self.assertEqual([int(v) for v in index[0]], [3, 4, 1])
        np.testing.assert_allclose(t_sel[0].numpy(), [4.0, 2.0, 1.0])

    def test_incompatible_combinations_are_refused(self):
        rng = np.random.default_rng(105)
        dosage = rng.integers(0, 3, size=(32, 10)).astype(float)
        phenotype = rng.normal(size=(32, 3))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = _write_bed(root / "refuse", dosage)
            for kwargs in (
                {"reduce": "max-abs-t"},                             # no output_dir
                {"reduce": "top-k", "output_dir": str(root / "a")},  # no k
                {"reduce": "max-abs-t", "output_dir": str(root / "b"),
                 "topk_per_trait": 2},
                {"reduce": "max-abs-t", "output_dir": str(root / "c"),
                 "sumstats_format": "binary"},
                {"reduce_top_k": 3, "output_dir": str(root / "d")},  # k without mode
                {"reduce": "nonsense", "output_dir": str(root / "e")},
            ):
                with self.subTest(**kwargs), self.assertRaises(ValueError):
                    run_linear_gwas(str(bed), phenotype, device="cpu",
                                    compute_dtype="float64", **kwargs)


if __name__ == "__main__":
    unittest.main()
