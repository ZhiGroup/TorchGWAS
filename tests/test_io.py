from __future__ import annotations

import threading
import time
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from torchgwas.bed import PlinkBedGenotype, resolve_plink_triplet
from torchgwas.io import load_bgen_genotype
from torchgwas.streaming import OrderedChunkLoader


def _write_plink_triplet(prefix: Path, dosage_a1: np.ndarray) -> tuple[Path, Path, Path]:
    dosage_a1 = np.asarray(dosage_a1, dtype=np.float32)
    n_samples, n_markers = dosage_a1.shape
    bed, bim, fam = resolve_plink_triplet(prefix)
    fam.write_text(
        "".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_samples)),
        encoding="utf-8",
    )
    bim.write_text(
        "".join(f"{1 + j} rs{j} 0 {100 + j} A{j % 4} C{j % 4}\n" for j in range(n_markers)),
        encoding="utf-8",
    )
    dosage_to_code = {0.0: 0b00, 1.0: 0b10, 2.0: 0b11}
    payload = bytearray(b"\x6c\x1b\x01")
    for marker in range(n_markers):
        for sample_start in range(0, n_samples, 4):
            byte = 0
            for offset in range(4):
                sample = sample_start + offset
                if sample >= n_samples or np.isnan(dosage_a1[sample, marker]):
                    code = 0b01
                else:
                    code = dosage_to_code[float(dosage_a1[sample, marker])]
                byte |= code << (2 * offset)
            payload.append(byte)
    bed.write_bytes(bytes(payload))
    return bed, bim, fam


class BedReaderTestCase(unittest.TestCase):
    def test_official_plink_documented_bytes_and_alleles(self):
        # Independent example: cog-genomics.org/plink/1.9/formats#bed.
        # Its PED has GG, AA, missing, AA, AA, AA for marker G/A.
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "official"
            Path(str(prefix)+".bed").write_bytes(bytes.fromhex("6c1b01dc0fe70f6b01"))
            Path(str(prefix)+".bim").write_text("1 snp1 0 1 G A\n1 snp2 0 2 1 2\n1 snp3 0 3 A C\n")
            Path(str(prefix)+".fam").write_text("".join(f"0 I{i} 0 0 0 -9\n" for i in range(6)))
            source = PlinkBedGenotype(prefix)
            np.testing.assert_allclose(source.read_chunk(0, 1)[:, 0],
                                       [0, 2, np.nan, 2, 2, 2], equal_nan=True)
            self.assertEqual(source.effect_alleles[0], "A")

    def test_decodes_bim_a2_dosage_and_metadata(self):
        expected = np.asarray(
            [
                [2, 0, 1],
                [1, 1, 0],
                [0, 2, np.nan],
                [2, 1, 2],
                [np.nan, 0, 1],
            ],
            dtype=np.float32,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "study.v1"
            bed, _, _ = _write_plink_triplet(prefix, expected)
            genotype = PlinkBedGenotype(bed, reader_workers=2, prefetch_chunks=2)
            observed = genotype.read_chunk(0, expected.shape[1])
            np.testing.assert_allclose(observed, expected, equal_nan=True)
            self.assertEqual(genotype.sample_ids.tolist(), ["I0", "I1", "I2", "I3", "I4"])
            self.assertEqual(genotype.marker_ids.tolist(), ["rs0", "rs1", "rs2"])
            self.assertEqual(genotype.effect_alleles.tolist(), ["C0", "C1", "C2"])

    def test_binary_bim_cache_avoids_reparsing_large_metadata(self):
        expected = np.asarray([[2, 0, 1], [1, 1, 0]], dtype=np.float32)
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bed, _, _ = _write_plink_triplet(root / "cached", expected)
            cache_dir = root / "metadata-cache"
            first = PlinkBedGenotype(bed, metadata_cache_dir=cache_dir)
            cache_path = first.metadata_cache_path
            self.assertTrue(cache_path.is_dir())
            self.assertTrue((cache_path / "manifest.json").is_file())

            original_read_csv = __import__("pandas").read_csv

            def reject_bim_parse(path, *args, **kwargs):
                if Path(path).suffix == ".bim":
                    raise AssertionError("BIM text should not be parsed when cache is valid")
                return original_read_csv(path, *args, **kwargs)

            with mock.patch("torchgwas.bed.pd.read_csv", side_effect=reject_bim_parse):
                second = PlinkBedGenotype(bed, metadata_cache_dir=cache_dir)

            self.assertEqual(second.marker_ids.tolist(), first.marker_ids.tolist())
            np.testing.assert_array_equal(second.positions, first.positions)

    def test_select_samples_decodes_in_requested_iid_order(self):
        expected = np.asarray(
            [
                [2, 0, 1],
                [1, 1, 0],
                [0, 2, 2],
                [2, 1, 2],
                [1, 0, 1],
                [0, 2, 0],
            ],
            dtype=np.float32,
        )
        selected_indices = np.asarray([5, 0, 3, 1])
        with tempfile.TemporaryDirectory() as tmpdir:
            bed, _, _ = _write_plink_triplet(Path(tmpdir) / "subset", expected)
            genotype = PlinkBedGenotype(bed).select_samples(
                [f"I{index}" for index in selected_indices]
            )
            observed = genotype.read_chunk(0, expected.shape[1])

        self.assertEqual(genotype.shape, (selected_indices.size, expected.shape[1]))
        self.assertEqual(genotype.sample_ids.tolist(), ["I5", "I0", "I3", "I1"])
        np.testing.assert_array_equal(observed, expected[selected_indices])

    def test_select_samples_rejects_missing_and_duplicate_iids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed, _, _ = _write_plink_triplet(
                Path(tmpdir) / "subset-errors",
                np.ones((6, 2), dtype=np.float32),
            )
            genotype = PlinkBedGenotype(bed)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                genotype.select_samples(["I0", "I0"])
            with self.assertRaisesRegex(ValueError, "absent"):
                genotype.select_samples(["I0", "not-present"])

    def test_dotted_prefix_is_not_truncated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "study.v1"
            bed, bim, fam = resolve_plink_triplet(prefix)
            self.assertEqual(bed.name, "study.v1.bed")
            self.assertEqual(bim.name, "study.v1.bim")
            self.assertEqual(fam.name, "study.v1.fam")

    def test_rejects_truncated_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "broken"
            bed, _, _ = _write_plink_triplet(prefix, np.ones((8, 4), dtype=np.float32))
            bed.write_bytes(bed.read_bytes()[:-1])
            with self.assertRaisesRegex(ValueError, "BED size mismatch"):
                PlinkBedGenotype(bed)

    def test_parallel_loader_preserves_order_and_propagates_errors(self):
        active = 0
        max_active = 0
        lock = threading.Lock()

        def reader(start: int, end: int, dtype: np.dtype) -> np.ndarray:
            nonlocal active, max_active
            with lock:
                active += 1
                max_active = max(max_active, active)
            time.sleep(0.02)
            with lock:
                active -= 1
            if start == 6:
                raise OSError("injected read failure")
            return np.full((3, end - start), start, dtype=dtype)

        loader = OrderedChunkLoader(
            n_markers=10,
            read_chunk=reader,
            chunk_size=2,
            dtype=np.float32,
            prefetch_chunks=4,
            reader_workers=4,
        )
        observed_starts = []
        with self.assertRaisesRegex(OSError, "injected read failure"):
            for start, _, _ in loader:
                observed_starts.append(start)
        self.assertEqual(observed_starts, [0, 2, 4])
        self.assertGreaterEqual(max_active, 2)

    def test_bgen_cache_quantizes_dosage_filters_variants_and_reuses_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bgen = root / "cohort.bgen"
            sample = root / "cohort.sample"
            cache = root / "cache"
            bgen.write_bytes(b"bgen")
            sample.write_text(
                "ID_1 ID_2 missing\n0 0 0\nF0 I0 0\nF1 I1 0\nF2 I2 0\n"
            )

            class FakeBgen:
                nvariants = 7
                nsamples = 3
                samples = np.asarray(["I0", "I1", "I2"])
                nalleles = np.asarray([2, 3, 2, 2, 2, 2, 2])
                phased = np.asarray([False, False, True, False, False, False, False])
                ids = np.asarray([f"id{j}" for j in range(7)])
                rsids = np.asarray(["rs0", ".", "rs2", "rs3", "rs4", "rs5", "rs6"])
                chromosomes = np.asarray(["1"] * 7)
                positions = np.arange(100, 107)
                allele_ids = np.asarray(["A,G", "A,C,G", "A,T", "C,T", "G,A", "T,C", "A,C"])

                def __init__(self):
                    self.closed = False
                    self.probabilities = np.zeros((3, 7, 3), dtype=np.float32)
                    self.probabilities[:, 0, :] = np.asarray(
                        [[0.8, 0.2, 0.0], [0.1, 0.4, 0.5], [0.0, 0.0, 1.0]]
                    )
                    self.probabilities[:, 3, :] = np.asarray([1.0, 0.0, 0.0])
                    self.probabilities[:, 4, :] = np.asarray([0.0, 1.0, 0.0])
                    self.probabilities[:, 5, :] = np.asarray(
                        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                    )
                    self.probabilities[:, 6, :] = np.asarray([0.2, 0.2, 0.2])
                    self.missing = np.zeros((3, 7), dtype=bool)
                    self.missing[1, 3] = True
                    self.ploidy = np.full((3, 7), 2, dtype=np.int8)
                    self.ploidy[2, 4] = 1

                def read(self, index, **kwargs):
                    indices = np.asarray(index[1])
                    return (
                        self.probabilities[:, indices, :].astype(kwargs["dtype"]),
                        self.missing[:, indices],
                        self.ploidy[:, indices],
                    )

                def close(self):
                    self.closed = True

            opened: list[FakeBgen] = []

            def fake_open(*args, **kwargs):
                instance = FakeBgen()
                opened.append(instance)
                return instance

            with mock.patch("torchgwas.bgen._open_bgen", side_effect=fake_open):
                genotype, _, _ = load_bgen_genotype(
                    bgen,
                    sample,
                    cache_dir=cache,
                    reader_workers=2,
                )
                # Variant 3 was excluded only because sample 1 is missing.
                # It is now masked and kept: its observed probabilities are
                # [1, 0, 0], so every dosage including the masked one is 0.
                self.assertEqual(genotype.shape, (3, 3))
                expected = np.asarray(
                    [[25, 0, 0], [178, 0, 128], [255, 0, 255]], dtype=np.uint8)
                np.testing.assert_array_equal(genotype.read_codes(0, 3), expected)
                np.testing.assert_allclose(
                    genotype.read_chunk(0, 3),
                    expected.astype(np.float32) / 127.5,
                )
                # rs3 is kept now: missingness was its only disqualifier.
                self.assertEqual(genotype.marker_ids.tolist(),
                                 ["rs0", "rs3", "rs5"])
                # rs3 has allele_ids "C,T", so its effect allele is T.
                self.assertEqual(genotype.effect_alleles.tolist(), ["G", "T", "C"])
                load_bgen_genotype(bgen, sample, cache_dir=cache)

            self.assertEqual(len(opened), 1)
            manifest_path = next(cache.glob("*.complete.json"))
            manifest = __import__("json").loads(manifest_path.read_text())
            self.assertEqual(manifest["zstd_level"], 15)
            self.assertEqual(manifest["chunk_size"], 2500)
            self.assertEqual(
                manifest["exclusion_counts"],
                {"invalid_probability": 1, "masked_missing": 1, "multiallelic": 1,
                 "non_diploid": 1, "phased": 1},
            )
            index = np.load(next(cache.glob("*.idx.npz")))
            self.assertTrue({"offs", "sizes", "nsnp", "nsamp", "chunk", "sub", "start"}.issubset(index.files))
            self.assertEqual(int(index["sub"]), 1)


if __name__ == "__main__":
    unittest.main()
