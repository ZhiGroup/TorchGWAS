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


def _write_plink_triplet(prefix: Path, dosage_a2: np.ndarray) -> tuple[Path, Path, Path]:
    dosage_a2 = np.asarray(dosage_a2, dtype=np.float32)
    n_samples, n_markers = dosage_a2.shape
    bed, bim, fam = resolve_plink_triplet(prefix)
    fam.write_text(
        "".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_samples)),
        encoding="utf-8",
    )
    bim.write_text(
        "".join(f"{1 + j} rs{j} 0 {100 + j} A{j % 4} C{j % 4}\n" for j in range(n_markers)),
        encoding="utf-8",
    )
    dosage_to_code = {2.0: 0b00, 1.0: 0b10, 0.0: 0b11}
    payload = bytearray(b"\x6c\x1b\x01")
    for marker in range(n_markers):
        for sample_start in range(0, n_samples, 4):
            byte = 0
            for offset in range(4):
                sample = sample_start + offset
                if sample >= n_samples or np.isnan(dosage_a2[sample, marker]):
                    code = 0b01
                else:
                    code = dosage_to_code[float(dosage_a2[sample, marker])]
                byte |= code << (2 * offset)
            payload.append(byte)
    bed.write_bytes(bytes(payload))
    return bed, bim, fam


class BedReaderTestCase(unittest.TestCase):
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
                self.assertEqual(genotype.shape, (3, 2))
                np.testing.assert_array_equal(
                    genotype.read_codes(0, 2),
                    np.asarray([[25, 0], [178, 128], [255, 255]], dtype=np.uint8),
                )
                np.testing.assert_allclose(
                    genotype.read_chunk(0, 2),
                    np.asarray([[25, 0], [178, 128], [255, 255]], dtype=np.float32) / 127.5,
                )
                self.assertEqual(genotype.marker_ids.tolist(), ["rs0", "rs5"])
                self.assertEqual(genotype.effect_alleles.tolist(), ["G", "C"])
                load_bgen_genotype(bgen, sample, cache_dir=cache)

            self.assertEqual(len(opened), 1)
            manifest_path = next(cache.glob("*.complete.json"))
            manifest = __import__("json").loads(manifest_path.read_text())
            self.assertEqual(manifest["zstd_level"], 15)
            self.assertEqual(manifest["chunk_size"], 2500)
            self.assertEqual(
                manifest["exclusion_counts"],
                {"invalid_probability": 1, "missing": 1, "multiallelic": 1, "non_diploid": 1, "phased": 1},
            )
            index = np.load(next(cache.glob("*.idx.npz")))
            self.assertTrue({"offs", "sizes", "nsnp", "nsamp", "chunk", "sub", "start"}.issubset(index.files))
            self.assertEqual(int(index["sub"]), 1)


if __name__ == "__main__":
    unittest.main()
