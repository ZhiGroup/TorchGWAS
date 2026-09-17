"""The store is a transport, so the only thing that matters is that it is exact.

It holds the same bytes a `.bed` holds and hands them back unchanged. Every
test here is therefore about identity, not about ratio -- a store that
compresses well and returns different genotypes is worthless, and the failure
would surface far downstream as a subtly wrong association.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from torchgwas.hardcall_store import (DEFAULT_FRAME_VARIANTS, HardcallStore,
                                      encode_hardcall_store)


def synthetic_packed(variants: int, samples: int, seed: int = 0) -> np.ndarray:
    """Packed two-bit rows with realistic structure.

    Not uniform random: real genotypes are dominated by one code and
    neighbouring variants are correlated, which is the structure the
    compressor exploits and the reason genome order beats MAF order. A
    uniform-random fixture would compress at 1.0x and so would not exercise
    the same code paths in zstd at all.
    """
    rng = np.random.default_rng(seed)
    bytes_per_variant = (samples + 3) // 4
    calls = np.zeros((variants, samples), dtype=np.uint8)
    for begin in range(0, variants, 64):
        # One frequency per block of neighbours, so a block is internally
        # correlated the way a linkage block is.
        freq = rng.uniform(0.001, 0.5)
        stop = min(begin + 64, variants)
        draws = rng.random((stop - begin, samples))
        block = np.where(draws < freq ** 2, 0,
                         np.where(draws < freq ** 2 + 2 * freq * (1 - freq),
                                  2, 3)).astype(np.uint8)
        # A sprinkling of the missing code, which is 1 in bed order.
        block[rng.random(block.shape) < 0.001] = 1
        calls[begin:stop] = block

    padded = bytes_per_variant * 4
    if padded != samples:
        calls = np.pad(calls, ((0, 0), (0, padded - samples)))
    quads = calls.reshape(variants, bytes_per_variant, 4)
    return (quads[:, :, 0] | (quads[:, :, 1] << 2)
            | (quads[:, :, 2] << 4) | (quads[:, :, 3] << 6)).astype(np.uint8)


class RoundTripTests(unittest.TestCase):
    VARIANTS = 5000
    SAMPLES = 1301           # deliberately not a multiple of 4

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.packed = synthetic_packed(self.VARIANTS, self.SAMPLES)
        self.prefix = Path(self.tmp.name) / "store"
        self.manifest = self.encode()

    def encode(self, **over):
        options = dict(frame_variants=512, workers=4)
        options.update(over)
        return encode_hardcall_store(
            lambda a, b: self.packed[a:b], self.VARIANTS, self.SAMPLES,
            self.prefix, **options)

    def test_every_byte_survives_the_round_trip(self):
        with HardcallStore(self.prefix) as store:
            got = store.read_packed(0, self.VARIANTS)
        np.testing.assert_array_equal(got, self.packed)

    def test_partial_reads_match_the_same_slice(self):
        """A chunk boundary need not align to a frame boundary.

        The reader decompresses the covering frames and slices, and getting
        the slice arithmetic wrong would shift genotypes by a few variants --
        which produces perfectly plausible association results for the wrong
        markers.
        """
        with HardcallStore(self.prefix) as store:
            for start, end in ((0, 1), (0, 512), (1, 513), (511, 1025),
                               (512, 1024), (1500, 1500), (4999, 5000),
                               (137, 4321), (0, self.VARIANTS)):
                np.testing.assert_array_equal(
                    store.read_packed(start, end), self.packed[start:end],
                    err_msg=f"[{start}, {end})")

    def test_reads_into_a_caller_owned_buffer(self):
        """The pinned-ring path writes into memory it already owns."""
        with HardcallStore(self.prefix) as store:
            target = np.zeros((600, store.bytes_per_variant), dtype=np.uint8)
            written = store.read_packed_into(target, 100, 700)
        self.assertEqual(written, 600 * ((self.SAMPLES + 3) // 4))
        np.testing.assert_array_equal(target, self.packed[100:700])

    def test_a_short_buffer_is_refused_rather_than_truncated(self):
        with HardcallStore(self.prefix) as store:
            target = np.zeros((10, store.bytes_per_variant), dtype=np.uint8)
            with self.assertRaises(ValueError):
                store.read_packed_into(target, 0, 600)

    def test_out_of_range_reads_raise(self):
        with HardcallStore(self.prefix) as store:
            for start, end in ((-1, 10), (0, self.VARIANTS + 1), (10, 5)):
                with self.assertRaises(IndexError):
                    store.read_packed(start, end)

    def test_the_file_does_not_depend_on_the_worker_count(self):
        """Frames are written in variant order however they were compressed.

        A store whose bytes depend on the thread count cannot be checksummed
        against a reference, which removes the cheapest correctness check the
        format has.
        """
        first = Path(f"{self.prefix}.zhc").read_bytes()
        for workers in (1, 2, 8):
            self.encode(workers=workers)
            self.assertEqual(Path(f"{self.prefix}.zhc").read_bytes(), first,
                             f"workers={workers} produced a different file")

    def test_a_trailing_partial_frame_round_trips(self):
        """5,000 variants in frames of 512 leaves a final frame of 392."""
        self.assertEqual(self.manifest["frames"], 10)
        self.assertEqual(self.manifest["frame_counts"][-1], 5000 - 9 * 512)
        with HardcallStore(self.prefix) as store:
            np.testing.assert_array_equal(
                store.read_packed(9 * 512, 5000), self.packed[9 * 512:])

    def test_samples_not_a_multiple_of_four_keep_their_stride(self):
        self.assertEqual(self.manifest["bytes_per_variant"],
                         (self.SAMPLES + 3) // 4)

    def test_the_manifest_reports_a_real_ratio(self):
        """Structured genotypes must actually compress.

        Real hard calls measured 6.4x at this level. The fixture is smaller and
        synthetic so it will not reach that, but anything at or below 1.0 means
        the fixture lost its structure and the other tests are checking a
        round trip through an effectively uncompressed file.
        """
        self.assertGreater(self.manifest["ratio"], 1.5)
        self.assertEqual(self.manifest["raw_bytes"],
                         self.VARIANTS * ((self.SAMPLES + 3) // 4))

    def test_a_foreign_manifest_is_refused(self):
        import json
        path = Path(f"{self.prefix}.zhc.json")
        manifest = json.loads(path.read_text())
        manifest["magic"] = "SOMETHINGELSE"
        path.write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            HardcallStore(self.prefix)


class FrameSizeTests(unittest.TestCase):
    """Frame size must not change what comes back, only how it is packed."""

    def test_every_frame_size_returns_identical_genotypes(self):
        packed = synthetic_packed(2048, 401, seed=7)
        with tempfile.TemporaryDirectory() as tmp:
            for frame in (64, 256, 1024, 2048, 4096, DEFAULT_FRAME_VARIANTS):
                prefix = Path(tmp) / f"store{frame}"
                encode_hardcall_store(lambda a, b: packed[a:b], 2048, 401,
                                      prefix, frame_variants=frame, workers=2)
                with HardcallStore(prefix) as store:
                    np.testing.assert_array_equal(
                        store.read_packed(0, 2048), packed,
                        err_msg=f"frame_variants={frame}")
                    # And a read that straddles frames at this size.
                    np.testing.assert_array_equal(
                        store.read_packed(33, 1777), packed[33:1777],
                        err_msg=f"frame_variants={frame} straddling")

    def test_a_frame_larger_than_the_file_is_one_frame(self):
        packed = synthetic_packed(100, 37, seed=3)
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "store"
            manifest = encode_hardcall_store(
                lambda a, b: packed[a:b], 100, 37, prefix,
                frame_variants=100_000, workers=2)
            self.assertEqual(manifest["frames"], 1)
            with HardcallStore(prefix) as store:
                np.testing.assert_array_equal(store.read_packed(0, 100), packed)


class BedSubstitutionTests(unittest.TestCase):
    """A store standing in for a `.bed` must be invisible to everything above.

    This is the property that makes the store safe to adopt: it is a transport
    swap, not a new genotype path. If decoded dosages differ by so much as one
    call, the substitution has introduced a correctness difference and the
    7.6x is worthless.
    """

    VARIANTS = 1500
    SAMPLES = 211

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        self.packed = synthetic_packed(self.VARIANTS, self.SAMPLES, seed=11)
        self.prefix = root / "cohort"

        # A real PLINK triplet, so the reader's own metadata path runs.
        (root / "cohort.fam").write_text("".join(
            f"FAM{i} IID{i} 0 0 1 -9\n" for i in range(self.SAMPLES)))
        (root / "cohort.bim").write_text("".join(
            f"1 rs{i} 0 {1000 + i} A G\n" for i in range(self.VARIANTS)))
        with (root / "cohort.bed").open("wb") as handle:
            handle.write(b"\x6c\x1b\x01")
            handle.write(self.packed.tobytes())

        encode_hardcall_store(lambda a, b: self.packed[a:b], self.VARIANTS,
                              self.SAMPLES, root / "cohort",
                              frame_variants=256, workers=3)

    def open_both(self):
        from torchgwas.bed import PlinkBedGenotype
        plain = PlinkBedGenotype(self.prefix)
        stored = PlinkBedGenotype(self.prefix, hardcall_store=self.prefix)
        self.addCleanup(plain.__del__)
        self.addCleanup(stored.__del__)
        return plain, stored

    def test_decoded_dosages_are_identical(self):
        plain, stored = self.open_both()
        for start, end in ((0, self.VARIANTS), (0, 256), (7, 1103),
                           (1499, 1500)):
            np.testing.assert_array_equal(
                stored.read_chunk(start, end), plain.read_chunk(start, end),
                err_msg=f"[{start}, {end})")

    def test_nan_for_missing_survives_the_store(self):
        """Missing calls decode to NaN, which array_equal would not catch."""
        plain, stored = self.open_both()
        a = plain.read_chunk(0, self.VARIANTS)
        b = stored.read_chunk(0, self.VARIANTS)
        self.assertTrue(np.isnan(a).any(), "fixture has no missing calls")
        np.testing.assert_array_equal(np.isnan(a), np.isnan(b))

    def test_metadata_comes_from_the_triplet_either_way(self):
        plain, stored = self.open_both()
        self.assertEqual(stored.shape, plain.shape)
        np.testing.assert_array_equal(stored.marker_ids, plain.marker_ids)
        np.testing.assert_array_equal(stored.sample_ids, plain.sample_ids)
        np.testing.assert_array_equal(stored.positions, plain.positions)

    def test_the_bed_need_not_exist_when_a_store_is_given(self):
        """Otherwise the store saves no disk, which is half its purpose."""
        from torchgwas.bed import PlinkBedGenotype
        Path(f"{self.prefix}.bed").unlink()
        genotype = PlinkBedGenotype(self.prefix, hardcall_store=self.prefix)
        self.addCleanup(genotype.__del__)
        np.testing.assert_array_equal(
            genotype.read_chunk(100, 400)[:, :],
            _decode(self.packed[100:400], self.SAMPLES))

    def test_a_store_from_a_different_cohort_is_refused(self):
        """Shapes that disagree would read plausible genotypes for wrong samples."""
        from torchgwas.bed import PlinkBedGenotype
        other = Path(self.tmp.name) / "other"
        encode_hardcall_store(
            lambda a, b: synthetic_packed(self.VARIANTS, self.SAMPLES)[a:b],
            self.VARIANTS - 10, self.SAMPLES, other,
            frame_variants=256, workers=2)
        with self.assertRaises(ValueError):
            PlinkBedGenotype(self.prefix, hardcall_store=other)


class LoaderRefusalTests(unittest.TestCase):
    """A store handed to a format that cannot use it must RAISE, not be dropped.

    Only the PLINK path substitutes a store for its genotype bytes. Silently
    ignoring the argument elsewhere would report a "store" measurement taken
    without the store -- a benchmark that lies in the flattering direction.
    This module has been bitten by exactly that before: `selected_sample_ids`
    was once accepted and dropped, and a sweep asking for 4,000 of 35,365
    samples got full-cohort scans at every requested size.
    """

    def test_a_store_on_a_pgen_load_is_refused(self):
        from torchgwas.io import load_genotype
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp) / "cohort.pgen"
            fake.write_bytes(b"\x6c\x1b\x10" + b"\0" * 64)
            with self.assertRaises(ValueError) as caught:
                load_genotype(fake, genotype_format="pgen",
                              hardcall_store=str(Path(tmp) / "store"))
            self.assertIn("PLINK", str(caught.exception))

    def test_the_message_names_the_offending_format(self):
        from torchgwas.io import load_genotype
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp) / "cohort.bgen"
            fake.write_bytes(b"\0" * 64)
            with self.assertRaises(ValueError) as caught:
                load_genotype(fake, genotype_format="bgen",
                              hardcall_store=str(Path(tmp) / "store"))
            self.assertIn("bgen", str(caught.exception))


def _decode(packed: np.ndarray, samples: int) -> np.ndarray:
    from torchgwas.bed import _A2_DOSAGE_LUT
    wide = _A2_DOSAGE_LUT[packed].reshape(packed.shape[0], -1)[:, :samples]
    return np.asarray(wide.T, dtype=np.float32, order="C")


if __name__ == "__main__":
    unittest.main()
