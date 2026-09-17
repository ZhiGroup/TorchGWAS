from concurrent.futures import ThreadPoolExecutor
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import torch
import zstandard as zstd
from torchgwas.zstd_native import ZstdInto, contiguous_frame_batches


class ZstdNativeTests(unittest.TestCase):
    def test_concurrent_direct_into_and_corruption(self):
        decoder = ZstdInto()
        expected = [bytes([i]) * (1000 + i) for i in range(12)]
        def one(raw):
            blob = bytearray(zstd.ZstdCompressor().compress(raw))
            out = bytearray(len(raw))
            decoder.decompress_into(blob, out)
            return bytes(out)
        with ThreadPoolExecutor(4) as pool:
            self.assertEqual(list(pool.map(one, expected)), expected)
        blob = zstd.ZstdCompressor().compress(expected[0])
        for invalid in (blob[:-1], blob + b'x', blob + blob):
            with self.assertRaises(ValueError):
                decoder.decompress_into(invalid, bytearray(1000))
        with self.assertRaises(ValueError):
            decoder.decompress_into(blob, bytearray(1001))
        with self.assertRaises(ValueError):
            decoder.decompress_into(blob, bytes(1000))

    def test_contiguous_read_order_short_read_and_early_close(self):
        blobs = [zstd.ZstdCompressor().compress(bytes([i]) * 1000) for i in range(9)]
        offsets = np.cumsum([0] + [len(x) for x in blobs[:-1]])
        sizes = [len(x) for x in blobs]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'data.zst'
            path.write_bytes(b''.join(blobs))
            batches = list(contiguous_frame_batches(path, offsets, sizes, target_bytes=50, prefetch_batches=4, read_workers=4))
            flat = [item for batch in batches for item in batch]
            self.assertEqual([i for i, _ in flat], list(range(9)))
            self.assertEqual([bytes(blob) for _, blob in flat], blobs)
            iterator = contiguous_frame_batches(path, offsets, sizes, target_bytes=50)
            next(iterator)
            iterator.close()
            path.write_bytes(path.read_bytes()[:-1])
            with self.assertRaises(OSError):
                list(contiguous_frame_batches(path, offsets, sizes))

    def test_public_reader_reblocks_once_and_returns_owned_chunks(self):
        from torchgwas.zstd_store import ZstdGenotype, encode_zstd_store
        raw = np.arange(35, dtype=np.uint8).reshape(5, 7)
        class Source:
            shape = raw.shape
            sample_ids = np.asarray(['a', 'b', 'c', 'd', 'e'])
            marker_ids = np.asarray([str(i) for i in range(7)])
            dosage_scale = 127.5
            variant_metadata = dict(chromosome=['1'] * 7, position=list(range(7)),
                                    effect_allele=['A'] * 7, other_allele=['G'] * 7)
            def iter_chunks(self, chunk_size, **kwargs):
                for start in range(0, 7, chunk_size):
                    end = min(start + chunk_size, 7)
                    yield start, end, raw[:, start:end]
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / 'public'
            encode_zstd_store(Source(), prefix, chunk_size=3, compression_workers=2)
            reader = ZstdGenotype(prefix, read_batch_bytes=30, read_workers=4, read_ahead_batches=4)
            try:
                for chunk_size in (1, 2, 3, 5, 20):
                    chunks = list(reader.iter_chunks(chunk_size, dtype=np.float32,
                                                     prefetch_chunks=2, reader_workers=2))
                    np.testing.assert_allclose(np.concatenate([x[2] for x in chunks], axis=1), raw / 127.5)
                self.assertTrue(reader.supports_fused_qc)
                self.assertEqual(reader.native_scale, 127.5)
                self.assertEqual(reader.native_dtype, np.uint8)
                aligned = list(reader.iter_native_chunks(3))
                self.assertTrue(all(x[2].T.flags.c_contiguous for x in aligned))
                np.testing.assert_array_equal(np.concatenate([x[2] for x in aligned], axis=1), raw)
                codes = list(reader.iter_native_chunks(2))
                self.assertTrue(all(x[2].T.flags.c_contiguous for x in codes))
                np.testing.assert_array_equal(np.concatenate([x[2] for x in codes], axis=1), raw)
                np.testing.assert_array_equal(reader.read_codes(1, 6), raw[:, 1:6])
            finally:
                reader.close()

    def test_scaled_reads_are_bit_identical_to_dividing_afterwards(self):
        """The fused widen-and-divide must not move a single bit.

        Widening the stored codes into the output block and scaling them used
        to be two separate sweeps; they are now one `np.divide(..., out=)`.
        That is only a safe change if the arithmetic is identical, which is why
        this asserts exact equality rather than a tolerance -- and why the
        implementation divides by `dosage_scale` instead of multiplying by its
        reciprocal, since the reciprocal of a scale like 127.5 or 254 is not
        exact in binary and would shift the last bits of every dosage.
        """
        from torchgwas.zstd_store import ZstdGenotype, encode_zstd_store

        rng = np.random.default_rng(17)
        raw = rng.integers(0, 256, size=(9, 20)).astype(np.uint8)
        for scale in (127.5, 254.0, 100.0, 1.0):
            class Source:
                shape = raw.shape
                sample_ids = np.asarray([f's{i}' for i in range(9)])
                marker_ids = np.asarray([str(i) for i in range(20)])
                dosage_scale = scale
                variant_metadata = dict(chromosome=['1'] * 20,
                                        position=list(range(20)),
                                        effect_allele=['A'] * 20,
                                        other_allele=['G'] * 20)

                def iter_chunks(self, chunk_size, **kwargs):
                    for start in range(0, 20, chunk_size):
                        end = min(start + chunk_size, 20)
                        yield start, end, raw[:, start:end]

            with tempfile.TemporaryDirectory() as tmp:
                prefix = Path(tmp) / 'scaled'
                encode_zstd_store(Source(), prefix, chunk_size=4,
                                  compression_workers=2)
                reader = ZstdGenotype(prefix, read_batch_bytes=64,
                                      read_workers=2, read_ahead_batches=2)
                try:
                    # Chunk sizes that do and do not divide the 4-row frames,
                    # so both the frame-aligned and the straddling copies run.
                    for chunk_size in (1, 3, 4, 6, 20):
                        chunks = list(reader.iter_chunks(
                            chunk_size, dtype=np.float64, prefetch_chunks=2,
                            reader_workers=2))
                        got = np.concatenate([x[2] for x in chunks], axis=1)
                        want = (raw.astype(np.float64) / scale if scale != 1.0
                                else raw.astype(np.float64))
                        with self.subTest(scale=scale, chunk_size=chunk_size):
                            np.testing.assert_array_equal(got, want)
                finally:
                    reader.close()

    def test_every_chunk_entry_point_accepts_a_variant_range(self):
        """A range must reach every loader, not just the ones already checked.

        `variant_range` has now been missed twice. It was added to `iter_chunks`
        on five sources and missed on `ZstdGenotype.iter_native_chunks`, which
        is the path `PinnedDosageLoader` actually prefers -- so the CUDA scan
        took the one route that lacked it. The failure read
        "ZstdGenotype cannot read a variant range, so it cannot be sharded
        across devices", which sounds like a considered limitation rather than a
        missing keyword.

        Checking the signatures directly is the only version of this test that
        catches the *next* one, since a path nobody exercises is exactly where
        this hides.
        """
        import inspect

        from torchgwas.bed import PlinkBedGenotype
        from torchgwas.bgen import BgenGenotype
        from torchgwas.pgen import PgenGenotype
        from torchgwas.zstd_store import ZstdGenotype

        wanted = ("iter_chunks", "iter_native_chunks", "iter_device_chunks",
                  "iter_packed_chunks")
        checked = 0
        for source in (PlinkBedGenotype, BgenGenotype, PgenGenotype, ZstdGenotype):
            for name in wanted:
                function = getattr(source, name, None)
                if function is None:
                    continue
                checked += 1
                parameters = inspect.signature(function).parameters
                with self.subTest(source=source.__name__, entry=name):
                    self.assertIn(
                        "variant_range", parameters,
                        f"{source.__name__}.{name} cannot be given a variant "
                        f"range, so a ranged or sharded scan through it fails")
        self.assertGreaterEqual(checked, 6, "expected at least six entry points")

    def test_decode_buffer_ring_survives_wrapping_under_concurrency(self):
        """The reused decompression buffers must wrap without aliasing.

        Buffer `n % depth` is rewritten by frame `n` while frames
        `n-depth+1 .. n-1` are still in flight, so the ring is only correct
        because the generator does not resume until its consumer has copied the
        value it was handed. That argument is worth a test with enough frames
        to wrap the ring many times and enough workers to overlap: the small
        fixtures elsewhere would pass even if the ring were one buffer.
        """
        from torchgwas.zstd_store import ZstdGenotype, encode_zstd_store

        rng = np.random.default_rng(31)
        n_samples, n_variants, scale = 23, 400, 254.0
        raw = rng.integers(0, 256, size=(n_samples, n_variants)).astype(np.uint8)

        class Source:
            shape = raw.shape
            sample_ids = np.asarray([f's{i}' for i in range(n_samples)])
            marker_ids = np.asarray([str(i) for i in range(n_variants)])
            dosage_scale = scale
            variant_metadata = dict(chromosome=['1'] * n_variants,
                                    position=list(range(n_variants)),
                                    effect_allele=['A'] * n_variants,
                                    other_allele=['G'] * n_variants)

            def iter_chunks(self, chunk_size, **kwargs):
                for start in range(0, n_variants, chunk_size):
                    end = min(start + chunk_size, n_variants)
                    yield start, end, raw[:, start:end]

        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / 'ring'
            # 5 rows per frame over 400 variants is 80 frames, so a depth of 3
            # wraps the ring 26 times.
            encode_zstd_store(Source(), prefix, chunk_size=5, compression_workers=2)
            reader = ZstdGenotype(prefix, read_batch_bytes=128, read_workers=3,
                                  read_ahead_batches=3)
            try:
                want = raw.astype(np.float64) / scale
                for chunk_size, depth, workers in ((7, 3, 4), (5, 2, 2),
                                                   (64, 4, 8), (13, 8, 3)):
                    chunks = list(reader.iter_chunks(
                        chunk_size, dtype=np.float64, prefetch_chunks=depth,
                        reader_workers=workers))
                    got = np.concatenate([c[2] for c in chunks], axis=1)
                    with self.subTest(chunk_size=chunk_size, depth=depth,
                                      workers=workers):
                        np.testing.assert_array_equal(got, want)
                # And a ranged read, where frames are trimmed at both ends.
                for lo, hi in ((3, 97), (0, 5), (395, 400), (128, 129)):
                    chunks = list(reader.iter_chunks(
                        11, dtype=np.float64, prefetch_chunks=3,
                        reader_workers=4, variant_range=(lo, hi)))
                    got = np.concatenate([c[2] for c in chunks], axis=1)
                    with self.subTest(variant_range=(lo, hi)):
                        np.testing.assert_array_equal(got, want[:, lo:hi])
            finally:
                reader.close()

    def test_lowrank_prefix_subframes_ring_and_close(self):
        file = Path(__file__).resolve().parents[1] / 'lowrank' / 'zstore.py'
        spec = importlib.util.spec_from_file_location('lowrank_zstore_test', file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        raw = np.arange(7 * 5, dtype=np.uint8).reshape(7, 5)
        # sub > final chunk rows exercises omitted empty subframes.
        chunk, sub = 3, 4
        frames = []
        for start in range(0, len(raw), chunk):
            block = raw[start:start + chunk]
            for k in range(sub):
                lo, hi = len(block) * k // sub, len(block) * (k + 1) // sub
                if hi > lo:
                    frames.append(zstd.ZstdCompressor().compress(block[lo:hi].tobytes()))
        original_empty = torch.empty
        def pageable(*args, **kwargs):
            kwargs.pop('pin_memory', None)
            return original_empty(*args, **kwargs)
        with tempfile.TemporaryDirectory() as tmp, patch.object(module.torch, 'empty', side_effect=pageable):
            prefix = str(Path(tmp) / 'store')
            Path(prefix + '.zst').write_bytes(b''.join(frames))
            np.savez(prefix + '.idx.npz', offs=np.cumsum([0] + [len(x) for x in frames[:-1]]),
                     sizes=[len(x) for x in frames], nsnp=7, nsamp=5, chunk=chunk, sub=sub)
            reader = module.ZReader(prefix, 5, 5, threads=3, depth=2, read_batch_bytes=30, read_workers=4, read_ahead_batches=4)
            got = []
            try:
                for bi, host, start, n in reader:
                    got.append(host.numpy().copy())
                    reader.release(bi)
                    with self.assertRaises(ValueError):
                        reader.release(bi)
            finally:
                reader.close()
            np.testing.assert_array_equal(np.concatenate(got), raw[:5])
            reader = module.ZReader(prefix, 7, 5, depth=1)
            iterator = iter(reader)
            next(iterator)
            reader.close()  # unreturned slot cannot deadlock producer shutdown
            iterator.close()
            self.assertFalse(reader.producer.is_alive())


if __name__ == '__main__':
    unittest.main()






class ZstdDirectFillTestCase(unittest.TestCase):
    """The pinned loader must fill from the store without an intermediate copy.

    Without `allows_direct_native_fill`, `PinnedDosageLoader` takes its general
    path, which is `np.copyto(buffer, array.T)` -- every decoded byte copied
    into pinned memory and transposed on the way. On the benchmark cohort that
    is 198.7 GB of strided host copying for one pass, and removing it measured
    **2.35x on the whole zstd scan** (33.6 -> 14.3 s, full file, K=128, cold,
    three rounds, with bitwise-identical t-statistics).

    These assert the two things that make it safe: the flag is actually set, so
    the fast path is the one taken, and the bytes the direct fill produces are
    the same ones `iter_chunks` produces.
    """

    def _store(self, tmpdir):
        from torchgwas.zstd_store import ZstdGenotype, encode_zstd_store

        rng = np.random.default_rng(20260913)
        samples, variants = 24, 700
        codes = rng.integers(0, 255, size=(samples, variants), dtype=np.uint8)

        class Source:
            shape = (samples, variants)
            sample_ids = np.array([f"s{i}" for i in range(samples)])
            family_ids = sample_ids
            marker_ids = np.array([f"rs{i}" for i in range(variants)])
            dosage_scale = 127.5
            variant_metadata = {
                "chromosome": np.full(variants, 1),
                "position": np.arange(variants),
                "effect_allele": np.full(variants, "A"),
                "other_allele": np.full(variants, "C"),
            }

            def iter_chunks(self, chunk_size, dtype=np.uint8, **kwargs):
                for start in range(0, variants, chunk_size):
                    end = min(variants, start + chunk_size)
                    yield start, end, codes[:, start:end]

        prefix = Path(tmpdir) / "store"
        encode_zstd_store(Source(), prefix, chunk_size=128,
                          compression_workers=2)
        return ZstdGenotype(prefix), codes

    def test_the_binary_variant_table_is_written_and_preferred(self):
        """Schema 2 reads variant metadata from arrays, not from the TSV.

        The TSV parse was **6.33 s of a 6.53 s store open** on the benchmark
        cohort -- 466 MB of text, 30% of a full scan's wall time, paid on every
        open. A reader-side cache hides it after the first open; writing arrays
        means a store never pays it at all.

        The check has to prove the arrays are actually *used*, not merely
        present, so it rewrites the TSV with wrong values: a reader still
        parsing text would hand those back.
        """
        import json

        from torchgwas.zstd_store import ZstdGenotype

        with tempfile.TemporaryDirectory() as tmpdir:
            store, _codes = self._store(tmpdir)
            prefix = Path(tmpdir) / "store"
            # Separate `.npy` files in a directory, not one `.npz`: the reader
            # memory-maps them, and `np.load` on a zip reads every array in
            # full.
            variants_dir = Path(tmpdir) / "store.variants"
            self.assertTrue(variants_dir.is_dir())
            self.assertEqual(
                sorted(path.stem for path in variants_dir.glob("*.npy")),
                ["chromosome", "effect_allele", "marker_id", "other_allele",
                 "position"])
            manifest = json.loads((Path(tmpdir) / "store.complete.json").read_text())
            self.assertEqual(manifest["cache_schema"], 2)
            expected = np.asarray(store.marker_ids, dtype=str)

            lines = (Path(tmpdir) / "store.variants.tsv").read_text().split("\n")
            poisoned = [lines[0]] + [
                line.replace("rs", "POISONED", 1) if line else line
                for line in lines[1:]]
            (Path(tmpdir) / "store.variants.tsv").write_text("\n".join(poisoned))

            got = np.asarray(ZstdGenotype(prefix).marker_ids, dtype=str)
            np.testing.assert_array_equal(got, expected)
            self.assertFalse(any(name.startswith("POISONED") for name in got),
                             "the reader parsed the TSV instead of the arrays")

    def test_a_schema_1_store_without_arrays_still_reads(self):
        """Stores written before the format change must keep working.

        Removing the `.variants` directory is what a schema 1 store looks like,
        and the reader has to fall back to parsing the TSV rather than fail to
        open.
        """
        import json
        import shutil

        from torchgwas.zstd_store import ZstdGenotype

        with tempfile.TemporaryDirectory() as tmpdir:
            store, _codes = self._store(tmpdir)
            prefix = Path(tmpdir) / "store"
            expected = np.asarray(store.marker_ids, dtype=str)
            positions = np.asarray(store.positions, dtype=np.int64)

            shutil.rmtree(Path(tmpdir) / "store.variants")
            manifest_path = Path(tmpdir) / "store.complete.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["cache_schema"] = 1
            manifest_path.write_text(json.dumps(manifest))

            legacy = ZstdGenotype(prefix)
            np.testing.assert_array_equal(
                np.asarray(legacy.marker_ids, dtype=str), expected)
            np.testing.assert_array_equal(
                np.asarray(legacy.positions, dtype=np.int64), positions)

    def test_the_store_declares_the_direct_fill(self):
        from torchgwas.zstd_store import ZstdGenotype

        self.assertTrue(
            getattr(ZstdGenotype, "allows_direct_native_fill", False),
            "PinnedDosageLoader will fall back to copying every decoded byte")

    def test_direct_fill_matches_iter_chunks_byte_for_byte(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            store, codes = self._store(tmpdir)
            owning = np.concatenate(
                [np.asarray(block) for _s, _e, block
                 in store.iter_chunks(chunk_size=200, dtype=np.uint8)], axis=1)
            filled = np.empty((codes.shape[1], codes.shape[0]), dtype=np.uint8)
            with store.native_reader_session() as read_into:
                for start in range(0, codes.shape[1], 200):
                    end = min(codes.shape[1], start + 200)
                    read_into(start, end, filled[start:end])
            self.assertTrue(np.array_equal(owning, filled.T))
            self.assertTrue(np.array_equal(codes, filled.T))

    def test_a_chunk_that_straddles_frames_is_still_exact(self):
        # Frames hold 128 variants; a 200-wide chunk never lands on a frame
        # boundary, which is the case the scratch branch exists for.
        with tempfile.TemporaryDirectory() as tmpdir:
            store, codes = self._store(tmpdir)
            for chunk in (128, 200, 300, 700):
                with self.subTest(chunk=chunk):
                    filled = np.empty((codes.shape[1], codes.shape[0]),
                                      dtype=np.uint8)
                    with store.native_reader_session() as read_into:
                        for start in range(0, codes.shape[1], chunk):
                            end = min(codes.shape[1], start + chunk)
                            read_into(start, end, filled[start:end])
                    self.assertTrue(np.array_equal(codes, filled.T))

    def test_reads_are_thread_safe(self):
        # The loader runs one fill per slot concurrently, so `read_into` is
        # called from several threads at once on disjoint ranges.
        from concurrent.futures import ThreadPoolExecutor

        with tempfile.TemporaryDirectory() as tmpdir:
            store, codes = self._store(tmpdir)
            filled = np.empty((codes.shape[1], codes.shape[0]), dtype=np.uint8)
            spans = [(s, min(codes.shape[1], s + 100))
                     for s in range(0, codes.shape[1], 100)]
            with store.native_reader_session() as read_into:
                with ThreadPoolExecutor(max_workers=8) as pool:
                    list(pool.map(
                        lambda span: read_into(span[0], span[1],
                                               filled[span[0]:span[1]]),
                        spans))
            self.assertTrue(np.array_equal(codes, filled.T))
