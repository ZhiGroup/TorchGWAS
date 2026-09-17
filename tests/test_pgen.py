from __future__ import annotations

import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from torchgwas.cli import _build_parser, _run_convert_pgen
from torchgwas.io import infer_genotype_format, load_genotype, load_pgen_genotype
from torchgwas.pgen import (
    DEFAULT_PGEN_COMPRESSION_WORKERS,
    DEFAULT_PGEN_DECODE_BATCH_SIZE,
    DEFAULT_PGEN_DECODE_WORKERS,
    PgenDosageSource,
    resolve_pgen_triplet,
)


HARDCALLS = np.asarray(
    [
        [0, 1, 2, 1],
        [2, -9, 0, 1],
        [1, 0, -9, 2],
        [4, 0, 1, 2],
        [2, 2, 1, 0],
    ],
    dtype=np.int8,
)


def _write_pgen_companions(root: Path, *, variant_count: int = 5) -> tuple[Path, Path, Path]:
    prefix = root / "study.v1"
    pgen, pvar, psam = (
        Path(f"{prefix}.pgen"),
        Path(f"{prefix}.pvar"),
        Path(f"{prefix}.psam"),
    )
    pgen.write_bytes(b"fake-pgen")
    psam.write_text("#FID IID SEX\nF0 I0 1\nF1 I1 2\nF2 I2 1\nF3 I3 2\n")
    alleles = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C"), ("A", "C")]
    pvar.write_text(
        "##fileformat=VCFv4.2\n#CHROM POS ID REF ALT\n"
        + "".join(
            f"{1 + index % 2} {100 + index} "
            f"{'.' if index == 1 else f'rs{index}'} "
            f"{alleles[index][0]} {alleles[index][1]}\n"
            for index in range(variant_count)
        )
    )
    return pgen, pvar, psam


class _FakePgenReader:
    matrix = HARDCALLS

    def __init__(self, filename, *, raw_sample_ct, variant_ct, sample_subset=None):
        del filename
        if raw_sample_ct != self.matrix.shape[1] or variant_ct != self.matrix.shape[0]:
            raise RuntimeError("dimension mismatch")
        self.subset = (
            np.arange(raw_sample_ct, dtype=np.uint32)
            if sample_subset is None
            else np.asarray(sample_subset, dtype=np.uint32)
        )
        self.closed = False

    def get_raw_sample_ct(self):
        return self.matrix.shape[1]

    def get_variant_ct(self):
        return self.matrix.shape[0]

    def read_range(self, start, end, output, *, allele_idx, sample_maj):
        if allele_idx != 1 or sample_maj:
            raise AssertionError("TorchGWAS must request variant-major ALT1 counts")
        output[:] = self.matrix[start:end, self.subset]

    def close(self):
        self.closed = True


class _FakeDosageReader(_FakePgenReader):
    dosages = np.asarray(
        [
            [0.0, 0.25, 1.5, 2.0],
            [2.0, 1.0, 0.5, 0.0],
            [1.0, -9.0, 1.0, 1.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.1, 0.2, 0.3, 0.4],
        ],
        dtype=np.float32,
    )

    def read_dosages_range(self, start, end, output, *, allele_idx, sample_maj):
        if allele_idx != 1 or sample_maj:
            raise AssertionError("TorchGWAS must request variant-major ALT1 dosage")
        output[:] = self.dosages[start:end, self.subset]


class PgenSourceTestCase(unittest.TestCase):
    def test_packed_transport_opt_in_contract_and_public_iteration(self):
        from torchgwas.pgen import PgenGenotype

        class PackedReader(_FakeDosageReader):
            matrix = np.where(HARDCALLS == 4, 0, HARDCALLS)

            def read_packed_range_into(self, start, end, output):
                output.fill(0)
                codes = np.where(self.matrix[start:end] == -9, 3, self.matrix[start:end]).astype(np.uint8)
                for sample in range(codes.shape[1]):
                    output[:, sample // 4] |= codes[:, sample] << (2 * (sample % 4))

        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, _, _ = _write_pgen_companions(Path(tmpdir))
            with mock.patch.dict(os.environ, {"TORCHGWAS_PGEN_PACKED": "1"}), mock.patch("torchgwas.pgen._open_pgen", side_effect=PackedReader):
                source = PgenGenotype(pgen, mode="hardcall")
                multiworker_source = PgenGenotype(pgen, mode="hardcall", reader_workers=2)
                self.assertEqual(multiworker_source.native_encoding, "pgen_2bit")
                self.assertEqual(source.native_encoding, "pgen_2bit")
                self.assertEqual(source.native_dtype, np.dtype(np.int8))
                self.assertEqual(source.native_transfer_dtype, np.dtype(np.uint8))
                self.assertEqual(source.native_row_width, 64)
                backing = np.empty(5 * 64 + 64, dtype=np.uint8)
                offset = (-backing.ctypes.data) % 64
                output = backing[offset:offset + 5 * 64].reshape(5, 64)
                with source.native_reader_session() as fill:
                    self.assertIs(fill(0, 5, output), output)
                    with self.assertRaisesRegex(ValueError, "aligned"):
                        fill(0, 1, backing[offset + 1:offset + 65].reshape(1, 64))
                    with self.assertRaisesRegex(ValueError, "exact shape"):
                        fill(0, 5, np.empty((5, 4), dtype=np.int8))
                with multiworker_source.native_reader_session() as fill:
                    fill(0, 5, output)
                decoded = np.stack([(output[:, 0] >> (2 * i)) & 3 for i in range(4)], axis=1).astype(np.int8)
                decoded[decoded == 3] = -9
                np.testing.assert_array_equal(decoded, PackedReader.matrix)
                self.assertFalse(output[:, 1:].any())
                chunks = list(source.iter_native_chunks(2))
                np.testing.assert_array_equal(np.concatenate([x[2] for x in chunks], axis=1).T, PackedReader.matrix)
                for kwargs in ({"mode": "auto"}, {"mode": "dosage"}, {"mode": "hardcall", "selected_sample_ids": ["I0", "I1"]}, {"mode": "hardcall", "selected_sample_ids": ["I1", "I0", "I2", "I3"]}):
                    with self.assertRaisesRegex(ValueError, "requires"):
                        PgenGenotype(pgen, **kwargs)
            with mock.patch.dict(os.environ, {"TORCHGWAS_PGEN_PACKED": "1"}), mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakePgenReader):
                # A reader with no `read_packed_range_into` still refuses, but
                # the reason is no longer "build the optional pgenlib bridge".
                # `NativePgenReader` provides that method itself, so the bridge
                # stopped being the only source of it; the refusal now names
                # the missing capability rather than a specific dependency, and
                # is a ValueError like the other two rejections above.
                with self.assertRaisesRegex(ValueError, "read_packed_range_into"):
                    PgenGenotype(pgen, mode="hardcall")
            with mock.patch.dict(os.environ, {"TORCHGWAS_PGEN_PACKED": "0"}), mock.patch("torchgwas.pgen._open_pgen", side_effect=PackedReader):
                source = PgenGenotype(pgen, mode="hardcall")
                self.assertEqual(source.native_encoding, "dosage")
                self.assertEqual(source.native_row_width, 4)
                self.assertEqual(source.native_transfer_dtype, np.dtype(np.int8))

    def test_pvar_cache_round_trips_and_is_rejected_when_stale(self):
        """A cached PVAR must equal a parsed one, and must not outlive its file.

        Worth caching at all because the parse is paid on every open whatever
        the scan then reads: 5.7-7.3 s on the 8.93M-variant cohort, which on
        the packed transport is 60% of the entire run.

        The cache lives in its own temporary directory here. A staleness test
        earlier in this project mutated a shared cache to prove rejection and
        never restored it, and every later run silently re-parsed -- which then
        read as "the cache gives no benefit".
        """
        from torchgwas.pgen import (_load_pgen_pvar_cache, _read_pvar,
                                    _read_pvar_cached)

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _pgen, pvar, _psam = _write_pgen_companions(root)
            cache_dir = root / "cache"
            direct = _read_pvar(pvar)

            first = _read_pvar_cached(pvar, cache_dir)
            self.assertIsNotNone(_load_pgen_pvar_cache(cache_dir, pvar)[0])
            second = _read_pvar_cached(pvar, cache_dir)
            for key, expected in direct.items():
                # Text columns come back as fixed-width unicode from the cache
                # and as object arrays from a fresh parse, exactly as the BIM
                # cache already behaves; compare values, not dtypes.
                np.testing.assert_array_equal(
                    np.asarray(first[key], dtype=str if key != "position" else np.int64),
                    np.asarray(expected, dtype=str if key != "position" else np.int64))
                np.testing.assert_array_equal(
                    np.asarray(second[key], dtype=str if key != "position" else np.int64),
                    np.asarray(expected, dtype=str if key != "position" else np.int64))

            # Touching the PVAR must invalidate: the cache keys on size and
            # mtime, so a rewritten file gets a different cache path and the
            # stale entry is simply never consulted again.
            cached, cache_path = _load_pgen_pvar_cache(cache_dir, pvar)
            self.assertIsNotNone(cached)
            text = pvar.read_text()
            self.assertIn("rs0", text)
            pvar.write_text(text.replace("rs0", "rs0x", 1))
            stale, new_path = _load_pgen_pvar_cache(cache_dir, pvar)
            self.assertIsNone(stale)
            self.assertNotEqual(new_path, cache_path)
            refreshed = _read_pvar_cached(pvar, cache_dir)
            self.assertIn("rs0x", [str(value) for value in refreshed["marker_id"]])

    def test_packed_transport_follows_the_native_statistics_backend(self):
        """With no explicit flag, packed engages exactly when it is readable.

        The packed row has one consumer -- the native fused statistics kernels
        -- and `native_scan.dosage_cuda_iterator` refuses `pgen_2bit` outright
        without them. So "on whenever eligible" is not safe: it would turn a
        working torch-backend scan into an error at scan time. It follows
        TORCHGWAS_NATIVE_STATS instead.

        This matters because the flag being off by default is what made every
        PGEN measurement in this project move 22,250 bytes per variant across
        PCIe where 5,568 would do -- four times the payload of the equivalent
        BED scan, on the same calls.
        """
        from torchgwas.pgen import PgenGenotype

        class PackedReader(_FakeDosageReader):
            matrix = np.where(HARDCALLS == 4, 0, HARDCALLS)

            def read_packed_range_into(self, start, end, output):
                output.fill(0)

        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, _, _ = _write_pgen_companions(Path(tmpdir))
            for native_stats, encoding, width, dtype in (
                    ("1", "pgen_2bit", 64, np.uint8),
                    ("0", "dosage", 4, np.int8)):
                environment = {"TORCHGWAS_NATIVE_STATS": native_stats}
                with self.subTest(native_stats=native_stats), \
                        mock.patch.dict(os.environ, environment), \
                        mock.patch("torchgwas.pgen._open_pgen",
                                   side_effect=PackedReader):
                    os.environ.pop("TORCHGWAS_PGEN_PACKED", None)
                    source = PgenGenotype(pgen, mode="hardcall")
                    self.assertEqual(source.native_encoding, encoding)
                    self.assertEqual(source.native_row_width, width)
                    self.assertEqual(source.native_transfer_dtype,
                                     np.dtype(dtype))

            # An explicit 0 still wins over the backend being on, which is what
            # an A/B of the two transports needs.
            with mock.patch.dict(os.environ, {"TORCHGWAS_NATIVE_STATS": "1",
                                              "TORCHGWAS_PGEN_PACKED": "0"}), \
                    mock.patch("torchgwas.pgen._open_pgen",
                               side_effect=PackedReader):
                source = PgenGenotype(pgen, mode="hardcall")
                self.assertEqual(source.native_encoding, "dosage")

    def test_metadata_alt_orientation_missing_policy_and_requested_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, pvar, psam = _write_pgen_companions(Path(tmpdir))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakePgenReader):
                genotype, sample_ids, marker_ids = load_pgen_genotype(
                    pgen,
                    pvar=pvar,
                    psam=psam,
                    selected_sample_ids=["I2", "I0", "I3"],
                    cache_dir=Path(tmpdir) / "cache",
                    pgen_mode="hardcall",
                    decode_workers=2,
                    decode_batch_size=2,
                    compression_workers=2,
                    zstd_chunk_size=2,
                )

            # The conversion writes into a staging directory and `io.py` moves
            # a fixed list of suffixes into the cache. `.variants.npz` -- the
            # schema 2 binary variant table -- was not on that list, so the
            # first real store was written correctly and then thrown away with
            # the temp directory, silently, after a 37-minute conversion. It
            # went unnoticed because every other test calls the writer directly
            # against a prefix and never reaches the publish step.
            published = [path for path in (Path(tmpdir) / "cache").glob("*.variants")
                         if path.is_dir()]
            self.assertTrue(published,
                            "the binary variant table was left in staging")
            self.assertEqual(
                sorted(path.stem for path in published[0].glob("*.npy")),
                ["chromosome", "effect_allele", "marker_id", "other_allele",
                 "position"])

            self.assertEqual(sample_ids.tolist(), ["I2", "I0", "I3"])
            # rs2 carries a missing call for the selected samples and is now
            # kept and masked, so it joins the marker list. rs3 still leaves:
            # its call of 4 is out of range, which masking does not excuse.
            self.assertEqual(marker_ids.tolist(), ["rs0", "2:101:C:T", "rs2", "rs4"])
            self.assertEqual(genotype.effect_alleles.tolist(), ["G", "T", "A", "C"])
            self.assertEqual(genotype.other_alleles.tolist(), ["A", "C", "G", "A"])
            # rs2 is [1, 0, -9, 2]; samples [I2, I0, I3] select [-9, 1, 2], so
            # the observed mean is 1.5 and the masked call becomes
            # rint(1.5 * 127.5) = 191.
            np.testing.assert_array_equal(
                genotype.read_codes(0, 4),
                np.asarray(
                    [[255, 0, 191, 128], [0, 255, 128, 255], [128, 128, 255, 0]],
                    dtype=np.uint8,
                ),
            )
            manifest = json.loads(genotype.manifest_path.read_text())
            self.assertEqual(
                manifest["missing_policy"],
                "mask_missing_calls_keep_variant",
            )
            self.assertEqual(manifest["exclusion_counts"],
                             {"invalid_call": 1, "masked_missing": 1})
            self.assertEqual(manifest["effect_allele"], "PVAR_ALT1_HARDCALL")
            genotype.close()

    def test_cache_reuse_invalidation_and_pgen_namespace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pgen, pvar, psam = _write_pgen_companions(root)
            cache = root / "cache"
            opened: list[_FakePgenReader] = []

            def open_fake(*args, **kwargs):
                reader = _FakePgenReader(*args, **kwargs)
                opened.append(reader)
                return reader

            with mock.patch("torchgwas.pgen._open_pgen", side_effect=open_fake):
                first, _, _ = load_pgen_genotype(
                    pgen,
                    cache_dir=cache,
                    pgen_mode="hardcall",
                    decode_workers=1,
                    compression_workers=1,
                )
                second, _, _ = load_pgen_genotype(
                    pgen,
                    cache_dir=cache,
                    pgen_mode="hardcall",
                    decode_workers=1,
                    compression_workers=1,
                )
                self.assertEqual(first.prefix, second.prefix)
                self.assertEqual(len(opened), 1)
                pvar.write_text(pvar.read_text() + "\n")
                third, _, _ = load_pgen_genotype(
                    pgen,
                    cache_dir=cache,
                    pgen_mode="hardcall",
                    decode_workers=1,
                    compression_workers=1,
                )

            self.assertEqual(len(opened), 2)
            self.assertNotEqual(first.prefix, third.prefix)
            self.assertIn("_pgen_", first.prefix.name)
            for path in cache.glob("*.complete.json"):
                self.assertEqual(json.loads(path.read_text())["source_format"], "pgen")
            first.close()
            second.close()
            third.close()

    def test_dosage_range_quantization_and_clear_unsupported_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, pvar, psam = _write_pgen_companions(Path(tmpdir))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakePgenReader):
                with self.assertRaisesRegex(RuntimeError, "read_dosages_range"):
                    PgenDosageSource(pgen, pvar=pvar, psam=psam, mode="dosage")
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakeDosageReader):
                genotype, _, _ = load_pgen_genotype(
                    pgen,
                    pvar=pvar,
                    psam=psam,
                    cache_dir=Path(tmpdir) / "dosage-cache",
                    pgen_mode="dosage",
                    decode_workers=1,
                    compression_workers=1,
                )
            # rs2 is [1.0, -9.0, 1.0, 1.0]: the observed calls are all 1.0, so
            # the masked one is too and the whole column is rint(1.0 * 127.5).
            np.testing.assert_array_equal(
                genotype.read_codes(0, 4),
                np.asarray(
                    [
                        [0, 255, 128, 13],
                        [32, 128, 128, 26],
                        [191, 64, 128, 38],
                        [255, 0, 128, 51],
                    ],
                    dtype=np.uint8,
                ),
            )
            self.assertEqual(genotype.marker_ids.tolist(),
                             ["rs0", "2:101:C:T", "rs2", "rs4"])
            genotype.close()

    def test_pvar_typed_positions_keep_lexical_metadata_and_ignore_annotations(self):
        from torchgwas.pgen import _read_pvar
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "lexical.pvar"
            for sep in ("\t", " "):
                path.write_text("##source=test\n" + sep.join(["#CHROM", "POS", "ID", "REF", "ALT", "INFO"]) + "\n"
                                + sep.join(["01", "42", "0007", "A", "G", "UNUSED=123"]) + "\n"
                                + sep.join(["NA", "43", "NA", "NA", "C", "UNUSED=456"]) + "\n"
                                + sep.join(["02", "44", ".", "G", "T", "."]) + "\n")
                result = _read_pvar(path)
                self.assertEqual(result["chromosome"].tolist(), ["01", "NA", "02"])
                self.assertEqual(result["marker_id"].tolist(), ["0007", "NA:43:NA:C", "02:44:G:T"])
                self.assertEqual(result["other_allele"].tolist(), ["A", "NA", "G"])
                self.assertEqual(result["position"].dtype, np.dtype(np.int64))
                np.testing.assert_array_equal(result["position"], [42, 43, 44])
                self.assertNotIn("INFO", result)

    def test_pvar_missing_columns_and_invalid_integer_positions(self):
        from torchgwas.pgen import _read_pvar
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "invalid.pvar"
            path.write_text("#CHROM POS ID REF\n1 42 rs1 A\n")
            with self.assertRaisesRegex(ValueError, "missing required columns.*ALT"):
                _read_pvar(path)
            for position in ("bad", "1.5", "NA"):
                path.write_text(f"#CHROM\tPOS\tID\tREF\tALT\n1\t{position}\trs1\tA\tG\n")
                with self.assertRaisesRegex(ValueError, "invalid integer POS"):
                    _read_pvar(path)

    def test_multiallelic_pvar_fails_before_opening_pgen(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, pvar, psam = _write_pgen_companions(Path(tmpdir))
            pvar.write_text(pvar.read_text().replace("rs0 A G", "rs0 A G,T"))
            with mock.patch("torchgwas.pgen._open_pgen") as opener:
                with self.assertRaisesRegex(ValueError, "biallelic"):
                    PgenDosageSource(pgen, pvar=pvar, psam=psam, mode="hardcall")
            opener.assert_not_called()

    def test_pgen_decode_workers_follows_reader_workers(self):
        """The reader count must actually reach PgenGenotype.

        Asserting only that the CLI default is None would stay green if
        someone restored a literal default, because the argument reaches the
        iterator either way; what matters is the count the reader is built
        with, since native_scan prefers source.decode_workers over its own
        argument and a reader session cannot be widened after construction.
        """
        import pgenlib

        from torchgwas.pgen import PgenGenotype

        with tempfile.TemporaryDirectory() as tmpdir:
            pgen, pvar, psam = _write_pgen_companions(Path(tmpdir))
            # The companions helper does not write genotype data, so put a
            # real five-variant PGEN behind it before opening it for real.
            with pgenlib.PgenWriter(os.fsencode(pgen), 4, variant_ct=5,
                                    nonref_flags=False) as writer:
                writer.append_biallelic_batch(HARDCALLS.copy())
            source, _, _, meta = load_genotype(
                pgen, genotype_format="pgen", pvar=pvar, psam=psam,
                reader_workers=6,
            )
            try:
                self.assertEqual(source.reader_workers, 6)
                self.assertEqual(source.decode_workers, 6)
                self.assertEqual(meta["reader_workers"], 6)
                self.assertIsNone(meta["pgen_decode_workers_requested"])
            finally:
                source.close()

            explicit, _, _, meta = load_genotype(
                pgen, genotype_format="pgen", pvar=pvar, psam=psam,
                reader_workers=6, pgen_decode_workers=2,
            )
            try:
                self.assertEqual(explicit.reader_workers, 2)
                self.assertEqual(meta["reader_workers"], 2)
                self.assertEqual(meta["pgen_decode_workers_requested"], 2)
            finally:
                explicit.close()

    def test_triplet_format_and_api_cli_dispatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "cohort.release.1"
            pgen, pvar, psam = _write_pgen_companions(Path(tmpdir))
            renamed = []
            for old, suffix in ((pgen, ".pgen"), (pvar, ".pvar"), (psam, ".psam")):
                new = Path(f"{prefix}{suffix}")
                old.replace(new)
                renamed.append(new)
            self.assertEqual(resolve_pgen_triplet(prefix), tuple(renamed))
            self.assertEqual(infer_genotype_format(prefix), "pgen")
            from torchgwas.pgen import PgenGenotype
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakeDosageReader), mock.patch(
                "torchgwas.io.load_pgen_genotype"
            ) as converter:
                observed, _, _, meta = load_genotype(
                    renamed[0], genotype_format="pgen", pvar=renamed[1], psam=renamed[2],
                    genotype_cache_dir=Path(tmpdir) / "unused-cache",
                )
            self.assertIsInstance(observed, PgenGenotype)
            self.assertEqual(meta["genotype_format"], "pgen")
            self.assertEqual(meta["genotype_backend"], "direct_pgen")
            self.assertFalse(meta["cache_conversion"])
            # The directory is no longer unused. It never meant "cache the
            # conversion" for the direct PGEN path -- that stays off, as
            # `cache_conversion` still asserts -- but it now holds the parsed
            # PVAR, which is what `genotype_cache_dir` already means for PLINK
            # (its BIM) and BGEN (its index). PGEN was the one format paying a
            # full metadata parse on every open: 5.7-7.3 s on the benchmark
            # cohort, 60% of a packed-transport run.
            self.assertTrue((Path(tmpdir) / "unused-cache").is_dir())
            converter.assert_not_called()
            args = _build_parser().parse_args(
                ["linear", "--genotype", str(renamed[0]), "--genotype-format", "pgen", "--output-dir", tmpdir]
            )
            self.assertEqual(args.genotype_format, "pgen")
            # The CLI must not pin a decode-worker count. PgenGenotype fixes
            # its reader count at construction and native_scan prefers
            # source.decode_workers over the argument it is handed, so a
            # literal default here silently caps PGEN at that many readers
            # however large --reader-workers is. None means "follow
            # --reader-workers"; an explicit value still wins.
            self.assertIsNone(args.pgen_decode_workers)
            self.assertEqual(args.pgen_decode_batch_size, DEFAULT_PGEN_DECODE_BATCH_SIZE)
            self.assertEqual(
                args.pgen_compression_workers, DEFAULT_PGEN_COMPRESSION_WORKERS
            )

            manifest_path = Path(tmpdir) / "converted.complete.json"
            manifest_path.write_text('{"source_format": "pgen", "n_variants": 5}\n')
            fake_genotype = SimpleNamespace(
                prefix=Path(tmpdir) / "converted", manifest_path=manifest_path
            )
            convert_args = _build_parser().parse_args(
                [
                    "convert-pgen",
                    "--genotype",
                    str(renamed[0]),
                    "--cache-dir",
                    str(Path(tmpdir) / "cache"),
                    "--output-json",
                    str(Path(tmpdir) / "conversion.json"),
                ]
            )
            with mock.patch(
                "torchgwas.cli.load_pgen_genotype",
                return_value=(fake_genotype, np.asarray(["I0"]), np.asarray(["rs0"])),
            ) as cli_loader:
                self.assertEqual(_run_convert_pgen(convert_args), 0)
            cli_loader.assert_called_once()
            self.assertEqual(
                cli_loader.call_args.kwargs["compression_workers"],
                DEFAULT_PGEN_COMPRESSION_WORKERS,
            )
            self.assertEqual(
                json.loads(Path(convert_args.output_json).read_text())["source_format"], "pgen"
            )

@unittest.skipUnless(importlib.util.find_spec("pgenlib"), "optional pgenlib is not installed")
class PgenIntegrationTestCase(unittest.TestCase):
    def test_real_pgenlib_hardcall_fixture_matches_trusted_matrix(self):
        import pgenlib

        hardcalls = np.asarray(
            [[0, 1, 2, 0], [2, 1, 0, 2], [1, -9, 1, 0], [2, 2, 1, 0]],
            dtype=np.int8,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            prefix = root / "trusted"
            pgen = Path(f"{prefix}.pgen")
            with pgenlib.PgenWriter(
                os.fsencode(pgen),
                hardcalls.shape[1],
                variant_ct=hardcalls.shape[0],
                nonref_flags=False,
            ) as writer:
                writer.append_biallelic_batch(hardcalls)
            Path(f"{prefix}.psam").write_text(
                "#IID\n" + "".join(f"S{index}\n" for index in range(hardcalls.shape[1]))
            )
            Path(f"{prefix}.pvar").write_text(
                "#CHROM POS ID REF ALT\n"
                + "".join(
                    f"1 {100 + index} rs{index} A G\n" for index in range(hardcalls.shape[0])
                )
            )
            genotype, _, _ = load_pgen_genotype(
                pgen,
                cache_dir=root / "cache",
                pgen_mode="hardcall",
                decode_workers=1,
                compression_workers=1,
                zstd_chunk_size=2,
            )
            # Variant 2 is [1, -9, 1, 0]. It used to be dropped; it is now
            # masked, so all four variants survive. Its observed calls are
            # 1, 1, 0, a mean of 2/3, and the masked call becomes
            # rint(2/3 * 127.5) = 85. The other three are unchanged.
            table = np.asarray([0, 128, 255], dtype=np.uint8)
            masked = hardcalls.copy()
            observed = masked[2] != -9
            expected_codes = np.empty((hardcalls.shape[1], hardcalls.shape[0]),
                                      dtype=np.uint8)
            expected_codes[:, [0, 1, 3]] = table[masked[[0, 1, 3]].T]
            column = table[np.where(observed, masked[2], 0)].astype(np.uint8)
            column[~observed] = np.rint(masked[2][observed].mean() * 127.5)
            expected_codes[:, 2] = column
            self.assertEqual(int(column[~observed][0]), 85)
            np.testing.assert_array_equal(genotype.read_codes(0, 4), expected_codes)
            np.testing.assert_allclose(
                genotype.read_chunk(0, 4), expected_codes.astype(np.float32) / 127.5
            )
            genotype.close()

    def test_real_pgenlib_dosage_range_fixture_preserves_quantized_dosage(self):
        import pgenlib

        dosages = np.asarray(
            [[0.0, 0.25, 1.5, 2.0], [2.0, 1.0, 0.5, 0.0], [1.0, -9.0, 1.0, 1.0]],
            dtype=np.float32,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            prefix = root / "trusted-dosage"
            pgen = Path(f"{prefix}.pgen")
            with pgenlib.PgenWriter(
                os.fsencode(pgen),
                dosages.shape[1],
                variant_ct=dosages.shape[0],
                nonref_flags=False,
                dosage_present=True,
            ) as writer:
                writer.append_dosages_batch(dosages)
            Path(f"{prefix}.psam").write_text("#IID\nS0\nS1\nS2\nS3\n")
            Path(f"{prefix}.pvar").write_text(
                "#CHROM POS ID REF ALT\n1 100 rs0 A G\n1 101 rs1 C T\n1 102 rs2 G A\n"
            )
            genotype, _, _ = load_pgen_genotype(
                pgen,
                cache_dir=root / "cache",
                decode_workers=1,
                compression_workers=1,
            )
            expected = np.rint(dosages[:2].T * 127.5).astype(np.uint8)
            np.testing.assert_array_equal(genotype.read_codes(0, 2), expected)
            manifest = json.loads(genotype.manifest_path.read_text())
            self.assertEqual(manifest["effect_allele"], "PVAR_ALT1_DOSAGE")
            self.assertEqual(manifest["pgen_mode"], "auto")
            self.assertEqual(manifest["resolved_pgen_mode"], "dosage")
            self.assertEqual(manifest["exclusion_counts"]["masked_missing"], 1)
            genotype.close()




class DirectPgenTestCase(unittest.TestCase):
    def test_direct_dosages_precision_order_missing_and_repeat(self):
        from torchgwas.pgen import PgenGenotype
        with tempfile.TemporaryDirectory() as tmp:
            pgen, _, _ = _write_pgen_companions(Path(tmp))
            # Exclude sample I2 with deliberately invalid mock dosage 3.0.
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakeDosageReader):
                source = PgenGenotype(pgen, selected_sample_ids=["I3", "I1", "I0"])
                self.assertEqual(source.shape, (3, 5))
                self.assertEqual(source.variant_metadata["effect_allele"][0], "G")
                expected = _FakeDosageReader.dosages[:, [3, 1, 0]].T.copy()
                expected[expected == -9] = np.nan
                for workers, chunk in [(1, 2), (2, 3), (1, 5)]:
                    blocks = list(source.iter_chunks(chunk, reader_workers=workers))
                    actual = np.concatenate([block for _, _, block in blocks], axis=1)
                    np.testing.assert_array_equal(actual, expected)
                    self.assertEqual(actual.dtype, np.float32)
                self.assertEqual(source.marker_ids.size, 5)

    def test_direct_hardcalls_are_unscaled_and_invalid_fails(self):
        from torchgwas.pgen import PgenGenotype
        with tempfile.TemporaryDirectory() as tmp:
            pgen, _, _ = _write_pgen_companions(Path(tmp))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=_FakePgenReader):
                source = PgenGenotype(pgen, mode="hardcall")
                native = source.read_chunk(0, 3, np.int8)
                np.testing.assert_array_equal(native, HARDCALLS[:3].T)
                actual = source.read_chunk(0, 3)
                expected = HARDCALLS[:3].T.astype(np.float32)
                expected[expected == -9] = np.nan
                np.testing.assert_array_equal(actual, expected)
                with self.assertRaisesRegex(ValueError, "invalid PGEN"):
                    source.read_chunk(3, 4)
                with self.assertRaisesRegex(ValueError, "float32 or float64"):
                    source.read_chunk(0, 1, np.uint8)

    def test_direct_early_close_releases_thread_readers(self):
        from torchgwas.pgen import PgenGenotype
        opened = []
        def create(*args, **kwargs):
            reader = _FakePgenReader(*args, **kwargs)
            opened.append(reader)
            return reader
        with tempfile.TemporaryDirectory() as tmp:
            pgen, _, _ = _write_pgen_companions(Path(tmp))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=create):
                source = PgenGenotype(pgen, mode="hardcall", selected_sample_ids=["I3"])
                iterator = source.iter_chunks(1, reader_workers=2)
                next(iterator)
                iterator.close()
                self.assertTrue(opened)
                self.assertTrue(all(reader.closed for reader in opened))


    @unittest.skipUnless(importlib.util.find_spec("pgenlib"), "optional pgenlib is not installed")
    def test_real_direct_round_trip_both_encodings(self):
        import pgenlib
        from torchgwas.pgen import PgenGenotype
        with tempfile.TemporaryDirectory() as tmp:
            pgen, pvar, psam = _write_pgen_companions(Path(tmp))
            for mode in ("hardcall", "dosage"):
                values = HARDCALLS.copy() if mode == "hardcall" else _FakeDosageReader.dosages.copy()
                values[3] = [0, 1, 2, 0]
                with pgenlib.PgenWriter(os.fsencode(pgen), 4, variant_ct=5,
                                       nonref_flags=False, dosage_present=mode == "dosage") as writer:
                    if mode == "hardcall":
                        writer.append_biallelic_batch(values)
                    else:
                        writer.append_dosages_batch(values)
                if mode == "hardcall" and hasattr(pgenlib.PgenReader, "read_packed_range_into"):
                    with mock.patch.dict(os.environ, {"TORCHGWAS_PGEN_PACKED": "1"}):
                        packed_source = PgenGenotype(pgen, mode="hardcall", reader_workers=2)
                    backing = np.empty(5 * 64 + 64, dtype=np.uint8)
                    offset = (-backing.ctypes.data) % 64
                    packed = backing[offset:offset + 5 * 64].reshape(5, 64)
                    with packed_source.native_reader_session() as fill:
                        fill(0, 5, packed)
                    calls = np.stack([(packed[:, 0] >> (2 * sample)) & 3 for sample in range(4)], axis=1).astype(np.int8)
                    calls[calls == 3] = -9
                    np.testing.assert_array_equal(calls, values)
                    self.assertFalse(packed[:, 1:].any())
                source = PgenGenotype(pgen, mode=mode, selected_sample_ids=["I3", "I0", "I1"])
                observed = np.concatenate([x for _, _, x in source.iter_chunks(2)], axis=1)
                expected = values[:, [3, 0, 1]].T.astype(np.float32)
                expected[expected == -9] = np.nan
                np.testing.assert_allclose(observed, expected, rtol=0, atol=1 / 32768)
                native_out = np.empty((5, 3), dtype=source.native_dtype)
                with source.native_reader_session() as fill:
                    fill(0, 5, native_out)
                np.testing.assert_allclose(native_out, values[:, [3, 0, 1]], rtol=0, atol=1 / 32768)
                if mode == "dosage":
                    self.assertGreater(abs(float(observed[1, 4]) - round(float(observed[1, 4]) * 127.5) / 127.5), 0.0001)


    def test_native_session_fills_owned_buffers_and_preserves_sample_order(self):
        from torchgwas.pgen import PgenGenotype
        for mode, reader_type, matrix in (("hardcall", _FakePgenReader, HARDCALLS),
                                          ("dosage", _FakeDosageReader, _FakeDosageReader.dosages)):
            for selected in (["I0", "I1", "I3"], ["I3", "I0", "I1"]):
                opened = []
                def create(*args, **kwargs):
                    reader = reader_type(*args, **kwargs)
                    opened.append(reader)
                    return reader
                with tempfile.TemporaryDirectory() as tmp:
                    pgen, _, _ = _write_pgen_companions(Path(tmp))
                    with mock.patch("torchgwas.pgen._open_pgen", side_effect=create):
                        source = PgenGenotype(pgen, mode=mode, selected_sample_ids=selected)
                        out = np.empty((3, 3), dtype=source.native_dtype)
                        with mock.patch.dict(os.environ, {"TORCHGWAS_SCAN_PROFILE": "1"}), source.native_reader_session() as fill:
                            self.assertIs(fill(0, 3, out), out)
                            np.testing.assert_array_equal(out, matrix[:3, [int(x[1:]) for x in selected]])
                            fill(1, 3, out[:2])
                            np.testing.assert_array_equal(out[:2], matrix[1:3, [int(x[1:]) for x in selected]])
                            self.assertEqual(len(opened), 2)  # probe + persistent fill reader
                            with self.assertRaisesRegex(ValueError, "variant-major"):
                                fill(0, 3, out.T)
                            with self.assertRaises(IndexError):
                                fill(-1, 2, out)
                        self.assertTrue(all(x.closed for x in opened))
                        with self.assertRaisesRegex(RuntimeError, "closed"):
                            fill(0, 3, out)
                        self.assertEqual(source.native_host_staging_copies, int(selected[0] == "I3"))
                        self.assertEqual(source.native_missing_value, -9)
                        self.assertEqual(source.last_native_reader_profile["fill_calls"], 2)
                        self.assertEqual(source.last_native_reader_profile["reader_count"], 1)
                        self.assertGreaterEqual(source.last_native_reader_profile["fill_worker_cpu_seconds"], 0)

    def test_native_session_parallel_readers_and_exception_cleanup(self):
        from concurrent.futures import ThreadPoolExecutor
        import threading
        from torchgwas.pgen import PgenGenotype
        opened = []
        barrier = threading.Barrier(2)
        class ParallelReader(_FakePgenReader):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                opened.append(self)
            def read_range(self, *args, **kwargs):
                barrier.wait(timeout=5)
                return super().read_range(*args, **kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            pgen, _, _ = _write_pgen_companions(Path(tmp))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=ParallelReader):
                source = PgenGenotype(pgen, mode="hardcall", selected_sample_ids=["I3"])
                with self.assertRaisesRegex(RuntimeError, "caller failed"):
                    with source.native_reader_session() as fill:
                        with ThreadPoolExecutor(max_workers=2) as pool:
                            outputs = [np.empty((2, 1), dtype=np.int8) for _ in range(2)]
                            jobs = [pool.submit(fill, i, i+2, out) for i, out in enumerate(outputs)]
                            for job in jobs:
                                job.result()
                            for i, out in enumerate(outputs):
                                np.testing.assert_array_equal(out, HARDCALLS[i:i+2, [3]])
                        raise RuntimeError("caller failed")
                self.assertEqual(len(opened), 3)
                self.assertTrue(all(x.closed for x in opened))

    def test_decode_tiles_are_independent_and_each_variant_read_once(self):
        from torchgwas.pgen import PgenGenotype
        calls = []
        class RecordingReader(_FakePgenReader):
            def read_range(self, start, end, output, **kwargs):
                calls.append((start, end))
                return super().read_range(start, end, output, **kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            pgen, _, _ = _write_pgen_companions(Path(tmp))
            with mock.patch("torchgwas.pgen._open_pgen", side_effect=RecordingReader):
                for decode_tile, compute_chunk in [(2, 3), (3, 2), (4, 3), (2, 1)]:
                    calls.clear()
                    source = PgenGenotype(pgen, mode="hardcall", selected_sample_ids=["I3"],
                                          decode_batch_size=decode_tile, reader_workers=2)
                    blocks = list(source.iter_native_chunks(compute_chunk))
                    np.testing.assert_array_equal(np.concatenate([b for _, _, b in blocks], axis=1),
                                                  HARDCALLS[:, [3]].T)
                    self.assertEqual(sorted(calls), [(i, min(5, i + decode_tile))
                                                    for i in range(0, 5, decode_tile)])
                    self.assertEqual([(start, end) for start, end, _ in blocks],
                                     [(i, min(5, i + compute_chunk)) for i in range(0, 5, compute_chunk)])

if __name__ == "__main__":
    unittest.main()
