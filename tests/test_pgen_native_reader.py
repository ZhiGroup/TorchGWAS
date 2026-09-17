"""The native PGEN reader must match the Python reference, and be chosen well.

`PgenDirectReader` is the oracle-validated definition of correct for this
format, so every decode here is compared against it rather than against a
recorded expectation. pgenlib is deliberately not required: it is installed in
none of the lab environments, which is the whole reason this backend exists.

The fixture is a real PGEN written by hand -- storage mode 0x10, uncompressed
records -- because building one with plink2 would make these tests depend on an
external binary that CI does not have, and the alternative (a mocked reader)
would test the mock.
"""

from __future__ import annotations

import struct
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import numpy as np

from torchgwas import pgen_native
from torchgwas.pgen import resolve_pgen_mode_and_backend
from torchgwas.pgen_reader import (
    PgenDirectReader,
    _bytes_per_sample_id,
    file_scope,
    pack_genovec,
    unpack_genovec,
)
from torchgwas.pgen_native_reader import NativePgenReader, native_backend_reason


VBLOCK_SIZE = 1 << 16


def write_pgen(path: Path, categories: np.ndarray) -> np.ndarray:
    """Write a storage-mode-0x10 PGEN of uncompressed records.

    `categories` is (variants, samples) with PGEN's own codes: 0 hom-ref,
    1 het, 2 hom-alt, 3 missing. Returns it unchanged, so a caller can keep
    one object as both the input and the expectation.
    """
    categories = np.ascontiguousarray(categories, dtype=np.uint8)
    variant_ct, sample_ct = categories.shape
    if variant_ct > VBLOCK_SIZE:
        raise ValueError("the fixture writer emits a single variant block")
    genovecs = [pack_genovec(row, sample_ct) for row in categories]
    record_bytes = genovecs[0].size
    # header_ctrl: two-byte record lengths (0x01), byte-wide record types
    # (0x04), so the index layout does not depend on guessing a nibble width.
    header_ctrl = 0x01 | 0x04
    vblock_ct = 1
    index_bytes = variant_ct * 1 + variant_ct * 2
    data_start = 12 + 8 * vblock_ct + index_bytes

    blob = bytearray()
    blob += b"\x6c\x1b"
    blob += bytes([0x10])
    blob += struct.pack("<II", variant_ct, sample_ct)
    blob += bytes([header_ctrl])
    blob += struct.pack("<Q", data_start)
    blob += bytes([0] * variant_ct)                       # vrtype 0: plain genovec
    blob += b"".join(struct.pack("<H", record_bytes) for _ in range(variant_ct))
    for genovec in genovecs:
        blob += genovec.tobytes()
    path.write_bytes(bytes(blob))
    return categories


def _uleb128(value: int) -> bytes:
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            out.append(byte | 0x80)
        else:
            out.append(byte)
            return bytes(out)


def _encode_difflist(entries, id_bytes: int) -> bytes:
    """A sparse (sample, category) list in the layout the decoder walks.

    Order is fixed by the format and is not the order the pieces are computed
    in: count, then one base sample id per group of 64, then a byte of run
    length for every group but the last, then the packed categories, then the
    deltas for every entry that does not open a group.
    """
    blob = bytearray(_uleb128(len(entries)))
    if not entries:
        return bytes(blob)
    groups = [entries[index:index + 64] for index in range(0, len(entries), 64)]
    deltas = []
    for group in groups:
        run = bytearray()
        previous = group[0][0]
        for sample, _category in group[1:]:
            run += _uleb128(sample - previous)
            previous = sample
        deltas.append(bytes(run))
    for group in groups:
        blob += int(group[0][0]).to_bytes(id_bytes, "little")
    blob += bytes(len(run) for run in deltas[:-1])
    packed = bytearray((len(entries) + 3) // 4)
    for index, (_sample, category) in enumerate(entries):
        packed[index >> 2] |= (int(category) & 3) << ((index & 3) * 2)
    blob += packed
    for run in deltas:
        blob += run
    return bytes(blob)


def _encode_onebit(row: np.ndarray, sample_ct: int, id_bytes: int,
                   low: int, high: int) -> bytes:
    """A one-bit record: two categories in a bit array, the rest sparse.

    The mode byte carries the unordered pair as ``3 * low + high``; the
    decoder inverts it as ``low = (mode - 1) // 3``.
    """
    bits = bytearray((sample_ct + 7) // 8)
    entries = []
    for sample in range(sample_ct):
        category = int(row[sample])
        if category == high:
            bits[sample >> 3] |= 1 << (sample & 7)
        elif category != low:
            entries.append((sample, category))
    return (bytes([3 * low + high]) + bytes(bits)
            + _encode_difflist(entries, id_bytes))


def write_pgen_onebit(path: Path, categories: np.ndarray, pairs) -> np.ndarray:
    """Write a PGEN whose every record is vrtype 1, the one-bit form.

    `write_pgen` emits only plain genovecs, so nothing in this file used to
    reach the one-bit decoder at all -- and that is the form the real cohort
    stores 22.6% of its variants in, and the one whose expansion was rewritten
    to work eight samples at a time. A synthetic file is the only way to cover
    it without a multi-gigabyte fixture.
    """
    categories = np.ascontiguousarray(categories, dtype=np.uint8)
    variant_ct, sample_ct = categories.shape
    if variant_ct > VBLOCK_SIZE:
        raise ValueError("the fixture writer emits a single variant block")
    id_bytes = _bytes_per_sample_id(sample_ct)
    records = [_encode_onebit(categories[index], sample_ct, id_bytes, *pairs[index])
               for index in range(variant_ct)]
    if any(len(record) > 0xFFFF for record in records):
        raise ValueError("the fixture writer stores two-byte record lengths")
    header_ctrl = 0x01 | 0x04
    index_bytes = variant_ct * 1 + variant_ct * 2
    data_start = 12 + 8 * 1 + index_bytes

    blob = bytearray()
    blob += b"\x6c\x1b"
    blob += bytes([0x10])
    blob += struct.pack("<II", variant_ct, sample_ct)
    blob += bytes([header_ctrl])
    blob += struct.pack("<Q", data_start)
    blob += bytes([1] * variant_ct)                       # vrtype 1: one-bit
    blob += b"".join(struct.pack("<H", len(record)) for record in records)
    for record in records:
        blob += record
    path.write_bytes(bytes(blob))
    return categories


def _encode_difflist_against(base: np.ndarray, target: np.ndarray,
                             id_bytes: int) -> bytes:
    """A difflist carrying every sample where `target` differs from `base`."""
    differing = np.flatnonzero(base != target)
    return _encode_difflist([(int(sample), int(target[sample]))
                             for sample in differing], id_bytes)


def write_pgen_mixed(path: Path, categories: np.ndarray, forms) -> np.ndarray:
    """Write a PGEN mixing record forms, including LD-compressed ones.

    Forms 2 and 3 are expressed against the most recent record that was not
    itself LD-compressed, which is the only thing that makes a read starting
    part-way into the file walk backwards and replay. Neither of the other
    fixture writers here produces one, so that replay -- and the decode path
    that handles it -- had no coverage at all.

    `forms` is one of 0 (plain genovec), 1 (one-bit), 2 (LD), 3 (LD inverted)
    or 4 (difflist against hom-ref) per variant. The first must not be an LD
    form, since there would be nothing to express it against.
    """
    categories = np.ascontiguousarray(categories, dtype=np.uint8)
    variant_ct, sample_ct = categories.shape
    id_bytes = _bytes_per_sample_id(sample_ct)
    if forms[0] in (2, 3):
        raise ValueError("the first record cannot be LD-compressed")

    records = []
    ld_base = None
    for index, form in enumerate(forms):
        row = categories[index]
        if form == 0:
            record = pack_genovec(row, sample_ct).tobytes()
        elif form == 1:
            counts = np.bincount(row, minlength=4)
            low, high = sorted(int(value) for value in np.argsort(counts)[-2:])
            record = _encode_onebit(row, sample_ct, id_bytes, low, high)
        elif form == 4:
            record = _encode_difflist(
                [(int(sample), int(row[sample]))
                 for sample in np.flatnonzero(row != 0)], id_bytes)
        elif form in (2, 3):
            if ld_base is None:
                raise ValueError("an LD record with no base before it")
            # Form 3 inverts categories 0 and 2 *after* the difflist is
            # applied, so the difflist has to describe the inverted target.
            wanted = row
            if form == 3:
                wanted = row.copy()
                wanted[row == 0] = 2
                wanted[row == 2] = 0
            record = _encode_difflist_against(ld_base, wanted, id_bytes)
        else:
            raise ValueError(f"the fixture writer does not emit form {form}")
        records.append(record)
        if form not in (2, 3):
            ld_base = row

    if any(len(record) > 0xFFFF for record in records):
        raise ValueError("the fixture writer stores two-byte record lengths")
    header_ctrl = 0x01 | 0x04
    index_bytes = variant_ct * 1 + variant_ct * 2
    data_start = 12 + 8 * 1 + index_bytes

    blob = bytearray()
    blob += b"\x6c\x1b"
    blob += bytes([0x10])
    blob += struct.pack("<II", variant_ct, sample_ct)
    blob += bytes([header_ctrl])
    blob += struct.pack("<Q", data_start)
    blob += bytes(int(form) for form in forms)
    blob += b"".join(struct.pack("<H", len(record)) for record in records)
    for record in records:
        blob += record
    path.write_bytes(bytes(blob))
    return categories


def _reference(path: Path, variant_ct: int) -> np.ndarray:
    """Categories per sample from the Python reference reader.

    `read_genovec` returns the packed two-bit vector, not one value per
    sample, so the comparison has to unpack it.
    """
    with PgenDirectReader(path) as reader:
        return np.stack([
            unpack_genovec(reader.read_genovec(index), reader.sample_ct)
            for index in range(variant_ct)
        ])


@unittest.skipUnless(pgen_native.available(),
                     "the native PGEN decoder is not built on this host")
class NativePgenReaderTestCase(unittest.TestCase):
    def _fixture(self, root: Path, variants=40, samples=37, seed=3):
        rng = np.random.default_rng(seed)
        categories = rng.integers(0, 4, size=(variants, samples)).astype(np.uint8)
        path = root / "fixture.pgen"
        write_pgen(path, categories)
        return path, categories

    def test_fixture_is_a_valid_pgen_the_reference_reader_accepts(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            np.testing.assert_array_equal(
                _reference(path, categories.shape[0]), categories)
            self.assertTrue(file_scope(path).supported)

    def test_read_range_matches_the_reference_including_missing_calls(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            want = categories.astype(np.int8)
            want[categories == 3] = -9
            with NativePgenReader(path) as reader:
                out = np.empty(categories.shape, dtype=np.int8)
                reader.read_range(0, categories.shape[0], out)
            np.testing.assert_array_equal(out, want)
            # Missing calls must survive as -9, not as the code 3 they are
            # stored under; the two differ by 12 in every downstream mean.
            self.assertTrue((out == -9).any())

    def test_the_one_pass_and_fallback_remaps_agree(self):
        """A stale library takes a different route; it must give the same bytes.

        `read_range` emits signed hard calls straight from C when the library
        has `torchgwas_pgen_expand_hardcall`, and otherwise expands to
        categories and remaps in numpy. The fast route is 10.4x on that step, so
        the two will not be exercised equally in practice -- which is exactly
        why the slow one needs pinning against it rather than being assumed.
        """
        from torchgwas import pgen_native

        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory), variants=33,
                                             samples=41, seed=9)
            fast = np.empty(categories.shape, dtype=np.int8)
            slow = np.empty(categories.shape, dtype=np.int8)
            with NativePgenReader(path) as reader:
                reader.read_range(0, categories.shape[0], fast)
                with unittest.mock.patch.object(
                        pgen_native, "expand_hardcall_available",
                        return_value=False):
                    reader.read_range(0, categories.shape[0], slow)
            np.testing.assert_array_equal(fast, slow)
            want = categories.astype(np.int8)
            want[categories == 3] = -9
            np.testing.assert_array_equal(fast, want)

    def test_partial_ranges_agree_with_the_whole(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            whole = categories.astype(np.int8)
            whole[categories == 3] = -9
            with NativePgenReader(path) as reader:
                for start, end in ((0, 1), (5, 9), (17, 40), (39, 40), (12, 12)):
                    out = np.empty((end - start, categories.shape[1]), dtype=np.int8)
                    reader.read_range(start, end, out)
                    with self.subTest(start=start, end=end):
                        np.testing.assert_array_equal(out, whole[start:end])

    def test_packed_rows_are_the_raw_two_bit_codes_and_padding_is_zeroed(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            variants, samples = categories.shape
            genovec_bytes = (samples + 3) // 4
            width = ((genovec_bytes + 63) // 64) * 64
            with NativePgenReader(path) as reader:
                out = np.full((variants, width), 0xFF, dtype=np.uint8)
                reader.read_packed_range_into(0, variants, out)
            expected = np.stack([pack_genovec(row, samples) for row in categories])
            np.testing.assert_array_equal(out[:, :genovec_bytes], expected)
            # Undefined padding would make two runs of one scan differ byte for
            # byte even though the genotypes are identical.
            np.testing.assert_array_equal(out[:, genovec_bytes:], 0)

    def test_one_bit_records_match_the_reference_at_every_sample_count(self):
        """The one-bit form, across sample counts that straddle its new stride.

        The expansion now writes two genovec bytes per byte of the bit array,
        so the samples past the last whole bit byte take a different path from
        the ones before it. 8 has no tail at all, 37 leaves five samples, 71
        leaves seven -- the widest tail there is -- and 300 crosses into
        two-byte difflist sample ids as well.
        """
        rng = np.random.default_rng(11)
        pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        for sample_ct in (8, 37, 64, 71, 300):
            with self.subTest(samples=sample_ct):
                variant_ct = len(pairs) * 2
                chosen = [pairs[index % len(pairs)] for index in range(variant_ct)]
                categories = np.empty((variant_ct, sample_ct), dtype=np.uint8)
                for index, (low, high) in enumerate(chosen):
                    row = rng.choice([low, high], size=sample_ct).astype(np.uint8)
                    # Every other variant carries exceptions, so the difflist
                    # that follows the dense fill is not always empty; 30% of
                    # 300 samples is 90 entries, which is two difflist groups.
                    if index % 2:
                        others = [c for c in range(4) if c not in (low, high)]
                        where = rng.random(sample_ct) < 0.3
                        row[where] = rng.choice(others, size=int(where.sum()))
                    categories[index] = row
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "onebit.pgen"
                    write_pgen_onebit(path, categories, chosen)
                    np.testing.assert_array_equal(
                        _reference(path, variant_ct), categories)
                    want = categories.astype(np.int8)
                    want[categories == 3] = -9
                    with NativePgenReader(path) as reader:
                        out = np.empty(categories.shape, dtype=np.int8)
                        reader.read_range(0, variant_ct, out)
                    np.testing.assert_array_equal(out, want)

    def _ld_fixture(self, root: Path, samples=71, seed=5):
        """A file whose middle is a long run of LD-compressed records.

        The run is what forces a read starting inside it to walk back to the
        last non-LD record and replay forward, so every start position below
        exercises a different replay length.
        """
        forms = [0, 4, 1] + [2, 2, 2, 3, 2, 2, 2, 2] + [0, 2, 2, 4, 2, 1, 2, 2]
        rng = np.random.default_rng(seed)
        categories = rng.integers(0, 4, size=(len(forms), samples)).astype(np.uint8)
        # LD records are only compact when they are close to their base, and a
        # base far from its dependants is also the case most likely to expose a
        # replay bug, so keep them near without making them identical.
        for index, form in enumerate(forms):
            if form in (2, 3) and index:
                keep = rng.random(samples) < 0.9
                categories[index] = np.where(keep, categories[index - 1],
                                             categories[index])
        path = root / "mixed.pgen"
        write_pgen_mixed(path, categories, forms)
        return path, categories, forms

    def test_ld_compressed_records_replay_from_any_start(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories, forms = self._ld_fixture(Path(directory))
            variant_ct, samples = categories.shape
            np.testing.assert_array_equal(
                _reference(path, variant_ct), categories)
            whole = categories.astype(np.int8)
            whole[categories == 3] = -9
            with NativePgenReader(path) as reader:
                for start in range(variant_ct):
                    end = min(start + 3, variant_ct)
                    out = np.empty((end - start, samples), dtype=np.int8)
                    reader.read_range(start, end, out)
                    with self.subTest(start=start, form=forms[start]):
                        np.testing.assert_array_equal(out, whole[start:end])

    def test_packed_transport_replays_ld_into_a_padded_buffer(self):
        """The packed path decodes into the caller's buffer, prefix and all.

        This is the case the direct fill had to get right: the buffer has room
        for the requested variants only, so the records replayed ahead of
        `start` cannot be written into it, and the padding past the genovec
        still has to come back zeroed.

        The comparison is on unpacked calls, not raw bytes, because the slots
        of the final genovec byte past the last real sample are *not* defined
        by the format -- a one-bit record leaves them holding its low category
        rather than zero, and has always done so. Nothing reads them: the
        preparation kernel iterates `sample < samples` with the logical count,
        so those slots never reach the statistics. The sibling test above can
        compare raw bytes only because its fixture is entirely plain genovecs,
        where the bytes come from `pack_genovec` in the first place.
        """
        with tempfile.TemporaryDirectory() as directory:
            path, categories, _forms = self._ld_fixture(Path(directory))
            variant_ct, samples = categories.shape
            genovec_bytes = (samples + 3) // 4
            width = ((genovec_bytes + 63) // 64) * 64
            with NativePgenReader(path) as reader:
                for start in range(variant_ct):
                    end = min(start + 4, variant_ct)
                    out = np.full((end - start, width), 0xFF, dtype=np.uint8)
                    reader.read_packed_range_into(start, end, out)
                    calls = np.stack([unpack_genovec(row[:genovec_bytes], samples)
                                      for row in out])
                    with self.subTest(start=start):
                        np.testing.assert_array_equal(calls, categories[start:end])
                        # Bytes past the genovec are a real contract: the
                        # transport row is padded to 64 bytes and undefined
                        # padding would make two runs differ byte for byte.
                        np.testing.assert_array_equal(out[:, genovec_bytes:], 0)

    def test_a_sample_subset_selects_the_right_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            subset = np.asarray([0, 3, 4, 11, 30, 36], dtype=np.uint32)
            whole = categories.astype(np.int8)
            whole[categories == 3] = -9
            with NativePgenReader(path, sample_subset=subset) as reader:
                out = np.empty((categories.shape[0], subset.size), dtype=np.int8)
                reader.read_range(0, categories.shape[0], out)
            np.testing.assert_array_equal(out, whole[:, subset])

    def test_packed_transport_refuses_a_sample_subset(self):
        # It is refused rather than silently re-packed: the packed path exists
        # to hand the GPU the file's own bytes, and a subset is not those bytes.
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            subset = np.asarray([1, 2, 3], dtype=np.uint32)
            with NativePgenReader(path, sample_subset=subset) as reader:
                out = np.zeros((4, 64), dtype=np.uint8)
                with self.assertRaises(ValueError):
                    reader.read_packed_range_into(0, 4, out)

    def test_declared_dimensions_are_checked_against_the_header(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            with self.assertRaises(ValueError):
                NativePgenReader(path, raw_sample_ct=categories.shape[1] + 1)
            with self.assertRaises(ValueError):
                NativePgenReader(path, variant_ct=categories.shape[0] + 1)

    def test_out_of_bounds_ranges_raise(self):
        with tempfile.TemporaryDirectory() as directory:
            path, categories = self._fixture(Path(directory))
            variants, samples = categories.shape
            with NativePgenReader(path) as reader:
                out = np.empty((2, samples), dtype=np.int8)
                with self.assertRaises(IndexError):
                    reader.read_range(variants - 1, variants + 1, out)

    def test_backend_choice_and_auto_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = self._fixture(Path(directory))
            self.assertIsNone(native_backend_reason(path))
            # A hard-call-only file under 'auto' must resolve to hardcall and
            # the native backend: resolving it to dosage is what made every
            # PGEN require pgenlib, including files holding no dosages at all.
            mode, backend, _ = resolve_pgen_mode_and_backend(path, "auto")
            self.assertEqual((mode, backend), ("hardcall", "native"))
            mode, backend, _ = resolve_pgen_mode_and_backend(path, "dosage")
            self.assertEqual((mode, backend), ("dosage", "pgenlib"))

    def test_a_file_the_index_cannot_parse_falls_back_rather_than_raising(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "not.pgen"
            path.write_bytes(b"fake-pgen")
            self.assertIsNotNone(native_backend_reason(path))
            mode, backend, reason = resolve_pgen_mode_and_backend(path, "auto")
            self.assertEqual(backend, "pgenlib")
            self.assertIn("did not parse", reason)


class BackendSelectionEnvTestCase(unittest.TestCase):
    """The env override must force the answer, and refuse loudly when it cannot."""

    def test_forcing_pgenlib_needs_no_native_library(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "whatever.pgen"
            path.write_bytes(b"fake-pgen")
            with unittest.mock.patch.dict(
                    "os.environ", {"TORCHGWAS_PGEN_BACKEND": "pgenlib"}):
                mode, backend, _ = resolve_pgen_mode_and_backend(path, "auto")
            self.assertEqual(backend, "pgenlib")
            # Forcing the backend must not change what 'auto' believes about
            # the data: this file's index does not parse, so nothing licenses
            # calling it hard-call-only.
            self.assertEqual(mode, "dosage")

    def test_forcing_native_on_an_unreadable_file_raises(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "whatever.pgen"
            path.write_bytes(b"fake-pgen")
            with unittest.mock.patch.dict(
                    "os.environ", {"TORCHGWAS_PGEN_BACKEND": "native"}):
                with self.assertRaises(ValueError):
                    resolve_pgen_mode_and_backend(path, "auto")

    def test_an_unknown_backend_name_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "whatever.pgen"
            path.write_bytes(b"fake-pgen")
            with unittest.mock.patch.dict(
                    "os.environ", {"TORCHGWAS_PGEN_BACKEND": "sqlite"}):
                with self.assertRaises(ValueError):
                    resolve_pgen_mode_and_backend(path, "auto")


if __name__ == "__main__":
    unittest.main()
