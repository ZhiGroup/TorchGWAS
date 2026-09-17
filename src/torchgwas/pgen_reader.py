"""Direct PGEN reader for biallelic hard calls.

Reads PLINK 2 ``.pgen`` storage mode 0x10/0x11 without pgenlib, and returns
the packed two-bit genotype vector the file already stores rather than one
byte per genotype. That is the representation the CUDA scan wants, so the
decode and the transport format agree and nothing has to be expanded and
repacked on the way through.

Scope is deliberate. Records carrying multiallelic, phase or dosage tracks
are refused, not guessed at; callers fall back to pgenlib for those. Every
benchmark file we target is biallelic hard calls, where the record types are:

====  ==========================================================
type  meaning
====  ==========================================================
0     genotype vector stored directly, no compression
1     "1-bit" form: a mode byte, a bitarray over the two most
      common categories, then a difflist for the remainder
2     LD-compressed: a difflist against the most recent
      non-LD-compressed variant
3     as 2, then categories 0 and 2 are swapped
4     difflist over the samples outside category 0
6     difflist over the samples outside category 2
7     difflist over the samples outside category 3
====  ==========================================================

Genotype categories are the on-disk two-bit codes: 0 homozygous reference,
1 heterozygous, 2 homozygous alternate, 3 missing.

The specification is at https://github.com/chrchang/plink-ng, pgen_spec.

Two details below were established by differential testing against pgenlib
rather than read from that specification, because neither appears in it.
They are recorded here because nobody will recover them from the spec:

Difflist group sizes are for skipping, not for advancing.
    A difflist stores one byte per group giving that group's delta-run size,
    so a reader looking for one sample can jump between groups. The runs are
    contiguous, so a reader decoding the whole list must ignore those sizes
    and simply keep reading. Honouring them as an advance produces sample
    ids past the end of the cohort.

The "1-bit" mode byte encodes an unordered category pair as 3*low + high.
    The record names two categories and its bitarray marks the
    higher-numbered one. Enumerating the six pairs with low < high as
    ``3 * low + high`` reproduces every mode byte observed, and is collision
    free over those six, so the pair inverts as ``low = (mode - 1) // 3``
    and ``high = mode - 3 * low``. Observed directly for the (0,1) and (1,2)
    pairs; the inversion is asserted rather than assumed, so a file using a
    pair this inference gets wrong raises instead of returning silent
    nonsense.
"""

from __future__ import annotations

import os
import struct
from dataclasses import dataclass

import numpy as np

PGEN_MAGIC = b"\x6c\x1b"
VBLOCK_SIZE = 1 << 16
DIFFLIST_GROUP = 64

# vrtype bits above the low three select optional tracks this reader refuses.
TRACK_MULTIALLELIC = 0x08
# Attribution established from fixtures rather than from the specification:
# a plink2 file written with dosage and no phasing carries record types
# 0x40, 0x41 and 0x60, so 0x20 belongs to the dosage encoding; a phased file
# with no dosage carries 0x10 alone. The union is unchanged either way, so
# this affects only which track a refusal names.
TRACK_PHASE = 0x10
TRACK_DOSAGE = 0x20 | 0x40 | 0x80
UNSUPPORTED_TRACKS = TRACK_MULTIALLELIC | TRACK_PHASE | TRACK_DOSAGE

# Two-bit code assumed by records that store only their exceptions.
BACKGROUND_FOR_TYPE = {4: 0, 6: 2, 7: 3}


class PgenFormatError(ValueError):
    """The file is not a PGEN this reader supports."""


def _bytes_per_sample_id(sample_ct: int) -> int:
    """Width of a difflist base sample id, as pgenlib sizes it."""
    if sample_ct <= 0:
        raise PgenFormatError("sample count must be positive")
    highest = sample_ct - 1
    if highest < (1 << 8):
        return 1
    if highest < (1 << 16):
        return 2
    if highest < (1 << 24):
        return 3
    return 4


def vrtype_index_layout(header_ctrl: int) -> tuple[int, int]:
    """Return (vrtype_bits, record_length_bytes) declared by ``header_ctrl``.

    Bits 0-1 give the record-length width minus one; bit 2 widens the record
    type from a nibble to a byte. Confirmed against four files whose layouts
    differ, each cross-checked by requiring the record lengths to sum to the
    bytes remaining after the header and index.
    """
    length_bytes = (header_ctrl & 0x03) + 1
    vrtype_bits = 8 if header_ctrl & 0x04 else 4
    return vrtype_bits, length_bytes


def _decode_lengths(raw: bytes, count: int, width: int) -> np.ndarray:
    if len(raw) != count * width:
        raise PgenFormatError("truncated record-length array")
    if width in (1, 2, 4):
        dtype = {1: np.uint8, 2: "<u2", 4: "<u4"}[width]
        return np.frombuffer(raw, dtype=dtype, count=count).astype(np.uint32)
    # Three-byte lengths: pad each to four and reinterpret.
    padded = np.zeros((count, 4), dtype=np.uint8)
    padded[:, :3] = np.frombuffer(raw, dtype=np.uint8).reshape(count, 3)
    return padded.view("<u4").reshape(count).astype(np.uint32)


def _uleb128(buffer: memoryview, position: int) -> tuple[int, int]:
    value = 0
    shift = 0
    while True:
        if position >= len(buffer):
            raise PgenFormatError("truncated variable-length integer")
        byte = buffer[position]
        position += 1
        value |= (byte & 0x7F) << shift
        if not byte & 0x80:
            return value, position
        shift += 7
        if shift > 63:
            raise PgenFormatError("variable-length integer too long")


def unpack_genovec(genovec: np.ndarray, sample_ct: int) -> np.ndarray:
    """Expand a packed two-bit vector to one uint8 category per sample."""
    expanded = np.empty(((genovec.size) * 4,), dtype=np.uint8)
    expanded[0::4] = genovec & 0x03
    expanded[1::4] = (genovec >> 2) & 0x03
    expanded[2::4] = (genovec >> 4) & 0x03
    expanded[3::4] = (genovec >> 6) & 0x03
    return expanded[:sample_ct]


def pack_genovec(categories: np.ndarray, sample_ct: int) -> np.ndarray:
    """Inverse of :func:`unpack_genovec`, used by tests and fixtures."""
    padded_ct = (sample_ct + 3) & ~3
    padded = np.zeros(padded_ct, dtype=np.uint8)
    padded[:sample_ct] = categories[:sample_ct]
    return (
        padded[0::4]
        | (padded[1::4] << 2)
        | (padded[2::4] << 4)
        | (padded[3::4] << 6)
    ).astype(np.uint8)


def ld_safe_start(vrtypes, index: int) -> int:
    """Earliest variant a decode must begin at to have an LD base for ``index``.

    Forms 2 and 3 are expressed against the most recent record that was not
    itself LD-compressed, so a decode entering the file at an arbitrary
    variant has to back up to that record and replay forward, discarding the
    rows before the one it wanted. Variant blocks never begin with an
    LD-compressed record, so the walk is bounded by the block; in practice it
    is a handful of variants, because the LD forms are a minority.
    """
    start = int(index)
    while start > 0 and (int(vrtypes[start]) & 0x07) in (2, 3):
        start -= 1
    return start


@dataclass
class PgenHeader:
    storage_mode: int
    variant_ct: int
    sample_ct: int
    header_ctrl: int
    vrtype_bits: int
    length_bytes: int
    vblock_offsets: np.ndarray
    vrtypes: np.ndarray
    record_offsets: np.ndarray
    record_lengths: np.ndarray


def read_header(path: str | os.PathLike) -> PgenHeader:
    """Parse the header, the variant-block index and every record locator.

    The per-block record lengths are validated against the distance between
    consecutive block offsets, so a misread index is detected here rather
    than surfacing later as a corrupt genotype.
    """
    with open(path, "rb") as handle:
        head = handle.read(12)
        if head[:2] != PGEN_MAGIC:
            raise PgenFormatError(f"{path} is not a PGEN file")
        storage_mode = head[2]
        if storage_mode == 0x01:
            raise PgenFormatError(
                "this PGEN is a plink1 BED layout (storage mode 0x01); its "
                "two-bit codes are BED's, not PGEN's, so read it with the BED "
                "reader rather than remapping it here"
            )
        if storage_mode not in (0x10, 0x11):
            raise PgenFormatError(
                f"unsupported PGEN storage mode {hex(storage_mode)}; this reader "
                "handles the variable-width modes 0x10 and 0x11, and plink2 has "
                "not been observed to emit the fixed-width modes"
            )
        variant_ct, sample_ct = struct.unpack("<II", head[3:11])
        header_ctrl = head[11]
        vblock_ct = (variant_ct + VBLOCK_SIZE - 1) // VBLOCK_SIZE
        vblock_offsets = np.frombuffer(
            handle.read(8 * vblock_ct), dtype="<u8", count=vblock_ct
        )

        # The record index geometry is declared, not fixed. Small files use
        # one-byte record lengths, and files whose record types need the high
        # bits (dosage, phase, multiallelic) widen the type field to a byte.
        # Assuming the benchmark file's 4-bit/2-byte layout made every other
        # file decode as truncated garbage.
        vrtype_bits, length_bytes = vrtype_index_layout(header_ctrl)
        vrtypes = np.empty(variant_ct, dtype=np.uint8)
        lengths = np.empty(variant_ct, dtype=np.uint32)
        offsets = np.empty(variant_ct, dtype=np.uint64)
        for block in range(vblock_ct):
            first = block * VBLOCK_SIZE
            in_block = min(VBLOCK_SIZE, variant_ct - first)
            raw = handle.read((in_block * vrtype_bits + 7) // 8)
            packed = np.frombuffer(raw, dtype=np.uint8)
            if vrtype_bits == 4:
                block_types = np.empty(len(packed) * 2, dtype=np.uint8)
                block_types[0::2] = packed & 0x0F
                block_types[1::2] = packed >> 4
            else:
                block_types = packed
            vrtypes[first:first + in_block] = block_types[:in_block]

            raw_lengths = handle.read(in_block * length_bytes)
            block_lengths = _decode_lengths(raw_lengths, in_block, length_bytes)
            lengths[first:first + in_block] = block_lengths
            starts = np.empty(in_block, dtype=np.uint64)
            starts[0] = vblock_offsets[block]
            np.cumsum(block_lengths[:-1], dtype=np.uint64, out=starts[1:])
            starts[1:] += vblock_offsets[block]
            offsets[first:first + in_block] = starts
            if block + 1 < vblock_ct:
                span = int(vblock_offsets[block + 1]) - int(vblock_offsets[block])
                if span != int(block_lengths.sum()):
                    raise PgenFormatError(
                        f"variant block {block} record lengths sum to "
                        f"{int(block_lengths.sum())} but the block spans {span} bytes"
                    )
    return PgenHeader(
        storage_mode=storage_mode,
        variant_ct=variant_ct,
        sample_ct=sample_ct,
        header_ctrl=header_ctrl,
        vrtype_bits=vrtype_bits,
        length_bytes=length_bytes,
        vblock_offsets=vblock_offsets,
        vrtypes=vrtypes,
        record_offsets=offsets,
        record_lengths=lengths,
    )


@dataclass
class ScopeReport:
    """Whether a whole file is within the direct reader's scope.

    Per-variant refusal is not a usable policy on its own. The real 135 GB
    dosage PGEN has 8,930,997 of 8,931,083 variants carrying a dosage track
    and 86 that do not, so a reader that refuses variant by variant would
    decode 86 of them and reject the rest. The backend has to be chosen once
    for the file, which is what this reports.
    """

    supported: bool
    variant_ct: int
    unsupported_variants: int
    reserved_form_variants: int
    multiallelic_variants: int
    phased_variants: int
    dosage_variants: int
    reason: str = ""


def scope_from_vrtypes(vrtypes) -> ScopeReport:
    """Classify every record type without touching the genotype data."""
    types = np.asarray(vrtypes, dtype=np.uint8)
    multiallelic = int(np.count_nonzero(types & TRACK_MULTIALLELIC))
    phased = int(np.count_nonzero(types & TRACK_PHASE))
    dosage = int(np.count_nonzero(types & TRACK_DOSAGE))
    tracked = int(np.count_nonzero(types & UNSUPPORTED_TRACKS))
    reserved = int(np.count_nonzero((types & 0x07) == 5))
    unsupported = int(
        np.count_nonzero((types & UNSUPPORTED_TRACKS) | ((types & 0x07) == 5))
    )
    reasons = []
    if dosage:
        reasons.append(f"{dosage} variants carry a dosage track")
    if phased:
        reasons.append(f"{phased} variants carry a phase track")
    if multiallelic:
        reasons.append(f"{multiallelic} variants are multiallelic")
    if reserved:
        reasons.append(f"{reserved} variants use a reserved record form")
    return ScopeReport(
        supported=unsupported == 0,
        variant_ct=int(types.size),
        unsupported_variants=unsupported,
        reserved_form_variants=reserved,
        multiallelic_variants=multiallelic,
        phased_variants=phased,
        dosage_variants=dosage,
        reason="; ".join(reasons),
    )


def file_scope(path) -> ScopeReport:
    """Report scope for a PGEN without decoding any genotypes.

    Reads only the header and record index, so it is cheap enough to run
    before choosing a backend even on a whole-cohort file.
    """
    return scope_from_vrtypes(read_header(path).vrtypes)


class PgenDirectReader:
    """Sequential two-bit genotype reader for biallelic hard-call PGEN."""

    def __init__(self, path: str | os.PathLike) -> None:
        self.path = str(path)
        self.header = read_header(path)
        self.sample_ct = self.header.sample_ct
        self.variant_ct = self.header.variant_ct
        self.genovec_bytes = (self.sample_ct + 3) // 4
        self._id_bytes = _bytes_per_sample_id(self.sample_ct)
        self._handle = open(self.path, "rb")
        # Records of type 2 and 3 are expressed against the most recent
        # record that was not itself LD-compressed.
        self._ld_base: np.ndarray | None = None
        self._ld_base_index: int | None = None
        self.scope = scope_from_vrtypes(self.header.vrtypes)

    def close(self) -> None:
        if getattr(self, "_handle", None) is not None:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> "PgenDirectReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _record(self, index: int) -> memoryview:
        offset = int(self.header.record_offsets[index])
        length = int(self.header.record_lengths[index])
        self._handle.seek(offset)
        payload = self._handle.read(length)
        if len(payload) != length:
            raise PgenFormatError(f"variant {index} record truncated")
        return memoryview(payload)

    def _parse_difflist(
        self, buffer: memoryview, position: int, with_genotypes: bool
    ) -> tuple[np.ndarray, np.ndarray | None, int]:
        """Decode a sparse (sample id, category) list.

        Layout: the entry count, then one base sample id per group of 64,
        then the byte size of every group's delta run except the last, then
        the packed categories when the record carries them, then the
        variable-length within-group id increments.
        """
        entry_ct, position = _uleb128(buffer, position)
        if entry_ct == 0:
            empty_ids = np.empty(0, dtype=np.int64)
            return empty_ids, (np.empty(0, dtype=np.uint8) if with_genotypes else None), position
        group_ct = (entry_ct + DIFFLIST_GROUP - 1) // DIFFLIST_GROUP

        bases = np.zeros(group_ct, dtype=np.int64)
        for group in range(group_ct):
            chunk = bytes(buffer[position:position + self._id_bytes])
            if len(chunk) != self._id_bytes:
                raise PgenFormatError("truncated difflist base id")
            bases[group] = int.from_bytes(chunk, "little")
            position += self._id_bytes

        group_sizes = np.zeros(group_ct, dtype=np.int64)
        if group_ct > 1:
            raw = buffer[position:position + group_ct - 1]
            if len(raw) != group_ct - 1:
                raise PgenFormatError("truncated difflist group sizes")
            group_sizes[:-1] = np.frombuffer(bytes(raw), dtype=np.uint8)
            position += group_ct - 1

        genotypes = None
        if with_genotypes:
            packed_ct = (entry_ct + 3) // 4
            raw = buffer[position:position + packed_ct]
            if len(raw) != packed_ct:
                raise PgenFormatError("truncated difflist genotypes")
            genotypes = unpack_genovec(
                np.frombuffer(bytes(raw), dtype=np.uint8), entry_ct
            )
            position += packed_ct

        ids = np.empty(entry_ct, dtype=np.int64)
        cursor = position
        written = 0
        for group in range(group_ct):
            remaining = min(DIFFLIST_GROUP, entry_ct - written)
            current = bases[group]
            ids[written] = current
            group_start = cursor  # retained for the skipping path below
            for _ in range(remaining - 1):
                delta, cursor = _uleb128(buffer, cursor)
                current += delta
                written += 1
                ids[written] = current
            written += 1
            # group_start is deliberately unused: see the module docstring on
            # difflist group sizes. The runs are contiguous, so advancing by
            # the recorded size would skip live delta bytes.
            del group_start
        return ids, genotypes, cursor

    def _decode(self, index: int) -> np.ndarray:
        vrtype = int(self.header.vrtypes[index])
        if vrtype & UNSUPPORTED_TRACKS:
            raise PgenFormatError(
                f"variant {index} carries a multiallelic, phase or dosage track "
                f"(vrtype {vrtype}); use the pgenlib backend for this file"
            )
        form = vrtype & 0x07
        buffer = self._record(index)

        if form == 0:
            if len(buffer) < self.genovec_bytes:
                raise PgenFormatError(f"variant {index} genotype vector truncated")
            return np.frombuffer(
                bytes(buffer[:self.genovec_bytes]), dtype=np.uint8
            ).copy()

        if form in BACKGROUND_FOR_TYPE:
            background = BACKGROUND_FOR_TYPE[form]
            categories = np.full(self.sample_ct, background, dtype=np.uint8)
            ids, genotypes, _ = self._parse_difflist(buffer, 0, with_genotypes=True)
            if ids.size:
                categories[ids] = genotypes
            return pack_genovec(categories, self.sample_ct)

        if form in (2, 3):
            if self._ld_base is None:
                raise PgenFormatError(
                    f"variant {index} is LD-compressed but no base variant precedes it"
                )
            categories = unpack_genovec(self._ld_base, self.sample_ct).copy()
            ids, genotypes, _ = self._parse_difflist(buffer, 0, with_genotypes=True)
            if ids.size:
                categories[ids] = genotypes
            if form == 3:
                swapped = categories.copy()
                swapped[categories == 0] = 2
                swapped[categories == 2] = 0
                categories = swapped
            return pack_genovec(categories, self.sample_ct)

        if form == 1:
            # The mode byte names an unordered pair of categories, and the
            # bitarray marks the higher-numbered one. Enumerating the pairs
            # with clear < set as 3*clear + set is collision free over the
            # six possibilities, so the pair inverts exactly.
            mode = buffer[0]
            low_category = (mode - 1) // 3
            high_category = mode - 3 * low_category
            if not 0 <= low_category < high_category <= 3:
                raise PgenFormatError(
                    f"variant {index} has unreadable 1-bit mode byte {mode}"
                )
            bitarray_bytes = (self.sample_ct + 7) // 8
            raw = buffer[1:1 + bitarray_bytes]
            if len(raw) != bitarray_bytes:
                raise PgenFormatError(f"variant {index} bitarray truncated")
            bits = np.unpackbits(
                np.frombuffer(bytes(raw), dtype=np.uint8), bitorder="little"
            )[:self.sample_ct]
            categories = np.where(bits, high_category, low_category).astype(np.uint8)
            ids, genotypes, _ = self._parse_difflist(
                buffer, 1 + bitarray_bytes, with_genotypes=True
            )
            if ids.size:
                categories[ids] = genotypes
            return pack_genovec(categories, self.sample_ct)

        raise PgenFormatError(f"variant {index} uses reserved record type {form}")

    def read_genovec(self, index: int) -> np.ndarray:
        """Return the packed two-bit vector for one variant.

        Reading strictly forward keeps the LD base in hand. A backward or
        skipping read re-establishes it by replaying from the most recent
        non-LD-compressed variant, which the format guarantees is close.
        """
        if not 0 <= index < self.variant_ct:
            raise IndexError(index)
        form = int(self.header.vrtypes[index]) & 0x07
        if form in (2, 3) and self._ld_base_index != index - 1:
            self._seek_ld_base(index)
        genovec = self._decode(index)
        if form not in (2, 3):
            self._ld_base = genovec
        self._ld_base_index = index
        return genovec

    def _seek_ld_base(self, index: int) -> None:
        start = index - 1
        while start >= 0 and (int(self.header.vrtypes[start]) & 0x07) in (2, 3):
            start -= 1
        if start < 0:
            raise PgenFormatError(f"no LD base precedes variant {index}")
        for position in range(start, index):
            form = int(self.header.vrtypes[position]) & 0x07
            genovec = self._decode(position)
            if form not in (2, 3):
                self._ld_base = genovec
            self._ld_base_index = position

    def iter_genovecs(self, start: int = 0, stop: int | None = None):
        stop = self.variant_ct if stop is None else stop
        for index in range(start, stop):
            yield index, self.read_genovec(index)
