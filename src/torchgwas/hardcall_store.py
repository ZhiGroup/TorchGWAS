"""A zstd-compressed store for hard calls, in the layout the scan already reads.

Every design choice here was made by measurement rather than taste, and the
measurements are in `benchmarks/direct_hardcall_compressibility.py` and
`benchmarks/direct_compression_decision.py`.

**What it is worth, and where.** Compression only pays where the READ binds, so
this is a LOW-K feature. MEASURED at full genome and full cohort, cold
everything, two interleaved rounds against the `.bed` it was built from:

    K = 1     bed 16.26 s   store  7.58 s    2.15x faster
    K = 512   bed 22.61 s   store 21.63 s    1.05x (inside the spread)

At K = 512 it is worth nothing, and claiming otherwise would be dishonest: both
formats converge on a floor set by the GEMM and the write. K = 1-8 is not a
corner case -- it is where every competing tool actually runs, since fastGWA's
`--mpheno` takes one column and plink2 was measured at K = 1.

The store is also far more STABLE than the bed (7.60/7.55 against 15.59/16.92
across rounds), because it reads 10.4 GB where the bed reads 79.0 GB. On a
shared machine that predictability can matter as much as the mean.

A first pass of the cost model predicted **7.63x** here. It was wrong: the
model has no decompression term at all, and expanding 79 GB of output per scan
is what eats the difference. Do not quote the model for this without the
measurement beside it.

**Why zstd and not deflate**, measured on 50,000 real variants:

    codec        ratio   decompress (1 thread)
    zstd 3       6.400x       1,342 MB/s
    zstd 10      6.740x       1,487 MB/s
    deflate 1    2.467x         224 MB/s
    deflate 6    2.880x         259 MB/s

Deflate was the interesting candidate because nvCOMP can decompress it ON THE
GPU, which would keep host cores out of the decode path entirely -- the bgen
reader already does this. On ratio and on CPU decode speed it loses on both
axes: 2.47x is barely better than the 2.20x pgen achieves for free, and it
decompresses 6x slower. zstd at level 3 gives 21.5 GB/s across 16 threads, far
above what the storage delivers, so CPU decode keeps up.

**That rejection used the wrong criterion, and it should be revisited.** Ratio
was never the binding constraint -- HOST CPU is. Measured with the sumstats
writer off, this store is 1.33x faster than the bed at K=512 (10.81 s against
14.43 s); with the writer on it is 1.4% slower, because zstd decode and the
writer fight for the same cores. bed decodes on the GPU and does not compete.
For the identical 36.68 GB write, bed pays +8.91 s and this store pays
+12.85 s. So a GPU-decoded deflate store would read 79/2.47 = 32 GB in ~5.8 s
against bed's ~14.2 s while taking nothing from the writer -- worse ratio,
better outcome. It is the one design that could win at high K, and it is not
built.

Level 3 and not 10: 6.400x against 6.740x is a 5% gain in size for 5.4x the
compression cost (433 MB/s against 80 MB/s), and encoding a full genome is
already the expensive part of adopting the format.

**Why genome order.** Sorting variants by MAF was measured and is WORSE --
86.8 MB against 67.6 MB for the same variants, a factor of 0.78. Neighbouring
variants are in linkage disequilibrium and therefore correlated, so a frame of
them shares structure the compressor exploits; MAF order destroys that locality
and offers nothing in its place. Genome order is also the order the scan reads
in, so there is no permutation to carry.

**Why the bed code convention.** Frames hold exactly the bytes a `.bed` holds
(0 = hom A1, 1 = missing, 2 = het, 3 = hom A2), so the existing CPU lookup
table and the GPU decode kernels work on them unchanged. The store is a
transport, not a new genotype encoding, and keeping it that way means it cannot
introduce a correctness difference.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

MAGIC = b"TGHC0001"

# Measured optimum, see the module docstring. Exposed because the encoder takes
# it as an argument and the reason for the value belongs next to the value.
DEFAULT_LEVEL = 3

# Variants per independently decodable frame. A read of any range decompresses
# only the frames it overlaps, so this trades edge waste (a chunk that straddles
# frames pays for the whole covering frames) against ratio (bigger frames give
# the compressor more context). 2,048 divides the 4,096 chunk the planner
# currently hands out; a caller running at the measured 1,024 optimum should
# pass 1,024 so frames and chunks align exactly.
DEFAULT_FRAME_VARIANTS = 2048


def _paths(prefix: str | Path) -> tuple[Path, Path]:
    prefix = Path(prefix)
    return Path(f"{prefix}.zhc"), Path(f"{prefix}.zhc.json")


def encode_hardcall_store(
    read_packed,
    variants: int,
    samples: int,
    output_prefix: str | Path,
    *,
    frame_variants: int = DEFAULT_FRAME_VARIANTS,
    level: int = DEFAULT_LEVEL,
    workers: int = 8,
    progress=None,
) -> dict:
    """Compress packed two-bit rows into independently decodable zstd frames.

    `read_packed(start, end)` must return a contiguous `(end - start,
    bytes_per_variant)` uint8 array in bed code order -- exactly what
    `PlinkBedGenotype` reads and what `PgenDirectReader.read_genovec` packs.

    Frames are compressed in parallel and written in variant order, so the file
    is byte-identical regardless of `workers`. That matters: a store whose
    contents depend on the thread count cannot be checksummed against a
    reference, and reproducibility is the cheapest correctness check available.
    """
    import zstandard as zstd

    if frame_variants <= 0 or workers <= 0:
        raise ValueError("frame_variants and workers must be positive")
    bytes_per_variant = (samples + 3) // 4
    data_path, manifest_path = _paths(output_prefix)
    data_path.parent.mkdir(parents=True, exist_ok=True)

    starts = list(range(0, variants, frame_variants))
    offsets = np.zeros(len(starts), dtype=np.int64)
    sizes = np.zeros(len(starts), dtype=np.int64)
    counts = np.zeros(len(starts), dtype=np.int32)

    def compress_one(index: int) -> tuple[int, bytes, int]:
        start = starts[index]
        stop = min(start + frame_variants, variants)
        block = np.ascontiguousarray(read_packed(start, stop), dtype=np.uint8)
        if block.shape != (stop - start, bytes_per_variant):
            raise ValueError(
                f"frame {index}: expected {(stop - start, bytes_per_variant)}, "
                f"got {block.shape}")
        blob = zstd.ZstdCompressor(level=level).compress(block.tobytes())
        return index, blob, stop - start

    position = 0
    raw_total = 0
    # A bounded window of in-flight frames, not a map over all of them: at full
    # genome this is 4,361 frames of 18 MB raw each, and submitting everything
    # at once would hold every compressed blob in memory waiting for the
    # writer. The deque is drained in submission order, which is what makes the
    # output byte-identical regardless of `workers`.
    from collections import deque

    window = max(workers * 2, 4)
    with data_path.open("wb", buffering=8 << 20) as handle, \
            ThreadPoolExecutor(max_workers=workers) as pool:
        inflight: deque = deque()

        def retire() -> None:
            nonlocal position, raw_total
            index, blob, count = inflight.popleft().result()
            handle.write(blob)
            offsets[index] = position
            sizes[index] = len(blob)
            counts[index] = count
            position += len(blob)
            raw_total += count * bytes_per_variant
            if progress is not None:
                progress(index + 1, len(starts))

        for index in range(len(starts)):
            inflight.append(pool.submit(compress_one, index))
            if len(inflight) >= window:
                retire()
        while inflight:
            retire()

    manifest = {
        "magic": MAGIC.decode("ascii"),
        "samples": int(samples),
        "variants": int(variants),
        "bytes_per_variant": int(bytes_per_variant),
        "frame_variants": int(frame_variants),
        "codec": "zstd",
        "level": int(level),
        "code_order": "bed",
        "frames": len(starts),
        "compressed_bytes": int(position),
        "raw_bytes": int(raw_total),
        "ratio": (raw_total / position) if position else 0.0,
        "frame_offsets": offsets.tolist(),
        "frame_sizes": sizes.tolist(),
        "frame_counts": counts.tolist(),
    }
    manifest_path.write_text(json.dumps(manifest))
    return manifest


class HardcallStore:
    """Read packed two-bit rows back out, decompressing only what is asked for.

    The public surface deliberately mirrors the packed side of
    `PlinkBedGenotype` -- `shape`, `bytes_per_variant`, `read_packed` and
    `read_packed_into` -- so the existing pinned-ring loader can be pointed at
    a store without the pipeline learning a third genotype path.
    """

    def __init__(self, prefix: str | Path) -> None:
        data_path, manifest_path = _paths(prefix)
        manifest = json.loads(Path(manifest_path).read_text())
        if manifest.get("magic") != MAGIC.decode("ascii"):
            raise ValueError(f"{manifest_path} is not a hard-call store")
        if manifest.get("codec") != "zstd":
            raise ValueError(f"unsupported codec {manifest.get('codec')!r}")
        self.manifest = manifest
        self.samples = int(manifest["samples"])
        self.variants = int(manifest["variants"])
        self.bytes_per_variant = int(manifest["bytes_per_variant"])
        self.frame_variants = int(manifest["frame_variants"])
        self._offsets = np.asarray(manifest["frame_offsets"], dtype=np.int64)
        self._sizes = np.asarray(manifest["frame_sizes"], dtype=np.int64)
        self._counts = np.asarray(manifest["frame_counts"], dtype=np.int64)
        self._starts = np.concatenate(
            ([0], np.cumsum(self._counts)[:-1])).astype(np.int64)
        self._fd = os.open(str(data_path), os.O_RDONLY)
        self.path = Path(data_path)

    @property
    def shape(self) -> tuple[int, int]:
        return self.samples, self.variants

    def close(self) -> None:
        if getattr(self, "_fd", None) is not None:
            os.close(self._fd)
            self._fd = None

    def __enter__(self) -> "HardcallStore":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # noqa: BLE001 - interpreter teardown
            pass

    def _frame_blob(self, index: int) -> memoryview:
        """Read one compressed frame, returning a view rather than a copy.

        This used to end with `bytes(out)`, which copied the whole frame for
        nothing -- `zstandard` accepts any buffer. At full genome that copy
        alone was 10.4 GB of pointless memory traffic per scan.
        """
        offset = int(self._offsets[index])
        want = int(self._sizes[index])
        out = bytearray(want)
        view = memoryview(out)
        got = 0
        while got < want:
            amount = os.preadv(self._fd, [view[got:]], offset + got)
            if amount == 0:
                raise OSError(f"unexpected EOF in {self.path} at {offset + got}")
            got += amount
        return view

    def read_packed_into(self, target, start: int, end: int) -> int:
        """Decompress variants [start, end) into `target`, a writable buffer.

        Writes straight into the caller's memory -- typically a slice of a
        pinned ring buffer -- so the compressed path allocates no more per
        chunk than the uncompressed one does. Returns the number of bytes.
        """
        import zstandard as zstd

        self._check_range(start, end)
        # An empty range is legitimate -- the last chunk of a range that
        # divides evenly asks for it -- and must not reach the memoryview cast,
        # which rejects a zero in the shape.
        if end == start:
            return 0
        view = memoryview(target).cast("B")
        stride = self.bytes_per_variant
        needed = (end - start) * stride
        if len(view) < needed:
            raise ValueError(f"target holds {len(view)} bytes, need {needed}")

        decompressor = zstd.ZstdDecompressor()
        first = int(np.searchsorted(self._starts, start, side="right")) - 1
        written = 0
        index = max(first, 0)
        while index < len(self._starts):
            frame_start = int(self._starts[index])
            frame_stop = frame_start + int(self._counts[index])
            if frame_start >= end:
                break
            blob = self._frame_blob(index)
            take_from = max(start, frame_start) - frame_start
            take_to = min(end, frame_stop) - frame_start
            span = (take_to - take_from) * stride
            if take_from == 0 and take_to == frame_stop - frame_start:
                # The whole frame is wanted, which is every frame but the two
                # at the ends. Decompress STRAIGHT into the caller's buffer --
                # typically a pinned ring slot. The obvious
                # `decompress()`-then-copy form allocates the full output and
                # then copies it, which at full genome is 79 GB allocated and
                # 79 GB copied per scan on top of the decompression itself.
                reader = decompressor.stream_reader(blob)
                try:
                    got = reader.readinto(view[written:written + span])
                finally:
                    reader.close()
                if got != span:
                    raise OSError(
                        f"frame {index} yielded {got} bytes, expected {span}")
            else:
                # An edge frame, at most two per read: decompress and slice.
                raw = decompressor.decompress(blob)
                view[written:written + span] = memoryview(raw)[
                    take_from * stride:take_to * stride]
            written += span
            index += 1
        if written != needed:
            raise OSError(
                f"store yielded {written} bytes for [{start}, {end}), "
                f"expected {needed}")
        return written

    def _check_range(self, start: int, end: int) -> None:
        if not (0 <= start <= end <= self.variants):
            raise IndexError(
                f"invalid variant range [{start}, {end}) for {self.variants}")

    def read_packed(self, start: int, end: int) -> np.ndarray:
        # Checked BEFORE allocating: a reversed range would otherwise reach
        # numpy as a negative dimension and raise ValueError, so a caller
        # catching IndexError for a bad range would miss it.
        self._check_range(start, end)
        out = np.empty((end - start, self.bytes_per_variant), dtype=np.uint8)
        self.read_packed_into(out, start, end)
        return out
