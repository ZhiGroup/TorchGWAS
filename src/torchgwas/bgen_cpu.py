"""ctypes binding for the C BGEN Layout 2 probability decoder.

Loading is optional and cached. If the library is missing or stale the caller
falls back to the numpy path in `bgen.py`, which stays the reference
implementation: the two are compared element for element in the tests, and the
C decoder refuses exactly the records the Python one refuses.

ctypes releases the GIL for the duration of a CDLL call, which is the point of
moving the decode here -- reader threads overlap instead of queuing.
"""
from __future__ import annotations

import ctypes as ct
import os
from functools import lru_cache
from pathlib import Path

import numpy as np

_ABI = 2


class BgenDecodeError(ValueError):
    """A record the C decoder refused, carrying the decoder's own message."""


@lru_cache(maxsize=1)
def load_cpu_decoder():
    path = Path(__file__).with_name('native') / 'libtorchgwas_bgen_cpu.so'
    if not path.is_file():
        raise ImportError('C BGEN decoder is not built; run bash build_bgen_cpu.sh')
    try:
        lib = ct.CDLL(str(path))
    except OSError as exc:
        raise ImportError(f'C BGEN decoder could not load: {exc}') from exc
    try:
        lib.tg_bgen_cpu_abi_version.restype = ct.c_int
        if lib.tg_bgen_cpu_abi_version() != _ABI:
            raise ImportError('C BGEN decoder ABI mismatch; rebuild with bash build_bgen_cpu.sh')
    except AttributeError as exc:
        raise ImportError('C BGEN decoder is stale; rebuild with bash build_bgen_cpu.sh') from exc
    lib.tg_bgen_cpu_message.argtypes = [ct.c_int]
    lib.tg_bgen_cpu_message.restype = ct.c_char_p
    lib.tg_bgen_cpu_decode.argtypes = [
        ct.c_void_p, ct.c_size_t, ct.c_int, ct.c_uint32, ct.c_uint32,
        ct.c_void_p, ct.c_uint32, ct.c_void_p, ct.c_size_t, ct.c_void_p,
        ct.c_int64,
    ]
    lib.tg_bgen_cpu_decode.restype = ct.c_int
    lib.tg_bgen_cpu_transpose.argtypes = [
        ct.c_void_p, ct.c_int64, ct.c_int64, ct.c_void_p,
    ]
    lib.tg_bgen_cpu_transpose.restype = None
    return lib


def transpose(source):
    """Transpose a 2-D float32 array, without holding the GIL.

    numpy's own transpose-and-copy is fine on its own but holds the GIL for the
    whole move, which serialises every reader thread behind it; this is the
    same blocked copy done where ctypes can release it.
    """
    if source.dtype != np.float32 or source.ndim != 2:
        raise TypeError('transpose expects a 2-D float32 array')
    source = np.ascontiguousarray(source)
    rows, columns = source.shape
    destination = np.empty((columns, rows), dtype=np.float32)
    load_cpu_decoder().tg_bgen_cpu_transpose(
        source.ctypes.data, rows, columns, destination.ctypes.data)
    return destination


def cpu_decoder_available(compression):
    """True when the C decoder is built and handles this file's codec.

    zstd is not wired through it, matching the GPU decoder; those files keep
    using the numpy path rather than silently decoding differently.
    """
    if os.environ.get('TORCHGWAS_BGEN_C_DECODE', '1') == '0':
        return False
    if compression not in (0, 1):
        return False
    try:
        load_cpu_decoder()
    except ImportError:
        return False
    return True


class CpuBgenDecoder:
    """Decodes one record at a time into a caller-owned column.

    The inflate scratch buffer is per-decoder, so each reader thread must hold
    its own; `read_chunk` creates one per call, which is once per chunk rather
    than once per variant.

    **Bind a decoder to a thread, not to a task index.** Sharing one between two
    concurrent decodes corrupts the scratch, and it does not announce itself as
    a race: it surfaces as `BgenDecodeError: zlib inflate failed`, which names
    the codec rather than the aliasing, and it appears only above one worker. A
    pool handing out decoders by `index % workers` looks right and is not --
    nothing stops task i and task i+workers from overlapping. The *source* is
    safe to share, by contrast: `BgenGenotype._records` opens its own descriptor
    and reads with `os.pread`, which carries its own offset.
    """

    def __init__(self, n_bgen_samples, sample_indices, compression):
        self.lib = load_cpu_decoder()
        self.n_bgen_samples = int(n_bgen_samples)
        self.compression = int(compression)
        self.indices = np.ascontiguousarray(sample_indices, dtype=np.uint32)
        self.scratch = np.empty(0, dtype=np.uint8)

    def _reserve(self, raw_length):
        if self.scratch.size < raw_length:
            self.scratch = np.empty(int(raw_length), dtype=np.uint8)

    def decode_into(self, payload, raw_length, out):
        """Write one variant's dosages into `out`, a float32 view of a column.

        `out` may be strided -- it is typically a column of a (samples,
        variants) chunk -- so the stride travels to C rather than forcing a
        contiguous temporary here.
        """
        if out.dtype != np.float32:
            raise TypeError('C BGEN decode writes float32')
        if out.shape != (len(self.indices),):
            raise ValueError('output length does not match the selected samples')
        stride = out.strides[0]
        if stride % out.itemsize:
            raise ValueError('output stride is not a whole number of floats')
        self._reserve(raw_length)
        # np.frombuffer takes a read-only memoryview without copying, so the
        # record reaches C as a pointer into the buffer pread already returned.
        source = np.frombuffer(payload, dtype=np.uint8)
        status = self.lib.tg_bgen_cpu_decode(
            source.ctypes.data, source.size, self.compression,
            int(raw_length), self.n_bgen_samples,
            self.indices.ctypes.data, len(self.indices),
            self.scratch.ctypes.data, self.scratch.size,
            out.ctypes.data, stride // out.itemsize)
        if status != 0:
            raise BgenDecodeError(
                self.lib.tg_bgen_cpu_message(status).decode('utf-8', errors='replace'))
