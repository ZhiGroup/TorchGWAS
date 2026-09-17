"""GIL-free single-frame zstd decode and bounded contiguous compressed reads."""
from __future__ import annotations
import ctypes as ct
import ctypes.util
import os
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from threading import local


class ZstdInto:
    """One native decompression context per calling thread.

    CDLL releases the GIL during native decompression. Input/output exporters
    remain strongly referenced until the synchronous native call finishes.
    Output must be writable and C-contiguous (including pinned NumPy views).
    """
    def __init__(self):
        name = ctypes.util.find_library('zstd')
        if not name:
            raise ImportError('libzstd is required for direct-to-buffer decompression')
        self.lib = ct.CDLL(name)
        self.lib.ZSTD_createDCtx.restype = ct.c_void_p
        self.lib.ZSTD_freeDCtx.argtypes = [ct.c_void_p]
        self.lib.ZSTD_freeDCtx.restype = ct.c_size_t
        self.lib.ZSTD_decompressDCtx.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_size_t, ct.c_void_p, ct.c_size_t]
        self.lib.ZSTD_decompressDCtx.restype = ct.c_size_t
        self.lib.ZSTD_findFrameCompressedSize.argtypes = [ct.c_void_p, ct.c_size_t]
        self.lib.ZSTD_findFrameCompressedSize.restype = ct.c_size_t
        self.lib.ZSTD_isError.argtypes = [ct.c_size_t]
        self.lib.ZSTD_isError.restype = ct.c_uint
        self.lib.ZSTD_getErrorName.argtypes = [ct.c_size_t]
        self.lib.ZSTD_getErrorName.restype = ct.c_char_p
        self.state = local()

    def _check(self, result):
        if self.lib.ZSTD_isError(result):
            raise ValueError('zstd: ' + self.lib.ZSTD_getErrorName(result).decode())
        return result

    def decompress_into(self, source, destination):
        source_view = memoryview(source).cast('B')
        output_view = memoryview(destination).cast('B')
        if output_view.readonly:
            raise ValueError('zstd destination must be writable')
        if not source_view.nbytes or not output_view.nbytes:
            raise ValueError('empty compressed input or output is not a genotype frame')
        # read-ahead supplies writable bytearrays, allowing zero-copy input.
        source_buffer = ((ct.c_ubyte * source_view.nbytes).from_buffer_copy(source_view)
                         if source_view.readonly else
                         (ct.c_ubyte * source_view.nbytes).from_buffer(source_view))
        output_buffer = (ct.c_ubyte * output_view.nbytes).from_buffer(output_view)
        frame_size = self._check(self.lib.ZSTD_findFrameCompressedSize(source_buffer, source_view.nbytes))
        if frame_size != source_view.nbytes:
            raise ValueError('zstd input must contain exactly one complete frame')
        if not hasattr(self.state, 'context'):
            self.state.context = _Context(self.lib)
        size = self._check(self.lib.ZSTD_decompressDCtx(self.state.context.ptr,
                          output_buffer, output_view.nbytes, source_buffer, source_view.nbytes))
        if size != output_view.nbytes:
            raise ValueError(f'zstd decoded {size} bytes; expected {output_view.nbytes}')
        return size


class _Context:
    def __init__(self, lib):
        self.lib = lib
        self.ptr = lib.ZSTD_createDCtx()
        if not self.ptr:
            raise MemoryError('ZSTD_createDCtx failed')

    def __del__(self):
        if getattr(self, 'ptr', None):
            self.lib.ZSTD_freeDCtx(self.ptr)


def contiguous_frame_batches(path, offsets, sizes, *, target_bytes=16 << 20, prefetch_batches=2, read_workers=1):
    """Yield ordered lists of (frame index, writable compressed memoryview).

    Independent I/O workers read large contiguous batches; submission and
    delivery remain ordered even if reads complete out of order.
    At most prefetch_batches read futures plus the yielded batch are retained
    here. Downstream consumers must independently bound decoder futures.
    A frame larger than target_bytes is read alone; the bound is therefore
    max(target_bytes, largest frame), not an unconditional byte ceiling.
    """
    if target_bytes <= 0 or prefetch_batches <= 0 or read_workers <= 0:
        raise ValueError('read batch size, depth and workers must be positive')
    if read_workers > 1 and not hasattr(os, 'pread'):
        raise ValueError('parallel batch reads require positional I/O support')
    if len(offsets) != len(sizes):
        raise ValueError('frame offset/size counts differ')
    previous_end = 0
    for offset, size in zip(offsets, sizes):
        if int(offset) < previous_end or int(size) <= 0:
            raise ValueError('frame index must be ordered, disjoint and nonempty')
        previous_end = int(offset) + int(size)
    fd = os.open(path, os.O_RDONLY | getattr(os, 'O_BINARY', 0))
    pending = deque()

    def read_group(first, end):
        start = int(offsets[first])
        length = int(offsets[end - 1]) + int(sizes[end - 1]) - start
        data = bytearray(length)
        view = memoryview(data)
        position = 0
        while position < length:
            if hasattr(os, 'preadv'):
                count = os.preadv(fd, [view[position:]], start + position)
            else:
                if hasattr(os, "pread"):
                    block = os.pread(fd, length - position, start + position)
                else:
                    os.lseek(fd, start + position, os.SEEK_SET)
                    block = os.read(fd, length - position)
                count = len(block)
                view[position:position + count] = block
            if not count:
                raise OSError(f'short zstd read at byte {start + position}')
            position += count
        return [(i, view[int(offsets[i]) - start:int(offsets[i]) - start + int(sizes[i])])
                for i in range(first, end)]

    try:
        with ThreadPoolExecutor(max_workers=read_workers, thread_name_prefix='zstd-read-ahead') as pool:
            first = 0
            while first < len(offsets) or pending:
                while first < len(offsets) and len(pending) < prefetch_batches:
                    end = first + 1
                    while end < len(offsets):
                        if int(offsets[end]) != int(offsets[end - 1]) + int(sizes[end - 1]):
                            break
                        if int(offsets[end]) + int(sizes[end]) - int(offsets[first]) > target_bytes:
                            break
                        end += 1
                    pending.append(pool.submit(read_group, first, end))
                    first = end
                yield pending.popleft().result()
    finally:
        os.close(fd)
