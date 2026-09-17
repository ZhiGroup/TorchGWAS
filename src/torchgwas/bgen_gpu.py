"""Optional nvCOMP decoder; only compact status arrays return to the CPU.

Each output is a fresh caller-owned torch tensor. Native decode synchronizes its
private stream before returning, so outputs are ready on every consumer stream.
The Python executor overlaps CPU reads and the next decode with previous GWAS
compute. No output or host input buffer is recycled while in flight.
"""
from __future__ import annotations
import ctypes as ct
import os
import time
from functools import lru_cache
from pathlib import Path
import numpy as np


@lru_cache(maxsize=1)
def load_decoder_library():
    # TORCHGWAS_BGEN_LIBRARY points at an alternative build so two builds
    # can be A/B'd interleaved in one run (as TORCHGWAS_PGEN_LIBRARY does
    # for the PGEN reader). Unset, the shipped library next to this file.
    override = os.environ.get('TORCHGWAS_BGEN_LIBRARY')
    path = Path(override) if override else Path(__file__).with_name('native') / 'libtorchgwas_bgen.so'
    if not path.is_file():
        raise ImportError('GPU BGEN decoder is not built; run bash build_direct_bgen.sh')
    try:
        lib = ct.CDLL(str(path))
    except OSError as exc:
        raise ImportError(f'GPU BGEN decoder could not load: {exc}') from exc
    try:
        version = lib.tg_bgen_abi_version
        version.restype = ct.c_int
        if version() != 4:
            raise ImportError("GPU BGEN decoder ABI mismatch; rebuild with bash build_direct_bgen.sh")
    except AttributeError as exc:
        raise ImportError("GPU BGEN decoder is stale; rebuild with bash build_direct_bgen.sh") from exc
    p = ct.c_void_p
    lib.tg_bgen_probe.argtypes = [ct.c_int]
    lib.tg_bgen_probe.restype = ct.c_int
    lib.tg_bgen_error.restype = ct.c_char_p
    lib.tg_bgen_create.argtypes = [ct.c_int, ct.c_size_t, p, ct.c_size_t]
    lib.tg_bgen_create.restype = p
    lib.tg_bgen_destroy.argtypes = [p]
    lib.tg_bgen_set_verify_adler.argtypes = [p, ct.c_int]
    lib.tg_bgen_read.argtypes = [ct.c_int, p, p, ct.c_size_t, ct.c_char_p, ct.c_size_t, p]
    lib.tg_bgen_read.restype = p
    lib.tg_bgen_free_batch.argtypes = [p]
    lib.tg_bgen_decode.argtypes = [p, p, p, p, p]
    lib.tg_bgen_decode.restype = ct.c_int
    if hasattr(lib,"tg_bgen_profile"):
        lib.tg_bgen_profile.argtypes=[p,p]
        lib.tg_bgen_profile.restype=None
    return lib


_LAUNCHABLE: dict[int, bool] = {}


def decoder_launchable(device) -> bool:
    """True when this library's kernels can actually launch on `device`.

    Loading is not the question. The library carries code only for the
    architectures it was built for, so on any other card it loads cleanly and
    then fails at *launch* -- `no kernel image is available for execution on
    the device` -- part way through a scan, where the caller has already
    committed to the GPU backend. `tg_bgen_probe` launches a kernel that does
    nothing so the question can be asked in advance and answered cheaply.

    Cached per device: the answer depends on the architecture and the library,
    and neither changes within a run.
    """
    import torch

    try:
        index = torch.device(device).index
        if index is None:
            index = torch.cuda.current_device()
    except Exception:  # noqa: BLE001 - a device we cannot name is not usable
        return False
    cached = _LAUNCHABLE.get(index)
    if cached is not None:
        return cached
    try:
        ok = bool(load_decoder_library().tg_bgen_probe(int(index)))
    except Exception:  # noqa: BLE001 - unbuilt, stale, or unlaunchable
        ok = False
    _LAUNCHABLE[index] = ok
    return ok


class GpuBgenDecoder:
    def __init__(self, n_samples, sample_indices, device):
        import torch
        self.lib = load_decoder_library()
        self.profile_enabled = os.environ.get("TORCHGWAS_BGEN_PROFILE", "0") != "0"
        self.host_timings = []
        self.device = torch.device(device)
        if self.device.type != 'cuda':
            raise ValueError('GPU BGEN decoding requires a CUDA device')
        self.device = torch.device('cuda', torch.cuda.current_device() if self.device.index is None else self.device.index)
        self.allocation_stream = torch.cuda.Stream(device=self.device)
        self.indices = np.ascontiguousarray(sample_indices, dtype=np.int32)
        self.handle = self.lib.tg_bgen_create(self.device.index, n_samples, self.indices.ctypes.data, len(self.indices))
        # Skipping checksum verification trades the only end-to-end integrity
        # check on this path for about an eighth of decode time. Opt in
        # explicitly; never inferred.
        self.verify_adler = os.environ.get("TORCHGWAS_BGEN_SKIP_ADLER", "0") != "1"
        if self.handle and not self.verify_adler:
            self.lib.tg_bgen_set_verify_adler(self.handle, 0)
        if not self.handle:
            self._error()

    def _error(self):
        raise RuntimeError(self.lib.tg_bgen_error().decode('utf-8', errors='replace'))

    def read(self, fd, offsets, lengths, metadata, positions):
        offsets = np.ascontiguousarray(offsets, dtype=np.uint64)
        lengths = np.ascontiguousarray(lengths, dtype=np.uint64)
        positions = np.ascontiguousarray(positions, dtype=np.int64)
        started = time.perf_counter()
        cpu_started = time.thread_time()
        handle = self.lib.tg_bgen_read(fd, offsets.ctypes.data, lengths.ctypes.data, len(offsets), metadata, len(metadata), positions.ctypes.data)
        if self.profile_enabled:
            self.host_timings.append(("read_native", time.perf_counter()-started, time.thread_time()-cpu_started))
        if not handle:
            self._error()
        return handle

    def decode(self, batch, count):
        import torch
        # A dedicated allocation stream permits overlap with GWAS compute. Native
        # decode waits for allocator dependencies and completes before return.
        try:
            with torch.cuda.device(self.device), torch.cuda.stream(self.allocation_stream):
                allocated_at = time.perf_counter()
                allocated_cpu = time.thread_time()
                out = torch.empty((count, len(self.indices)), dtype=torch.float32, device=self.device)
                statuses = np.empty(count, dtype=np.int32)
                native_at = time.perf_counter()
                native_cpu = time.thread_time()
                if self.profile_enabled:
                    self.host_timings.append(("output_allocation", native_at-allocated_at, native_cpu-allocated_cpu))
                result = self.lib.tg_bgen_decode(self.handle, batch, out.data_ptr(), statuses.ctypes.data, torch.cuda.current_stream(self.device).cuda_stream)
                if self.profile_enabled:
                    self.host_timings.append(("decode_native", time.perf_counter()-native_at, time.thread_time()-native_cpu))
                if result:
                    self._error()
                bad = np.flatnonzero(statuses)
                if bad.size:
                    names = {1:'DEFLATE error',2:'inflated size mismatch',3:'sample/allele count mismatch',4:'non-diploid',6:'phased',7:'invalid probability bit width',8:'probability payload length mismatch',9:'invalid probabilities',10:'Adler32 checksum mismatch'}
                    i = int(bad[0])
                    raise ValueError(f'BGEN variant {i} in batch: {names.get(int(statuses[i]), str(statuses[i]))}; no variants were silently removed')
                return out
        finally:
            # Even an output allocation failure must release the prepared batch.
            self.lib.tg_bgen_free_batch(batch)

    def profile(self):
        if not hasattr(self.lib,'tg_bgen_profile'):
            return {'available':False,'enabled':False}
        values=np.empty(16,dtype=np.float64)
        self.lib.tg_bgen_profile(self.handle,values.ctypes.data)
        keys=('enabled','batches','compressed_bytes','inflated_bytes','dosage_bytes',
              'h2d_and_descriptors_ms','inflate_ms','probability_decode_ms','adler32_ms','status_d2h_ms',
              # Device capacities, high-water marks. `nvcomp_temp_capacity_bytes`
              # is the term the memory model cannot derive: nvCOMP sizes its own
              # scratch at runtime and it varies by version and GPU, so a shipped
              # calculator has to measure it on the host it lands on.
              'nvcomp_temp_capacity_bytes','compressed_capacity_bytes',
              'raw_capacity_bytes','nvcomp_temp_last_bytes','records_last')
        result=dict(zip(keys,values.tolist()))
        result['available']=True
        result['completion_wait']='blocking_event' if os.environ.get('TORCHGWAS_BGEN_BLOCKING_SYNC')=='1' else 'stream'
        result['enabled']=bool(result['enabled'])
        result['host_stages'] = {}
        for stage in {row[0] for row in self.host_timings}:
            rows=[row for row in self.host_timings if row[0]==stage]
            result['host_stages'][stage]={'calls':len(rows), 'sum_wall_seconds':sum(row[1] for row in rows), 'sum_thread_cpu_seconds':sum(row[2] for row in rows), 'max_wall_seconds':max(row[1] for row in rows)}
        return result

    def close(self):
        if getattr(self, 'handle', None):
            self.lib.tg_bgen_destroy(self.handle)
            self.handle = None
