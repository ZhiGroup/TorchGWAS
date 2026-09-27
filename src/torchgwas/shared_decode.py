"""Decode each genotype chunk once and feed it to several concurrent GPU scans.

Phenotype tiles scan the same variants on different GPUs. Without sharing,
every tile opens its own reader and decodes the whole genotype again, so N
tiles cost N decode passes of host CPU; on lab-2080ti that made 4 tiles 50%
slower than 2. The hub runs one PinnedDosageLoader and hands every decoded
pinned chunk to each subscriber (one per tile scan). Each GPU copies the chunk
from pinned host memory over its own PCIe link; a buffer returns to the
decoder only after every live subscriber has released it (after its H2D).

A subscriber has the PinnedDosageLoader consumer interface (iterate, release,
close), so dosage_cuda_iterator uses it unchanged in place of its own loader.
Every subscriber sees the same chunk sequence, including sizes chosen by an
adaptive selector, so all tiles stay aligned. The slowest tile paces the
decoder; the ring bounds how far the fastest can run ahead.
"""
from __future__ import annotations

import queue
import threading

from .streaming import PinnedDosageLoader

_DONE = object()


def root_slots(host_slots):
    """Chunk buffers the fan-out ring keeps on the root GPU."""
    return max(2, min(int(host_slots), 4))


class SharedDecodeHub:
    """fanout_device: copy each chunk once over PCIe into a ring on this GPU;
    subscribers then pull it GPU-to-GPU (NVLink peer copies where present)
    instead of each copying from pinned host memory over its own PCIe link."""

    def __init__(self, source, chunk_size, depth, reader_workers, *, subscribers,
                 variant_range=None, chunk_size_selector=None, record_timing=False,
                 fanout_device=None):
        if type(subscribers) is not int or subscribers < 1:
            raise ValueError('At least one subscriber required')
        self.capacity = int(chunk_size)
        self.subscribers = subscribers
        self._loader = PinnedDosageLoader(source, chunk_size, depth, reader_workers,
                                          variant_range=variant_range,
                                          chunk_size_selector=chunk_size_selector,
                                          record_timing=record_timing)
        self.fanout_device = None if fanout_device is None else str(fanout_device)
        if self.fanout_device is not None:
            import torch
            shape, dtype = self._loader.slot_shape, self._loader.slot_dtype
            with torch.cuda.device(torch.device(self.fanout_device)):
                self._fanout_stream = torch.cuda.Stream(device=self.fanout_device)
                with torch.cuda.stream(self._fanout_stream):
                    # A root slot is held only from its PCIe copy until every
                    # peer copy lands; decode lookahead lives in the host ring.
                    self._root = [torch.empty(shape, dtype=dtype, device=self.fanout_device)
                                  for _ in range(root_slots(self._loader.depth))]
                self._root_ready = [torch.cuda.Event() for _ in self._root]
            self._free_root = queue.Queue()
            for slot in range(len(self._root)):
                self._free_root.put(slot)
            self._peer_streams = {}
        self._lock = threading.Lock()
        self._queues = [queue.Queue() for _ in range(subscribers)]
        self._refs = {}
        self._outstanding = [set() for _ in range(subscribers)]
        self._closed = [False]*subscribers
        self._handed = [False]*subscribers
        self._loader_closed = False
        self._stopping = threading.Event()
        self.chunks = 0
        self._thread = threading.Thread(target=self._dispatch, daemon=True, name='torchgwas-shared-decode')
        self._thread.start()

    def _to_root(self, item):
        """One PCIe copy into the root GPU ring; returns the fan-out item."""
        import torch
        index, host, start, end = item[:4]
        while True:
            try:
                slot = self._free_root.get(timeout=0.1)
                break
            except queue.Empty:
                if self._stopping.is_set():
                    return None
        count = end-start
        with torch.cuda.device(torch.device(self.fanout_device)), torch.cuda.stream(self._fanout_stream):
            self._root[slot][:count].copy_(host, non_blocking=True)
            self._root_ready[slot].record(self._fanout_stream)
        self._root_ready[slot].synchronize()  # the pinned host buffer is free again
        self._loader.release(index)
        timing = item[4] if len(item) > 4 else None
        return (slot, self._root[slot][:count], start, end, timing, self._root_ready[slot])

    def peer_stream(self, device):
        """A dedicated stream on the root GPU for one destination's peer copies.

        PyTorch runs a cross-device copy on the source GPU's current stream;
        without this it would queue behind the root tile's own compute."""
        import torch
        with self._lock:
            stream = self._peer_streams.get(str(device))
            if stream is None:
                stream = torch.cuda.Stream(device=self.fanout_device)
                self._peer_streams[str(device)] = stream
            return stream

    def _dispatch(self):
        end = _DONE
        try:
            for item in self._loader:
                if self.fanout_device is not None:
                    item = self._to_root(item)
                    if item is None:
                        break
                index = item[0]
                with self._lock:
                    live = [i for i in range(self.subscribers) if not self._closed[i]]
                    if not live:
                        self._give_back(index)
                        continue
                    self._refs[index] = len(live)
                    self.chunks += 1
                    for i in live:
                        self._outstanding[i].add(index)
                        self._queues[i].put(item)
        except BaseException as error:  # propagate decode errors to every scan
            end = error
        for q in self._queues:
            q.put(end)

    def _give_back(self, index):
        """A buffer every live subscriber is done with: root slot or host slot."""
        if self.fanout_device is not None:
            self._free_root.put(index)
        elif not self._loader_closed:
            self._loader.release(index)

    def abandon_untaken(self):
        """Close subscribers no scan ever took, so the decoder can shut down."""
        for index in range(self.subscribers):
            if not self._handed[index]:
                self._handed[index] = True
                self._close(index)

    def audit(self):
        return dict(subscribers=self.subscribers, capacity=self.capacity, chunks=self.chunks,
                    decode_workers=getattr(self._loader, 'decode_workers_effective', None),
                    transfer=('pcie_to_root_then_gpu_peer' if self.fanout_device else 'pcie_per_gpu'),
                    fanout_device=self.fanout_device)

    def subscriber(self, index):
        if not 0 <= index < self.subscribers or self._handed[index]:
            raise ValueError('Each hub subscriber can be taken once')
        self._handed[index] = True
        return _Subscriber(self, index)

    def _release(self, subscriber, index):
        with self._lock:
            if index not in self._outstanding[subscriber]:
                return  # already released (e.g. while closing)
            self._outstanding[subscriber].discard(index)
            self._refs[index] -= 1
            if self._refs[index] == 0:
                del self._refs[index]
                self._give_back(index)

    def _close(self, subscriber):
        with self._lock:
            if self._closed[subscriber]:
                return
            self._closed[subscriber] = True
            held = list(self._outstanding[subscriber])
        # A scan that stops early (error, cancellation) must not stall the rest.
        for index in held:
            self._release(subscriber, index)
        while True:
            try:
                item = self._queues[subscriber].get_nowait()
            except queue.Empty:
                break
            if isinstance(item, tuple):
                self._release(subscriber, item[0])
        with self._lock:
            finished = all(self._closed)
            if finished:
                self._loader_closed = True
        if finished:
            self._stopping.set()
            self._loader.close()
            self._loader.ready.put(None)  # unblock the dispatcher if it waits
            self._thread.join(timeout=10)


class _Subscriber:
    """PinnedDosageLoader-compatible view of one scan's share of the hub."""

    def __init__(self, hub, index):
        self.hub, self.index = hub, index
        self.capacity = hub.capacity
        self.fanout = hub.fanout_device is not None

    def peer_stream(self, device):
        return self.hub.peer_stream(device)

    def __iter__(self):
        while True:
            item = self.hub._queues[self.index].get()
            if item is _DONE:
                return
            if isinstance(item, BaseException):
                raise item
            yield item

    def release(self, index):
        self.hub._release(self.index, index)

    def close(self):
        self.hub._close(self.index)


def _smi_gpus(devices):
    """(nvidia-smi index, PCI bus id) per CUDA device, matched by UUID, or None.

    Matching by UUID handles CUDA_VISIBLE_DEVICES renumbering."""
    import subprocess
    import torch
    try:
        out = subprocess.check_output(['nvidia-smi', '--query-gpu=index,uuid,pci.bus_id', '--format=csv,noheader'],
                                      text=True, timeout=20)
        table = {}
        for line in out.splitlines():
            index, uuid, bus = (x.strip() for x in line.split(','))
            table[uuid.lower().removeprefix('gpu-')] = (int(index), bus)
        return [table[str(torch.cuda.get_device_properties(torch.device(d).index).uuid).lower().removeprefix('gpu-')]
                for d in devices]
    except (OSError, subprocess.SubprocessError, ValueError, KeyError, AttributeError):
        return None  # AttributeError: torch without .uuid


def nvlink_root(devices):
    """devices[0] when NVLink connects it to every other device, else None.

    Reads `nvidia-smi topo -m` (NV# entries) and also requires CUDA peer
    access. PCIe peer paths (PIX/PXB) do not qualify: they share the same
    links a per-GPU host copy would use.
    """
    import re
    import subprocess
    import torch
    if len(devices) < 2:
        return None
    gpus = _smi_gpus(devices)
    try:
        topo = subprocess.check_output(['nvidia-smi', 'topo', '-m'], text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return None
    rows = {}
    for line in re.sub(r'\x1b\[[0-9;]*m', '', topo).splitlines():
        tokens = line.split()
        if tokens and re.fullmatch(r'GPU\d+', tokens[0]):
            rows[int(tokens[0][3:])] = tokens[1:]
    try:
        physical = [index for index, _ in gpus]
        linked = all(rows[physical[0]][p].startswith('NV') for p in physical[1:])
    except (TypeError, KeyError, IndexError):
        return None
    cuda = [torch.device(d).index for d in devices]
    peers = all(torch.cuda.can_device_access_peer(i, cuda[0]) for i in cuda[1:])
    return str(devices[0]) if linked and peers else None


def pcie_root_port(bus_id, sysfs='/sys/bus/pci/devices'):
    """'pciDDDD:BB/<root port>' above a PCI device, from the sysfs tree, or None.

    Devices below the same root port share its link to the CPU: lab-a100
    GPUs 0 and 1 hang off one switch under root port 0000:00:01.1, while every
    lab-h100 GPU has its own."""
    import os
    try:
        domain, rest = bus_id.split(':', 1)
        name = f'{int(domain, 16):04x}:{rest.lower()}'
    except ValueError:
        return None
    parts = os.path.realpath(os.path.join(sysfs, name)).split('/')
    if 'devices' not in parts:
        return None
    tail = parts[parts.index('devices')+1:]
    if len(tail) < 3 or not tail[0].startswith('pci'):
        return None  # not found, or directly on the host bridge
    return '/'.join(tail[:2])


def shared_uplink(devices):
    """True when two of the devices sit below one PCIe root port."""
    gpus = _smi_gpus(devices)
    if gpus is None:
        return False
    ports = [pcie_root_port(bus) for _, bus in gpus]
    known = [p for p in ports if p is not None]
    return len(known) != len(set(known))


def fanout_root(devices):
    """gpu_fanout='uplink': devices[0] when NVLink connects the devices and at
    least two of them share a PCIe uplink, so per-GPU host copies split its
    bandwidth (lab-a100 probe: 3.3 GB/s each versus 6.7 alone). Opt-in: the
    end-to-end runs so far were not transfer-bound, and fan-out did not pay."""
    if len(devices) < 2 or not shared_uplink(devices):
        return None
    return nvlink_root(devices)


def shared_decode_eligible(genotype, *, devices, tiles, compute_dtype):
    """True when every tile can share one host decode pass.

    Requires the native dosage pipeline with host decoding (PGEN, zstd store,
    BGEN on the CPU decoder), CUDA float32, more than one tile and no more
    tiles than GPUs, so all tiles scan concurrently.
    """
    import torch
    if tiles < 2 or tiles > len(devices) or compute_dtype != 'float32':
        return False
    if not getattr(genotype, 'supports_fused_qc', False) or hasattr(genotype, 'iter_packed_chunks'):
        return False
    if hasattr(genotype, 'iter_device_chunks'):
        backend = (genotype.resolve_decode_backend(torch.device(devices[0]))
                   if hasattr(genotype, 'resolve_decode_backend') else 'gpu')
        if backend != 'cpu':
            return False
    return all(torch.device(d).type == 'cuda' for d in devices)
