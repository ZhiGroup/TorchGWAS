"""Independent pinned H2D/D2H service for declared GPU groups.

Each round measures every single GPU, every declared shared link and the whole
active set, in a seeded random order, from pinned host buffers whose NUMA page
residency is checked. Rounds alternate between calibration and held-out halves.
A price is published only when every held-out median falls inside the
(slightly widened) calibration range and near its median, the selected GPUs
were idle before any of our traffic and carried no foreign process in any
round, and a group never outran the sum of its members. Activity on other GPUs,
which may share a PCIe switch or root complex, is recorded but not excluded. The low/median/high scenarios are
observed spans of repeat medians, not confidence intervals or hardware ceilings.
Use `high` for necessary resource floors and `low` for conservative ceilings.
No link topology is inferred: shared links are declared by the caller.
"""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import ctypes
import os
import random
import statistics
import subprocess
import threading
import time

from .gpu_identity import canonical_cuda_device

PROTOCOL = 'torchgwas.pinned_transfer_groups.v1'
KIND = 'transfer_capacity'
NAME = 'pinned_transfer_groups.v1'
SCENARIOS = ('low', 'median', 'high')
DIRECTIONS = ('h2d', 'd2h')


def _index(device):
    return int(canonical_cuda_device(device)[5:])


def declared_groups(devices, links=()):
    """Singles, declared links and the whole set, each exactly once, in order."""
    devices = [canonical_cuda_device(device) for device in devices]
    if not devices or len(set(devices)) != len(devices):
        raise ValueError('Unique explicit transfer devices required')
    groups = [(device,) for device in devices]
    for link in links:
        link = tuple(canonical_cuda_device(device) for device in link)
        if (len(link) < 2 or len(set(link)) != len(link) or not set(link) <= set(devices)):
            raise ValueError('Declared links need two or more unique active GPUs')
        if any(set(link) == set(group) for group in groups):
            raise ValueError('Duplicate declared transfer link')
        groups.append(link)
    if len(devices) > 1 and not any(set(group) == set(devices) for group in groups):
        groups.append(tuple(devices))
    return groups


def pinned_page_nodes(address, nbytes, *, max_pages=256):
    """NUMA node of sampled resident pages via move_pages(2) query mode.

    Passing no target nodes only reports placement; nothing is migrated.
    Negative statuses (for example -ENOENT for an unpopulated page) are kept.
    """
    if type(address) is not int or address <= 0 or type(nbytes) is not int or nbytes <= 0:
        raise ValueError('Positive host address and size required')
    page = os.sysconf('SC_PAGE_SIZE')
    first = address - address % page
    count = (address + nbytes - first + page - 1) // page
    step = max(1, count // max_pages)
    pages = [first + i*page for i in range(0, count, step)]
    library = ctypes.CDLL('libnuma.so.1', use_errno=True)
    move = library.move_pages
    move.argtypes = [ctypes.c_int, ctypes.c_ulong, ctypes.POINTER(ctypes.c_void_p),
                     ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    move.restype = ctypes.c_long
    request = (ctypes.c_void_p*len(pages))(*pages)
    status = (ctypes.c_int*len(pages))()
    if move(0, len(pages), request, None, status, 0) != 0:
        raise OSError(ctypes.get_errno(), 'move_pages query failed')
    counts = {}
    for value in status:
        counts[str(value)] = counts.get(str(value), 0) + 1
    return dict(sampled_pages=len(pages), total_pages=count, status_counts=counts)


def gpu_load(uuids):
    """Utilization, memory and foreign compute processes on every physical GPU.

    nvidia-smi utilization averages a recent window, so after our own copies
    it reports our traffic; only foreign processes identify competing work.
    Unselected GPUs are neighbour context: they may share a PCIe switch.
    """
    rows = subprocess.check_output(['nvidia-smi', '--query-gpu=uuid,utilization.gpu,memory.used',
                                    '--format=csv,noheader,nounits'], text=True)
    apps = subprocess.check_output(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid',
                                    '--format=csv,noheader'], text=True)
    wanted = {uuid.lower().removeprefix('gpu-') for uuid in uuids}
    load = {}
    for line in rows.splitlines():
        uuid, util, memory = (field.strip() for field in line.split(','))
        load[uuid.lower().removeprefix('gpu-')] = dict(
            utilization_percent=int(util), memory_used_mib=int(memory), foreign_pids=[])
    for line in apps.splitlines():
        if not line.strip():
            continue
        uuid, pid = (field.strip() for field in line.split(','))
        uuid = uuid.lower().removeprefix('gpu-')
        if uuid in load and int(pid) != os.getpid():
            load[uuid]['foreign_pids'].append(int(pid))
    if not wanted <= set(load):
        raise ValueError('Selected GPU UUID missing from nvidia-smi')
    return dict(selected={uuid: load[uuid] for uuid in sorted(wanted)},
                neighbours={uuid: row for uuid, row in sorted(load.items()) if uuid not in wanted})


def host_load():
    with open('/proc/loadavg', encoding='ascii') as stream:
        one, five, fifteen = (float(x) for x in stream.read().split()[:3])
    with open('/proc/stat', encoding='ascii') as stream:
        fields = [int(x) for x in stream.readline().split()[1:]]
    return dict(loadavg=[one, five, fifteen], cpu_jiffies=fields, cpu_count=os.cpu_count())


def _busy_fraction(before, after):
    delta = [b-a for a, b in zip(before['cpu_jiffies'], after['cpu_jiffies'])]
    total = sum(delta)
    idle = delta[3] + (delta[4] if len(delta) > 4 else 0)
    return None if total <= 0 else 1. - idle/total


class _Buffers:
    """One pinned host and one device buffer per GPU, allocated by this thread."""

    def __init__(self, devices, nbytes):
        import torch
        self.torch = torch
        self.rows = {}
        for device in devices:
            with torch.cuda.device(_index(device)):
                host = torch.empty(nbytes, dtype=torch.uint8, pin_memory=True)
                host.fill_(_index(device) % 251)  # populate every page before placement query
                gpu = torch.empty(nbytes, dtype=torch.uint8, device=device)
                gpu.fill_(_index(device) % 251)
                stream = torch.cuda.Stream(device=device)
                torch.cuda.synchronize(device)
            self.rows[device] = (host, gpu, stream)

    def placement(self):
        return {device: pinned_page_nodes(host.data_ptr(), host.numel())
                for device, (host, _, _) in self.rows.items()}


def _observe(group, buffers, *, direction, copies, pool):
    torch = buffers.torch
    gate = threading.Barrier(len(group)+1)

    def copy(device):
        host, gpu, stream = buffers.rows[device]
        source, target = (host, gpu) if direction == 'h2d' else (gpu, host)
        with torch.cuda.device(_index(device)), torch.cuda.stream(stream):
            begin = torch.cuda.Event(enable_timing=True)
            end = torch.cuda.Event(enable_timing=True)
            gate.wait()
            wall_begin = time.perf_counter()
            begin.record(stream)
            for _ in range(copies):
                target.copy_(source, non_blocking=True)
            end.record(stream)
            end.synchronize()
            return device, wall_begin, time.perf_counter(), begin.elapsed_time(end)/1000.

    pending = [pool.submit(copy, device) for device in group]
    gate.wait()
    spans = [future.result() for future in pending]
    wall = max(row[2] for row in spans) - min(row[1] for row in spans)
    if wall <= 0:
        raise ValueError('Positive simultaneous copy span required')
    nbytes = buffers.rows[group[0]][0].numel()
    return dict(wall_seconds=wall, aggregate_bytes_per_second=nbytes*copies*len(group)/wall,
                event_seconds={row[0]: row[3] for row in spans},
                start_skew_seconds=max(row[1] for row in spans)-min(row[1] for row in spans))


def measure_transfer_groups(devices, links=(), *, size_bytes=32 << 20, copies=8, samples=3,
                            warmups=2, rounds=10, seed=0, idle_utilization_limit=10,
                            uuids=None):
    """Raw per-round observations; no summary, qualification or publication.

    Each round runs every (direction, group) cell `samples` times in a fresh
    seeded order and records GPU/host load before and after the round.
    """
    groups = declared_groups(devices, links)
    devices = [group[0] for group in groups if len(group) == 1]
    if (type(size_bytes) is not int or not 1 << 20 <= size_bytes <= 1 << 30 or
            type(copies) is not int or not 1 <= copies <= 64 or
            type(samples) is not int or not 1 <= samples <= 16 or
            type(warmups) is not int or not 0 <= warmups <= 16 or
            type(rounds) is not int or not 4 <= rounds <= 64 or rounds % 2 or
            type(idle_utilization_limit) is not int or not 0 <= idle_utilization_limit <= 100):
        raise ValueError('Bounded explicit transfer calibration arguments required')
    import torch
    if uuids is None:
        uuids = {device: str(torch.cuda.get_device_properties(device).uuid).lower().removeprefix('gpu-')
                 for device in devices}
    initial = gpu_load(uuids.values())  # before any of our own traffic
    initial_busy = sorted(uuid for uuid, row in initial['selected'].items()
                          if row['utilization_percent'] > idle_utilization_limit or row['foreign_pids'])
    buffers = _Buffers(devices, size_bytes)
    placement = buffers.placement()
    cells = [(direction, group) for direction in DIRECTIONS for group in groups]
    generator = random.Random(seed)
    records = []
    with ThreadPoolExecutor(max_workers=len(devices)) as pool:
        for direction, group in cells:
            for _ in range(warmups):
                _observe(group, buffers, direction=direction, copies=copies, pool=pool)
        for round_index in range(rounds):
            order = list(cells)
            generator.shuffle(order)
            gpu_before, host_before = gpu_load(uuids.values()), host_load()
            started = time.time()
            observations = []
            for direction, group in order:
                rows = [_observe(group, buffers, direction=direction, copies=copies, pool=pool)
                        for _ in range(samples)]
                observations.append(dict(direction=direction, devices=list(group), samples=rows,
                    bytes_per_second=statistics.median(row['aggregate_bytes_per_second']
                                                       for row in rows)))
            gpu_after, host_after = gpu_load(uuids.values()), host_load()
            busy = sorted({uuid for snapshot in (gpu_before, gpu_after)
                           for uuid, row in snapshot['selected'].items() if row['foreign_pids']})
            active_neighbours = sorted({uuid for snapshot in (gpu_before, gpu_after)
                                        for uuid, row in snapshot['neighbours'].items()
                                        if row['utilization_percent'] > idle_utilization_limit})
            records.append(dict(round=round_index, half='calibration' if round_index % 2 == 0
                                else 'holdout', started_unix_seconds=started,
                                order=[[d, list(g)] for d, g in order], observations=observations,
                                gpu_before=gpu_before, gpu_after=gpu_after,
                                host_before=host_before, host_after=host_after,
                                host_busy_fraction=_busy_fraction(host_before, host_after),
                                busy_selected_gpus=busy, active_neighbour_gpus=active_neighbours))
    placement_after = buffers.placement()
    return dict(devices=devices, groups=[list(group) for group in groups], uuids=uuids,
                initial_gpu_load=initial, initial_busy_selected_gpus=initial_busy,
                size_bytes=size_bytes, copies=copies, samples=samples, warmups=warmups,
                rounds=rounds, seed=seed, idle_utilization_limit=idle_utilization_limit,
                placement_before=placement, placement_after=placement_after, records=records)


def _cell_key(direction, group):
    return direction + ':' + ','.join(group)


def summarize_transfer_groups(raw, *, tolerance=0.10, expected_nodes=None):
    """Qualify calibration scenarios against the held-out half.

    `expected_nodes` optionally maps each device to the NUMA node its pinned
    buffer must occupy. Returns `qualified` with a transfer value, or the
    failing checks. Nothing is dropped: a busy or inconsistent round fails.
    """
    if (isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or
            not 0 < tolerance < 1):
        raise ValueError('Held-out tolerance must lie in (0, 1)')
    groups = [tuple(group) for group in raw['groups']]
    devices = list(raw['devices'])
    failures = []
    if raw.get('initial_busy_selected_gpus'):
        failures.append(dict(check='selected_gpu_idle_at_start',
                             uuids=raw['initial_busy_selected_gpus']))
    for record in raw['records']:
        if record['busy_selected_gpus']:
            failures.append(dict(check='selected_gpu_idle', round=record['round'],
                                 uuids=record['busy_selected_gpus']))
    for stage in ('placement_before', 'placement_after'):
        for device, placement in raw[stage].items():
            nodes = {key for key, count in placement['status_counts'].items() if count}
            if any(int(key) < 0 for key in nodes) or len(nodes) != 1:
                failures.append(dict(check='pinned_single_numa_node', stage=stage, device=device,
                                     status_counts=placement['status_counts']))
            elif expected_nodes is not None and nodes != {str(expected_nodes[device])}:
                failures.append(dict(check='pinned_expected_numa_node', stage=stage,
                                     device=device, observed=sorted(nodes),
                                     expected=expected_nodes[device]))
    halves = {}
    for record in raw['records']:
        for row in record['observations']:
            key = _cell_key(row['direction'], row['devices'])
            halves.setdefault(key, dict(calibration=[], holdout=[]))[record['half']].append(
                row['bytes_per_second'])
    cells = {}
    for direction in DIRECTIONS:
        for group in groups:
            key = _cell_key(direction, group)
            calibration, holdout = halves[key]['calibration'], halves[key]['holdout']
            low, middle, high = min(calibration), statistics.median(calibration), max(calibration)
            held = statistics.median(holdout)
            cell = dict(direction=direction, devices=list(group), low=low, median=middle,
                        high=high, holdout_median=held, holdout_ratio=held/middle,
                        calibration=calibration, holdout=holdout)
            cells[key] = cell
            # A tight calibration span must not reject a held-out median a
            # fraction of a percent outside it, so the span is widened by half
            # the tolerance; the median drift test carries the main check.
            margin = tolerance/2
            if not low*(1-margin) <= held <= high*(1+margin) or abs(held/middle-1) > tolerance:
                failures.append(dict(check='holdout', cell=key, low=low, median=middle, high=high,
                                     holdout_median=held))
            if len(group) > 1:
                ceiling = sum(cells[_cell_key(direction, (device,))]['high'] for device in group)
                if high > ceiling*(1+tolerance):
                    failures.append(dict(check='group_not_superadditive', cell=key,
                                         group_high=high, member_high_sum=ceiling))
    summary = dict(cells=cells, failures=failures, tolerance=tolerance,
                   qualified=not failures)
    if not failures:
        summary['value'] = transfer_value(cells, devices, groups, raw)
    return summary


def transfer_value(cells, devices, groups, raw):
    whole = tuple(devices)
    shared = [group for group in groups if len(group) > 1 and set(group) != set(whole)]
    scenarios = {}
    for name in SCENARIOS:
        def rate(direction, group):
            return cells[_cell_key(direction, group)][name]
        full = next(group for group in groups if set(group) == set(whole))
        scenarios[name] = dict(
            per_device={device: dict(h2d=rate('h2d', (device,)), d2h=rate('d2h', (device,)))
                        for device in devices},
            shared_transfer_capacities=dict(h2d=rate('h2d', full), d2h=rate('d2h', full)),
            shared_links=[dict(devices=list(group), h2d_bytes_per_second=rate('h2d', group),
                               d2h_bytes_per_second=rate('d2h', group)) for group in shared])
    return dict(protocol=PROTOCOL, devices=devices, uuids=raw['uuids'], scenarios=scenarios,
                pinned_placement=raw['placement_before'],
                scope='Pinned uint8 copies on explicit streams; observed repeat-median spans under recorded load, not guaranteed ceilings.')


def measurement_protocol(raw, *, tolerance, implementation_sha256):
    return dict(operation=PROTOCOL, size_bytes=raw['size_bytes'], copies=raw['copies'],
                samples=raw['samples'], warmups=raw['warmups'], rounds=raw['rounds'],
                seed=raw['seed'], groups=raw['groups'],
                idle_utilization_limit=raw['idle_utilization_limit'], tolerance=tolerance,
                split='alternate rounds: even calibration, odd held-out',
                aggregation='median of samples per round; low/median/high over calibration rounds',
                implementation_sha256=implementation_sha256)


def apply_transfer_prices(context, value, scenario):
    """Copy one scenario into a detailed-profile context (returns a new dict)."""
    if scenario not in SCENARIOS:
        raise ValueError('Unknown transfer scenario')
    chosen = value['scenarios'][scenario]
    context = deepcopy(context)
    if set(context.get('devices', ())) != set(value['devices']):
        raise ValueError('Transfer prices must cover exactly the context devices')
    for device, rates in chosen['per_device'].items():
        context['profiles'][device]['h2d_bytes_per_second'] = rates['h2d']
        context['profiles'][device]['d2h_bytes_per_second'] = rates['d2h']
    context['shared_transfer_capacities'] = deepcopy(chosen['shared_transfer_capacities'])
    context['shared_links'] = deepcopy(chosen['shared_links'])
    return context


def transfer_price_targets(context_index, value, scenario):
    """Binding targets that tie every applied transfer value to the record."""
    if type(context_index) is not int or context_index < 0 or scenario not in SCENARIOS:
        raise ValueError('Explicit context index and scenario required')
    base = ['scenarios', scenario]
    targets = []
    for device in value['devices']:
        for direction in DIRECTIONS:
            targets.append(dict(context_path=[context_index, 'profiles', device,
                                              direction+'_bytes_per_second'],
                                value_path=base+['per_device', device, direction]))
    for direction in DIRECTIONS:
        targets.append(dict(context_path=[context_index, 'shared_transfer_capacities', direction],
                            value_path=base+['shared_transfer_capacities', direction]))
    for index, _ in enumerate(value['scenarios'][scenario]['shared_links']):
        for field in ('h2d_bytes_per_second', 'd2h_bytes_per_second'):
            targets.append(dict(context_path=[context_index, 'shared_links', index, field],
                                value_path=base+['shared_links', index, field]))
    return targets


def earliest_observation(raw):
    return min(record['started_unix_seconds'] for record in raw['records'])


def attach_transfer_prices(profile, record_path, *, context_name, scenario, max_age_seconds=None):
    """Rebind a detailed profile with one fresh transfer record applied.

    The record must have been produced under the profile's exact source and
    execution context (including input/output storage identity), so run the
    producer with the job's own paths and CPU affinity. Existing price bindings
    and artifacts are kept; overlapping targets are refused by the binder, so
    this cannot silently replace an already-bound transfer price. Rebinding
    never renews the record's original observation age.
    """
    import json
    from pathlib import Path
    from .calibration_cache import read_calibration_record
    from .detailed_calibration import bind_detailed_profile, sha256_file
    path = Path(record_path).expanduser().resolve(strict=True)
    dependencies = json.loads(path.read_text(encoding='utf-8'))['dependencies']
    if (dependencies.get('source_sha256') != profile['source_sha256'] or
            dependencies.get('execution_context') != profile['execution_context']):
        raise ValueError('Transfer record source or execution context differs from profile')
    record = read_calibration_record(path, kind=KIND, name=NAME, dependencies=dependencies,
                                     max_age_seconds=max_age_seconds)['record']
    found = [index for index, row in enumerate(profile['contexts']) if row.get('name') == context_name]
    if len(found) != 1:
        raise ValueError('One named detailed context required')
    index = found[0]
    contexts = deepcopy(profile['contexts'])
    contexts[index] = apply_transfer_prices(contexts[index], record['value'], scenario)
    binding = dict(artifact=str(path), kind=KIND, name=NAME, dependencies=dependencies,
                   max_age_seconds=max_age_seconds,
                   targets=transfer_price_targets(index, record['value'], scenario))
    artifacts = dict(profile['component_artifacts'])
    artifacts[str(path)] = sha256_file(path)
    return bind_detailed_profile(contexts, profile['execution_context'],
        sources=profile['source_sha256'], component_artifacts=artifacts,
        limitations=list(profile['limitations'])+[
            f'Transfer prices from {NAME} scenario {scenario!r}: observed pinned-copy spans under recorded load, not guaranteed link ceilings.'],
        price_bindings=list(profile.get('price_bindings') or [])+[binding])
