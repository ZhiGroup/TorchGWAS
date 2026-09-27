"""Declared shared transfer links for compact multi-GPU work floors.

A transfer may traverse several nested links. Each link constrains the sum of
bytes from its device set; the necessary time is the maximum link load, not
the sum of link times. Topology and service ceilings are supplied externally.
"""
import math


def transfer_link_loads(work_by_device, links, *, direction, max_links=32):
    """Return per-link byte loads and a necessary service floor.

    `work_by_device` maps every active GPU to mandatory bytes. Links may
    overlap to represent nested shared buses; no topology is inferred.
    """
    if not isinstance(work_by_device, dict) or not work_by_device or any(
            not isinstance(key, str) or not key or type(value) is not int or value < 0
            for key, value in work_by_device.items()):
        raise ValueError('Explicit nonnegative device transfer bytes required')
    if direction not in ('h2d', 'd2h'):
        raise ValueError('Explicit transfer direction required')
    if type(max_links) is not int or max_links < 1:
        raise ValueError('Positive transfer-link budget required')
    if links is None:
        links = ()
    if not isinstance(links, (list, tuple)) or len(links) > max_links:
        raise ValueError('Bounded transfer-link list required')
    known = set(work_by_device)
    rows = []
    normalized = []
    for index, link in enumerate(links):
        if (not isinstance(link, dict) or
                set(link) != {'devices', 'h2d_bytes_per_second', 'd2h_bytes_per_second'}):
            raise ValueError('Existing execution-graph shared-link fields required')
        devices = link['devices']
        if (not isinstance(devices, (list, tuple)) or not devices or
                any(not isinstance(device, str) for device in devices) or
                len(set(devices)) != len(devices) or not set(devices) <= known):
            raise ValueError('Link devices must be unique active GPUs')
        for name in ('h2d_bytes_per_second', 'd2h_bytes_per_second'):
            value = link[name]
            if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(value) or value <= 0):
                raise ValueError('Positive finite transfer-link ceiling required')
        normalized.append(dict(devices=list(devices),
                               h2d_bytes_per_second=link['h2d_bytes_per_second'],
                               d2h_bytes_per_second=link['d2h_bytes_per_second']))
        capacity = link[direction + '_bytes_per_second']
        amount = sum(work_by_device[device] for device in devices)
        seconds = amount / capacity
        if not math.isfinite(seconds):
            raise ValueError('Transfer-link load overflow')
        rows.append(dict(index=index, direction=direction, devices=list(devices),
                         bytes=amount, bytes_per_second=capacity,
                         floor_seconds=seconds))
    return dict(declarations=normalized, links=rows,
                floor_seconds=max((row['floor_seconds'] for row in rows), default=0.))
