"""Read-only bounded observation of a multi-GPU indexed result queue."""
import queue
import time


def snapshot_indexed_result_queue(result_queue, bounds, devices, finished, *, max_items=64):
    """Copy queue references under its mutex, then count host-array payload."""
    if (not isinstance(result_queue, queue.Queue) or
            type(result_queue.maxsize) is not int or result_queue.maxsize < 1 or
            type(max_items) is not int or max_items < 1 or
            result_queue.maxsize > max_items or
            not isinstance(bounds, (list, tuple)) or
            not isinstance(devices, (list, tuple)) or
            len(bounds) != len(devices) or len(bounds) < 2 or
            len(set(devices)) != len(devices)):
        raise ValueError('Bounded indexed multi-GPU result queue required')
    spans = []
    for span, device in zip(bounds, devices):
        if (not isinstance(span, (list, tuple)) or len(span) != 2 or
                any(type(value) is not int for value in span) or
                not 0 <= span[0] < span[1] or
                not isinstance(device, str) or not device):
            raise ValueError('Explicit indexed queue source partitions required')
        spans.append((span[0], span[1], device))
    if any(a < d and c < b for i,(a,b,_) in enumerate(spans)
           for c,d,_ in spans[i+1:]):
        raise ValueError('Indexed queue source partitions overlap')
    began = time.perf_counter()
    with result_queue.mutex:
        if len(result_queue.queue) > max_items:
            raise ValueError('Indexed queue observation exceeds item budget')
        items = list(result_queue.queue)
        anchor = time.perf_counter()
    rows = []
    sentinels = 0
    payload = 0
    for item in items:
        if item is finished:
            sentinels += 1
            continue
        if (not isinstance(item, tuple) or len(item) < 4 or
                type(item[0]) is not int or type(item[1]) is not int or
                item[0] >= item[1]):
            raise ValueError('Indexed queue has an unbound result')
        first, last = item[:2]
        owners = [device for lo,hi,device in spans if lo <= first < last <= hi]
        if len(owners) != 1:
            raise ValueError('Indexed queue result crosses source partitions')
        bytes_here = 0
        for value in item[2:]:
            if value is None:
                continue
            nbytes = getattr(value, 'nbytes', None)
            if type(nbytes) is not int or nbytes < 0:
                raise ValueError('Indexed queue result has unknown host payload')
            bytes_here += nbytes
        rows.append(dict(device=owners[0], variant_range=[first,last],
                         resident_array_bytes=bytes_here))
        payload += bytes_here
    return dict(kind='torchgwas.indexed_result_queue_observation.v1',
        capture_started_seconds=began,capture_anchor_seconds=anchor,
        capture_finished_seconds=time.perf_counter(),
        capacity=result_queue.maxsize,queued_items=len(items),
        queued_results=len(rows),queued_sentinels=sentinels,
        resident_array_bytes=payload,results=rows,
        scope='One atomic snapshot of queued owned host results in the multi-GPU indexed producer queue. A result already delivered to the indexed writer, producer/GPU work, part encoding/fsync and final manifest are outside this queue.')
