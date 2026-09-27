"""Source-derived sync_file_range requests after each binary write.

This ledger counts requests and dependencies, not Linux writeback latency.
The last submitted interval and unsent tail still require closing fsync.
"""


def binary_writeback_work(events, interval_bytes=64 << 20, enabled=True):
    if isinstance(interval_bytes, bool) or not isinstance(interval_bytes, int) or interval_bytes < 0:
        raise ValueError('Writeback interval must be a nonnegative integer')
    active = bool(enabled and interval_bytes)
    offset = started = waited = 0
    rows = []
    for index, event in enumerate(events):
        length = event['bytes']
        if isinstance(length, bool) or not isinstance(length, int) or length < 0:
            raise ValueError('Write size must be a nonnegative integer')
        offset += length
        actions = []
        while active and offset - started >= interval_bytes:
            actions.append(dict(kind='submit',offset=started,bytes=interval_bytes))
            started += interval_bytes
            if started - waited >= 2 * interval_bytes:
                actions.append(dict(kind='wait',offset=waited,bytes=interval_bytes))
                actions.append(dict(kind='drop_cache',offset=waited,bytes=interval_bytes))
                waited += interval_bytes
        rows.append(dict(write_index=index,offset_after=offset,actions=actions))
    return dict(enabled=active,interval_bytes=interval_bytes,events=rows,
                submitted_bytes=started,waited_bytes=waited,
                submit_calls=started // interval_bytes if active else 0,
                wait_calls=waited // interval_bytes if active else 0,
                fadvise_calls=waited // interval_bytes if active else 0,
                submitted_not_waited_bytes=started-waited,
                unsubmitted_tail_bytes=offset-started,
                not_explicitly_waited_before_fsync_bytes=offset-waited,
                scope='Source requests per array. Unwaited bytes may already be durable; '
                      'these are not inferred dirty bytes or additional physical traffic.')
