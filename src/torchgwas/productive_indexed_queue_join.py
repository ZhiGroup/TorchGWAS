"""Bind queued JAGWAS results to an exact productive source/output frontier."""

def bind_jagwas_result_queue(boundary, observation, *, active_writer=None):
    if (not isinstance(boundary, dict) or
            boundary.get('kind') != 'torchgwas.productive_output_boundary.v1' or
            boundary.get('reduction') != 'jagwas' or
            boundary.get('valid') is not True or
            not isinstance(boundary.get('partitions'), list) or
            not isinstance(observation, dict) or
            observation.get('kind') != 'torchgwas.indexed_result_queue_observation.v1' or
            not isinstance(observation.get('results'), list) or
            observation.get('queued_results') != len(observation['results'])):
        raise ValueError('Bound JAGWAS output and live result queue required')
    pending = {}
    for part in boundary['partitions']:
        if (not isinstance(part, dict) or
                not isinstance(part.get('id'), str) or
                not isinstance(part.get('device'), str) or
                not isinstance(part.get('issued_not_indexed_written'), list)):
            raise ValueError('Bound indexed source partition required')
        for span in part['issued_not_indexed_written']:
            if (not isinstance(span, list) or len(span) != 2 or
                    any(type(value) is not int for value in span) or
                    not 0 <= span[0] < span[1] or
                    tuple(span) in pending):
                raise ValueError('Unique pending JAGWAS source ranges required')
            pending[tuple(span)] = (part['id'], part['device'],part['trait_range'])
    queued = {}
    payload = 0
    for row in observation['results']:
        if (not isinstance(row, dict) or
                not isinstance(row.get('variant_range'), list) or
                len(row['variant_range']) != 2 or
                any(type(value) is not int for value in row['variant_range']) or
                not isinstance(row.get('device'), str) or
                type(row.get('resident_array_bytes')) is not int or
                row['resident_array_bytes'] < 0):
            raise ValueError('Bound queued JAGWAS result required')
        span = tuple(row['variant_range'])
        if span in queued or pending.get(span, (None,None))[1] != row['device']:
            raise ValueError('Queued JAGWAS result was not one unfinished issued chunk')
        queued[span] = dict(partition_id=pending[span][0],
                            device=row['device'],variant_range=list(span),
                            resident_array_bytes=row['resident_array_bytes'])
        payload += row['resident_array_bytes']
    if payload != observation.get('resident_array_bytes'):
        raise ValueError('Queued JAGWAS host payload differs from snapshot')
    active=None
    if active_writer is not None:
        if (not isinstance(active_writer,dict) or
                not isinstance(active_writer.get('variant_range'),list) or
                len(active_writer['variant_range'])!=2 or
                any(type(value) is not int for value in active_writer['variant_range']) or
                not isinstance(active_writer.get('trait_range'),list) or
                len(active_writer['trait_range'])!=2 or
                any(type(value) is not int for value in active_writer['trait_range']) or
                not isinstance(active_writer.get('device'),str)):
            raise ValueError('Bound active JAGWAS indexed writer range required')
        span=tuple(active_writer['variant_range'])
        owner=pending.get(span)
        if (span in queued or owner is None or
                owner[1]!=active_writer['device'] or
                owner[2]!=active_writer['trait_range']):
            raise ValueError('Active JAGWAS writer was not one unfinished issued chunk')
        active=dict(partition_id=owner[0],device=owner[1],
                    variant_range=list(span),trait_range=list(owner[2]))
    outside = [dict(partition_id=part_id,device=device,
                    variant_range=list(span))
               for span,(part_id,device,_) in pending.items() if span not in queued]
    unresolved=[row for row in outside if active is None or
                row['variant_range']!=active['variant_range']]
    return dict(kind='torchgwas.bound_jagwas_result_queue.v1',
        issued_revision=boundary['issued_revision'],
        written_events=boundary['written_events'],
        capture_anchor_seconds=observation['capture_anchor_seconds'],
        queued_source_chunks=list(queued.values()),
        queued_chunks=len(queued),queued_resident_array_bytes=payload,
        active_writer_source_chunk=active,
        completed_producer_chunks=len(queued)+int(active is not None),
        issued_not_part_written_chunks=len(pending),
        issued_not_part_written_outside_queue=outside,
        issued_not_part_written_upstream_or_unresolved=unresolved,
        prediction_complete=False,selection_validated=False,
        scope='Exact issued JAGWAS chunks currently in the shared result queue at one observed instant. Queued chunks have produced host results but still need indexed part consumption. An active writer chunk is producer-complete only when stable across the queue anchor. Other unfinished issued chunks may be upstream or between stages; their full-chunk workload remains conservative. Final metadata durability is not included.')
