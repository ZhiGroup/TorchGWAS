"""Remove completed JAGWAS producer stages visible in the result queue."""
from copy import deepcopy


_FIELDS = ('read_bytes', 'decode_input_bytes',
           'native_ld_base_update_bytes', 'native_ld_replay_packed_bytes',
           'h2d_bytes', 'fp32_gemm_flops', 'fp64_projection_flops')


def refine_issued_jagwas_with_queue(issued, joined):
    """Partition the issued full-chunk upper workload by observed queue state."""
    if (not isinstance(issued, dict) or
            issued.get('kind') != 'torchgwas.productive_issued_work.v1' or
            issued.get('reduction') != 'jagwas' or
            not isinstance(issued.get('chunks'), list) or
            issued.get('pending_chunks') != len(issued['chunks']) or
            not isinstance(joined, dict) or
            joined.get('kind') != 'torchgwas.bound_jagwas_result_queue.v1' or
            any(issued.get(name) != joined.get(name)
                for name in ('issued_revision', 'written_events')) or
            not isinstance(joined.get('queued_source_chunks'), list) or
            not isinstance(joined.get('issued_not_part_written_outside_queue'), list)):
        raise ValueError('Same-revision issued JAGWAS work and bound queue required')
    by_key = {}
    for chunk in issued['chunks']:
        key=(chunk['id'],chunk['device'],tuple(chunk['variant_range']))
        if key in by_key:
            raise ValueError('Duplicate issued JAGWAS chunk')
        by_key[key]=chunk
    queued = set()
    for row in joined['queued_source_chunks']:
        key=(row['partition_id'],row['device'],tuple(row['variant_range']))
        if key not in by_key or key in queued:
            raise ValueError('Queued result differs from issued full chunk')
        queued.add(key)
    if (joined.get('queued_chunks') != len(queued) or
            joined.get('issued_not_part_written_chunks') != len(by_key)):
        raise ValueError('Issued and queued JAGWAS chunk counts differ')
    active_row=joined.get('active_writer_source_chunk')
    active=None
    if active_row is not None:
        active=(active_row['partition_id'],active_row['device'],
                tuple(active_row['variant_range']))
        if active not in by_key or active in queued:
            raise ValueError('Active writer differs from issued full chunk')
    completed_keys=queued | ({active} if active is not None else set())
    if joined.get('completed_producer_chunks')!=len(completed_keys):
        raise ValueError('Bound JAGWAS producer-complete count differs')
    expected_outside = {key for key in by_key if key not in queued}
    reported_outside = {(row['partition_id'], row['device'],
                         tuple(row['variant_range']))
                        for row in joined['issued_not_part_written_outside_queue']}
    if reported_outside != expected_outside or len(reported_outside) != len(
            joined['issued_not_part_written_outside_queue']):
        raise ValueError('Bound queued and outside source chunks do not conserve')
    remaining = [row for key,row in by_key.items() if key not in completed_keys]
    completed = [row for key,row in by_key.items() if key in completed_keys]
    def sum_work(rows):
        return {name:sum(row[name] for row in rows) for name in _FIELDS}
    upper=sum_work(remaining)
    done=sum_work(completed)
    queued_done=sum_work([by_key[key] for key in queued])
    active_done=sum_work([] if active is None else [by_key[active]])
    if any(upper[name]+done[name]!=issued['total_work'][name]
           for name in _FIELDS):
        raise ValueError('Issued producer work does not conserve under queue split')
    per_device={}
    for row in remaining:
        device=per_device.setdefault(row['device'],dict(h2d_bytes=0,
            fp32_gemm_flops=0,fp64_projection_flops=0))
        for name in device:device[name]+=row[name]
    return dict(kind='torchgwas.issued_jagwas_queue_refinement.v1',
        input_identity=deepcopy(issued['input_identity']),
        issued_revision=issued['issued_revision'],
        written_events=issued['written_events'],
        queued_chunks_with_completed_producer=len(queued),
        active_writer_chunks_with_completed_producer=int(active is not None),
        queued_producer_work_completed=queued_done,
        active_writer_producer_work_completed=active_done,
        observed_producer_work_completed=done,
        upstream_or_unresolved_full_chunk_work_upper=upper,
        upstream_or_unresolved_chunks=len(remaining),
        per_device_upstream_or_unresolved_upper=per_device,
        prediction_complete=False,selection_validated=False,
        scope='Source read/decode, H2D and GPU work is already complete for queued and stably active JAGWAS results at the bound snapshot. Other issued chunks retain full-chunk upper workload because their producer state is unresolved. Output selection/encoding/fsync, future work, capacity contention and final publication are outside this refinement; no elapsed completion bound or switch authorization.')
