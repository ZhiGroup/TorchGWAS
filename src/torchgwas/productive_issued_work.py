"""Bounded full-chunk workload for issued source awaiting matrix output.

An output event does not locate a chunk within the reader, GPU streams or
result queue. Replaying the entire issued chunk is therefore a conservative
nominal workload ceiling for those stages at the held checkpoint. Counts are
not remaining service times: prices, contention and unobserved queue state
still need a finite conditional completion schedule.
"""

from copy import deepcopy

from .analytical_plan_cache import input_identity
from .layout_frontier import KIND as FRONTIER_KIND
from .pgen_work_bounds import PgenHeaderWork


def _span(value, name):
    if (not isinstance(value, (list, tuple)) or len(value) != 2 or
            any(type(x) is not int for x in value) or
            not 0 <= value[0] < value[1]):
        raise ValueError('Nonempty ' + name + ' required')
    return tuple(value)


def productive_issued_work(boundary, frontier, header, *, covariate_rank,
                           max_partitions=16, max_pending_chunks=64,
                           max_source_records=1_000_000,
                           max_signatures=256):
    """Inspect only source chunks issued but lacking a matrix/part completion.

    Dense matrix progress can split an issued chunk. Its whole source chunk is
    retained, so this report may deliberately count work already performed.
    Indexed progress identifies completed chunks exactly, including empty
    output. Every pending chunk keeps its original device and phenotype panel.
    """
    for name, value in (('max_partitions', max_partitions),
                        ('max_pending_chunks', max_pending_chunks),
                        ('max_source_records', max_source_records),
                        ('max_signatures', max_signatures)):
        if type(value) is not int or value < 1:
            raise ValueError('Positive bounded ' + name + ' required')
    if (not isinstance(header, PgenHeaderWork) or
            not isinstance(frontier, dict) or
            frontier.get('kind') != FRONTIER_KIND or
            not isinstance(boundary, dict) or
            boundary.get('kind') != 'torchgwas.productive_output_boundary.v1' or
            boundary.get('valid') is not True or
            any(boundary.get(name) != frontier.get(name)
                for name in ('issued_revision', 'written_events', 'reduction')) or
            frontier.get('input_identity') != header.input_identity or
            input_identity(header.input_identity['path']) != header.input_identity or
            not isinstance(boundary.get('partitions'), list) or
            not 0 < len(boundary['partitions']) <= max_partitions):
        raise ValueError('Current bound source and output frontier required')
    samples = header._header.sample_ct
    if (type(covariate_rank) is not int or
            not 0 <= covariate_rank < samples - 2):
        raise ValueError('Bound residual covariate rank required')
    mode = frontier['reduction']
    traits = frontier['total_traits']
    job_lo, job_hi = _span(frontier['job_variant_range'], 'job variant range')
    original = {}
    pending = []
    records = 0
    for part in boundary['partitions']:
        if (not isinstance(part, dict) or
                not isinstance(part.get('id'), str) or not part['id'] or
                part['id'] in original or
                not isinstance(part.get('device'), str) or
                not part['device'].startswith('cuda:')):
            raise ValueError('Unique bound issued partition required')
        lo, hi = _span(part.get('variant_range'), 'issued variant range')
        a, b = _span(part.get('trait_range'), 'issued phenotype range')
        if (not job_lo <= lo < hi <= job_hi or b > traits or
                (mode == 'jagwas' and (a, b) != (0, traits))):
            raise ValueError('Issued partition differs from original job')
        cursor = part.get('issued_to')
        spans = part.get('issued_ranges')
        if (type(cursor) is not int or not lo <= cursor <= hi or
                not isinstance(spans, list)):
            raise ValueError('Bound issued source prefix required')
        expected = lo
        reserved = []
        for span in spans:
            first, last = _span(span, 'issued chunk range')
            if first != expected or last > cursor:
                raise ValueError('Issued chunks are not a contiguous prefix')
            expected = last
            reserved.append((first, last))
        if expected != cursor:
            raise ValueError('Issued cursor differs from reserved source')
        if mode is None:
            matrix_to = part.get('matrix_written_to')
            if type(matrix_to) is not int or not lo <= matrix_to <= cursor:
                raise ValueError('Bound dense matrix prefix required')
            incomplete = [span for span in reserved if span[1] > matrix_to]
        else:
            unfinished = part.get('issued_not_indexed_written')
            if not isinstance(unfinished, list):
                raise ValueError('Bound indexed unfinished chunks required')
            incomplete = [_span(span, 'indexed unfinished chunk') for span in unfinished]
            if len(set(incomplete)) != len(incomplete) or any(
                    span not in reserved for span in incomplete):
                raise ValueError('Indexed unfinished chunk was not issued once')
        original[part['id']] = dict(device=part['device'],
                                    variant_range=[lo, hi],
                                    trait_range=[a, b], issued_to=cursor)
        for first, last in incomplete:
            pending.append(dict(id=part['id'], device=part['device'],
                                variant_range=[first, last],
                                trait_range=[a, b]))
            records += last - first
    if len(pending) > max_pending_chunks or records > max_source_records:
        raise ValueError('Issued full-chunk workload exceeds bounded budget')
    for rectangle in frontier['rectangles']:
        spec = original.get(rectangle['id'])
        if (spec is None or spec['device'] != rectangle['device'] or
                spec['trait_range'] != rectangle['trait_range'] or
                [spec['issued_to'], spec['variant_range'][1]] !=
                rectangle['variant_range']):
            raise ValueError('Unissued frontier differs from issued partition')
    chunks = []
    total = dict(read_bytes=0, decode_input_bytes=0,
                 native_ld_base_update_bytes=0,
                 native_ld_replay_packed_bytes=0,
                 h2d_bytes=0, fp32_gemm_flops=0,
                 fp64_projection_flops=0)
    units = {}
    per_device = {}
    for part in pending:
        first, last = part['variant_range']
        work = header.bounds(first, last, max_records=max_source_records,
                             max_signatures=max_signatures)
        if (work['input_identity'] != header.input_identity or
                tuple(work['variant_range']) != (first, last)):
            raise ValueError('Issued source chunk changed during inspection')
        width = part['trait_range'][1] - part['trait_range'][0]
        markers = last - first
        row = dict(part, read_bytes=work['read_bytes'],
                   decode_input_bytes=work['decode_input_bytes'],
                   native_ld_base_update_bytes=work['native_ld_base_update_bytes'],
                   native_ld_replay_packed_bytes=work['native_ld_replay_packed_bytes'],
                   source_units=deepcopy(work['source_units']),
                   h2d_bytes=markers * samples,
                   fp32_gemm_flops=2 * markers * samples *
                       (width + covariate_rank + 1),
                   fp64_projection_flops=(2 * markers * width * width
                       if mode == 'jagwas' else 0))
        chunks.append(row)
        device_work = per_device.setdefault(part['device'],
            dict(h2d_bytes=0, fp32_gemm_flops=0,
                 fp64_projection_flops=0))
        for name in total:
            total[name] += row[name]
            if name in device_work:
                device_work[name] += row[name]
        for name, (low, high) in row['source_units'].items():
            pair = units.setdefault(name, [0, 0])
            pair[0] += low
            pair[1] += high
    if input_identity(header.input_identity['path']) != header.input_identity:
        raise ValueError('PGEN input changed during issued-work inspection')
    return dict(kind='torchgwas.productive_issued_work.v1',
                input_identity=deepcopy(header.input_identity),
                issued_revision=boundary['issued_revision'],
                written_events=boundary['written_events'],
                reduction=mode, samples=samples,
                covariate_rank=covariate_rank,
                pending_chunks=len(chunks), pending_source_records=records,
                original_partitions=deepcopy(original),
                chunks=chunks, total_work=total,
                total_source_units=units, per_device_work=per_device,
                prediction_complete=False, selection_validated=False,
                scope='Bounded full-chunk nominal read/decode/H2D/GEMM workload for issued chunks without matrix or indexed-part completion. All may already be in flight or partly complete; this is a conservative conditional workload ceiling, not remaining service time, a completion upper, or switch authorization.')
