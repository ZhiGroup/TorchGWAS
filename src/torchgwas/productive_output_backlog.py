"""Conditional array payload still associated with issued, unwritten output.

This counts a conservative *maximum* payload for the held output checkpoint.
An indexed part may already be partially written, and a dense statistic array
may lead the reported common prefix. No positive remaining-byte lower bound is
claimed, and framing, fsync, queues and service time are outside this ledger.
"""

from copy import deepcopy

from .productive_occupancy import whole_layout_survivors


def productive_output_backlog(boundary, *, total_traits, store_beta=True,
                              store_variant_df=False, occupancy_scenario=None,
                              max_pending_ranges=1024,
                              max_cluster_visits=100_000):
    """Count scenario-bound array payload for issued output without completion.

    Reduced-output scenarios use the same global pair/variant coordinate as
    unissued candidate layouts, so retiling cannot invent a different set of
    survivors. The result is an upper *workload* input, not an elapsed-time
    bound or permission to change the running layout.
    """
    if (not isinstance(boundary, dict) or
            boundary.get('kind') != 'torchgwas.productive_output_boundary.v1' or
            boundary.get('valid') is not True or
            not isinstance(boundary.get('partitions'), list) or
            type(total_traits) is not int or total_traits < 1 or
            type(store_beta) is not bool or
            type(store_variant_df) is not bool or
            type(max_pending_ranges) is not int or max_pending_ranges < 1 or
            type(max_cluster_visits) is not int or max_cluster_visits < 1):
        raise ValueError('Valid bounded productive output checkpoint required')
    mode = boundary.get('reduction')
    if mode not in (None, 'significant', 'jagwas'):
        raise ValueError('Unknown productive output mode')
    if ((mode is None) != (occupancy_scenario is None)):
        raise ValueError('Reduced output requires one explicit occupancy scenario')
    if mode == 'jagwas' and not store_beta:
        raise ValueError('JAGWAS has no configurable beta output')
    if mode is not None and store_variant_df:
        raise ValueError('Dense df sidecar applies only to dense output')
    pending = []
    rows = []
    for part in boundary['partitions']:
        if (not isinstance(part, dict) or
                not isinstance(part.get('id'), str) or not part['id'] or
                not isinstance(part.get('trait_range'), (list, tuple)) or
                len(part['trait_range']) != 2 or
                any(type(value) is not int for value in part['trait_range']) or
                not 0 <= part['trait_range'][0] < part['trait_range'][1] <= total_traits):
            raise ValueError('Bounded output partition geometry required')
        trait = part['trait_range']
        if mode == 'jagwas' and list(trait) != [0, total_traits]:
            raise ValueError('JAGWAS backlog requires the full phenotype panel')
        if mode is None:
            pairs = part.get('issued_not_matrix_written_pairs')
            df_markers = part.get('issued_not_df_written_markers')
            if (type(pairs) is not int or pairs < 0 or
                    type(df_markers) is not int or df_markers < 0):
                raise ValueError('Bound dense pending pair/df counts required')
            rows.append(dict(id=part['id'], pending_pairs=pairs,
                             pending_df_markers=df_markers,
                             selected_rows=None,
                             array_payload_bytes_upper=(
                                 4 * (1 + store_beta) * pairs +
                                 (4 * df_markers if store_variant_df else 0))))
        else:
            spans = part.get('issued_not_indexed_written')
            if not isinstance(spans, list):
                raise ValueError('Bound indexed pending ranges required')
            keys = []
            for span in spans:
                if (not isinstance(span, (list, tuple)) or len(span) != 2 or
                        any(type(value) is not int for value in span) or
                        not 0 <= span[0] < span[1]):
                    raise ValueError('Nonempty indexed pending source range required')
                key = str(len(pending))
                pending.append(dict(id=key, variant_range=list(span),
                                    trait_range=list(trait)))
                keys.append(key)
            rows.append(dict(id=part['id'], pending_ranges=len(keys),
                             pending_range_keys=keys))
    if len(pending) > max_pending_ranges:
        raise ValueError('Issued output exceeds bounded pending-range budget')
    scenario_report = None
    if mode is not None:
        retained = {}
        if pending:
            scenario_report = whole_layout_survivors(
                dict(kind='torchgwas.pgen_layout_source_floor.v1',
                     reduction=mode, total_traits=total_traits,
                     partitions=pending), occupancy_scenario,
                max_partitions=max_pending_ranges,
                max_cluster_visits=max_cluster_visits)
            retained = scenario_report['retained_ranges']
        unit = (28 if store_beta else 24) if mode == 'significant' else 16
        for row in rows:
            selected = sum(retained[key][0] for key in row.pop('pending_range_keys'))
            row['selected_rows'] = selected
            row['array_payload_bytes_upper'] = unit * selected
    return dict(kind='torchgwas.productive_output_backlog.v1',
                issued_revision=boundary['issued_revision'],
                written_events=boundary['written_events'],
                reduction=mode, store_beta=store_beta,
                store_variant_df=store_variant_df,
                occupancy_scenario=deepcopy(occupancy_scenario),
                partitions=rows,
                array_payload_bytes_upper=sum(row['array_payload_bytes_upper']
                                              for row in rows),
                scenario_report=scenario_report,
                prediction_complete=False, selection_validated=False,
                scope='Conditional maximum pending array payload at a held issued/output checkpoint. Dense array prefixes and indexed part writes can lead their completion events, so no positive remaining-byte lower bound is claimed. Excludes framing, queue/stream state, writer service, fsync and manifest/directory publication; not elapsed time or switch authorization.')
