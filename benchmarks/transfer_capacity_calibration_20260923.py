"""Produce a held-out-checked pinned transfer_capacity record for one GPU set.

The record's dependencies are this process's source identity, detailed
execution context (CPU affinity, NUMA policy, GPU/PCI identity, environment)
and measurement protocol, so it binds only to a job run in the same context.
The report is always written; the record only when every check passes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import (bind_detailed_profile, execution_context,
    sha256_file, source_identity, validate_detailed_profile)
from torchgwas.transfer_calibration import (KIND, NAME, attach_transfer_prices,
    earliest_observation, measure_transfer_groups, measurement_protocol,
    summarize_transfer_groups)


def _links(values):
    return [tuple(item.split(',')) for item in values]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--links', nargs='*', default=[],
                        help='comma-separated device groups sharing a link, e.g. cuda:0,cuda:1')
    parser.add_argument('--cpus', type=int, nargs='*')
    parser.add_argument('--expect-local-node', action='store_true',
                        help='require each pinned buffer on its GPU PCI NUMA node')
    parser.add_argument('--size-mib', type=int, default=32)
    parser.add_argument('--copies', type=int, default=8)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--warmups', type=int, default=2)
    parser.add_argument('--rounds', type=int, default=10)
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--tolerance', type=float, default=0.10)
    parser.add_argument('--idle-utilization-limit', type=int, default=10)
    parser.add_argument('--max-age-seconds', type=float, default=6*3600.)
    parser.add_argument('--input-path', type=Path, required=True)
    parser.add_argument('--output-path', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.cpus:
        os.sched_setaffinity(0, set(args.cpus))
    devices = list(args.devices)
    implementation = hashlib.sha256(Path(__file__).read_bytes()+
        Path(measure_transfer_groups.__code__.co_filename).read_bytes()).hexdigest()
    before = execution_context(devices, input_path=args.input_path, output_path=args.output_path)
    source = source_identity()
    expected = None
    if args.expect_local_node:
        expected = {device: int(before['devices'][device]['numa_node']) for device in devices}
        if any(node < 0 for node in expected.values()):
            raise ValueError('GPU PCI device reports no NUMA node')
    raw = measure_transfer_groups(devices, _links(args.links), size_bytes=args.size_mib << 20,
        copies=args.copies, samples=args.samples, warmups=args.warmups, rounds=args.rounds,
        seed=args.seed, idle_utilization_limit=args.idle_utilization_limit,
        uuids={device: before['devices'][device]['uuid'] for device in devices})
    after = execution_context(devices, input_path=args.input_path, output_path=args.output_path)
    summary = summarize_transfer_groups(raw, tolerance=args.tolerance, expected_nodes=expected)
    if after != before:
        summary['failures'].append(dict(check='execution_context_unchanged'))
        summary['qualified'] = False
        summary.pop('value', None)
    protocol = measurement_protocol(raw, tolerance=args.tolerance,
                                    implementation_sha256=implementation)
    report = dict(kind='torchgwas.transfer_capacity_calibration_report.v1',
                  created_unix_seconds=time.time(), protocol=protocol, raw=raw,
                  summary=summary, published=None)
    if summary['qualified']:
        deps = dict(source_sha256=source, execution_context=before, measurement_protocol=protocol)
        report['published'] = CalibrationParameterCache(args.cache).store(
            KIND, NAME, summary['value'], dependencies=deps,
            provenance=dict(report=str(args.report.resolve()), host=os.uname().nodename),
            max_age_seconds=args.max_age_seconds, observed_unix_seconds=earliest_observation(raw))
        # Real-context check: attach to a skeleton profile, then validate it
        # against a freshly captured context as a later job would.
        skeleton = bind_detailed_profile([dict(name='transfer', devices=devices,
            profiles={device: {} for device in devices})], before, sources=source,
            component_artifacts={report['published']['path']: sha256_file(report['published']['path'])},
            limitations=['Skeleton profile: transfer prices only.'])
        report['binding_check'] = {}
        for scenario in ('low', 'high'):
            profile = attach_transfer_prices(skeleton, report['published']['path'],
                                             context_name='transfer', scenario=scenario)
            current = execution_context(devices, input_path=args.input_path,
                                        output_path=args.output_path)
            evidence = validate_detailed_profile(profile, current, sources=source)['price_evidence']
            report['binding_check'][scenario] = dict(status=evidence['status'],
                verified_targets=evidence['verified_targets'],
                age_seconds=evidence['bindings'][0]['age_seconds'])
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=1, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    cells = {key: dict(low=round(c['low']/1e9, 3), median=round(c['median']/1e9, 3),
                       high=round(c['high']/1e9, 3), holdout=round(c['holdout_median']/1e9, 3))
             for key, c in summary['cells'].items()}
    print(json.dumps(dict(qualified=summary['qualified'], failures=summary['failures'],
                          GBps=cells, published=report['published'],
                          binding_check=report.get('binding_check')), indent=1))


if __name__ == '__main__':
    main()
