"""Compare audited API CPU with a matched unhooked-before/after bracket.

A successful bracket supports total CPU transfer only in the observed context.
It cannot certify GIL partition, scheduling accuracy, or transfer to other shapes.
"""
import math
import statistics
from .gil_service import gil_service_prices


def _positive(value, name, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError('Invalid ' + name)
    return value


def _workers(probe):
    result = {}
    for context in probe['results']:
        for worker in context['workers']:
            key = (context['mode'], worker['device'])
            if key in result:
                raise ValueError('Duplicate worker context')
            result[key] = worker
    if not result:
        raise ValueError('Worker contexts required')
    return result


def _unhooked_rows(worker):
    if worker.get('complete') is not True:
        raise ValueError('Complete unhooked worker required')
    grouped = {}
    for row in worker['rows']:
        count, repeat, phase = row['sample_count'], row['repeat'], row['phase']
        if any(isinstance(v, bool) or not isinstance(v, int) for v in (count, repeat, phase)) or count < 1 or repeat < 0 or phase not in (0, 1):
            raise ValueError('Invalid unhooked sample context')
        name = row['primitive']
        if not isinstance(name, str) or not name:
            raise ValueError('Primitive name required')
        raw, total, empty, cpu = [_positive(row[key], key, zero=True) for key in (
            'raw_cpu_seconds_per_call', 'raw_cpu_sum_seconds', 'empty_cpu_seconds_per_call', 'cpu_seconds_per_call')]
        if not math.isclose(raw * count, total, rel_tol=1e-10, abs_tol=1e-15) or not math.isclose(cpu, max(0., raw-empty), rel_tol=1e-10, abs_tol=1e-15):
            raise ValueError('Unhooked CPU observations do not conserve')
        pair = grouped.setdefault(name, {}).setdefault(repeat, {})
        if phase in pair:
            raise ValueError('Duplicate unhooked phase')
        pair[phase] = (cpu, count)
    if not grouped or any(len(repeats) < 3 or any(set(pair) != {0, 1} for pair in repeats.values()) for repeats in grouped.values()):
        raise ValueError('At least three complete unhooked repeats required')
    return grouped


def compare_gil_probe_bracket(before, audited, after, *, maximum_relative_difference):
    """Use an explicit diagnostic tolerance, never infer an autotune gate."""
    tolerance = _positive(maximum_relative_difference, 'relative comparison tolerance', zero=True)
    for probe in (before, after):
        if 'ld_audit' not in probe or probe['ld_audit'] is not None:
            raise ValueError('Unhooked LD_AUDIT control required')
        for key in ('host', 'torch_version', 'affinity', 'devices', 'primitive_bank', 'fixed_shape'):
            if key not in audited or probe.get(key) != audited[key]:
                raise ValueError('Probe context mismatch: ' + key)
        # These are the common executable sources. The two Python drivers
        # intentionally differ; the bank and clock/hook helper must not.
        expected = audited['source_sha256']
        common = set(probe['source_sha256']) & set(expected)
        required_suffixes = ('direct_jagwas_host_primitives.py', 'direct_calculator_gil_audit.c', 'libtorchgwas_gil_audit.so')
        for suffix in required_suffixes:
            keys = [key for key in common if key.endswith('/' + suffix)]
            if len(keys) != 1 or probe['source_sha256'][keys[0]] != expected[keys[0]]:
                raise ValueError('Shared measurement source mismatch: ' + suffix)
    if before['source_sha256'] != after['source_sha256']:
        raise ValueError('Unhooked driver changed during bracket')
    workers = [_workers(probe) for probe in (before, audited, after)]
    if not set(workers[0]) == set(workers[1]) == set(workers[2]):
        raise ValueError('Worker context coverage differs')
    contexts = {}
    for mode, device in sorted(workers[1]):
        key = (mode, device)
        price = gil_service_prices(audited, mode=mode, device=device,
            torch_version=audited['torch_version'], cpu_affinity=audited['affinity'])
        ends = [_unhooked_rows(worker[key]) for worker in (workers[0], workers[2])]
        rows = workers[1][key]['rows']
        coverage = {}
        for row in rows:
            coverage.setdefault(row['primitive'], {}).setdefault(row['repeat'], {})[int(row['recording'])] = row['sample_count']
        if any(set(end) != set(coverage) for end in ends):
            raise ValueError('Primitive coverage differs')
        observations = {}
        for name, repeats in coverage.items():
            for end in ends:
                if set(end[name]) != set(repeats):
                    raise ValueError('Repeat coverage differs')
                for repeat, pair in repeats.items():
                    if set(pair) != {0, 1} or len(set(pair.values())) != 1 or any(value[1] != pair[0] for value in end[name][repeat].values()):
                        raise ValueError('Sample coverage differs')
            repeat_means = [[statistics.fmean(v[0] for v in pair.values()) for pair in end[name].values()] for end in ends]
            before_cpu, after_cpu = [statistics.median(values) for values in repeat_means]
            center = (before_cpu + after_cpu) / 2
            reference_low, reference_high = min(before_cpu, after_cpu), max(before_cpu, after_cpu)
            audited_cpu = price['cpu_primitives'][name]
            disabled_cpu = statistics.median(row['cpu_seconds_per_call'] for row in rows if row['primitive'] == name and not row['recording'])
            resolved = reference_low > 0
            drift = abs(after_cpu-before_cpu) / center if center else None
            # Conservative comparison with BOTH endpoint medians. An overlap
            # with a broad drift envelope alone cannot pass this criterion.
            error = max(abs(audited_cpu/value-1) for value in (before_cpu, after_cpu)) if resolved else None
            stable = resolved and drift <= tolerance
            compatible = stable and error <= tolerance
            observations[name] = dict(before_cpu_seconds=before_cpu, after_cpu_seconds=after_cpu,
                unhooked_cpu_seconds=center, unhooked_repeat_mean_range=[min(sum(repeat_means, [])), max(sum(repeat_means, []))],
                corrected_audited_cpu_seconds=audited_cpu, recording_disabled_cpu_seconds=disabled_cpu,
                audit_relative_delta_interval=[audited_cpu/reference_high-1, audited_cpu/reference_low-1] if resolved else None,
                disabled_relative_delta_interval=[disabled_cpu/reference_high-1, disabled_cpu/reference_low-1] if resolved else None,
                resolution_limited=not resolved,
                bracket_relative_drift=drift, maximum_endpoint_relative_error=error,
                bracket_stable=stable, total_cpu_compatible=compatible,
                serial_transfer_qualified=False)
        contexts[mode+':'+str(device)] = dict(primitives=observations,
            total_cpu_compatible=all(row['total_cpu_compatible'] for row in observations.values()))
    return dict(contexts=contexts, maximum_relative_difference=tolerance,
        total_cpu_compatible=all(row['total_cpu_compatible'] for row in contexts.values()),
        serial_transfer_qualified=False, prediction_complete=False,
        scope='Matched fixed-API CPU before/audited/after comparison. Endpoint medians and repeat ranges describe observations, not confidence bounds. Changing server load remains a confounder. Passing an explicit CPU tolerance does not establish GIL occupancy transfer, kernel costs, shape transfer, end-to-end prediction, or autotune readiness.')
