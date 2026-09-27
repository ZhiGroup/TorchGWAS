"""Reconcile the metadata selector trace to the recorded actual CUDA calls."""
import argparse
import hashlib
import json
from pathlib import Path
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record

parser = argparse.ArgumentParser()
parser.add_argument('--census', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()
source = source_identity()
capture = json.loads(Path(args.census).read_text())
for name in ('reduce.py', 'selection_geometry.py'):
    assert capture['source_sha256'][name] == source[name]
assert capture['durations_recorded'] is False
fields = ['shape', 'stride', 'dtype', 'device_type', 'bytes']
def operations(steps):
    return [dict(op=row['op'], **{key:[{field:v[field] for field in fields} for v in row[key]]
        for key in ['inputs', 'outputs']}) for row in steps
        if row['op'] != 'aten.detach.default']
checks = []
for row in capture['rows']:
    counts = [b['retained'] for b in row['blocks']]
    predicted = device_significant_tensor_work(40, row['B'], row['K'], counts, max_cells=row['max_cells'])
    assert predicted['torch_version'] == capture['torch_version']
    expected = operations(row['observed_tensor_steps']); actual = operations(predicted['steps'])
    if expected != actual:
        first = next((i for i, (a,b) in enumerate(zip(expected,actual)) if a != b), min(len(expected),len(actual)))
        raise AssertionError(dict(shape=[row['B'],row['K']], mode=row['mode'], first=first,
            expected=expected[first:first+1], actual=actual[first:first+1], lengths=[len(expected),len(actual)]))
    assert row['independent_cpu_arrays_equal'] is True
    assert predicted['selected_payload_d2h_bytes'] == sum(b['payload_bytes'] for b in row['blocks'])
    nonzero_copies = [v for v in row['transfer_events'] if 'Device -> Pinned' in v['name']]
    selected_copies = [v for v in row['transfer_events'] if 'Device -> Pageable' in v['name']]
    assert len(nonzero_copies) == predicted['nonzero_calls'] and all(v['bytes'] == 4 for v in nonzero_copies)
    assert len(selected_copies) == predicted['selected_copy_calls']
    assert sum(v['bytes'] for v in selected_copies) == predicted['selected_payload_d2h_bytes']
    assert len(row['transfer_events']) == len(nonzero_copies) + len(selected_copies)
    assert predicted['nonzero_count_d2h_bytes'] == sum(v['bytes'] for v in nonzero_copies)
    assert predicted['selector_d2h_bytes'] == sum(v['bytes'] for v in row['transfer_events'])
    checks.append(dict(B=row['B'], K=row['K'], max_cells=row['max_cells'], mode=row['mode'],
        tensor_operations=len(actual), nonzero_calls=predicted['nonzero_calls'],
        nonzero_count_bytes=4*predicted['nonzero_calls'], selected_copy_calls=predicted['selected_copy_calls'],
        selected_payload_d2h_bytes=predicted['selected_payload_d2h_bytes'], operations_and_transfers_exact=True))
assert source == source_identity()
write_record(Path(args.out) / 'audit.json', dict(source_sha256=source,
    census_sha256=hashlib.sha256(Path(args.census).read_bytes()).hexdigest(),
    torch_version=capture['torch_version'], device=capture['device'], checks=checks,
    scope='Exact source operations/shapes/dtypes/strides and separately captured CUDA transfer counts/bytes. '
          'CPU NumPy detach aliases are excluded from operation comparison. No duration-based prices or autotune qualification.'))
print(json.dumps(dict(cases=len(checks), all_exact=True)), flush=True)
