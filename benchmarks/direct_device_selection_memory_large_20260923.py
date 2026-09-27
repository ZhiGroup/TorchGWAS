"""Check exact installed nonzero requests for the huge-panel selector shape."""
import argparse
import json
from pathlib import Path

from torchgwas.detailed_calibration import sha256_file, source_identity
from torchgwas.geometry_collection import write_record
from torchgwas.nonzero_memory import device_selection_memory
from torchgwas.selection_geometry import device_selection_shape


parser = argparse.ArgumentParser()
parser.add_argument('--census', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()
census = json.loads(Path(args.census).read_text())
source = source_identity()
for name in ('reduce.py', 'selection_geometry.py', 'nonzero_memory.py'):
    if census['source_sha256'].get(name) != source[name]:
        raise ValueError('Workspace census predates the current selector: ' + name)
if device_selection_shape(512, 600000, 1 << 20) != (2048, 512, 293):
    raise ValueError('Huge-panel production selector shape changed')
work = device_selection_memory(22250, 512, 600000,
    census=census, context=census['context'])
shapes = {(row['rows'], row['traits'], row['cells']) for row in work['block_shapes']}
if shapes != {(512, 2048, 1048576), (512, 1984, 1015808)}:
    raise ValueError('Memory admission omitted a production mask extent')
report = dict(census_sha256=sha256_file(args.census),
    source_sha256=work['source_sha256'], context=census['context'],
    selection_shape=[2048, 512, 293], block_shapes=work['block_shapes'],
    selection_gpu_bytes=work['selection_gpu_bytes'],
    maximum_selection_cells=work['maximum_selection_cells'],
    selected_payload_bytes_per_block=work['selected_payload_bytes_per_block'],
    scope='Exact installed CUB requests at both mask extents plus a conservative selector-local tensor bound for one 512-by-600000 source chunk. No whole-job memory or runtime claim.')
write_record(args.out, report)
print(json.dumps(dict(selection_shape=report['selection_shape'],
    mask_extents=sorted(row['cells'] for row in work['block_shapes']),
    selection_gpu_bytes=work['selection_gpu_bytes'])), flush=True)
