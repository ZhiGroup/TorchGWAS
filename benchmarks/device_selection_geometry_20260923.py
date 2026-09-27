"""Matched empty-selector control for source-equivalent block geometries."""
import hashlib
import argparse
import json
import os
from pathlib import Path
import statistics
import time

import torch

import torchgwas.reduce as reduction
from torchgwas.selection_geometry import device_selection_shape


ROWS = 512
TRAITS = 600000
MAX_CELLS = 1 << 20


def legacy_shape(rows, traits, max_cells):
    width = min(traits, max_cells)
    height = max(1, max_cells // width)
    blocks = ((rows + height - 1) // height) * ((traits + width - 1) // width)
    return width, height, blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scenario', choices=('empty', 'sparse'), default='empty')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError('A CUDA device is required')
    device = torch.device('cuda:0')
    t_stat = torch.zeros((ROWS, TRAITS), device=device, dtype=torch.float32)
    expected_retained = 0
    if args.scenario == 'sparse':
        row_index = torch.arange(ROWS, device=device)
        t_stat[row_index, (row_index * 1171) % TRAITS] = 10.
        expected_retained = ROWS
    beta = torch.zeros_like(t_stat)
    status = torch.zeros(ROWS, device=device, dtype=torch.uint8)
    df = torch.full((ROWS,), 38., device=device, dtype=torch.float32)
    critical = torch.ones(100, device=device, dtype=torch.float32)
    geometries = {'legacy': legacy_shape(ROWS, TRAITS, MAX_CELLS),
                  'optimized': device_selection_shape(ROWS, TRAITS, MAX_CELLS)}
    if geometries['optimized'][2] >= geometries['legacy'][2]:
        raise AssertionError('Control shape does not reduce nonzero calls')
    original = reduction.device_selection_shape
    rows = []
    try:
        for name in ('legacy', 'optimized'):
            reduction.device_selection_shape = (
                legacy_shape if name == 'legacy' else original)
            parts = list(reduction.device_significant_pairs(
                beta, t_stat, status, df, critical, max_cells=MAX_CELLS))
            torch.cuda.synchronize(device)
            if (len(parts) != geometries[name][2] or
                    sum(part[2].size for part in parts) != expected_retained):
                raise AssertionError('Selector warmup differs from source control')
        for index, name in enumerate(['legacy', 'optimized', 'optimized', 'legacy'] * 2):
            reduction.device_selection_shape = (
                legacy_shape if name == 'legacy' else original)
            torch.cuda.synchronize(device)
            began = time.perf_counter()
            cpu = time.process_time()
            parts = list(reduction.device_significant_pairs(
                beta, t_stat, status, df, critical, max_cells=MAX_CELLS))
            torch.cuda.synchronize(device)
            elapsed = time.perf_counter() - began
            cpu_seconds = time.process_time() - cpu
            retained = sum(part[2].size for part in parts)
            if len(parts) != geometries[name][2] or retained != expected_retained:
                raise AssertionError('Selector output differs from source control')
            rows.append(dict(index=index, geometry=name, wall_seconds=elapsed,
                             process_cpu_seconds=cpu_seconds, parts=len(parts),
                             retained_pairs=retained,
                             nonempty_parts=sum(bool(part[2].size) for part in parts)))
    finally:
        reduction.device_selection_shape = original
    medians = {name: statistics.median(row['wall_seconds'] for row in rows
                                      if row['geometry'] == name)
               for name in geometries}
    output = dict(shape=dict(rows=ROWS, traits=TRAITS, max_cells=MAX_CELLS),
                  scenario=args.scenario, expected_retained_pairs=expected_retained,
                  geometries={name:dict(width=shape[0], height=shape[1], blocks=shape[2])
                              for name, shape in geometries.items()},
                  device=torch.cuda.get_device_name(device), torch_version=torch.__version__,
                  cuda_visible_devices=os.getenv('CUDA_VISIBLE_DEVICES'),
                  selector_sha256=hashlib.sha256(Path(reduction.__file__).read_bytes()).hexdigest(),
                  observations=rows, median_wall_seconds=medians,
                  scope='Alternating single-process selector calls with identical resident tensors and no writer or GWAS pipeline. Shared GPU/CPU load is uncontrolled; these times are not a throughput guarantee or JIT price.')
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=False)
    with path.open('x') as stream:
        json.dump(output, stream, indent=2)
    print(json.dumps(dict(path=str(path), geometries=output['geometries'],
                          median_wall_seconds=medians), sort_keys=True))


if __name__ == '__main__':
    main()
