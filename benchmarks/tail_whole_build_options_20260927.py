"""Inductor options for the whole-tail graph: GPU ns per cell against the block build.

The one-graph export (tails._tail_whole) needs one host call per strip but
ran 0.82-0.85 ns per cell against 0.36-0.42 for the 12-call block build
(benchmarks/tail_graph_forms_20260927.py), which forces the (c, d, h) state
to memory every four iterations. Builds the whole graph under several
Inductor config patches into a scratch directory and times each on a
(1024, 8192) and a (4096, 512) FP64 input.

    python benchmarks/tail_whole_build_options_20260927.py cuda:4 /data/zxie3/tail_build_scratch
"""
import json
from pathlib import Path
import sys
import time

import torch
from torch.export import Dim

from torchgwas import tails

# First sweep (H100, ns per cell at (1024, 8192) / (4096, 512)): default
# 0.818 / 0.855, realize_opcount_threshold 8 0.557 / 0.607, max_fusion_size 16
# 0.736 / 0.776, both 1.442 / 1.522; the block build 0.358 / 0.424.
VARIANTS = {
    'opcount4': {'realize_opcount_threshold': 4},
    'opcount2': {'realize_opcount_threshold': 2},
    'opcount4_reads2': {'realize_opcount_threshold': 4, 'realize_reads_threshold': 2},
    'opcount8_reads2': {'realize_opcount_threshold': 8, 'realize_reads_threshold': 2},
}


def main():
    device = torch.device(sys.argv[1])
    torch.cuda.set_device(device)
    scratch = Path(sys.argv[2]); scratch.mkdir(parents=True, exist_ok=True)
    rows, traits = Dim('rows', min=2, max=1 << 24), Dim('traits', min=2, max=1 << 24)

    class Whole(torch.nn.Module):
        def forward(self, t, df):
            return tails._tail_whole(t, df)

    example = torch.linspace(0.0, 40.0, 64 * 33, dtype=torch.float64, device=device).reshape(64, 33)
    for name, options in VARIANTS.items():
        path = scratch/f'{name}.so'
        started = time.perf_counter()
        try:
            torch._export.aot_compile(Whole(), (example, torch.full_like(example, 30.0)),
                                      dynamic_shapes=({0: rows, 1: traits}, {0: rows, 1: traits}),
                                      options={'aot_inductor.output_path': str(path), **options})
        except Exception as error:  # noqa: BLE001 - report and continue
            print(json.dumps(dict(variant=name, error=repr(error)[:300])), flush=True)
            continue
        build = time.perf_counter() - started
        call = torch._export.aot_load(str(path), device=str(device))
        row = dict(variant=name, options=options, build_seconds=round(build, 1))
        for shape in ((1024, 8192), (4096, 512)):
            t = torch.randn(shape, dtype=torch.float64, device=device) * 4
            df = torch.full(shape, 22_238.0, dtype=torch.float64, device=device)
            got = call(t, df)
            got = got[0] if isinstance(got, (list, tuple)) else got
            want = tails._evaluate(tails._eager_stages(), t, df, tails._starts(device))
            torch.cuda.synchronize()
            started = time.perf_counter()
            for _ in range(20):
                call(t, df)
            torch.cuda.synchronize()
            seconds = (time.perf_counter() - started) / 20
            row[str(shape)] = dict(ns_per_cell=round(1e9 * seconds / (shape[0] * shape[1]), 3),
                                   max_rel=float(((got - want).abs() / want.abs().clamp_min(1e-300)).max()))
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
