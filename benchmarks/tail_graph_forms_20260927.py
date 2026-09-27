"""Device tail: one exported graph against the 12-call block build, per call and per cell.

Loads both AOTInductor builds for one GPU (the block build from its own
.build-libs directory, given on the command line) and times, on one thread,
calls on shapes from a min-p winner column to a K = 8,192 strip:
milliseconds per call (host-bound when small) and ns per cell (GPU-bound when
large). Also times the whole graph through the low-level runner, which skips
aot_load's Python wrapper.

    python benchmarks/tail_graph_forms_20260927.py cuda:4 .build-libs/device_tail_<block key>
"""
import json
from pathlib import Path
import sys
import time

import torch

from torchgwas import tails


def timed(call, t, df, repeats):
    call(t, df)
    torch.cuda.synchronize()
    started = time.perf_counter()
    for _ in range(repeats):
        call(t, df)
    torch.cuda.synchronize()
    return (time.perf_counter() - started) / repeats


def main():
    device = torch.device(sys.argv[1])
    torch.cuda.set_device(device)
    blocks = Path(sys.argv[2])
    pro, blk, epi = [torch._export.aot_load(str(blocks/f'{name}.so'), device=str(device))
                     for name in ('prologue', 'block', 'epilogue')]
    starts = tails._starts(device)

    def block_form(t, df):
        first, second, argument, c, d, h, x, y, reflect = pro(t, df)
        for start in starts:
            c, d, h = blk(first, second, argument, c, d, h, start)
        return epi(df, x, y, reflect, h)

    tails.prepare_device_tail(device)
    whole = tails._STAGES[device][0]
    directory = tails.device_tail_directory(torch.cuda.get_device_capability(device))
    runner = torch._C._aoti.AOTIModelContainerRunnerCuda(str(directory/'whole.so'), 1, str(device))

    def runner_form(t, df):
        return runner.run([t, df])

    for shape in ((2048, 2), (4096, 64), (4096, 512), (1024, 8192), (512, 8192)):
        t = torch.randn(shape, dtype=torch.float64, device=device) * 4
        df = torch.full(shape, 22_238.0, dtype=torch.float64, device=device)
        cells = shape[0] * shape[1]
        repeats = 200 if cells < 1 << 20 else 20
        row = dict(shape=shape)
        for name, call in (('blocks', block_form), ('whole', whole), ('whole_runner', runner_form)):
            seconds = timed(call, t, df, repeats)
            row[name] = dict(ms_per_call=round(1e3 * seconds, 3), ns_per_cell=round(1e9 * seconds / cells, 3))
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
