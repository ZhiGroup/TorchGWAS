"""What the indexed store's publication step costs after the scan.

Times writing and fsyncing variant_ids.npy (the whole variant-ID column, as
the writer stores it) and a small manifest, in a scratch directory on the
target filesystem. Formats: fixed-width unicode (current), fixed-width bytes.

    python benchmarks/publication_cost_20260926.py --pvar input.pvar --dir /data/zxie3/tmp --repeats 3
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np


def ids_from_pvar(path):
    ids = []
    with open(path) as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            ids.append(line.split('\t', 3)[2])
    return ids


def timed_save(target, array):
    started = time.perf_counter()
    with target.open('wb') as handle:
        np.save(handle, array, allow_pickle=False)
        handle.flush()
        written = time.perf_counter()
        os.fsync(handle.fileno())
    finished = time.perf_counter()
    return dict(bytes=target.stat().st_size, write_seconds=written - started, fsync_seconds=finished - written)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pvar', type=Path, required=True)
    parser.add_argument('--dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--out', type=Path, default=None)
    args = parser.parse_args()
    started = time.perf_counter()
    ids = ids_from_pvar(args.pvar)
    parse = time.perf_counter() - started
    started = time.perf_counter()
    unicode_ids = np.asarray(ids, dtype=str)
    unicode_convert = time.perf_counter() - started
    started = time.perf_counter()
    byte_ids = np.char.encode(unicode_ids, 'ascii')
    byte_convert = time.perf_counter() - started
    rows = []
    with open('/proc/pressure/io') as stream:
        pressure = stream.readline().strip()
    for repeat in range(args.repeats):
        scratch = Path(tempfile.mkdtemp(prefix='.publication-probe-', dir=args.dir))
        try:
            for name, array in (('unicode', unicode_ids), ('bytes', byte_ids)):
                row = dict(repeat=repeat, format=name, dtype=str(array.dtype), **timed_save(scratch / f'{name}.npy', array))
                rows.append(row)
                print(json.dumps(row), flush=True)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
    report = dict(variants=len(ids), parse_seconds=parse, unicode_convert_seconds=unicode_convert,
                  bytes_convert_seconds=byte_convert, directory=str(args.dir), io_pressure=pressure, rows=rows)
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}))
    if args.out:
        args.out.write_text(json.dumps(report, indent=1) + '\n')


if __name__ == '__main__':
    main()
