"""PGEN packed decode alone: variants per second against decoder threads.

At K = 512 on several GPUs the real-data scan is now bounded by host decode
(docs/autotune_design_20260924.md): 16 decoder threads did no better than 8.
Here N threads each fill disjoint 4,096-variant chunks through the source's
native_reader_session (the scan's own fill, packed two-bit rows) for
`--seconds`, from a warm page cache, and the aggregate rate is reported. Linear
scaling means decode is CPU-bound per thread; a plateau means a shared limit
(memory bandwidth, the page cache, or cores the host leaves free).

    python benchmarks/decode_thread_scaling_20260928.py --pgen /data/zxie3/torchgwas_bench/full_scale_k512_20260926/input.pgen
"""
import argparse
import json
import os
import threading
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pgen', required=True)
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 2, 4, 8, 16, 24])
    parser.add_argument('--seconds', type=float, default=6.0)
    parser.add_argument('--chunk', type=int, default=4096)
    args = parser.parse_args()
    os.environ.setdefault('TORCHGWAS_PGEN_BACKEND', 'native')
    os.environ.setdefault('TORCHGWAS_STATS_BACKEND', 'triton')  # packed rows
    from torchgwas.pgen import PgenGenotype
    source = PgenGenotype(args.pgen, mode='hardcall', reader_workers=max(args.threads),
                          metadata_cache_dir=os.path.join(os.path.dirname(args.pgen), 'metadata_cache'))
    variants = source.shape[1]
    chunks = variants // args.chunk
    for count in args.threads:
        done = [0] * count
        stop = threading.Event()
        with source.native_reader_session() as read_into:
            errors = []

            def worker(index):
                # 64-byte aligned, as the scan's pinned rows are.
                shape = (args.chunk, source.native_row_width)
                size = int(np.prod(shape)) * np.dtype(source.native_transfer_dtype).itemsize
                raw = np.empty(size + 64, dtype=np.uint8)
                offset = (-raw.ctypes.data) % 64
                out = raw[offset:offset + size].view(source.native_transfer_dtype).reshape(shape)
                position = index
                try:
                    run(index, out, position)
                except Exception as error:  # noqa: BLE001 - reported below
                    errors.append(repr(error))

            def run(index, out, position):
                while not stop.is_set():
                    start = (position % chunks) * args.chunk
                    read_into(start, start + args.chunk, out)
                    done[index] += 1
                    position += count
            threads = [threading.Thread(target=worker, args=(i,)) for i in range(count)]
            began = time.perf_counter()
            for thread in threads:
                thread.start()
            time.sleep(args.seconds)
            stop.set()
            for thread in threads:
                thread.join()
            seconds = time.perf_counter() - began
        if errors:
            raise SystemExit(errors[0])
        rate = sum(done) * args.chunk / seconds
        print(json.dumps(dict(threads=count, variants_per_second=round(rate), us_per_variant=round(1e6 / rate, 3),
                              per_thread_us=round(1e6 * count / rate, 3), encoding=source.native_encoding)),
              flush=True)


if __name__ == '__main__':
    main()
