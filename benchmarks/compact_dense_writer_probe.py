"""One-process, metadata-only comparison of fixed-chunk dense writer ledgers."""
import argparse
import json
import resource
import time

from torchgwas.binary_output_work import binary_output_work, compact_binary_output_work


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kind', choices=('compact', 'expanded'), required=True)
    parser.add_argument('--markers', type=int, default=8_086_101)
    parser.add_argument('--traits', type=int, default=128)
    parser.add_argument('--chunk-markers', type=int, default=128)
    args = parser.parse_args()
    settings = dict(markers=args.markers, traits=args.traits,
                    chunk_markers=args.chunk_markers, block_bytes=16 << 20,
                    queue_depth=3, borrow_chunks=False, store_beta=True,
                    fsync=True, writeback_bytes=64 << 20,
                    sync_file_range=True, store_variant_df=True)
    started = time.process_time()
    if args.kind == 'compact':
        work = compact_binary_output_work(**settings)
        counts = dict(payload_bytes=work['payload_bytes'],
                      staging_copy_calls=work['staging_copy_calls'],
                      write_calls_minimum=work['write_calls_minimum'])
    else:
        work = binary_output_work(**settings)
        counts = dict(payload_bytes=work['binary_payload_bytes'],
                      staging_copy_calls=work['staging_copy_calls'],
                      write_calls_minimum=work['binary_write_calls_minimum'])
    print(json.dumps(dict(kind=args.kind, markers=args.markers,
                          traits=args.traits, chunk_markers=args.chunk_markers,
                          process_cpu_seconds=time.process_time()-started,
                          peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          **counts), sort_keys=True))


if __name__ == '__main__':
    main()
