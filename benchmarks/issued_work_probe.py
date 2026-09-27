"""Read-only, post-output issued-work cost probe on a native PGEN index."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import resource
import sys
import time

from torchgwas.layout_frontier import unissued_frontier
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.productive_boundary import ProductiveBoundaryProgress
from torchgwas.productive_issued_work import productive_issued_work
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.sumstats import DenseWriteProgress


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--chunk', type=int, default=128)
    parser.add_argument('--issued-chunks', type=int, default=64)
    parser.add_argument('--traits', type=int, default=128)
    args = parser.parse_args()
    if args.chunk < 1 or not 2 <= args.issued_chunks <= 65 or args.traits < 1:
        parser.error('Positive chunk/traits and 2..65 issued chunks required')
    started = time.perf_counter()
    cpu_started = time.process_time()
    header = PgenHeaderWork(args.input)
    header_wall = time.perf_counter() - started
    header_cpu = time.process_time() - cpu_started
    markers = header._header.variant_ct
    if args.chunk * args.issued_chunks >= markers:
        parser.error('Probe needs an unissued source suffix')
    partitions = [dict(id='original', device='cuda:0',
                       variant_range=[0, markers], trait_range=[0, args.traits])]
    run = ProductiveTuningRun(partitions,
        chunk_sizes=[args.chunk, args.chunk * 2], initial=args.chunk)
    boundary = ProductiveBoundaryProgress(partitions, None)
    reserve = run.for_partition('original')
    for index in range(args.issued_chunks):
        first = index * args.chunk
        if reserve(first, markers, run.capacity) != args.chunk:
            raise RuntimeError('Synthetic issue frontier changed')
    event = DenseWriteProgress(0, args.chunk, (0, args.traits),
        args.chunk * args.traits * 8, args.chunk,
        time.perf_counter(), 'read-only-probe', 'cuda:0')
    run.output_written(event)
    boundary.observe(event)
    snapshot = run.snapshot()
    frontier = unissued_frontier(snapshot, source_identity=header.input_identity,
        reduction=None, total_traits=args.traits,
        job_variant_range=[0, markers])
    held = boundary.bind(snapshot)
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    started = time.perf_counter()
    cpu_started = time.process_time()
    report = productive_issued_work(
        held, frontier, header, covariate_rank=3,
        max_pending_chunks=args.issued_chunks - 1,
        max_source_records=args.chunk * (args.issued_chunks - 1))
    ledger_wall = time.perf_counter() - started
    ledger_cpu = time.process_time() - cpu_started
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    run.finish(successful=False)
    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    code_sha256 = dict(
        probe=digest(__file__),
        issued_work=digest(inspect.getfile(productive_issued_work)),
        boundary=digest(inspect.getfile(ProductiveBoundaryProgress)),
        source_header=digest(inspect.getfile(PgenHeaderWork)))
    result = dict(kind='torchgwas.issued_work_probe.v1',
        input_path=str(Path(args.input).resolve()),
        input_size_bytes=Path(args.input).stat().st_size,
        input_identity=header.input_identity,
        python_executable=sys.executable,
        code_sha256=code_sha256,
        samples=header._header.sample_ct, variants=markers,
        traits=args.traits, chunk_markers=args.chunk,
        issued_chunks=args.issued_chunks,
        header_wall_seconds=header_wall, header_cpu_seconds=header_cpu,
        ledger_wall_seconds=ledger_wall, ledger_cpu_seconds=ledger_cpu,
        peak_rss_before_ledger_kib=rss_before,
        peak_rss_after_ledger_kib=rss_after,
        pending_chunks=report['pending_chunks'],
        pending_source_records=report['pending_source_records'],
        total_work=report['total_work'],
        source_units=report['total_source_units'],
        scope='Read-only index probe, synthetic first-output checkpoint. Header parse and issued-work ledger are timed separately. No genotype payload, GPU or writer execution; not a scan throughput or JIT-switch result.')
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
