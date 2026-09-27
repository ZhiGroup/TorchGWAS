"""Summarize layout_profile_20260924.py JSON lines: medians per layout and selection backend.

Usage: python layout_profile_summary_20260924.py RESULTS.jsonl [...]
"""
import json
import statistics
import sys
from collections import defaultdict


def main(paths):
    for path in paths:
        groups = defaultdict(list)
        for line in open(path):
            record = json.loads(line)
            groups[record['layout'], record['backend']].append(record)
        print(f'== {path}')
        print(f"{'layout':9s} {'select':6s} {'n':>2s} {'exec_s':>8s} {'api_s':>7s} {'setup_s':>8s} "
              f"{'gpu_s/scan':>10s} {'fetch_s':>8s} {'user_s':>7s} {'sys_s':>6s}  exec runs")
        for (layout, backend), rows in sorted(groups.items()):
            median = lambda values: statistics.median(values) if values else float('nan')
            execs = [r['executor_seconds'] for r in rows if r['executor_seconds'] is not None]
            setup = median([max(s['setup_seconds'] for s in r['scans']) for r in rows])
            gpu = median([statistics.mean(s['gpu_compute_milliseconds'] for s in r['scans'])/1e3 for r in rows])
            fetch = median([max(s['fetch_seconds'] for s in r['scans']) for r in rows])
            print(f"{layout:9s} {backend:6s} {len(rows):2d} {median(execs):8.1f} {median([r['api_seconds'] for r in rows]):7.1f} "
                  f"{setup:8.2f} {gpu:10.2f} {fetch:8.2f} {median([r['process_cpu'][0] for r in rows]):7.1f} "
                  f"{median([r['process_cpu'][1] for r in rows]):6.1f}  {[round(x, 1) for x in execs]}")


if __name__ == '__main__':
    main(sys.argv[1:])
