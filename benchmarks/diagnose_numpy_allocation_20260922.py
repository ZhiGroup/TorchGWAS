"""Counterbalanced fresh-process NumPy huge-page-setting diagnostics.

No price publication. Named concurrent controls prevent an isolation claim.
All observations include separate user/system CPU and faults.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

ROOT = Path('results/numpy_allocation_setting_20260922')


def child(label):
    import faulthandler
    import numpy as np
    faulthandler.dump_traceback_later(30, repeat=True)
    records = []
    def stage(name, function):
        print('BEGIN', name, flush=True)
        before = resource.getrusage(resource.RUSAGE_THREAD)
        start = time.perf_counter()
        result = function()
        wall = time.perf_counter()-start
        after = resource.getrusage(resource.RUSAGE_THREAD)
        record = dict(name=name, wall_seconds=wall,
            user_cpu_seconds=after.ru_utime-before.ru_utime,
            system_cpu_seconds=after.ru_stime-before.ru_stime,
            minor_faults=after.ru_minflt-before.ru_minflt,
            major_faults=after.ru_majflt-before.ru_majflt)
        records.append(record); print(json.dumps(record), flush=True)
        return result
    shape = (1024, 8193)
    a = stage('ones', lambda: np.ones(shape, np.float32))
    flat = stage('flat_index', lambda: np.arange(a.size, dtype=np.int64))
    rows = stage('row_index', lambda: flat//shape[1])
    columns = stage('column_index', lambda: flat%shape[1])
    result = stage('matrix_gather', lambda: a[rows,columns])
    copied = stage('copy', lambda: a.copy())
    assert result.size == a.size and copied.shape == a.shape
    hashes = [hashlib.sha256(memoryview(v).cast('B')).hexdigest() for v in (a,result,copied)]
    assert len(set(hashes)) == 1
    faulthandler.cancel_dump_traceback_later()
    report = dict(label=label, records=records, numpy=np.__version__,
        numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
        declared_madvise=os.getenv('NUMPY_MADVISE_HUGEPAGE'), affinity=sorted(os.sched_getaffinity(0)),
        thp_enabled=Path('/sys/kernel/mm/transparent_hugepage/enabled').read_text().strip(),
        thp_defrag=Path('/sys/kernel/mm/transparent_hugepage/defrag').read_text().strip(),
        output_sha256=hashes[0], script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope=__doc__)
    (ROOT/(label+'.json')).write_text(json.dumps(report, indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--child'); args=parser.parse_args()
    if args.child: return child(args.child)
    ROOT.mkdir(parents=True, exist_ok=False)
    schedule=[('disabled_1','0'),('enabled_1','1'),('enabled_2','1'),('disabled_2','0')]
    (ROOT/'schedule.json').write_text(json.dumps(dict(schedule=schedule,
        concurrent_controls=[2295665,2481000], purpose=__doc__), indent=2)+'\n')
    for label,flag in schedule:
        env=dict(os.environ,NUMPY_MADVISE_HUGEPAGE=flag)
        with (ROOT/(label+'.log')).open('x') as log:
            subprocess.run([sys.executable,__file__,'--child',label],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        record=json.loads((ROOT/(label+'.json')).read_text())
        print(json.dumps(dict(label=label,records=record['records'])),flush=True)
    (ROOT/'complete.json').write_text(json.dumps(dict(complete=True,labels=[x[0] for x in schedule]))+'\n')


if __name__=='__main__': main()
