"""Stage diagnostics for unexpectedly slow large dense generic selection.

Not independent price publication. Tracebacks and per-stage clocks identify
where work is spent; this process may overlap the named running controls.
"""
import faulthandler
import hashlib
import json
import os
from pathlib import Path
import resource
import time
import numpy as np
from torchgwas.host_significance import ceil_float32, fill_predicate_mask


def main():
    root = Path('results/diagnose_host_selector_20260922')
    root.mkdir(parents=True, exist_ok=False)
    records = []
    def stage(name, function):
        print('BEGIN', name, flush=True)
        faults = resource.getrusage(resource.RUSAGE_THREAD)
        wall = time.perf_counter(); cpu = time.thread_time()
        result = function()
        cpu = time.thread_time()-cpu; wall = time.perf_counter()-wall
        after = resource.getrusage(resource.RUSAGE_THREAD)
        record = dict(name=name, wall_seconds=wall, cpu_seconds=cpu,
            minor_faults=after.ru_minflt-faults.ru_minflt, major_faults=after.ru_majflt-faults.ru_majflt)
        records.append(record)
        print(json.dumps(record), flush=True)
        (root/'progress.json').write_text(json.dumps(records, indent=2)+'\n')
        return result
    faulthandler.dump_traceback_later(20, repeat=True)
    shape = (1024, 8193)
    values = stage('create_values', lambda: np.ones(shape, np.float32))
    df = np.full((shape[0],1), 31., np.float32)
    critical = np.zeros((shape[0],1), np.float64)
    limits = stage('critical_round', lambda: np.broadcast_to(ceil_float32(critical), shape))
    keep = stage('mask_allocate', lambda: np.empty(shape, bool))
    stage('predicate', lambda: fill_predicate_mask(values, limits, keep))
    rows = stage('flatnonzero', lambda: np.flatnonzero(keep).astype(np.int64, copy=False))
    columns = stage('columns_allocate', lambda: np.empty_like(rows))
    stage('coordinate_divmod', lambda: np.divmod(rows, shape[1], out=(rows,columns)))
    selected_df = stage('broadcast_df_gather', lambda: np.broadcast_to(df,shape)[rows,columns])
    selected_t = stage('matrix_gather', lambda: values[rows,columns])
    stage('inplace_rebase', lambda: np.add(rows, 123, out=rows))
    assert len(rows) == shape[0]*shape[1] and rows[0] == 123 and rows[-1] == 1146
    assert columns[0] == 0 and columns[-1] == 8192
    assert selected_t[0] == selected_t[-1] == 1. and selected_df[0] == selected_df[-1] == 31.
    faulthandler.cancel_dump_traceback_later()
    (root/'report.json').write_text(json.dumps(dict(records=records, shape=shape,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        numpy=np.__version__, affinity=sorted(os.sched_getaffinity(0)),
        backend=os.getenv('TORCHGWAS_HOST_PREDICATE'), scope=__doc__), indent=2)+'\n')


if __name__=='__main__': main()
