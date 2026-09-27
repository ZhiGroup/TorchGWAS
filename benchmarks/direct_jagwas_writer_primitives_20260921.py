"""Fixed NumPy operations and two-array NPZ costs; no association inputs."""
import argparse
import hashlib
import io
import json
import os
import platform
import statistics
import time
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='1',NUMPY_MADVISE_HUGEPAGE='0')
import numpy as np
from torchgwas.detailed_calibration import source_identity
from torchgwas.reduced_output_work import jagwas_indexed_part_work


class CountingFile:
    def __init__(self):self.position=self.size=self.written=0
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        if whence not in (0,1,2):raise ValueError('Unknown seek origin')
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        return self.position
    def write(self,value):
        count=len(value);self.position+=count;self.written+=count;self.size=max(self.size,self.position)
        return count
    def flush(self):pass
    def read(self,*args):raise io.UnsupportedOperation('write-only counting sink')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    args=parser.parse_args();root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,range(12,20))
    observations=[];started_load=list(os.getloadavg())
    def measure(name,units,function,pattern,batch=1):
        for _ in range(3):value=function()
        rows=[]
        for repeat in range(9):
            begin=time.perf_counter();cpu=time.thread_time()
            for _ in range(batch):value=function()
            rows.append(dict(repeat=repeat,cpu_seconds=(time.thread_time()-cpu)/batch,
                wall_seconds=(time.perf_counter()-begin)/batch))
        record=dict(primitive=name,units=units,pattern=pattern,batch=batch,observations=rows)
        observations.append(record)
        (root/'observations.json').write_text(json.dumps(observations,indent=2))
        return statistics.median(row['cpu_seconds'] for row in rows)
    extent=1<<20
    single=np.linspace(1.,10.,extent,dtype=np.float32)
    double=single.astype(np.float64)
    full=np.ones(extent,bool);empty=np.zeros(extent,bool);sparse=empty.copy();sparse[::64]=True
    indices=np.flatnonzero(full);sparse_indices=np.flatnonzero(sparse);none=indices[:0]
    # CPU-owned contiguous input memory is a declared context, not evidence
    # about buffers freshly written by DMA or shared LLC/NUMA behavior.
    banks={
        'fp32_to_fp64_view':[(0,lambda:np.asarray(single[:0],dtype=np.float64).reshape(-1),'empty'),
            (extent,lambda:np.asarray(single,dtype=np.float64).reshape(-1),'contiguous')],
        'finite_fp64':[(0,lambda:np.isfinite(double[:0]),'empty'),
            (extent,lambda:np.isfinite(double),'finite')],
        'flatnonzero_empty':[(0,lambda:np.flatnonzero(empty[:0]),'empty'),
            (extent,lambda:np.flatnonzero(empty),'all_false')],
        'flatnonzero_nonempty':[(0,lambda:np.flatnonzero(empty[:0]),'empty'),
            (extent,lambda:np.flatnonzero(full),'all_true'),
            (extent,lambda:np.flatnonzero(sparse),'stride64')],
        'index_add':[(0,lambda:8192+none,'empty'),(extent,lambda:8192+indices,'dense')],
        'fp64_gather':[(0,lambda:double[none],'empty'),
            (extent,lambda:double[indices],'dense'),
            (len(sparse_indices),lambda:double[sparse_indices],'stride64')],
    }
    prices={}
    for name,controls in banks.items():
        fixed=measure(name,0,controls[0][1],controls[0][2],batch=64)
        rates=[]
        for units,function,pattern in controls[1:]:
            rates.append((measure(name,units,function,pattern)-fixed)/units)
        prices[name]=dict(call_cpu_seconds=fixed,unit_cpu_seconds=max(0.,max(rates)))
    controls=[];archives=[]
    for count in [0,65536,extent]:
        arrays=dict(variant_index=np.arange(count,dtype=np.int64),chi2=np.ones(count,np.float64))
        amount=sum(value.nbytes for value in arrays.values())
        def serialize(arrays=arrays):
            sink=CountingFile();np.savez(sink,**arrays);return sink
        sink=serialize()
        if count:
            work=jagwas_indexed_part_work(count)
            assert sink.size==work['file_bytes']
            assert sink.written==work['file_bytes']+sum(row['local_header_bytes'] for row in work['arrays'])
        archives.append(dict(rows=count,payload_bytes=amount,file_bytes=sink.size,submitted_bytes=sink.written))
        controls.append((amount,measure('npz_jagwas',amount,serialize,'discarding_seekable_file',batch=8 if not count else 1)))
    fixed=controls[0][1]
    archive=dict(field_schema=[['variant_index','<i8'],['chi2','<f8']],call_cpu_seconds=fixed,
        byte_cpu_seconds=max(0.,max((value-fixed)/amount for amount,value in controls[1:])))
    report=dict(prices=prices,archive=archive,archive_work=archives,observations=observations,
        source_sha256=source_identity(),benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        context=dict(host=os.uname().nodename,python=platform.python_version(),numpy=np.__version__,
            affinity=sorted(os.sched_getaffinity(0)),numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
            load_before=started_load,load_after=list(os.getloadavg()),input_memory='CPU-owned contiguous arrays'),
        prediction_complete=False,
        scope='Independent fixed empty/1M-element NumPy controls and two-array seekable NPZ serialization. '
            'Median repeat batch means; maximum of the declared bulk patterns is a supplied scenario, not a bound. '
            'The sink counts final extent and header rewrites but excludes OS writes, fsync and storage. '
            'GIL share, fresh-DMA/cache/NUMA transfer, setup and full runtime accuracy are not established.')
    (root/'prices.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(dict(prices=prices,archive=archive,archive_work=archives)),flush=True)


if __name__=='__main__':main()
