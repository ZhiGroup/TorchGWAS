"""Paired source-selector control: pageable copies versus one owned pinned payload."""
import gc
import hashlib
import inspect
import json
import os
from pathlib import Path
import random
import statistics
import textwrap
import time
from direct_device_significance_primitive_measure_20260921 import telemetry

import numpy as np
import torch
from torchgwas import reduce
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record

ROOT = Path('results/device_significance_pinned_control_20260921')
ROOT.mkdir(parents=True, exist_ok=False)
(ROOT/'harness.py').write_bytes(Path(__file__).read_bytes())
os.sched_setaffinity(0, range(12, 20))
torch.set_num_threads(2); torch.set_num_interop_threads(1)
DEVICE = 'cuda:0'; LIMIT = 512 << 20
torch.cuda.set_device(DEVICE)
torch.cuda.set_per_process_memory_fraction(LIMIT/torch.cuda.get_device_properties(DEVICE).total_memory, DEVICE)
if torch.cuda.mem_get_info(DEVICE)[0] < LIMIT: raise ValueError('Insufficient free memory')
source = source_identity()
original = textwrap.dedent(inspect.getsource(reduce.device_significant_pairs))
old = '''            vi=(ri+first+start).cpu().numpy()
            trait=(ti+left).cpu().numpy()
            selected_beta=beta[first:last,left:right][ri,ti].cpu().numpy()
            selected_t=values[ri,ti].cpu().numpy()
            selected_df=df[ri].cpu().numpy()
            yield (start+first,start+last,vi,trait,selected_beta,selected_t,selected_df)'''
new = '''            count=indices.shape[0]
            payload=torch.empty(28*count,dtype=torch.uint8,device="cpu",pin_memory=True)
            vi=payload[:8*count].view(torch.int64)
            trait=payload[8*count:16*count].view(torch.int64)
            selected_beta=payload[16*count:20*count].view(torch.float32)
            selected_t=payload[20*count:24*count].view(torch.float32)
            selected_df=payload[24*count:].view(torch.float32)
            vi.copy_(ri+first+start,non_blocking=True)
            trait.copy_(ti+left,non_blocking=True)
            selected_beta.copy_(beta[first:last,left:right][ri,ti],non_blocking=True)
            selected_t.copy_(values[ri,ti],non_blocking=True)
            selected_df.copy_(df[ri],non_blocking=True)
            torch.cuda.current_stream(beta.device).synchronize()
            yield (start+first,start+last,vi.numpy(),trait.numpy(),selected_beta.numpy(),selected_t.numpy(),selected_df.numpy())'''
assert original.count(old) == 1
proposal = original.replace(old, new)
namespace = dict(reduce.device_significant_pairs.__globals__)
exec(compile(proposal, 'pinned_payload_candidate.py', 'exec'), namespace)
functions = {'baseline': reduce.device_significant_pairs, 'pinned_payload': namespace['device_significant_pairs']}
(ROOT/'original.py').write_text(original)
(ROOT/'proposal.py').write_text(proposal)
rng = random.Random(9219201)
rows = []
states = []
write_record(ROOT/'protocol.json',dict(repetitions=9,calls_per_observation=4,seed=9219201,device=DEVICE,
    shapes=[[1024,1024],[1024,4096]],occupancies=['empty','one_per_block','sparse','dense'],
    source_sha256=source,pinned_selected_bytes_per_block_upper=28*(1<<20),
    boundary='Source selector iteration, output destruction and final device synchronization after warm-up; no preload shim, no profiler.',
    scope='Frozen reversible implementation control, no published speed or autotune qualification.'))


def arrays(function, tensors):
    blocks = list(function(*tensors, start=17))
    return [np.concatenate([v[i] for v in blocks]) for i in range(2, 7)]


def consume(function, tensors):
    count = 0
    for block in function(*tensors, start=17):
        count += len(block[2])
        del block
    return count


for b, k in [(1024, 1024), (1024, 4096)]:
    for occupancy in ['empty', 'one_per_block', 'sparse', 'dense']:
        beta = torch.arange(b*k, dtype=torch.float32, device=DEVICE).reshape(b, k)
        values = torch.zeros_like(beta)
        if occupancy == 'one_per_block': values.view(-1)[::(1 << 20)] = 2.
        elif occupancy == 'sparse': values.view(-1)[::1024] = 2.
        elif occupancy == 'dense': values.fill_(2.)
        status = torch.zeros(b, dtype=torch.uint8, device=DEVICE)
        df = torch.full((b,), 38., device=DEVICE)
        critical = torch.ones(41, device=DEVICE); critical[0] = torch.inf
        tensors = (beta, values, status, df, critical)
        expected = arrays(functions['baseline'], tensors)
        observed = arrays(functions['pinned_payload'], tensors)
        for left, right in zip(expected, observed): np.testing.assert_array_equal(left, right)
        retained = len(expected[0]); del expected, observed
        held=list(functions['pinned_payload'](*tensors,start=17))
        saved=[[value.copy() for value in block[2:]] for block in held]
        changed=(beta+11,values,status,df,critical)
        for _ in range(3): consume(functions['pinned_payload'],changed)
        torch.cuda.synchronize(DEVICE); gc.collect()
        for block,truth in zip(held,saved):
            for value,expected_value in zip(block[2:],truth):np.testing.assert_array_equal(value,expected_value)
        del held,saved,changed
        peaks = {}
        for mode, function in functions.items():
            for _ in range(3): assert consume(function, tensors) == retained
            torch.cuda.synchronize(DEVICE); gc.collect()
            base = torch.cuda.memory_allocated(DEVICE); torch.cuda.reset_peak_memory_stats(DEVICE)
            assert consume(function, tensors) == retained
            peaks[mode] = torch.cuda.max_memory_allocated(DEVICE)-base
        observations = []
        for repeat in range(9):
            states.append(dict(B=b,K=k,occupancy=occupancy,repeat=repeat,telemetry=telemetry()))
            order = list(functions); rng.shuffle(order)
            for mode in order:
                torch.cuda.synchronize(DEVICE)
                cpu = time.thread_time(); wall = time.perf_counter()
                for _ in range(4): count = consume(functions[mode], tensors)
                torch.cuda.synchronize(DEVICE)
                observations.append(dict(repeat=repeat, mode=mode,
                    wall_seconds=(time.perf_counter()-wall)/4,
                    caller_thread_cpu_seconds=(time.thread_time()-cpu)/4, retained=count))
                assert count == retained
        medians = {mode:statistics.median(v['wall_seconds'] for v in observations if v['mode'] == mode) for mode in functions}
        row = dict(B=b, K=k, occupancy=occupancy, retained=retained, arrays_equal=True,retained_buffers_survive_later_calls=True,
            observations=observations, median_wall_seconds=medians,
            speedup=medians['baseline']/medians['pinned_payload'], peak_extra_allocated_bytes=peaks)
        rows.append(row); write_record(ROOT/f'{b}_{k}_{occupancy}.json', row)
        print(json.dumps({key:row[key] for key in ['B','K','occupancy','retained','median_wall_seconds','speedup','peak_extra_allocated_bytes']}), flush=True)
        del tensors, beta, values, status, df, critical
        gc.collect(); torch.cuda.empty_cache()
assert source == source_identity()
write_record(ROOT/'telemetry.json',states)
write_record(ROOT/'report.json', dict(source_sha256=source, observations=rows,
    context=dict(torch_version=torch.__version__, device=torch.cuda.get_device_properties(DEVICE).name,
        affinity=sorted(os.sched_getaffinity(0)), torch_threads=torch.get_num_threads(),
        torch_interop_threads=torch.get_num_interop_threads(), memory_limit_bytes=LIMIT),
    boundary='Complete source selector iteration and output destruction, final device synchronize; fixed prepared synthetic statistics.',
    scope='Implementation A/B only, not component prices or GWAS ranking qualification. Original and proposed sources are frozen alongside all randomized repetitions.'))
