# Pinned-transfer diagnostic for JIT multi-GPU pricing

The first-chunk candidate factory requires explicit shared H2D and D2H
capacities, but existing detailed profiles do not contain an independent
measurement producer for those fields. The bounded
[`shared_transfer_capacity_probe_20260923.py`](../benchmarks/shared_transfer_capacity_probe_20260923.py)
records raw pinned-host copy service on each selected GPU and on the pair
concurrently. It checks the PyTorch-to-physical GPU UUID mapping, records
`nvidia-smi` load snapshots, process CPU affinity, wall and CUDA-event spans,
and publishes its JSON with exclusive creation. It performs no association
scan and does not run in the cold or post-output JIT planner.

On lab-a100, GPUs 0 and 7 are A100-SXM4-80GB cards attached to CPU NUMA nodes
0 and 1 respectively (`nvidia-smi topo -m`). Each observation copied a 16 MiB
pinned buffer eight times per GPU. There were two warmups and eight measured
repeats for each single-GPU and concurrent group. The selected GPUs showed at
most 20% utilization before each run. The raw reports are
[`a100_0_7_v1.json`](../results/shared_transfer_diagnostic_20260923/a100_0_7_v1.json),
[`a100_0_7_cpu_node0_v1.json`](../results/shared_transfer_diagnostic_20260923/a100_0_7_cpu_node0_v1.json),
[`a100_0_7_cpu_node1_v1.json`](../results/shared_transfer_diagnostic_20260923/a100_0_7_cpu_node1_v1.json),
and the final UUID-checked
[`a100_0_7_cpu_node1_uuid_checked_v2.json`](../results/shared_transfer_diagnostic_20260923/a100_0_7_cpu_node1_uuid_checked_v2.json).
The final report's script SHA-256 is
`519efdb88074a36f8f732bdc56fb336b3539ccfbb8e7a2fde51de635ee31f032`,
matching the local source after the UUID check. PyTorch omits the `GPU-`
prefix returned by `nvidia-smi`; the check normalizes that prefix. Direct
UUID inspection confirmed that CUDA indices 0 and 7 matched physical GPUs 0
and 7 in these runs.

Median aggregate decimal GB/s:

| CPU binding | Direction | GPU 0 alone | GPU 7 alone | GPUs 0+7 together |
|---|---|---:|---:|---:|
| CPU 0 (node 0) | H2D | 6.61 | 9.46 | 13.21 |
| CPU 0 (node 0) | D2H | 6.47 | 5.26 | 9.01 |
| CPU 24 (node 1), UUID-checked | H2D | 6.62 | 22.08 | 13.10 |
| CPU 24 (node 1), UUID-checked | D2H | 6.50 | 13.85 | 12.99 |

The earlier CPU-24 run measured GPU 7 at 24.59 GB/s H2D and 19.07 GB/s D2H,
so even the nominally matched binding did not yield one stable D2H price.
Changing CPU binding strongly changed GPU 7 service; GPU 0 remained near
6.6 GB/s in these observations. The script did not verify the actual NUMA
residency of every pinned allocation, so CPU binding is a controlled input,
not proof of its memory placement. Several other GPUs on the shared server
were busy, and the selected GPUs already held substantial device memory.
The two load snapshots do not exclude transient competing transfers.

These results demonstrate that a GPU-agnostic or device-count-scaled
transfer rate is unsafe. They are **diagnostic observations, not calibrated
capacity ceilings** and must not populate `shared_transfer_capacities` or
authorize a JIT chunk/GPU switch. The detailed profile already binds GPU UUID,
PCI placement, CPU affinity and NUMA policy; future transfer artifacts also
need independently checked pinned allocation placement and a controlled
load condition. A producer should measure per-device and simultaneous-group
service for the active device set, preserve raw repeats and provenance, and
validate the resulting scenario against held-out output-inclusive native
PGEN jobs. The finite continuation must still include issued GPU work,
producer queues, active output, final writeback/fsync and publication before
any first-chunk selection can be applied.
