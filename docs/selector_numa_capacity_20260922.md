# Independent selector service under NUMA policies

Eighteen fixed controls do not reproduce the large loaded-scan selector penalty.
With physically written resident inputs, the frozen selector used about
6.4-7.0 ms of CPU per call with one worker and 7.2-7.5 ms with two. Changing
private NUMA policy changed fault behavior but did not produce a consistent
throughput improvement across the three repeats. No new rate or multiplier is
installed, and the production allocation policy is unchanged.

## Experiment

H100 job `20260922-155341-1105236` completed three randomized complete repeats
of default, node-1 bind, and node-1 bind with balancing, crossed with one or two
workers. Each observation used a fresh process, CPUs 12-19 on node 1, the frozen
`host_significance.py`, NumPy's huge-page advice disabled, two warmup calls and
1,024 measured calls per worker. GPU 1 was occupied by another job when inspected,
so this was a CPU component experiment, not a new matched GPU scan.

Policy was set before importing NumPy or creating workers. Every worker checked
inheritance, and the original calling-thread policy was restored afterward.
No global kernel setting was changed. Node-local free memory had to exceed
1 GiB before starting each child. Peak process RSS was about 116-196 MiB.

Each worker allocated private beta and t arrays and filled every element with
nonzero values before measurement. Shape was 1,024 x 8,193 for worker zero and
1,024 x 8,192 for worker one, matching the frozen two-tile layout. Empty output
was checked on every call. A small selected case checked the count and
coordinates before timing; unchanged inputs were checked afterward. These are
component controls, not new association-reference validation.

The experiment made 27,648 measured calls across 27 worker lifetimes. All 918
sampled-page queries (17 pages before and after each worker) returned node 1.
This verifies sampled input placement, not every input, scratch or output page.

## Results

| Policy | Workers | Worker CPU ms/call, range | Aggregate calls/s, median | Aggregate calls/s, range |
| --- | ---: | ---: | ---: | ---: |
| Default | 1 | 6.507-6.966 | 152.1 | 137.0-152.3 |
| Bind node 1 | 1 | 6.487-6.538 | 152.4 | 150.5-154.2 |
| Bind node 1 with balancing | 1 | 6.432-6.778 | 154.2 | 128.7-155.5 |
| Default | 2 | 7.198-7.468 | 262.2 | 261.8-275.9 |
| Bind node 1 | 2 | 7.194-7.351 | 271.7 | 270.4-271.8 |
| Bind node 1 with balancing | 2 | 7.261-7.515 | 263.1 | 229.1-271.0 |

Throughput uses the earliest worker loop start through the latest loop end.
It excludes imports, input allocation, warmup, final input verification and
thread joins. Loop CPU includes Python dispatch, empty-result checks and fixed
sampling overhead. Whole-process CPU and join spans remain separately recorded.
The apparent difference between policy medians is not a general speedup claim;
the default/bind direction reverses in the third two-worker repeat.

Strict bind workers incurred at most two sampled minor faults after their first
sample across the remaining 63 sampled calls. Default and balanced-bind workers
sometimes accumulated thousands, yet their mean CPU costs remained close. Raw
NUMA scheduler endpoints also changed with policy. Those fields can already be
nonzero before this work begins; their absolute values are not migration work
performed by this experiment. They are not used as per-stage capacity prices.

The Linux [memory-policy interface](https://man7.org/linux/man-pages/man2/set_mempolicy.2.html)
allows the private binding/balancing comparison, and the
[kernel documentation](https://docs.kernel.org/admin-guide/sysctl/kernel.html#numa-balancing)
describes why balancing can generate faults. The measurements here do not
attribute each fault to a particular kernel mechanism.

## Consequence for the calculator

These sustained controls are consistent with the original roughly 7 ms
independent whole-selector estimate. They do not explain the larger loaded
CPU costs in the [frozen scan diagnosis](frozen_runtime_diagnosis_20260922.md).
The next investigation must retain concurrent pipeline work and buffer
lifetimes; a uniform NUMA penalty inferred from loaded spans would not follow
from this evidence. The existing [NUMA context binding](numa_context_binding_20260922.md)
remains necessary for compatibility even though binding did not establish a
throughput improvement here.

The controls reuse resident CPU-written inputs and emit no selected pairs.
They exclude DMA, full-pipeline memory pressure, decode work, completion-event
waiters and cross-thread allocation ownership. External CPU and memory load was
recorded, not controlled. Three repeats of one shape do not qualify general
capacity transfer, nonempty output or the calculator's large-job predictions.

## Artifacts

The separate diagnostic project contains `selector_numa_capacity.py`,
`summarize_selector_numa_capacity.py`, and
`results/selector_numa_capacity_20260922/` with the frozen schedule, all 18 raw
records, logs, completion record and derived summary. The preceding three-case
smoke control is retained under `results/selector_numa_capacity_smoke_20260922/`.

The executed harness hash is
`aa09e80a7813848fae14c2ced2b6a37e91cb599cd04d122aadd676f556b4bc0b`;
the summarizer hash is
`07edbc28a58e234c2a576c2875510e235d50b0babb437cd9b63ca02ebbf6d5a6`.
Both remote script hashes and the frozen selector hash matched local files after
pulling. The summarizer verified all run identities, policy inheritance, NumPy
binary identity, fixed call counts and the complete schedule.
