# Continuing the calculator from an in-flight pipeline

Starting a new graph for the remaining variants repeats setup, forgets occupied
slots and pending output, and can overstate the benefit of switching. The
existing analytical scheduler now supports checkpoints and continuation.
This is model state, not an observation of the actual executor.

## State and compatible future work

`ExecutionGraph.checkpoint(at_seconds)` advances the Python event solver to a
model decision time and retains completed/start timestamps, nominal remaining
service, FIFO order, free tokens, conditional wakeups and resource waits.
`resume(checkpoint)` continues from that state and returns `remaining_seconds`
in addition to the absolute modeled completion time. A checkpoint can also be
advanced with `checkpoint(later, initial=checkpoint)`.

Remaining service is work at the node's nominal capacity. For two ten-second
CPU tasks sharing one core, six elapsed seconds leave seven nominal seconds
per task. If two cores then become available, completion takes seven more
seconds. Subtracting six from each task or starting both tasks again gives the
wrong answer. Running demands and acquired tokens remain active across the
decision boundary.

Started node durations, demands, dependencies, token actions and queue contracts
must match. Token capacities and resource waits already entered are preserved.
Future unstarted nodes and resource availability may change. Queued results
cannot lose their consumer or change arrival order. A completed checkpoint
performs no additional work. Inputs are copied; resuming does not mutate the
earlier state. JSON round trips preserve the contract.

The generic DAG cannot identify physical source ranges or allocations.
`adaptive_candidate.future_chunk_candidate` constructs one source-aware
counterfactual from a candidate, a fine census, a permitted size set and an
explicit issued-chunk count for every tile/shard. It keeps every issued range
and its complete encoded/LD replay work, changes only the future ranges, and
keeps devices, phenotype partition, readers, depth, output and ring capacity.
It validates the source identity and enforces finite census/graph expansion
budgets. It does not enumerate alternatives or read the genotype file.

Issued counts must come from the executor's issue frontier, including reads
whose results have not arrived. Completed-output counts are insufficient.
Each candidate still requires shared-allocation admission, valid independent
component evidence and live resource checks. JAGWAS retains its complete panel
and factor on every active GPU. No output mode or significance threshold is
treated as an optimization action.

## Verification and timing scope

Tests resume at interior and exact event times, including repeated checkpoints,
shared resources, blocked FIFO consumers, occupied buffer slots, held Python
critical sections, conditional wakeups, running services and completed graphs.
They compare all start/end times and wait resource work with uninterrupted
schedules. Invalid histories and changes to started contracts are refused.
Existing native C++ and Python event solver comparisons remain in the suite;
continuation currently uses the Python solver.

An actual-shape two-GPU JAGWAS graph is checkpointed at its first modeled durable
part completion. Alternatives keep the issued source prefix and change future
sizes among 128, 256 and 512, retaining a fixed 512-row ring. Every started and
completed prefix timestamp remains unchanged. The unchanged choice reproduces
the full original completion time. Tests also cover successive future revisions,
source tampering, admission-size mismatch, exceeded budgets and record aliasing.

The separate `direct_candidate_continuation_20260922.py` audit uses the retained
native PGEN input from the earlier real two-GPU execution. It reports calculator
construction and continuation CPU/wall time separately from forecast gains.
Component prices, setup services and the decision checkpoint are model controls;
neither that checkpoint nor those forecasts are measured live pipeline progress.

A100 job `20260922-001037-814423` passed all 277 tests with no failures, errors
or skips in 67.22 seconds. The optional CPU extension was built in the preceding
validation job, so the native/Python comparisons ran rather than being skipped.
The final test artifacts are in `results/execution_checkpoint_v3_20260922`.

The same job produced `results/candidate_continuation_v1_20260922`. At the modeled
decision, 857 nodes had completed, three still had service remaining, and four
plus three chunks had been issued on the two GPUs. Their source ranges and paid
timestamps were preserved in all alternatives. The unchanged 512-row choice
reproduced the original schedule, with 0.037858 seconds of modeled work left.

| Future size | Source chunks by GPU | Modeled remaining seconds | Calculator build wall seconds | Resume wall seconds |
|---|---|---:|---:|---:|
| 128 | 8, 4 | 0.051308 | 0.279429 | 0.090342 |
| 256 | 6, 4 | 0.043301 | 0.941320 | 0.070943 |
| 512 | 5, 4 | 0.037858 | 0.221097 | 0.056458 |

These measured calculator costs are substantial relative to the remaining time
in this synthetic-service case. They motivate a cheap cost gate and reuse of
already-constructed work before evaluating another alternative. They do not
establish a production planning-cost bound, total JIT overhead or a GWAS speedup.
The report binds the earlier execution report, five unchanged input/metadata
hashes and the current calculator sources. After pulling, all 119 package source
hashes, the harness and its three helper hashes matched the local checkout.

## Connection to JIT tuning

These calls permit a deferred planner to score unpaid work while carrying its
modeled pipeline state. They do not yet observe live queue/kernel progress,
provide a Bayesian posterior, impose a hard planning CPU deadline or choose an
execution action. The public API still needs its nonblocking startup and bounded
incremental decision loop. Full graph construction also needs to be amortized
or made incremental; moving it to a thread does not make its cost disappear.

Only reusable structural/component evidence may transfer across jobs under the
cache validity rules. A pipeline checkpoint is job-specific progress and must
not be reused as a later job's execution state or relabeled as fresh capacity
calibration. Future Bayesian updates must preserve contributing observation
times and remain separate immutable records. The requirements for time to first
useful output, matched end-to-end overhead and cost-aware exploration remain in
initial_chunk_calibration_20260922.md.
