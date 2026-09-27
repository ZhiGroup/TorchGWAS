# Calculator proposals from a productive source prefix

`analytical_chunk_proposal` now supplies the three-field proposal consumed by
`ProductiveTuningRun.planning_step`. It preserves the actual reserved source
intervals, including prefetch and shifted LD restarts, instead of assuming that
the old regular-size candidate describes the executed prefix. It evaluates
one alternative on fixed admitted partitions; it does not search a grid before
useful output. JAGWAS remains full-panel on every active GPU. Dense and host
significant-pairs phenotype tiles use the same source-prefix mechanism.

The callback uses the existing graph factories and component prices. No
association duration is fitted. A missing decoder/kernel/writer price still
fails explicitly; the proposer does not invent a fallback value. Memory
admission, parameter-record dependency/age validation and live resource checks
remain caller responsibilities. The cheap cost gate runs before this callback.

## A sufficient condition with unknown in-flight progress

Reserved source intervals do not identify the exact state of CUDA streams,
reader workers, result queues or the writer. Creating a model checkpoint from
the last written row would invent that state. The initial policy instead uses
conditional bounds within the existing fluid-resource graph.

For the current size, all descendants of unissued source-submission nodes are
known unstarted. Mandatory resource work divided by capacity and the
capacity-adjusted dependency path through these nodes provide a remaining-work
lower bound L. Optional handoff delays and active waits are omitted from this
floor; they can only add work or delay. Already-issued chunks do not contribute
mandatory unfinished work merely because they were reserved.

For the alternative, the ceiling U conservatively includes the complete graph's
service, including prefix/setup work that may already have finished. It bounds
any feasible continuation under the declared model without guessing which
prefix services remain. Let s_i be a node's remaining nominal service,
a_ir = demand_ir / capacity_r, and let W bound the sum of normalized active
wait demands. Define the service potential

    P = sum_i s_i * (1 + sum_r a_ir).

The graph's progress rate is

    v_i = min(1, min_r capacity_r / total_active_demand_r).

With n positive-service nodes active and A = sum_i,r a_ir, all active rates
are at least 1/(1+A+W). Thus the potential decreases at a rate at least
(n+A)/(1+A+W), which is at least 1/(1+W). A valid ceiling is therefore

    U = (1+W) * sum_i full_nominal_service_i * (1 + sum_r a_ir).

Optional delays are included at full service in this ceiling. The bound
assumes a feasible reachable execution, so it is not a validator for arbitrary
token/FIFO deadlocks. It is also not a hardware-time or empirical uncertainty
bound: missing costs or stale/wrong component prices remain model limitations.

Wait concurrency is bounded with explicit chains. For each consecutive pair,
the previous wait's completion endpoint must be an ancestor of the next wait's
start endpoint. The code proves that relation in the dependency DAG under a
traversal budget. Each chain contributes at most one wait demand per resource;
unassigned waits get independent slots. Input-release and result-worker chains
continue across sequential phenotype tiles on the same GPU, preventing an
incorrect assumption that every tile's workers are simultaneously active.

The proposal uses baseline_seconds=L and candidate_seconds=U. A change is
eligible only when L-U exceeds actual planning wall time plus the declared
switching cost, and the existing time/CPU/count limits still allow it. The
issue frontier is held during the step, so newly unissued work cannot quietly
start while the proposal is computed. Failing this sufficient condition does
not prove the current size optimal. The ceiling can be too loose to act.

The bounds refer to the supplied service/capacity and output-occupancy scenario.
They are not a posterior over changing shared-machine conditions or future
survivor counts. A robust or Bayesian controller must explicitly represent
those additional uncertainties, as well as measurement and exploration cost.

## Verification and current limitation

Tests compare bounds with actual model continuations on random shared-resource
DAGs and on source-generated graphs, including active waits, optional delays,
one/two devices and sequential phenotype tiles. They reject incorrect wait
chains, incomplete prefixes, changed partition identities, unadmitted sizes
and missing prices. A bridge integration test invokes the source calculator
only after a writer-completion signal. Its deliberately generous test horizon
checks wiring, not planner profitability.

The retained two-GPU native-PGEN execution supplies a real initial prefix for
`direct_jit_proposal_20260922.py`. Its N=2,049, M=4,097, K=512, C=2 input and
source reservations are real. Component prices remain explicit synthetic test
controls, including six named decoder primitives needed by shifted LD starts.
The first audit correctly refused these missing primitives; the audit bank was
then completed explicitly, leaving production refusal unchanged.

The source/prefix integration suite passed 110 tests
(`results/jit_proposal_v2_20260922`). After adding cross-tile wait chains and
bridge/output-mode coverage, all 56 final targeted tests passed
(`results/jit_proposal_v4_20260922`), with no failures, errors or skips. These
suites overlap; their counts should not be added as distinct tests.

The final actual-prefix audit is
`results/jit_proposal_actual_prefix_v3_20260922/report.json`. All 123 package
source hashes, the benchmark and three helper hashes, and the execution-report
hash match the delivered source. Input hashes were rechecked on the server.
The current-size unissued-work floor was 0.14103816 seconds under the synthetic
prices. The final sequential callback measurements were:

| Proposed size | Conditional continuation ceiling | Planner CPU | Planner wall |
| --- | ---: | ---: | ---: |
| 512 | 1.79950 s | 1.35178 s | 1.37439 s |
| 256 | 1.97427 s | 0.38391 s | 0.38941 s |
| 128 (unchanged control) | 2.17014 s | 0.18392 s | 0.19923 s |

Six proven wait chains cover both GPUs, with no unassigned waits. Every gain
floor is negative, so none authorizes a switch. These are individual callback
timings, not repeated benchmark medians. The first callback includes uncached
calculation; later callbacks share a bounded computational cache. The bounds
are deliberately loose, and neither their values nor these synthetic prices
establish the actual remaining hardware runtime.

This is the first calculator-driven proposal path, not a completed public
autotuner. The demonstrated small-job bounds do not justify a change, and
planning is expensive relative to that job's remaining work. Public startup,
cheap cost/horizon forecasts, automatic candidate choice, identified parameter
updates and safe phenotype/GPU reassignment remain open. Existing immutable
evidence and observation-age rules continue to apply.
