# Reusable work and cost gates for incremental planning

The continuation audit found that constructing and solving one candidate could
cost more than the modeled work left in a short job. Two mechanisms now address
that problem without changing the calculator's resource equations: a bounded
session cache and a cheap admission gate for individual planning steps.

## Exact computational reuse

`PlanningWorkCache` is explicitly activated around one planning step and can
survive across later steps in the same job. It stores the existing source-derived
tensor component calculation for an exact N/B/K/C geometry. The binding includes
the compiled statistics/joint kernels, resource and host prices, CPU fraction,
reduction, range validation, default dtype and implementation identities.
Changed bindings miss. Consumers receive independent copies; changing a returned
graph cannot poison another prediction. Both entry count and the estimated size
of retained Python ledgers are bounded. This estimate is not an RSS guarantee.

The expensive statistics trace, analytical cache/traffic calculation and joint
projection can therefore be reused across chunks, GPU shards and future-size
evaluations with identical component inputs. Source-range decoding, reader
initialization, output selection, queues and writer work remain range-specific
and are reconstructed. No empirical whole-GWAS timing table is introduced.

The existing dense, significant-host and JAGWAS bounded planners also receive
this reuse within their ordinary calls. An enclosing incremental session is
preserved rather than replaced by a nested whole-plan call. This improves the
existing opt-in planner; it does not make that blocking search the JIT startup
path.

This cache is an optimization of arithmetic, not a calibration authority. It
does not publish evidence or reset observation times. Immutable empirical
records must still pass dependency and age checks; live resource admission
remains mandatory. In particular, an unchanged numeric price that has expired
does not become eligible because its arithmetic is cached. There is no automatic
reuse of these process-local entries across jobs.

## One step at a time

`IncrementalPlanningBudget` refuses to invoke its evaluator until the caller
signals useful output delivery. Repeating that signal never renews the early
window. The caller supplies forecasts for remaining completion time, the next
step's CPU/wall cost and switching time; it may also supply an expected gain.
A negative or insufficient gain is declined. Evaluation is also refused when
the remaining horizon cannot cover planning plus switching, when the next step
would exceed the CPU/window/count budgets, or while another step is in flight.

An accepted callback performs one bounded unit of planning. Its actual planner
thread CPU and wall duration are charged, including failed calls. A call that
overruns the budget, returns after the useful horizon, or completes after the
job closes cannot authorize a switch. Its calculated values may still serve as
arithmetic work, but they are not an accepted execution decision.

These limits are cooperative: an already-running Python/native calculation is
not preempted. Forecast costs must be independently justified and refreshed when
their context changes. The default 50 ms CPU limit is a development budget,
separate from observation-callback accounting, not a validated production
optimum. Moving work to a thread does not remove interpreter, CPU or memory
contention. This gate is not a Bayesian posterior or a knowledge-gradient rule;
those must account for uncertainty and information value in addition to cost.

## Verification and remaining integration

Tests compare cached and uncached component results, mutate returned values,
change prices/geometries/dtype/implementation identity, exercise bounded eviction
and nested sessions, and check exceptions and concurrent access. Clock-controlled
tests verify that rejected steps never call their evaluators, startup is deferred,
deadlines do not reset, switching fits the early window, actual costs accumulate,
and late/failed results cannot authorize an action.

The matched audit in `direct_planning_reuse_20260922.py` alternates cached and
uncached evaluation order on the retained real native-PGEN source. It compares
entire graph hashes and continuation scores, reports construction/resume cost
separately, and verifies that the short-horizon gate performs zero evaluations.
Its component prices and modeled horizon are synthetic controls. It does not
measure production JIT overhead, time to first useful association, Bayesian
filter quality or a GWAS throughput gain.

The A100 integration run passed 209 tests with no failures or skips
(`results/planning_session_v2_20260922`). The final signed-gain and switching-
window amendments passed all 25 focused tests
(`results/planning_session_v3_20260922`). The final matched audit is saved in
`results/planning_reuse_v2_20260922/report.json`; its 120 package-source hashes,
benchmark hash and prior execution-report hash match the delivered source.

For N=2,049, M=4,097, K=512 and C=2, the three alternating-order repetitions
gave the following median calculator costs. Total includes constructing the
candidate graph and resuming the analytical checkpoint.

| Future chunk size | Uncached total (ms) | Reused total (ms) | Reused construction (ms) |
| --- | ---: | ---: | ---: |
| 128 | 250.06 | 112.05 | 51.07 |
| 256 | 201.44 | 97.05 | 43.77 |
| 512 | 184.74 | 89.03 | 38.47 |

All graph hashes and modeled remaining times matched exactly. The bounded
cache retained four entries (605,692 estimated bytes), with 29 hits and four
misses. Planning was never invoked before the first-output signal; the
short-horizon check also performed zero evaluations. Even the reused costs
are substantial relative to this small fixture's modeled remaining work.
These timings establish arithmetic reuse on this fixture, not production
autotuner profitability.

The public JIT lifecycle still needs to connect first useful output, the issue
frontier, observations, live admission and safe chunk control to these single
planning steps. The internal bridge and real writer/source boundary audit are
now described in productive_run_20260922.md. Source-generated single-option
proposals are described in jit_proposal_20260922.md; automatic candidate choice
and public integration remain open. Future phenotype-tile/GPU reassignment also requires the proper
safe boundaries. JAGWAS retains the full phenotype panel on each GPU. The full
objective remains open; neither a passing cache test nor a cheaper calculator
establishes the requested end-to-end behavior.
