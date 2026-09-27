# Remaining-work scenarios and accountable tuning cost

`window_forecast.py` estimates remaining work from three bounded analytical
window comparisons. It keeps pipeline fill/drain separate from marginal work
and makes extrapolation assumptions visible. It does not fit measured GWAS
durations or authorize an automatic change. The three inputs are evaluations
of the existing source/resource graph with independent component prices.

## Binding and scope

Each comparison carries a JSON-stable contract covering caller-validated
source/runtime identity, input identity and dimensions, fixed data settings,
complete component profiles, reader ordinal, devices, phenotype ranges,
starting variant, chunk sizes, output mode/threshold, writer settings and
shared capacities. The caller supplies an already validated model identity;
comparisons do not repeatedly read and hash the package themselves. A loaded
report and a fresh report have the same binding after JSON serialization.

Exactly three increasing horizons must use identical contracts. Their actual
reported rectangles must equal their comparison coverage, and each partition
must grow in the same proportion. This prevents a larger horizon from quietly
changing the workload balance, scientific output, price context or prepared
reader state. A pair of layouts can differ from each other, but each layout
stays fixed across its own three horizons.

Work is counted in variant–phenotype pairs, including full-panel JAGWAS.
JAGWAS still permits only variant sharding, with the complete phenotype panel
on each GPU. Dense and significant output retain their existing partition
rules. For significant output the comparison uses explicit survivor bins;
it does not infer the selectivity of future variants from a single mean.

## Calculation

Let the three modeled work amounts be W1 < W2 < W3 and their model times T1,
T2 and T3. For each layout:

    s1 = (T2 - T1) / (W2 - W1)
    s2 = (T3 - T2) / (W3 - W2)
    change = abs(s2 - s1) / max(s1, s2)

Both marginal costs must be positive. The caller declares an allowed relative
slope change, model-error allowance e, maximum R/W3 extrapolation ratio and
signed continuation adjustments [a_low, a_high]. For remaining work R >= W3:

    low  = max(0, (T3 + (R-W3)*min(s1,s2))*(1-e) + a_low)
    high = max(0, (T3 + (R-W3)*max(s1,s2))*(1+e) + a_high)

Anchoring at T3 avoids paying isolated-window fill/drain once per repeated
window. These intervals express declared scenarios and tolerance. They are
not statistical confidence intervals or proven hardware-time bounds. The
caller must state assumptions about future source work, output occupancy,
resource capacity and partition balance. A small extrapolation limit must not
be silently enlarged to make a large job eligible.

The signed boundary adjustments represent work missing from an isolated
prepared window: already-issued reads, active kernels, writer staging,
preparation, tails and final API metadata/directory publication. Supplying
zero is an explicit empty-state scenario, not a measurement of a live job.
The caller must separately establish the actual unissued extent, current
admission/freshness and applicable boundary adjustments.

Dense continuation requires the already-resolved writer block size. An
automatic block size would incorrectly resolve again from the new first chunk.
The forecast also flags an unseen periodic-writeback regime if extrapolation
would introduce submissions beyond a horizon that has not yet modeled both
submission and the previous-range wait/drop cycle. This detects one concrete
nonlinearity; it is not a proof that every future service regime is stationary.

Statuses distinguish `stable_scenario`, `unstable_marginal_cost` and
`unmodeled_writer_regime`. Unstable or unmodeled cases cannot pass the payback
prerequisite. Even a stable scenario retains `selection_validated=False`.

## Costs during productive tuning

The numerical payback prerequisite is:

    baseline.low - candidate.high > planning + switching + publication + reserve

All four costs are explicit. Cache publication is part of output-inclusive
execution even when deferred until successful completion. A reused structural
record is not republished, but checking that no publication is required can
still cost time. Timing three model comparisons excludes caller binding/header
work and the final forecast arithmetic; an execution controller must charge
those as well.

`IncrementalPlanningBudget` and `ProductiveTuningRun` now charge cumulative
planning wall time, including earlier unsuccessful steps, when checking the
next decision. The cheap precheck includes publication/reserve forecasts, and
the final decision uses actual accumulated planning cost. Thus spreading work
over early chunks cannot erase its cost. Each decision conservatively repays
all recorded planning costs plus its proposed switch and deferred costs;
previous forecast savings are not banked as credit.

When structural persistence is enabled, the controller requires an explicit
publication-cost forecast before opening the cache or evaluating a proposal.
Missing cost evidence skips that step and leaves the scan running. An explicit
zero is permitted; responsibility for that forecast remains with the caller.
The actual publication cost remains in the completed run report. Structural
cache construction/loading is lazy, after written useful output, and is
included in the planning step. The default planning budget is unchanged.

Persisted evidence keeps its existing semantics: structural work is reusable
only while its dependencies match; empirical prices retain their original
measurement time and expiry; available memory/contention must be read again.
Neither an early-chunk validation nor a cache hit renews an observation's age.

## Remote verification

Job `20260922-060517-938919` completed with **238 tests passing in 62.86 s**.
The log is `results/window_forecast_v5_20260922/tests.log`. Tests include
contract/coverage mismatches, reader-state changes, proportional partition
growth, JSON round trips, writer-regime changes, extrapolation limits,
continuation adjustments, stable/unstable slopes, and actual controller gates
for prior planning and publication cost. Existing header, window, cache and
productive-execution regressions also passed.

The two final analytical audits are:

- `results/window_forecast_arithmetic_short_v3_20260922/report.json`
- `results/window_forecast_arithmetic_long_v3_20260922/report.json`

Both use the same immutable plain-genovec PGEN fixture: N=2,049 and M=8,193,
on XFS `/dev/md0` mounted at `/data`. They compare B=128 against B=512 on two
modeled GPUs. Dense/significant output uses two 512-trait partitions; JAGWAS
retains the complete 512-trait panel on each variant shard. Captured GPU
geometry is real; all component prices and survivor occupancy are synthetic
controls. Offline fixture profile construction includes exact censuses and
is outside the comparison timers. These are not productive startup runs.

The short audit models 512/1,024/1,536 variants per device and checks a fourth
horizon of 2,048. The long audit models 1,024/2,048/3,072 and checks 4,096.
Both declare 5% model-error allowance, a 10% slope-change limit, zero boundary
adjustments and maximum extrapolation ratio 2. The JAGWAS shard starts differ
between the two audits to keep the respective horizons disjoint; the plain
record form/work per variant is homogeneous. No threshold was relaxed after
seeing a refusal.

| Output | Short maximum slope change | Long maximum slope change | Long three-comparison CPU cost |
| --- | ---: | ---: | ---: |
| Dense | 6.43% | 0.000893% | 1,455.94 ms |
| Significant, empty | 12.08% | 0.648% | 1,019.74 ms |
| Significant, sparse | 12.08% | 1.544% | 404.92 ms |
| Significant, dense | 14.54% | 5.824% | 248.84 ms |
| JAGWAS | 21.40% | 0.445% | 1,141.84 ms |

Four of five short comparisons refuse extrapolation for unstable marginal
cost; all five long comparisons pass the declared stability test. All 20
held-out model times (two layouts, five modes, two audits) fall inside the
declared scenario intervals. For the long audit, absolute midpoint error is
at most 0.577%. This establishes arithmetic consistency on this fixture, not
accuracy against measured GWAS elapsed time or guaranteed interval coverage.

No comparison passes payback: the candidate is slower in the held-out model
and has a negative scenario gain floor. The explicit payback scenarios add
60 ms for publication and 10 ms reserve, with zero switching cost. These are
declared costs, not new measurements. The long three-comparison wall costs
range from 250.79 to 1,494.47 ms and exceed the default planning allowance;
no default-budget compliance or calculator speedup is established.

Both reports' 131 package-source hashes, six helper/geometry hashes and
benchmark hash match the tested local source. Their input hashes agree, and
the script checks source/input identities again after each audit. Reports
remain in the shared project results area; the input is on local `/data`.

## Remaining integration

This patch does not enable public automatic switching. The subsequent
[productive forecast bridge](productive_forecast_20260922.md) binds bounded
comparisons to exact held unissued ranges and explicit unequal-partition
scenarios. [Indexed output evidence](indexed_output_evidence_20260922.md) binds
observed survivors to their producing partitions. Remaining work must resolve
actual writer state and future occupancy, and validate component prices
and prediction error against real output-inclusive runs. Scalable phenotype
tile/device reassignment remains separate. A homogeneous small fixture cannot
establish production-scale prediction accuracy or profitable tuning.
