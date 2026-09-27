# Matched scenarios for productive chunk decisions

The public initial-chunk tuner now compares baseline and candidate within the
same declared host-sharing and output-occupancy scenario. Its previous
aggregation subtracted the largest candidate upper forecast from the smallest
baseline lower forecast, even when those forecasts came from different
scenarios. This could reject a candidate that saved time in every scenario.

For example, paired forecasts of 10 to 8 seconds and 100 to 80 seconds imply
savings of 2 and 20 seconds. The previous crossed comparison was 10 minus 80,
or minus 70 seconds. The corrected minimum paired saving is 2 seconds.

For scenario `s`, the bounded-window calculator continues to produce a baseline
lower forecast `L_s` and candidate upper forecast `U_s`. The productive decision
now uses:

```
conditional_gain = min_s(L_s - U_s)
switch only if conditional_gain > cumulative_planning_wall_seconds
                                  + switching_seconds
                                  + publication_seconds
                                  + reserve_seconds
```

The pair with the smallest saving is returned to the existing payback controller.
Every paired baseline forecast must also remain within the caller's declared
remaining-time horizon; checking only the limiting-gain pair could overlook a
longer scenario. The audit retains all scenarios, their paired gains, the
limiting scenario index and the largest baseline lower forecast.

These are matched counterfactual conditions for identical remaining source and
scientific output. Different chunk sizes still have their own modeled resource
demands and dependency schedules. This change does not assert that future load
is stationary or that the declared scenario set covers the real machine.
Independent model-error and continuation-boundary intervals remain inside
each paired comparison. An unstable scenario still prevents a switch, and all
existing price freshness, input identity, memory and planning-budget checks
remain in force. No measured coefficient or scientific threshold changes.

## Verification scope

The new regression controls exercise dense, significant-pair and full-panel
JAGWAS routes, both scenario orders, a gain in every scenario, a loss in one
scenario, exact payback equality, previously accumulated planning cost, an
unstable scenario and a horizon violation outside the limiting-gain scenario.
They verify that issued source ranges stay intact and subsequent ranges use
the chosen size, or the original size after rejection. Scenario times in these
tests are explicit synthetic controls, not independent hardware measurements.

Existing tests separately cover actual GPU output parity across chunk changes,
stale evidence, writer failure and cleanup, as well as construction of real
dense/significant calculator forecasts. The remote artifacts for this change
are in `results/matched_productive_scenarios_20260922/`.

A100 job `20260922-142826-1079293` passed **255 tests in 86.72 seconds**,
without failures or skips. This includes all 36 new matched-condition cases
and the existing public GPU execution, forecast, productive-budget,
price-binding, adaptive-chunk and digest-lifecycle tests. The pulled source
manifest matches the final implementation and test file. No measured
component-rate or whole-job performance claim follows from this test batch.

The earlier public audit in
`results/public_initial_chunks_20260922/execution_v2/report.json` used only one
host/output scenario. Its baseline/candidate marginal relative changes were
0.1645/0.3279 and 0.2166/0.2840 in the deferred/reuse jobs, exceeding the declared
0.1 tolerance. Both remain rejected for unstable extrapolation; this fix does
not reinterpret them as profitable decisions. Those component prices were
largely synthetic. Calculator accuracy, coverage of realistic resource
conditions and profitable tuning on production jobs remain unqualified.
