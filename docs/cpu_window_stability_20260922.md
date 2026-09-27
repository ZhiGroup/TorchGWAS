# Rejecting drift within a CPU-cost measurement window

`CpuServiceRefresh` now checks temporal consistency as well as overall spread.
Previously, a seven-sample window of 10, 10, 10, 20, 30, 30, 30 milliseconds
passed the default fourfold spread limit and could publish a 20 ms median,
although its recent cost had reached 30 ms. The reverse transition had the
same problem.

The calculator compares the medians of the first and last nonoverlapping
halves, using the existing drift-ratio and absolute CPU tolerances. For seven
samples these are the first three and last three; the middle sample still
contributes to the final median and spread. A changing window is rejected as
`unstable_measurement`, even when its overall spread is acceptable.

The same check applies to a short validation window before reuse. Two samples
whose median matches an old coefficient cannot authorize reuse when the two
samples disagree beyond the drift tolerance. They remain available for the
existing longer measurement window. Cached records are also checked from
their original ordered samples. Original values, hashes and observation times
are preserved; an accepted replacement is a new record. Snapshots and new
record provenance include the temporal-check result.

This is a bounded heuristic, not a confidence interval or a prediction of
future capacity. It does not change the measured coefficient, thresholds,
public API, probes or scheduling.

Remote A100 job `20260922-130821-1042151` passed **107 tests** in 22.84 seconds.
Tests cover transitions in both directions, a misleading short-check median,
legacy records with drifting samples, immutable reuse, expiry, budgeted
measurement, existing price binding and the absolute/relative tolerances.

The same job audited the first seven repeat means from each independent
primitive/control in four existing selector banks: **116 windows**. The old
spread rule accepted 104, and the new temporal check accepted the same 104.
It found no additional rejections in those saved windows. All four source
artifacts remained byte-identical and no measurements were published or given
new observation dates. This audit does not explain the remaining dense-output
timing error or qualify the calculator for automatic use.

Evidence:

- `results/cpu_window_stability_checks_v2_20260922/pytest.txt`
- `results/cpu_window_stability_audit_20260922/report.json`
- `benchmarks/audit_cpu_window_stability_20260922.py`
