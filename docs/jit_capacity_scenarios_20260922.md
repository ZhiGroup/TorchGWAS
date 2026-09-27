# Conditional shared capacity in bounded JIT comparisons

The calculator now keeps calibrated per-window/per-tile service profiles fixed
when a caller supplies less shared CPU, DRAM, input or output capacity. It rejects
capacity above any participating profile's calibrated value. Bounded comparison
contracts retain both the unchanged profile digests and the conditional shared
capacities. This permits a resource-availability sensitivity calculation without
publishing a new component price or mutating immutable evidence.

`initial_chunks.capacity_scenarios` optionally names fractions for all four
shared resources. It requires an all-ones nominal case; each fraction is in
`(0, 1]`. Together with host-sharing and output-occupancy scenarios, at most four
combinations are admitted. Every combination compares baseline and candidate
under the same capacity and contributes to the minimum paired gain and remaining
time check. The extra model work occurs after useful output and is charged to
the existing CPU/wall planning budget. The default remains the single nominal
case, so startup and callback work do not increase by default. These fractions
are assumptions supplied by the caller, not a live-capacity estimator.

The read-only frozen H100 significant-pair diagnostic
`/home/x/work/torchGWAS-calculator-diagnostics-h100/results/frozen_cpu_capacity_scenarios_v2_20260922/report.json`
tested one already selected `N=35,365`, `M=1,048,576`, `K=16,385` candidate.
The frozen full model required changing the tile-profile CPU capacity together
with shared capacity; that diagnostic is therefore a conditional sensitivity,
not an immutable-profile comparison. It left component prices, source, plan and
profile artifacts untouched and used no observed runtime to fit a parameter.

| Assumed CPU capacity | Predicted executor | CPU work floor |
| ---: | ---: | ---: |
| 6.865 cores (100%) | 15.862 s | 14.413 s |
| 5.149 cores (75%) | 15.879 s | 14.636 s |
| 3.433 cores (50%) | 16.914 s | 16.050 s |
| 1.716 cores (25%) | 27.136 s | 27.118 s |
| 0.858 cores (12.5%) | 48.964 s | 48.964 s |

The prior observed executor median was 55.769 s. The capacity sensitivity does
not identify the cause: matching that median by CPU capacity alone would require
an even smaller unverified allowance than 0.858 core. Modeled active-wait CPU
work also changes as the graph is solved, so it cannot be treated as an
independent fixed demand. The full-candidate preparation itself took 52.292 s
for this diagnostic; it is not the deferred startup path. The bounded JIT window
model still needs validation against real loaded stage and durable-output
measurements before a profitable switch is trusted.

On lab-a100, job `20260922-204359-1205016` passed 215 focused tests across
window composition, forecasts, JIT lifecycle, price binding and tiled modeling.
After adding integrated conditional-capacity cases, job
`20260922-204559-1205351` passed three targeted tests: the first-chunk
controller priced both nominal and half-CPU scenarios for dense and significant
output using the real bounded model, and the real JAGWAS window retained the
full phenotype panel and identical payload under lower shared CPU capacity.
These fixture prices are synthetic; the tests establish contract and execution
behavior, not runtime prediction accuracy.
