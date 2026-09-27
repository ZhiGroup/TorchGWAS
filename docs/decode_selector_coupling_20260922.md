# Independent decode–selector coupling control on H100

The frozen H100 candidate predicted 15.862 s and had a 55.769 s measured
executor median. Earlier complete scans showed excess loaded CPU work in PGEN
decode and significant-pair selection, but generic selectors alone did not
reproduce the worst loaded penalty. This control measures whether those two
unchanged production primitives slow each other when they run together. It
does not fit a whole-GWAS duration or replace an independent capacity price.

The diagnostic project is `/home/x/work/torchGWAS-calculator-diagnostics-h100`.
Job `20260922-200008-1171545` ran
`qualify_decode_selector_coupling.py` on `lab-h100`, using the frozen
`torchGWAS-calculator-h100` source, the real local `/data` PGEN prefix, and
physically written empty-output selector buffers. It used the original
eight-core affinity, two-second measurement windows, and three randomized
repeats of selector-only (two workers), decoder-only (two workers), and mixed
(two plus two) conditions. Reader construction, array preparation and the
first decode are outside each timed window. Source and input identities were
checked; every repeated first-chunk decode matched its initial byte digest.
The completed schedule, nine run records and independently checked summary
are under `results/decode_selector_coupling_20260922/v1/` in that project.

| Matched repeat | Mixed / selector-only median CPU per call | Mixed / decoder-only median CPU per call | Selector throughput ratio | Decoder throughput ratio |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.170 | 1.053 | 0.845 | 0.642 |
| 1 | 1.050 | 1.006 | 0.529 | 0.715 |
| 2 | 1.064 | 1.018 | 0.726 | 0.744 |

Selector-only median CPU service ranged from 6.81 to 7.29 ms/call; mixed
service ranged from 7.25 to 8.31 ms/call. Decoder-only median CPU service was
8.37–8.45 ms/call; mixed was 8.50–8.84 ms/call. Mixed-worker runnable scheduler
wait was substantial (selector totals 0.99–1.91 s and decoder totals
0.88–1.63 s per two-second condition). Selector-only wait itself ranged from
0.40 to 1.84 s, showing that external shared-host pressure varied within the
schedule. CPU service, elapsed throughput and runnable wait are separate
observations; the smaller throughput is not a measured increase in primitive
CPU work of the same size.

The control covers a hot bounded PGEN prefix, empty significant output, and
CPU decode/selection without GPU transfers, result workers, or the full
executor's memory lifetimes. Two-second randomized repeats do not certify
steady long-job service. GPU 1 was occupied by another job, so the original
GPU 1+2 frozen scan was not rerun. The observed 5–17% selector CPU increase
cannot explain a loaded worker mean near 28 ms/call or the full 15.862-to-
55.769 s gap on its own. It does establish that runnable delay and current
CPU availability need to be observed in a full, matched scan before treating
the isolated 6.865-core model capacity as available during tuning. Do not
install these loaded spans as component rates or infer a CPU-capacity
multiplier from them.
