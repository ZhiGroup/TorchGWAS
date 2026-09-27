# Output feedback in the first-chunk controller

The productive controller now retains at most 32 completed indexed-output
observations while counting later events in constant space. A sample binds the
producer device, admitted variant and phenotype partition, source interval,
reduction mode, survivor count, part bytes, part fsync flag and writer-active
wall interval. Significant-pair capacity is variants times phenotype width;
JAGWAS capacity is one statistic per variant and still requires the complete
phenotype panel. Dense writer progress remains a byte/row aggregate because
its blocks may split or combine scan chunks.

Before building a reduced-output forecast, the controller compares these
actual completed counts with the declared occupancy scenarios. A nonempty
chunk refutes an empty-only scenario; an empty or partially retained chunk
refutes a full-only scenario. In either case optional planning stops and the
scientific scan continues with its admitted chunk size. Mixed empty/full
scenarios remain eligible. Unbound or invalid output observations also stop
optional planning. This is a falsification check, not an estimator of future
selectivity. The first written event can already exercise it before any
source-window/model work begins.

At that event, a cheap issue-frontier check also stops optional planning when
an active partition has less unissued source than the largest forecast horizon,
or more than eight partitions exceed the model limit. It runs before price
revalidation, cache lookup, source-window construction and incremental probes.
The admitted scan continues unchanged. The proposer repeats the check on its
held frontier, so concurrent source reservations cannot make a stale cheap
check authorize a forecast.

Writer-active intervals are audited but never installed as independent writer
capacities. They include formatting, optional part fsync and shared-machine
effects. The observation does not change an immutable component price, its
original timestamp, or a cached source/work record. It adds no upfront
benchmark or candidate search.

The final focused remote A100 batch passed 148 tests across productive control,
significant-pairs and JAGWAS paths. It covers reduced-output counts, bounded
retention, unbound/duplicate refusal, dense progress exclusion, a contradicted
scenario rejected before model construction, and a short remaining-source job
that skips price validation and planning.

A fresh public two-GPU JAGWAS lifecycle run on the final source wrote equal
control/deferred/reuse outputs (maximum absolute difference zero) and completed
33 indexed chunks in each deferred pass. The audit retained the first 32 and
counted all 33; each had its full-panel producer identity, source range and
128/128 retained variants, with zero invalid or unbound events. Both deferred
decisions remained unapplied. The pulled artifact is
`results/jit_output_feedback_20260922/public_jagwas_v2/report.json`. Its
synthetic component rates and varying shared-server API times do not establish
throughput gain or calculator accuracy.

This is a safety and measurement step toward JIT tuning, not a profitable
configuration change. Next work is to bound independent decoder/selection/
writer component measurement across productive callbacks, preserve immutable
records for later jobs, and make a real chunk-size switch repay all observed
planning and switching costs. Phenotype tiling and device reassignment are
still separate, admitted moves; JAGWAS never partitions phenotypes.
