# First-chunk staged screen integration

The public productive controller can opt in to an evidence-only staged
screen with `register_staged_screen(candidate_factory, options)` or the
native `initial_chunks.staged_screen` configuration described in
[native first-chunk candidate construction](jit_native_staged_chunk_candidates_20260923.md).
Registration is permitted once, before first written output, and does no
PGEN scan or calculator work. The explicit hook's caller supplies source
capacities, occupancy scenario, finite budgets and independently priced
layouts; the native configuration derives same-partition chunk candidates
from the active context after useful output.

Written output continues to earn one bounded source-metadata step per event.
When the staged ledger completes, its callback starts a separate screen worker
outside the writer callback and source-stage lock. The worker binds actual
written output to one issued-source snapshot, validates the current profile
and live capacity, builds the unissued frontier, and passes that frontier and
the completed source stage to `productive_staged_partial_screen`. The whole
worker CPU and wall time are retained with the screen. It checks the issue and
output token and the profile digest again afterward, recording
`current_evidence`, `stale_frontier`, `stale_binding`, `screen_budget`, or
`screen_error`. A job finishing before stage completion cannot launch a new
screen. Shutdown joins an active screen before finalizing its audit.

The screen's mode-specific partial resource floors, real output boundary and
staged source cost remain evidence. A result marked `current_evidence` means
the issue/output token and active profile stayed unchanged and the supplied
candidate source/work prices passed their declared-target audits. It does not
include queued/in-flight output or predict finite completion. No chunk,
phenotype tile, or GPU change can be applied by this path.

The optional `max_rebases` setting permits at most two extra background
screens if source issue or written output advances during a screen. Each retry
captures a new atomic source snapshot, binds the current written-output
boundary, and asks the candidate factory for layouts covering the new
unissued frontier. The completed metadata ledger is reused. The original
CPU/wall budget is shared across attempts; each screen cooperatively stops
after a candidate, and the audit retains bounded attempt summaries and the
latest full report. A pending profile refresh, job finish, exhausted source
or budget prevents another retry. The default is zero retries. A fast scan
can still outrun the bounded attempts and leave `stale_frontier` evidence.

The [source-price binding audit](jit_staged_source_price_binding_20260923.md)
and [work-price binding audit](jit_staged_work_price_binding_20260923.md)
reject supplied rates that disagree with the active profile and report
matching but undeclared measurements separately. Shared transfer capacity,
GPU shape service and some reduced-output prices are still unbound.

The A100 regression exercises a real PGEN header and four synthetic typed
written-output events, then holds the screen while the source revision or
profile changes. It checks the off-writer thread, exact remaining range,
bound output, retry budget and status of each case. The latest focused A100
source-stage/controller/screen suite passed 120 tests. These are lifecycle
and accounting checks, not a loaded GWAS throughput or switching validation.
