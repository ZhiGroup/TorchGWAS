# Measurement freshness after reading a saved cache

`CalibrationParameterCache.lookup` now rechecks measurement age after reading
and validating the cache directory. Previously it compared every record with
the time at lookup entry. For example, an observation at time 90 with a
30-second lifetime could be returned with age 10 when a lookup started at time
100 but finished reading at time 120. It is now rejected at the expiry boundary.

The effective lifetime remains the minimum of the immutable producer lifetime
and any shorter caller limit. A valid hit reports its age at the post-read
check. Neither its original observation time, publication time, lifetime nor
record bytes change. Records already unusable at entry are discarded before
retaining their payloads; surviving records are checked again after I/O.

Selection still prefers the newest valid observation. If that record expires
during lookup, an older record with a longer, still-valid lifetime can be used.
Directory traversal order does not affect that decision. A backwards wall-clock
step between lookup entry and the final check returns a cache miss with reason
`clock_moved_backwards`. This is a guard for that lookup, not a solution to clock
changes between separate jobs.

Structural records continue to depend on exact source, implementation, device
and other declared dependencies, without empirical expiry. Live memory and
contention records remain audit-only. CPU/GPU/transfer/storage capacities and
stage observations retain their finite lifetimes. Stage observations are still
not independent component prices.

Consumers must recheck exact bound evidence before using it for a tuning
decision: a cache hit is not a reservation or a guarantee that the record stays
fresh indefinitely. The public price-binding and CPU-refresh paths already
perform these later checks. This change fixes the generic lookup contract and
its reported age; it does not qualify measured rates or calculator accuracy.

## Validation

Regression tests inject time advancing inside actual cache-file reads, including
expiry at the producer and caller boundaries, updated hit ages, both traversal
orders with an older valid fallback, backwards clock movement, unchanged
record bytes, and structural reuse after a long read. The remote test artifact
is `results/cache_read_freshness_20260922/`.

The final A100 job `20260922-142135-1078815` passed **267 tests in 72.76
seconds**, without failures or skips. These include exact bound-record and
price-binding checks, CPU refresh/reuse, digest and structural caches, initial
observation collection, and public initial-chunk execution and cleanup. The
pulled `final_source.sha256` matches the final implementation and regression
test file. The preceding run also passed 267 tests before the small change that
retains early rejection of already-unusable record payloads; its log and source
hashes are preserved separately.

The in-job host-memory check was also reviewed. It compares the full admitted
host envelope against current available memory, so allocations already owned by
the running job can make it conservatively stop tuning. It remains unchanged:
crediting all process RSS would also credit unrelated allocations and does not
establish ownership or the residency of the modeled buffers. A future correction
needs explicit owned-allocation accounting. This limitation is separate from
calibration freshness and runtime-prediction qualification.
