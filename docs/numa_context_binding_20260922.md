# NUMA policy in calibration compatibility

Detailed calibration contexts now include the calling thread's memory policy,
policy node mask, allowed memory nodes and the host's automatic NUMA balancing
setting. CPU affinity alone does not identify these settings. Existing exact
context comparisons reject a saved profile if any of these fields differ or
the new field is missing. They do not rewrite the original record, observation
date or expiry date.

`src/torchgwas/numa_context.py` reads Linux policy through `get_mempolicy` in
`libnuma.so.1`. It preserves raw policy flags, uses a mask sized for possible
nodes, and separately queries allowed nodes. Only the library handle is cached;
policy values are read afresh. Repeated policy and balancing reads detect changes
during capture. Unsupported, failed, inconsistent or malformed queries fail
closed instead of treating an unknown policy as compatible. Linux detailed
context collection therefore requires this library and readable policy state.

This collector is read-only. It does not set memory policy, alter CPU affinity
or change the system's automatic balancing setting.

## Scope

The context describes the calling thread's default policy and allowed nodes.
New threads inherit that policy, but existing threads can have other policies,
and individual mappings can override the default. It does not certify actual
page residency or available capacity. Those are separate from profile identity.
See [set_mempolicy](https://man7.org/linux/man-pages/man2/set_mempolicy.2.html)
and the read-only [get_mempolicy](https://man7.org/linux/man-pages/man2/get_mempolicy.2.html)
interface.

Recording automatic balancing matters because it can introduce page faults
and migration while assessing memory placement. A minor-fault count alone does
not distinguish this from other causes. The
[kernel documentation](https://docs.kernel.org/admin-guide/sysctl/kernel.html#numa-balancing)
describes the mechanism; the [selector diagnostic](selector_stage_diagnostic_20260922.md)
records the observed ambiguity. No measured throughput gain is claimed here.

## Verification

Final A100 job `20260922-153744-1098179` passed 143 tests in 25.94 seconds across
NUMA context, detailed calibration, NumPy context, cache, price binding,
resident-copy refresh and productive digest suites. The 15 NUMA tests cover
high node indices, raw flags, fresh reads, legacy rejection, query errors,
changes during capture, and real private-thread policy changes and inheritance.

The same job ran `benchmarks/direct_numa_context_binding_20260922.py` against
the actual collector. A profile saved under the original policy was rejected
after privately setting mode 2 (`MPOL_BIND`) on node 0, then mode 8194
(`MPOL_BIND | MPOL_F_NUMA_BALANCING`). Restoring the original thread policy
restored the complete original context and allowed reuse of the unchanged
profile. The harness restored policy in a `finally` block and changed no global
NUMA settings. Its synthetic coefficient is solely a compatibility control,
not a measured capacity or a qualification of calculator accuracy.

Pulled evidence is in `results/numa_context_binding_20260922/`, including
`tests_final.log`, `final_source.sha256`, and `audit/report.json`. All 141
executed package source hashes, the audit harness hash, and the final test/source
hashes matched local files. An earlier 143-test run in 34.46 seconds is preserved
separately; the final run followed the C return-type correction and extension
of the real policy test to cover the balancing flag.
