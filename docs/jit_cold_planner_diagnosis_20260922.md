# Cold JIT validation and PGEN header diagnosis

This turn separated source opening from the productive planner callback. All
diagnostics ran on lab-a100 against the current local-WSL-pushed source; the
read-only evidence is in `results/jit_context_profile_20260922/evidence_v1/`.
The five pulled files matched their remote SHA-256 values.

Five repeated execution-context captures took 0.483 seconds under `cProfile`.
File hashing accounted for 0.253 seconds cumulative, per-CPU topology reads
about 0.080 seconds, and two mount lookups per call about 0.048 seconds. These
are loaded validation costs, not throughput capacities. In the prior public
JAGWAS audit without digest reuse, each of the two productive context captures
took roughly 0.125–0.141 seconds. Price comparison validation was only about
0.006 seconds per scenario, so batching it would barely help that case.

The existing `BindingDigestCache` already reuses file bytes while checking
fresh filesystem identity and retaining original price ages. With a structural
cache directory, the public audit's cold productive checks took 0.151 and
0.032 seconds for execution-context capture; the next job took 0.042 and
0.029 seconds with four disk digest hits. A second public audit now exercises
the code change: when planning-cost history has a cache directory and no
separate digest, refresh, or structural directory is supplied, it becomes the
digest-cache parent. Explicit `binding_digest_cache_dir=None` still disables
that fallback. The first job hashed 73,729,398 file bytes and used four
in-memory hits; the next had four disk hits and hashed zero file bytes. Its
context captures were 0.144/0.032 seconds cold and 0.035/0.035 seconds on
reuse. Control, deferred and reuse JAGWAS outputs matched exactly in both
audits. Shared-server API times are not a causal speed comparison. The final
fallback regression passed 103 focused tests.

The same 8,086,101-variant PGEN whose header parse previously took 11.608 CPU
seconds parsed in 0.079 CPU seconds in the current header-only profile.
Eight alternating warm-cache header parses with NumPy huge-page advice on/off
produced identical marker count, final record offset and mean record length.
The four on timings were 0.059–0.078 CPU seconds; the four off timings were
0.095–0.123 seconds. Huge-page advice also reduced minor faults in this
control. The earlier multi-second result is therefore not an established
algorithmic lower bound or evidence to disable huge-page advice. No parser or
reader implementation was changed on this evidence.

The remaining productive callback cost in the current small public case is
roughly 0.12–0.18 seconds for scenario graph comparison and, on the first
proposal, another approximately 0.18 seconds for source windows. The large-job
H100 absolute prediction gap remains unresolved; its original GPU pair was
still occupied when checked. A production JIT switch still needs independently
qualified prices, output-inclusive repayment, and a matched full-executor
validation rather than a loaded header or stage duration used as a price.

The later [current-code large-source compact-admission probe](../results/large_compact_admission_20260923/report.json)
(A100 job `20260922-232838-1248971`) used the server-local 8,086,101-marker
PGEN without reading genotype payload. In one ordered process, header parsing
took 1.114 CPU seconds and the compact 128-marker memory layout took 0.139
CPU seconds. Fixed and shifted memory envelope calculations for 128, 1,024
and 4,096 markers each took under 0.003 CPU seconds. Peak process RSS was
137,748 KiB. The script and input/model identities match the pushed source.
This isolates a currently cheap compact envelope; it excludes profile checks,
tensor memory, public startup, GPU and writer, and its loaded/cached state is
uncontrolled. The variable header times across these diagnostics do not
establish a stable cold-start cost or justify a parser change yet.
