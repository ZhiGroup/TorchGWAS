# Held-out pinned transfer calibration for JIT multi-GPU pricing

The [transfer diagnostic](jit_shared_transfer_diagnostic_20260923.md) showed
that one GPU-agnostic or device-count-scaled copy rate is unsafe, but it had
no producer that could populate `shared_transfer_capacities`, per-device
H2D/D2H prices or `shared_links`. This change adds that producer and a
binding path. It does not enable a JIT layout switch.

## Producer and qualification

[`transfer_calibration.py`](../src/torchgwas/transfer_calibration.py) measures
pinned uint8 copies on explicit streams for every single GPU, every
caller-declared shared link and the whole active set, in both directions.
Each round runs all cells in a fresh seeded order; a round's cell value is the
median of its samples. Even rounds are calibration, odd rounds held-out.

A `transfer_capacity` record (`pinned_transfer_groups.v1`) is published only
if every check passes; otherwise the report lists each failure and nothing is
stored:

- selected GPUs idle before any of our allocation or traffic, and no foreign
  compute process on them before or after any round;
- every sampled pinned-buffer page on one NUMA node, before and after, and
  with `--expect-local-node` on the GPU's PCI NUMA node (`move_pages(2)`
  query mode; nothing migrates);
- each held-out median within 10% of the calibration median and inside the
  calibration span widened by half that tolerance;
- no group's high value above the sum of its members' high values;
- an unchanged detailed execution context before and after measurement.

The value holds `low`, `median` and `high` scenarios: minimum, median and
maximum of calibration-round medians. They are observed spans, not confidence
intervals or hardware ceilings. `high` is appropriate for necessary resource
floors; `low` for a conservative completion ceiling (handoff step 2).
Topology is never inferred: shared links are declared with `--links`.

Record dependencies are the process's source identity, full
`execution_context` (CPU affinity, NUMA policy, GPU UUID/PCI, environment,
input/output storage) and the measurement protocol. A record therefore binds
only to a job with the same context: **run the producer with the job's own
input/output paths and CPU affinity.** `attach_transfer_prices(profile,
record, context_name=..., scenario=...)` applies one scenario to a named
context, adds exact price-binding targets for every per-device, shared and
link value, and keeps existing bindings. An overlapping target (an already
bound transfer price) is refused, and rebinding never renews the original
observation age. The driver
[`transfer_capacity_calibration_20260923.py`](../benchmarks/transfer_capacity_calibration_20260923.py)
also attaches the new record to a skeleton profile and validates it against a
freshly captured context, as a later job would.

## Observations

Protocol: 32 MiB × 8 copies, 3 samples per cell, 2 warmups, 10 rounds, seed
20260923. Records expire six hours after their first observation.
Decimal GB/s, calibration median (held-out median in parentheses):

| Host, GPUs, CPUs | Cell | H2D | D2H |
|---|---|---:|---:|
| 2080 Ti 1,2,5; CPUs 0,5,11 (node 0) | each GPU alone | 11.99–12.03 (12.00–12.03) | 11.73–12.23 (11.62–12.05) |
| | 1+2, same PCIe switch (declared link) | 12.13–12.14 (12.13–12.15) | 11.55–12.02 (11.56–11.96) |
| | 1+2+5 | 18.13–18.19 (18.18–18.18) | 12.65–13.03 (12.95–13.15) |
| H100 0,3; CPUs 7–10 (node 0) | each GPU alone | 54.62–54.72 (54.68–54.74) | 54.46–54.79 (54.41–54.68) |
| | 0+3 | 106.08 (106.85) | 95.39 (95.66) |

2080 Ti ranges span the qualified runs v2 and v3. Two GPUs behind one 2080 Ti
PCIe switch receive no more bandwidth than one. All three together stay near
18 GB/s H2D, below the declared-link-plus-single sum of about 24 GB/s, so a
further shared limit exists upstream; it is observed, not identified. The
2080 Ti whole-set cell had one slow H2D round (low 14.7 GB/s in v2) and fast
D2H rounds (high 19.6 GB/s in v2, 15.7 in v3); these stay in the scenarios.
H100 H2D scales almost linearly across two GPUs on one node; D2H does not.

The first 2080 Ti attempt (v1) did not publish. Its idle gate read
`nvidia-smi` utilization right after each round, which averages our own copy
traffic. The gate now checks utilization once before allocation and foreign
processes per round, recording utilization on other GPUs as neighbour context.
Its only held-out failure was a 0.5% miss of a very tight span, which led to
the half-tolerance widening. v1's values agree with v2 and v3 and remain in
its report.

No A100 record exists. Every A100 GPU carried another user's compute process
throughout this work, which the gate refuses by design.

| Report | Qualified | Record |
|---|---|---|
| `2080ti_1_2_5_v1.json` | no (gate design, above) | — |
| `2080ti_1_2_5_v2.json` | yes (earlier source revision) | `024573f8…c273b9` |
| `2080ti_1_2_5_v3.json` | yes; real binding check 10/10 targets | `697d4886…09493a` |
| `h100_0_3_v1.json` | yes; real binding check 6/6 targets | `c192242c…11cf28` |

Artifacts: `results/transfer_capacity_20260923/` (reports and content-addressed
cache records), pulled locally. Final source SHA-256:
`transfer_calibration.py` `d052ad40…ef05`, driver `e3f34244…12fd`.

## Tests

`tests/test_transfer_calibration.py` (13 tests) covers group declaration,
stable qualification, each gate, held-out drift, the widened span, expected
NUMA node, exact binding of all targets with tamper and expiry refusal,
link-load consumption, attachment preserving existing bindings and refusing
rebinding or context mismatch, and a real-GPU measurement with page placement.
It passed on the 2080 Ti with the price-binding and detailed-calibration
suites (56 tests before the attach helper; 13 after).

A full-suite run on the 2080 Ti gave 3,957 passed, 67 failed, 54 skipped and
11 collection errors. No existing file was modified. The collection errors
import benchmark or script files absent from this checkout. The failures
examined are host-bound: for example
`test_device_significant_screen_charges_empty_block_barriers` refuses a
"stale installed device-selector launch profile" bound to A100. The broad
regression should be repeated on A100 when its GPUs allow it.

## Limits

- Pinned placement is verified for the producer's own buffers. A job's pinned
  rings are allocated by its own threads; a job whose pages land on another
  node is not detected by binding.
- Traffic from unselected GPUs sharing a switch or root complex is recorded,
  not excluded. Host CPU load is recorded (2080 Ti whole-machine busy fraction
  0.71–0.90 across rounds; H100 about 0.26).
- These are copy-service observations. Handoff step 5 still has to show that
  the bound scenarios bracket transfer service inside held-out native jobs.
