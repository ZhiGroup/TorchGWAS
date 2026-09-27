# Memory admission for chunk changes during a job

Changing the next chunk size does not resize the scan's allocated rings, and it
can move a native PGEN read to an LD restart absent from a regular fixed-chunk
schedule. The new adaptive_candidate module accounts for that shared allocation
before an online policy can choose among chunk sizes. It reuses the existing
tensor, decoder, pinned-buffer and output-memory ledgers; association runtimes
are not fitted.

This is a prerequisite for public first-chunk tuning. It does not yet connect
the measurement controller to candidate scoring or automatically choose an
execution configuration.

## Finite execution shapes

AlignedChunkSizeControl requires an explicit finite size set whose smallest
member divides every choice. It changes only reads not yet issued. The largest
size remains the ring capacity; devices, reader budgets, prefetch depth,
phenotype tiles and variant intervals stay fixed.

If fewer than the requested rows remain, the control issues the largest
admitted size that fits, followed eventually by the final remainder below the
smallest size. For choices 128/256/512, 385 remaining rows become 256+128+1;
they do not introduce an uncalibrated 385-row kernel.

Thus each fixed source interval requires its admitted full shapes and at most
one remainder shape. Existing general-purpose private chunk control continues
to support arbitrary sizes; it does not receive this restricted geometry
contract automatically. The memory report lists required and missing geometry.
Coverage alone does not validate a compiled kernel family or timing profile.

## Shared tensor and decoder allocations

The device calculation retains depth times N times capacity int8 input bytes.
It takes separate maxima for current temporaries, retained output, the previous
converted genotype, and any previous reduced result. This permits an old result
from one shape to coexist with the current expression from another. A maximum
of whole isolated-plan peaks can miss that combination. Setup, library
workspace, full JAGWAS factor and explicit reserves retain their existing roles.
Pinned allocations remain those of the largest fixed ring.

A complete exact encoded census on the smallest chunk grid covers every
possible read start. For each start, the native input buffer bound includes
the next capacity-sized payload and that start's LD prefix. Interior grid
prefixes are not added: the native reader restarts only at the actual read
start. The separate LD scratch bound covers packed base data and relative
offset storage. These reader buffers are grow-only.

The concurrent-reader bound uses the smallest admitted chunk size. A short
source may have one chunk at the allocated capacity but several active reads
at a smaller size. Genotype rings still retain their original capacity.

The candidate calculation starts from the existing dense, significant-host or
JAGWAS memory ledger, preserving writer, queue, predicate and selection
reserves. It adds the largest per-device increase in decoder workspace and
applies the mixed-shape tensor bound. Sequential phenotype tiles retain the
existing pinned-cache accounting. JAGWAS always retains the full phenotype
panel on each active GPU.

The caller must check supplied/live memory budgets and missing geometry before
using this report to execute. It remains conservative source accounting plus
explicit reserves, with unresolved allocator, driver, library and host-memory
terms. It is not a complete RSS or CUDA-reservation guarantee.

## Verification

A100 job 20260921-233111-796860 passed all 109 tests in 35.35 seconds, with no
failures or skips. Artifacts: results/adaptive_candidate_v1_20260922.

The new tests enumerate all reachable transitions for several aligned size
sets and tails, check a mixed-lifetime counterexample, and compare the LD
buffer bound against exact independently collected censuses at every possible
small-grid start. A source whose coarse starts are non-LD records requires a
larger read buffer when a changed chunk starts inside an LD block; the test
requires that increase to appear in host admission. Malformed source coverage,
wrong capacity/device binding and exceeded budgets fail before tensor work.

A real two-GPU test changes future chunk sizes while retaining fixed buffers
and complete JAGWAS factors. It checks exact coverage, allowed work shapes,
independent FP64 statistics and worker cleanup. Existing adaptive execution,
measurement-controller, tensor-memory and public reduction-autotune tests also
passed.

The separate native-PGEN execution audit in A100 job 20260921-233536-798146
completed on N=2049, M=4097, K=512 and two covariates. Inputs and metadata were
on /data (/dev/md0, XFS). It applied the memory report before execution, retained
a 512-row ring on both GPUs, and requested 128 -> 256 -> 512 after two and six
completed initial measurements. The transition sequence is an explicit test
control, not an automatic tuning policy. All early association results were
retained.

The emitted shapes were exactly 1, 128, 256 and 512. All 4097 variants appeared
once, 16 observations completed with none pending, and both GPUs kept the full
512-trait factor. The maximum absolute difference from the fixed-512 run was
0.0001220703125. Thirty-three independent FP64 OLS/correlation checks had maximum
absolute error 0.0001234262185789703.

The decoder correction added 1,350,570 and 675,432 host bytes for cuda:0 and
cuda:1. Modeled device budgets were 386,061,340 bytes each, including the explicit
256 MiB reserve. Observed peak allocated tensor bytes were 56,906,752 and
39,017,984; both devices also had a physical 2 GiB allocator cap. This one
fixture verifies the stated requests and execution, not a general reserved-VRAM
or RSS ceiling.

Artifacts are results/adaptive_candidate_execution_v1_20260922, generated by
benchmarks/direct_adaptive_candidate_20260922.py. The final verifier checked
unchanged input files. After retrieval, all 118 package-source hashes and the
harness hash matched the local checkout.

The remaining online work is to evaluate the actual remaining source schedule
under this fixed allocation, connect independent component measurement/reuse to
the shared analytical objective, and apply bounded decisions at safe boundaries.
Loaded CUDA stage spans remain drift observations, not hardware capacities.
