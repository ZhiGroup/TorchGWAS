# Bounded device selection geometry

The device significant-pair selector and the calculator now call the same
`device_selection_shape(rows, traits, max_cells)` helper. It chooses the block
height and trait width that minimize the number of blocking CUDA nonzero calls
while keeping every block within the existing cell cap. For a fixed height,
the widest allowed trait strip cannot increase the number of blocks. The
remaining objective is `ceil(rows / height) * ceil(traits / width)`; both
factors are constant between quotient breakpoints, so the helper skips those
intervals. Ties favor the old wide-strip geometry. Full source chunks and a
short tail can choose different shapes, matching the actual per-call selector.
The compact D2H, indexed-archive and count-latency floors use the same block
census. No CUDA or Triton kernel was added.

For 512 variants by 600,000 traits with the production one-million-cell cap,
the old 1-by-600,000 shape made 512 selection blocks. The chosen 512-by-2,048
shape makes 293. Exhaustive small-shape tests check the minimizing count;
source-traced selector, real device-output and calculator regressions check
the changed geometry and numerical output. The focused A100 batch passed 80
tests in `20260923-012009-1284623` after updating count-floor expectations.

Matched selector-only A100 controls used identical resident tensors and
alternating old/new calls. The [empty-output report](../results/device_selection_geometry_20260923_v1/report.json)
has median wall times 0.672/0.180 seconds across four calls each. The
[sparse report](../results/device_selection_geometry_sparse_20260923_v1/report.json)
retains the same 512 planted pairs in both arms; nonempty parts fall from
512 to 293 and median wall times are 0.307/0.131 seconds. Every adjacent
old/new pair favored the new shape. The reports were pulled from A100 and
SHA-256 matched remote files and the executed selector source. These controls
exclude genotype decode, the statistics kernel, writer durability, queueing
and a full GWAS job; shared-host load varied, especially in the empty run.
They demonstrate a selector-level gain for this shape, not an end-to-end
speedup or an independent JIT timing price.

The selector's GPU memory admission now calls the same shape helper. Before
this correction it could request workspace evidence for the old mask extent
while executing a larger new block. Admission now requires exact installed
nonzero allocations for every full/tail mask extent; an old-shape census
fails. The source-traced selector also hashes `selection_geometry.py`, so a
saved kernel census cannot silently outlive a geometry change. A fresh
[A100 launch census](../results/device_significance_geometry_20260923_v3/census.json)
passed source-operation and launch audits for all 12 empty, sparse, dense and
invalid cases. The immutable test fixture is a copy of that measured census,
not a changed hash on the old one. The [large-panel memory report](../results/device_selection_memory_large_20260923_v1/report.json)
uses separate A100 allocation captures for 1,048,576 and 1,015,808 mask
cells. It computes a conservative 138,489,856-byte selector-local GPU budget
for one 512-by-600,000 source chunk. The focused current-code regression
passed 76 tests in A100 job `20260923-013929-1307377`; the final cross-mode
regression passed 342 tests in `20260923-014111-1307765`, including native
multi-GPU significant output.

The optimizer still needs a finite whole-job continuation, loaded calibration
and safe layout transition. Fewer nonzero calls lower one service term but do
not prove the best phenotype tile, chunk width or GPU assignment.
