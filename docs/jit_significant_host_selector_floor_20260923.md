# Compact host selector for significant pairs

`native_layout_significant_host_selection_floor` counts the public
host-selected significant-pair predicate, NumPy nonzero, gathers and index
rebasing over a fixed unissued phenotype tile. It requires explicit
retained-pair intervals, the active host-selector/NumPy protocol, a bound
threshold-one setting and the independently priced primitives already used
by `host_selection_service`. The public significant scan returns beta for
selection even if the indexed writer omits beta from t-only output.

NumPy nonzero has empty, sparse and dense regimes. The compact calculation
checks which regimes remain possible for each full or short-tail chunk from
the tile's retained-count interval. Every other primitive has fixed
per-chunk dispatch plus a common per-retained-pair slope. This yields
conditional CPU and logical host-memory work intervals in constant work per
tile; no full source chunk graph is expanded. The partial envelope charges
the selector alongside source, DMA and the indexed archive under shared CPU
and host-memory capacities. A GPU's sequential phenotype tiles also retain
their selector CPU chain; host selection occurs in the producers, while NPZ
writing uses the separate single consumer.

These are partial service bounds. First-use allocation, queueing, device
selection, current in-flight buffers and final metadata are not priced. The
upper endpoint is not a completion ceiling, and neither endpoint authorizes
a JIT chunk, phenotype tile or GPU switch.

The A100 layout/significant/JAGWAS batch passed 154 tests in job
`20260923-004814-1263921`, including empty, sparse, dense and short-tail
regimes. A final 32-test batch in `20260923-004947-1264219` checked selector
and archive resource work together. These are component-equation tests;
loaded whole-job calibration and a finite completion ceiling remain open.
