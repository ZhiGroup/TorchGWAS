# Sparse output for whole-source JIT layouts

The compact source floor can compare alternate tile and GPU ownership over
millions of unissued variants, but reduced output also needs the same
survivors in both layouts. `whole_layout_survivors` supplies an explicit
conditional scenario on a global coordinate: significant pairs use
`variant * total_traits + trait`, and JAGWAS uses `variant`. Empty and
dense scenarios are exact endpoints; sparse scenarios use the declared
rational retained fraction and `spread` or `clustered` placement.
Candidate rectangles count the same global pattern after retiling or
variant sharding. JAGWAS still requires the complete phenotype panel in
every partition.

Spread counts use integer floor sums, so a huge phenotype/variant rectangle
does not expand into individual pairs or markers. Clustered placement uses
its periodic variant pattern and an explicit visit budget. The report gives
one exact conditional retained count per partition, suitable for
`native_layout_output_floor`; it does not infer future significance from
the first chunks or renew component prices.

The short-window first-chunk comparator now uses the same global-coordinate
ledger in its fine bins. A full finite continuation must carry those counts
through every candidate and output-service component. It still needs
independently priced selector and
archive service, in-flight/writer state, a completion ceiling, and a
qualified calibration before it can authorize a public JIT switch.

The A100 global-ledger and output regression passed 262 tests in
`20260923-021817-1322713`. The short-window implementation then reused
one occupancy ledger across host/capacity scenarios; 201 focused JIT and
model tests passed in `20260923-022013-1325839`. These validate
conditional accounting and exact rebinning, not future selectivity.
