# Shared transfer links and host-memory load in the compact JIT floor

The compact H2D and D2H floors now accept the same `shared_links` entries as
the finite execution graph: each entry declares `devices`,
`h2d_bytes_per_second` and `d2h_bytes_per_second`. Multiple entries can
overlap to represent a switch and its upstream root. For each direction, the
floor divides the sum of mandatory bytes for devices on each link by that
link's declared capacity. The largest link constraint joins the existing
global and per-device transfer constraints. Link times are not added: a byte
can traverse nested links while other devices transfer concurrently. The
caller must supply independently supported service ceilings; no topology or
live availability is inferred from GPU names.

`native_layout_partial_envelope` requires matching link declarations in the
compute and output reports. It also adds a single shared host-memory
constraint:

`(source decoder/read DRAM work + genotype H2D bytes + result D2H bytes)
 / declared shared DRAM capacity`.

The native decoder's logical traffic includes writing the unpacked host
genotype buffer; the H2D transfer later reads that buffer, and D2H writes a
separate host result buffer. These are distinct service demands on the same
memory system. The source traffic is an existing model approximation, so the
combined value remains a conditional necessary floor, not a physical-memory
counter or an elapsed-time upper bound. Output occupancy gives an interval
for the D2H term in significant-pair mode.

An illustrative two-GPU fixture with equal per-device work gives an H2D
floor of 5.16 time units when both GPUs share a 300-byte/s link and 2.58
when each has its own such link. Its D2H floor likewise changes from 27 to
13.5 with a 20-byte/s shared versus separate link. These rates are fixture
values, not A100/H100 measurements or a GPU recommendation. A separate
host-memory test shows the combined floor exceeds the source-only floor even
when GPU links and output capacity are loose.

The final 71-test A100 source/frontier/compute/output batch passed as job
`20260923-000737-1255007`, after correcting test arithmetic and accepting the
source composer's floating-point logical-byte sum. The two earlier jobs
`20260923-000207-1254111` and `20260923-000451-1254592` were development
checks; their failures were in the new assertions and type validation. The
current code and tests were pushed with `proj` before the passing run.

This improves candidate floors for multi-GPU allocation. It does not price
kernel setup, selector or writer service, full-duplex versus half-duplex
contention beyond the declared directional links, live device memory, queued
work or final drain. The frozen H100 absolute prediction gap remains. No
layout or GPU switch uses this floor yet.
