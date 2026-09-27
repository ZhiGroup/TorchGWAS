# Compact dense-writer service for a post-output JIT pass

`native_layout_dense_writer_service_floor` prices the exact compact dense
writer counts already attached to a fixed native PGEN output layout. It uses
the finite writer graph's independent CPU copy, zeroing, executor dispatch,
page-cache, writeback-request, eviction and fsync service fields. The existing
legacy NumPy-copy fallback is accepted but records its missing per-call price;
this is not a newly calibrated service.

The report has two necessary, conditional constraints: total CPU-equivalent
writer work divided by the declared shared CPU capacity, and serial
write/close chains for tiles assigned to each GPU. The separate output floor
already counts durable payload against the shared storage service and D2H
payload against the GPU and shared-link capacities. The partial envelope takes
the maximum of these simultaneous loads; it never adds their times or treats
the result as a completion ceiling. The report binds the current PGEN identity,
each source partition, writer geometry and a digest of the compact output
work. It is computed only after useful output in a productive JIT step.

In the differential fixture, the compact CPU-equivalent count equals the sum
of the expanded `BinaryWriterSchedule` graph's CPU demands, and its service
floor is no greater than that graph's solved completion time. The latest A100
targeted batch, job `20260923-001603-1256153`, passed 416 tests. This tests
equation consistency, not source-matched price quality or end-to-end speed.

The compact report lacks queue stalls, dirty-page coupling, short writes,
current in-flight writer state, final metadata publication and a qualified
conditional upper completion bound. It currently applies only to dense
output. Significant pairs and JAGWAS still need their mode-specific selector
and indexed-writer service in a finite continuation. A measured full-job
control must qualify any switch before the public tuner uses this floor to
select a chunk, tile or GPU layout.
