# NumPy identity for reusable calculator measurements

Saved measurements remain immutable evidence. A new job can reuse them only
when their dependencies still match and their original observation age is
acceptable. Reading a measurement, checking it against new observations, or
deriving another profile does not renew that age. The existing independent
[nonzero branch correction](nonzero_branch_pricing_20260922.md) demonstrates
reusing observations while retaining their artifact hashes and timestamps.

This change closes a compatibility gap in `detailed_calibration.py`: a NumPy
version string does not uniquely identify its compiled implementation or
enabled CPU features. Execution contexts now include the SHA-256 of the loaded
core module's installed binary, its runtime CPU feature map, and its compiled
baseline and dispatch feature lists. The context also captures `NPY_*` and
`NUMPY_*` environment variables. Actual runtime metadata and environment
declarations are separate because feature controls are read at import time.
Changing an environment variable afterward does not change the imported CPU
feature map; the changed declaration still invalidates the previous context.

Only the binary digest can use `BindingDigestCache`. Its existing host/boot,
path, inode, size, mtime and ctime checks apply, including checks before
publication. Runtime CPU features are read again on each context check.
This is ordinary filesystem-change detection, not protection against hostile
metadata manipulation. Without a digest cache the binary bytes are read every
time. Missing or unrecognized runtime CPU metadata rejects profile binding.

The reuse rules remain distinct:

| Parameter class | Reuse condition | What happens when it changes |
| --- | --- | --- |
| Structural work and implementation identity | Matching source, binary, hardware and relevant settings | A different dependency key; preserve the old evidence |
| Measured CPU/GPU/transfer/storage costs | Matching dependencies, original age limit, and any required drift check | Collect independent observations and publish a new record |
| Free memory and current resource contention | Fresh observation for this job | Re-evaluate admission and decisions |

Existing profiles without the NumPy core identity fail context validation.
Their immutable artifacts remain readable but are not automatically certified
under the stronger binding. A newly bound profile needs evidence establishing
the appropriate measurement context; adding today's identity to an old
observation is not a valid migration.

Startup still requires valid supplied evidence. If an optional initial-chunk
planning check detects a change, the existing public policy keeps the current
chunk size and continues the scientific scan. This patch changes context
validation only; it does not change statistical kernels, output, thresholds,
chunk policies or the distinction between diagnostics and independent prices.

Remote A100 job `20260922-122854-1030922` passed **242 targeted tests in 53.37
seconds**. These cover same-version binary replacement, missing metadata,
feature changes, rejection of legacy profiles, cross-job digest reuse,
immutable observation ages, productive tuning checks and resident-copy refresh.
The log is `results/numpy_context_binding_checks_20260922/pytest.txt`.

The same job completed three fresh-process checks using the real NumPy 2.2.6
core and installed libraries. The first process hashed 145 package/library
files in four groups. Both later processes reused all four groups without
reading those file contents, and the four original digest records stayed
byte-identical. The unchanged context passed validation. Disabling AVX2 before
import changed the reported runtime feature and rejected the previous profile.
A post-import environment change left the runtime feature map unchanged but
was rejected as a changed declaration. Reports, hashes and contexts are in
`results/numpy_context_reuse_20260922/`; the script is
`benchmarks/numpy_context_reuse_20260922.py`.

This audit performs no association scan and measures no resource capacity.
Broader calculator runtime accuracy and profitable large-job tuning remain
unqualified.
