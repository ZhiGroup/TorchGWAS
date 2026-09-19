# TorchGWAS runtime predictor and pipeline planner

## Purpose and scope

This document explains the performance models developed alongside TorchGWAS from first principles. It is intended for readers who know basic programming and algebra but do not yet have a background in performance modeling.

The project contains two related tools:

1. **The cross-tool runtime predictor** estimates the runtime of TorchGWAS, PLINK 2.0, or fastGWA-lr for a stated workload and hardware profile. It can also search for a workload at which the predicted ordering of two tools changes. The top-level entry point is `runtime_calculator.py`.
2. **The TorchGWAS pipeline planner** estimates resource use for TorchGWAS itself and searches for a feasible chunk size, queue depth, and worker configuration. Its implementation is `src/torchgwas/pipeline_model.py`.

These tools answer different questions. The cross-tool predictor asks, “Under these assumptions, which program is predicted to finish first?” The pipeline planner asks, “How should TorchGWAS divide and schedule this workload on this machine?” A memory-feasible plan is not automatically the fastest plan, and a predicted crossover is not a universal performance guarantee.

The public TorchGWAS package currently ships the pipeline planner. The complete cross-tool predictor and its frozen hardware profiles remain research-analysis components used for the manuscript and are not yet exposed as a supported public command-line interface. File names for those components are retained below so the analysis can be audited against the research archive; their appearance here does not mean that they are installed by the public package.

The calculators are research tools. They expose their assumptions and return conditional estimates. They must not be described as exact simulators, fitted timing curves, confidence intervals, or guarantees that TorchGWAS will be faster for every dataset above one workload threshold.

## 1. The central idea

A program uses resources to perform work. Runtime prediction therefore has two parts:

1. Count the work required by the source code.
2. Estimate how quickly the target hardware can service that work.

For example, reading 20 GB from storage at a sustained rate of 2 GB/s requires at least

```text
20 GB / 2 GB/s = 10 s
```

Similarly, performing 10 trillion floating-point operations at an effective rate of 5 trillion operations/s requires at least 2 s of GPU service.

Real programs are more complicated because resources can operate at the same time. While the GPU processes one genotype block, CPU workers can decode a later block and a writer can save an earlier result block. Consequently, stage times usually cannot be added directly. The calculator represents the finite sequence of blocks and their dependencies so that overlap, startup, queue limits, backpressure, and final drain are included.

An accessible analogy is an assembly line:

- storage reads a block;
- CPU workers decode it;
- the block is transferred to the GPU;
- the GPU calculates associations;
- results return to the host;
- a writer saves the results.

The throughput of a long assembly line is governed mainly by its slowest repeatedly used station. A short run also pays for filling and draining the line. If a downstream station is full, an upstream station must wait. This is why the calculator needs a finite schedule rather than only a table of independent stage times.

## 2. Notation

The calculator uses the following dimensions.

| Symbol | Meaning |
| --- | --- |
| `N` | Number of analyzed samples |
| `M` | Number of tested variants |
| `K` | Number of analyzed quantitative phenotypes |
| `C` | Number of non-intercept covariates |
| `B` | Number of variants in one compute chunk |
| `D` | Pipeline depth, or the number of in-flight buffer slots |

The manuscript writes products such as `MN` and `KMN`. Multiplication is commutative, so `MN = NM`, but using one order consistently makes the text easier to read. In this project, `MN` denotes the number of sample-variant values and `KMN` denotes the number of phenotype-variant-sample combinations.

In code, `C` is normally a count. In the statistical methods, the same letter may denote a covariate matrix. A document that includes both topics should explicitly distinguish “`C`, the covariate matrix” from “`C`, the number of covariates,” or use a different symbol for one of them.

## 3. What is measured, counted, modeled, and validated

The calculator keeps four kinds of evidence separate.

### 3.1 Source-derived counts

These are quantities obtained from source code, file structure, and workload dimensions. Examples include:

- numbers of floating-point operations;
- bytes read from a genotype file;
- bytes transferred between host and GPU;
- result bytes written;
- numbers of genotype records in different PGEN encodings;
- numbers of missing-genotype branches;
- queue slots and chunk counts.

These are not timing measurements.

### 3.2 Independently measured resource rates

These characterize the machine without timing a complete GWAS. Examples include:

- sustained storage read and write bandwidth;
- host-memory bandwidth;
- CPU parsing and genotype-decoding throughput;
- host-to-device and device-to-host transfer rates;
- effective GPU throughput for the relevant kernels;
- kernel launch and host-dispatch service.

The unit is part of the value. Storage and memory rates are in bytes/s, compute rates are in operations/s, and small services may be recorded in seconds/call or CPU-seconds/record.

### 3.3 Modeled runtime

The calculator combines work counts, resource rates, and scheduling dependencies. This produces a conditional estimate. It is conditional because the answer changes if CPU availability, storage caching, output format, software version, or data encoding changes.

### 3.4 End-to-end validation measurements

Complete program timings are held out from model construction and used to evaluate the prediction. They are not silently converted into correction factors for the general calculator. A model can have the correct qualitative crossover while still making substantial numerical errors in the exact runtime.

The environment-ready validation bundle contains 72 runs across three servers, four workloads, three methods, and two repetitions. Depending on server and method, mean absolute percentage error ranged from 16.6% to 59.4%. These errors are model limitations and must remain visible.

## 4. Repository map

The most important files are:

| File | Role | Availability |
| --- | --- | --- |
| `runtime_calculator.py` | Top-level command dispatcher | Research archive |
| `src/torchgwas/first_principles.py` | Auditable source-work ledger | Research archive |
| `src/torchgwas/first_principles_cli.py` | `work` command-line interface | Research archive |
| `src/torchgwas/mechanistic_torch.py` | TorchGWAS runtime candidate | Research archive |
| `src/torchgwas/mechanistic_plink.py` | PLINK 2.0 runtime candidate | Research archive |
| `src/torchgwas/mechanistic_cpu.py` | fastGWA-lr runtime candidate | Research archive |
| `src/torchgwas/execution_graph.py` | Finite dependency scheduling | Research archive |
| `src/torchgwas/binary_schedule.py` | Bounded binary-writer scheduling | Research archive |
| `src/torchgwas/crossover_scenario.py` | Explicit synthetic workload census | Research archive |
| `src/torchgwas/crossover.py` | Discrete crossover search | Research archive |
| `src/torchgwas/crossover_cli.py` | `crossing` command-line interface | Research archive |
| `src/torchgwas/pipeline_model.py` | TorchGWAS resource and plan search | Public package |

The principal frozen validation bundle is:

```text
paper/calculator_environment_ready_20260917/
```

Its `inputs.json` file contains the frozen resource profiles and source-count inputs used for the three-server comparison. Its `report/` directory contains the predicted-versus-measured accuracy tables and crossover report.

Several older files document intermediate calculators. Their retained results are useful audit history, but old commands such as `ties`, `marker-ties`, `nm-ties`, or `resource` are not the current top-level interface. The current entry-point modes are `work`, `candidate`, and `crossing`.

## 5. Exact time decomposition

This section gives the equations implemented by the code. It is intentionally more detailed than a normal methods description. A reader should be able to start from a workload census and a hardware profile, follow the equations, and reconstruct the reported runtime.

For every composite equation, the explanation should identify each symbol, decompose each grouped term into concrete arrays or operations, state the units, and explain why terms are added, maximized, or multiplied. A compact expression is not a substitute for this accounting.

There are two generations of runtime engine that must not be confused:

1. **The current source/resource candidate engine** in `src/torchgwas/mechanistic_*.py` is called by `runtime_calculator.py candidate` and `crossing`. The regenerated Figure 2 boundary uses this engine.
2. **The retained component predictor** in `benchmarks/direct_marker_ties.py`, `benchmarks/direct_nm_model.py`, and `_scratch/plot_panel_b_plink48_projected.py` generated an earlier Figure 2 draft. It remains audit history but is no longer the figure source.

The two engines share the same scientific intent, but their inputs and scheduling abstractions differ. Sections 5.3–5.5 give the exact current equations, and Section 5.7 states the figure-specific resource scenario and boundary search.

### 5.1 The universal work-to-time rule

For work amount `W` and available service rate `R`, the code uses

```text
service(W, R) = 0,       if W = 0
                W / R,  if W > 0 and R > 0
                infinity, if W > 0 and R = 0.
```

`W` is a non-negative amount of required work and `R` is the available capacity in matching work units per second. Zero work requires zero time. Positive work at positive capacity uses the usual work/rate relation. Positive work with zero available capacity is impossible and is represented by infinity instead of silently returning zero or dividing by zero.

The units must cancel. Examples are:

```text
bytes / (bytes/s) = s
FLOPs / (FLOPs/s) = s
calls × (s/call) = s
records × (CPU-s/record) = CPU-s.
```

Let `q` be the scheduling fraction available to one serial CPU task, with `0 < q <= 1`. A task requiring `U` CPU-seconds has conditional wall service

```text
T_serial = U / q.
```

`U` is CPU service measured in CPU-seconds at one fully scheduled logical CPU. `q` is the fraction of that logical CPU available to this serial task. The task cannot usefully spread over multiple cores, so reduced scheduling share lengthens wall time inversely. `q=1` leaves the service unchanged; `q=0.5` doubles it.

Let `C_cpu` be the total CPU capacity available to the whole process, expressed in cores. A collection of independent CPU tasks requiring total CPU work `U_total` cannot complete faster than

```text
T_CPU_floor = U_total / C_cpu.
```

`U_total` is the sum of CPU-seconds required by all CPU nodes, including nodes that may run concurrently. `C_cpu` is the total core-equivalent capacity available to the process. Even with perfect parallel scheduling, the run cannot complete faster than total CPU work divided by total CPU capacity. This is a whole-run lower bound and does not replace serial dependencies between particular nodes.

The serial-task term and whole-process CPU floor serve different purposes. The first represents one task being scheduled intermittently. The second prevents many simultaneous tasks from consuming more CPU capacity than the process actually has.

### 5.2 Sum versus maximum

The most important modeling decision is whether two services occur sequentially or can overlap.

Use a **sum** when the same execution path performs operations one after another:

```text
T_sequential = T_1 + T_2 + ... + T_j.
```

Every `T_h` is the service time of one required phase. The sum is used only when phase `h+1` cannot start before phase `h` finishes. No overlap credit is permitted in this expression.

Examples are TorchGWAS covariate-basis construction followed by pinned-buffer allocation, or the fill, matrix operation, and regression steps executed sequentially within one PLINK worker.

Use a **maximum** for independent capacity constraints that must all be satisfied during the same interval:

```text
T_capacity = max(T_CPU, T_DRAM, T_storage, T_GPU, ...).
```

Each term is the time the same piece of work would require if constrained by one resource. When those demands can be served concurrently, completion must wait for the slowest resource, hence the maximum. Using a maximum does not assert overlap unless the execution graph permits it.

This does not mean that all programs are represented by one global `max`. Exact program ordering is represented by a finite dependency schedule. The `max` is used only at a point where the model explicitly permits overlap or imposes a whole-run resource floor.

For a CPU or GPU micro-stage that can be limited by instructions, arithmetic, or data movement, the code uses

```text
T_stage = max(T_instruction, T_arithmetic, T_memory).
```

The instruction, arithmetic, and memory terms are alternative capacity limits on the same micro-stage. They are not three sequential operations. The largest service determines the roofline estimate; adding them would assume complete serialization and double-count work that modern execution overlaps.

For a cache-modeled CPU stage,

```text
T_memory = max(
    bytes_L1   / rate_L1,
    bytes_L2   / rate_L2,
    bytes_L3   / rate_L3,
    bytes_DRAM / rate_DRAM
).
```

`bytes_L1`, `bytes_L2`, `bytes_L3`, and `bytes_DRAM` are the traffic assigned to each cache or memory level by the declared cache scenario. Each is divided by the corresponding sustained rate. These levels participate in servicing the same stage, so the slowest level supplies the memory bound. The equation is conditional on the cache-traffic census; it is not valid if all logical bytes are charged independently at every level.

The stages within one variant are then summed because that worker executes them sequentially.

For one modeled GPU kernel group,

```text
active_SM_fraction = min(1, grid_blocks / SM_count)

T_math = issued_FLOPs /
         (FP32_rate × active_SM_fraction × GPU_fraction)

T_memory = max(
    HBM_bytes / (HBM_rate × GPU_fraction),
    L2_bytes  / (L2_rate  × GPU_fraction)
)

T_kernel = kernel_count × kernel_launch_seconds / GPU_fraction
           + max(T_math, T_memory).
```

`grid_blocks/SM_count` estimates the fraction of streaming multiprocessors that can be occupied by the compiled grid, capped at one. `issued_FLOPs` counts padded instructions actually launched, not only useful matrix entries. `FP32_rate × active_SM_fraction × GPU_fraction` is the available arithmetic capacity. HBM and L2 byte counts are divided by their available bandwidths, and the slower supplies the memory limit. Kernel launches are serial host/device overheads, so `kernel_count × kernel_launch_seconds` is added to the arithmetic-or-memory roofline. `GPU_fraction` represents the fraction of the GPU available to this process.

Kernel groups are serialized in the eager Torch stream. Host dispatch may overlap earlier GPU execution, so the start of a kernel is

```text
kernel_start = max(previous_GPU_finish, corresponding_host_dispatch_finish)
kernel_finish = kernel_start + T_kernel.
```

The kernel cannot begin before both the preceding GPU work and its own host dispatch have finished, so its start is their maximum. Once started, its modeled service `T_kernel` is added to obtain completion. Applying these recurrences to every kernel preserves stream ordering while allowing host dispatch to overlap earlier device execution.

### 5.3 Exact TorchGWAS top-level decomposition

The current source/resource candidate is `torch_runtime()` in `src/torchgwas/mechanistic_torch.py`. Its supported validation branch is `K = 1`, eight covariates, complete phenotype data, and matching sample order.

The final equation is a top-level **sum**:

```text
T_TorchGWAS = T_CUDA_context
            + T_metadata_QC_basis
            + T_pinned_allocation
            + T_GPU_design
            + T_first_use
            + T_scan_binary_close
            + T_sidecars_exit.
```

The terms are sequential at the outermost level. `T_CUDA_context` initializes CUDA for a fresh analytical run. `T_metadata_QC_basis` parses and validates tables and constructs the covariate basis. `T_pinned_allocation` allocates page-locked host rings. `T_GPU_design` residualizes and standardizes phenotypes and constructs the resident design. `T_first_use` accounts for lazy numerical-library initialization. `T_scan_binary_close` is the complete finite overlapped genotype-read, decode, transfer, association, significance, and binary-write graph. `T_sidecars_exit` publishes identifiers and small metadata after the dense streams. The internal nodes of `T_scan_binary_close` must not be added again at this level.

These terms are sequential phases at the process level. `T_scan_binary_close` is itself a finite overlapped schedule and must not be replaced by a sum of its internal stage clocks.

#### 5.3.1 Cold-start GPU initialization

The paper-facing calculator represents the analytical work of a fresh, cold-input run. Python imports, shell or process bootstrap, thread-pool creation, and unrelated environment setup are not assigned runtime terms. CUDA-context creation and lazy GPU-library initialization remain inside the modeled run:

```text
T_CUDA_context = cuda_context_seconds.
```

Here, “cold” means a fresh run whose genotype input pages are not assumed to have been warmed by a previous association run. Genotype-file reading and decoding therefore remain part of the scan pipeline even though general environment-startup cost is omitted.

#### 5.3.2 Metadata reads

Metadata-file reads are small relative to the genotype payload and are not assigned a separate runtime term in the paper-facing calculator:

```text
T_metadata_IO = 0.
```

This omission applies only to the metadata read itself. Metadata parsing, sample matching, phenotype and covariate checks, and covariate-basis construction remain modeled below.

#### 5.3.3 Metadata parsing, quality control, and covariate basis

The CPU work is

```text
U_setup = u_pvar_row × M
        + u_psam_row × N
        + u_pheno_QC_cell × N × K
        + u_covar_QC_cell × N × C
        + u_basis × (N C^2 + C^3)
        + u_pgen_index × M × initial_index_parses
        + u_numpy_copy_byte × [12N(K+C) + 12NC].

T_metadata_QC_basis = U_setup / q.
```

Each `u_*` is an independently supplied CPU service per row, cell, record, or copied byte.

The dimensions are `N` samples, `M` variants, `K` analyzed phenotype columns, and `C` retained covariate columns. `initial_index_parses` is the number of complete PGEN-index traversals before the decoder-worker pipeline begins; the frozen profiles use two, one while establishing the source scope and one for the initial reader.

Every `u_*` coefficient has units of **CPU-seconds per unit**. Each is the median of independent microbenchmarks of the corresponding source routine on the target host, not a coefficient fitted from a complete GWAS. Multiplication by a row, cell, record, or byte count therefore produces CPU-seconds.

For example, the frozen H100-host profile used five repetitions per primitive and retained these medians:

| Coefficient | Frozen value | Unit and measured primitive |
| --- | ---: | --- |
| `u_pvar_row` | `9.26409e-7` | CPU-seconds per `.pvar` row parsed |
| `u_psam_row` | `3.54621e-7` | CPU-seconds per `.psam` row parsed |
| `u_pheno_QC_cell` | `1.73981e-8` | CPU-seconds per phenotype cell checked |
| `u_covar_QC_cell` | `1.27679e-8` | CPU-seconds per covariate cell checked |
| `u_basis` | `2.92666e-9` | CPU-seconds per `NC^2+C^3` basis-work unit |
| `u_pgen_index` | `9.81067e-9` | CPU-seconds per PGEN index record parsed |
| `u_numpy_copy_byte` | `3.13271e-11` | CPU-seconds per logical NumPy copy byte |

These numbers are machine- and software-build-specific. The equation structure and source counts can be reused, but another host must remeasure the primitive rates. The full setup expression is therefore **source-informed and mechanistic**, but not purely derivable from a CPU's advertised clock rate: Python, NumPy, allocation, branching, and library calls require measured effective service unless a complete instruction/cache model with guaranteed bounds is available.

1. **Variant metadata parsing**

   ```text
   U_pvar = u_pvar_row × M.
   ```

   `u_pvar_row` is the measured CPU time per `.pvar` row for `_read_pvar()`. The routine locates and validates the header, reads `CHROM`, `POS`, `ID`, `REF`, and `ALT`, converts `POS` to an integer, preserves the lexical fields, and constructs the returned arrays. There is one row per variant, hence `M`. This is parsing and object-construction work after the bytes are available; metadata storage-read time is excluded. The generic microbenchmark contains ordinary well-formed biallelic rows, so unusually long identifiers and malformed-row branches are not separately priced.

2. **Sample metadata parsing**

   ```text
   U_psam = u_psam_row × N.
   ```

   `u_psam_row` is the measured CPU time per `.psam` row for `_read_psam()`. The routine splits fields, extracts `IID` and optional `FID`, checks row width, verifies that identifiers are present and unique, and constructs the ID arrays. There is one row per sample, hence `N`. Metadata storage-read time is again excluded.

3. **Phenotype quality control**

   ```text
   U_pheno_QC = u_pheno_QC_cell × N × K.
   ```

   The phenotype matrix contains `NK` cells. `_phenotype_column_mask()` checks for invalid infinities, identifies observed and missing values, counts observations, computes each phenotype's observed-value mean and centered sum of squares, and retains columns with sufficient observations and non-zero variance. `u_pheno_QC_cell` is CPU-seconds per examined cell. The `NK` count normalizes the complete column-wise routine; it does not imply one machine instruction per cell.

4. **Covariate quality control**

   ```text
   U_covar_QC = u_covar_QC_cell × N × C.
   ```

   The covariate matrix contains `NC` cells. `column_std_mask()` determines which covariate columns have non-zero variance; the surrounding preprocessing also requires finite values and aligned sample rows. `u_covar_QC_cell` is CPU-seconds per examined covariate cell. Constant columns are removed before the basis is used.

5. **Covariate-basis construction**

   ```text
   U_basis = u_basis × (N C^2 + C^3).
   ```

   `_covariate_basis()` centers and scales the retained `N × C` matrix, performs a thin singular-value decomposition, determines its numerical rank, and returns an orthonormal basis for the covariate column space. For a tall matrix, dense basis construction has leading work proportional to `NC^2`; factorization and rank work on the small `C × C` problem contributes a term proportional to `C^3`. The sum is a **work normalizer**. `u_basis` is obtained by timing the complete basis routine and dividing by this normalizer, so `NC^2 + C^3` is not claimed to be an exact NumPy/LAPACK FLOP count.

6. **Initial PGEN-index parsing**

   ```text
   U_initial_index = u_pgen_index × M × initial_index_parses.
   ```

   The PGEN index contains information for each of the `M` variant records. `u_pgen_index` is measured CPU-seconds per index record for `read_header()`. One complete traversal costs `u_pgen_index × M`; the final multiplier counts the complete traversals performed before streaming. Additional decoder-worker index initialization is charged separately to each worker's first block in Section 5.3.7.

7. **Explicit dtype-conversion copies**

   ```text
   U_array_copy = u_numpy_copy_byte × [12N(K+C) + 12NC].
   ```

   `u_numpy_copy_byte` is measured CPU-seconds per logical byte copied by `numpy.copyto()`. The two byte counts are:

   - `12N(K+C)` for converting the loaded phenotype and covariate arrays from float64 to the float32 working representation. Each of the `N(K+C)` elements entails an 8-byte read and a 4-byte write, giving `8 + 4 = 12` logical bytes.
   - `12NC` for promoting the float32 covariate matrix back to float64 for the stable rank decision inside `_covariate_basis()`. Each covariate element entails a 4-byte read and an 8-byte write, again giving 12 logical bytes.

   Covariates occur in both terms: first in the common float64-to-float32 conversion and then in the float32-to-float64 promotion. These are logical memory-copy bytes, not bytes reread from the genotype file. Internal basis arithmetic and workspace are represented by `U_basis` rather than enumerated as copies here.

The seven components are sequential setup work and are added:

```text
U_setup = U_pvar + U_psam + U_pheno_QC + U_covar_QC
        + U_basis + U_initial_index + U_array_copy.
```

`U_setup` is CPU service, not wall time. The path is serial, so the calculator divides by the scheduling fraction `q` available to this task:

```text
T_metadata_QC_basis = U_setup / q.
```

At `q = 1`, wall time equals CPU service. At `q = 0.5`, the task receives half of one logical CPU on average and its modeled wall time doubles.

#### 5.3.4 Pinned allocation

For chunk size `B` and pipeline depth `D`, the current candidate allocates pinned input and result rings of

```text
pinned_bytes = D × [N B + (16K + 5)B].
```

This is the allocated size of the **pinned host-memory rings**. It is not GPU memory and is not the process's total host-memory use. Here `D` is the pipeline depth, or number of reusable chunk slots that may be in flight; `B` is the allocated variant capacity of one slot; `N` is the sample count; and `K` is the number of phenotypes analyzed together.

Each input slot holds one decoded genotype block with shape `N × B` in signed 8-bit form. One byte is stored for each sample-variant value:

```text
input_bytes_per_slot = N × B × 1 byte = NB.
```

Each result slot holds five outputs for the same `B` variants:

```text
coefficient array: B × K float32 values = 4BK bytes
t-statistic array: B × K float32 values = 4BK bytes
-log10(P) array:   B × K float64 values = 8BK bytes
status array:      B uint8 values       = B bytes
residual-df array: B float32 values     = 4B bytes.
```

The coefficient and t-statistic matrices contribute `8BK` bytes. The underflow-safe significance calculation is performed in float64 on the GPU and remains float64 during device-to-host transfer, contributing another `8BK` bytes. It is cast to float32 only when the binary writer constructs `neglog10p.f32`. The per-variant status and residual degrees of freedom contribute `B + 4B = 5B` bytes. Therefore,

```text
result_bytes_per_slot
    = 4BK + 4BK + 8BK + B + 4B
    = (16K + 5)B.
```

One complete slot consequently occupies

```text
bytes_per_slot = NB + (16K + 5)B.
```

There are `D` simultaneously allocated slots, giving the stated expression. The formula uses the configured capacity `B`, rather than the number of variants in the final partial chunk, because every reusable slot is allocated at the maximum chunk size.

The operating system pins memory by page, so the next step rounds the allocation upward to a whole number of pages:

Let

```text
pages = ceil(pinned_bytes / 4096).
```

Linux pins complete 4,096-byte pages. Division gives the exact page-equivalent count and `ceil` rounds any partial final page upward; even one byte on that page requires the whole page to be pinned.

Then

```text
T_pinned_allocation = pages × pin_CPU_seconds_per_page / q
                    + pages × pin_driver_seconds_per_page.
```

`pin_CPU_seconds_per_page` is the CPU service for bookkeeping one newly pinned page and is divided by the serial CPU scheduling fraction `q`. `pin_driver_seconds_per_page` is the remaining elapsed driver service per page and is not divided by `q`. Both occur for every page and are sequential contributions, so they are added.

#### 5.3.5 GPU design setup

This term covers phenotype residualization and standardization on the GPU, construction of the resident phenotype-plus-covariate design, and the transfers needed to build that design. The source primitive was measured once at the small reference shape

```text
N_0 = 32, C_0 = 8, K_0 = 1.
```

For the supported candidate, `C = 8` and `K = 1`. The model avoids charging the first 32 rows twice by defining

```text
dN = N - N_0 = N - 32.
```

Because the candidate refuses `N < 32`, `dN` is non-negative. The fixed reference-shape service is split into the CPU time consumed while issuing and preparing the work and the remaining elapsed service:

```text
T_tiny_CPU     = tiny_setup_CPU_seconds / q
T_tiny_non_CPU = tiny_setup_non_CPU_seconds / GPU_fraction.
```

`tiny_setup_CPU_seconds` and `tiny_setup_non_CPU_seconds` are empirical primitive measurements of the complete `N=32, C=8, K=1` setup call. They are not derived from an instruction census, so this part of the candidate is hybrid rather than strictly first-principles. Division by `q` prices reduced CPU scheduling availability; division by `GPU_fraction` prices reduced GPU availability.

The additional projection arithmetic beyond the 32-row reference is

```text
F_projection = 4 dN C K FLOPs.
```

Residualizing `K` phenotypes against a rank-`C` basis performs two dense matrix products. The product `Q^T Y` costs approximately `2dNCK` FLOPs and `Q(Q^T Y)` costs another `2dNCK`, producing the factor four.

The modeled additional device-memory traffic is

```text
H_design = 4 dN (6C + 18K + 4) bytes.
```

The leading `4` is the size of a float32 value. The ledger assigns six basis/covariate-related value movements per added sample and covariate, eighteen phenotype/intermediate value movements per added sample and phenotype, and four additional values per sample for intercept and small design auxiliaries. These coefficients aggregate logical reads and writes across centering, projection, standardization, concatenation, design copying, and sum-of-squares construction. They are a source-informed logical-traffic approximation, not a claim about exact HBM transactions after caching or kernel fusion. The current source does not preserve a finer line-by-line derivation of the `6`, `18`, and `4`; therefore the guide must not present them as exact hardware traffic.

Arithmetic and HBM traffic constrain the same GPU setup work. The roofline service is consequently the slower of the two:

```text
T_extra_GPU = max(
    F_projection / GPU_FP32_rate,
    H_design / HBM_rate
) / GPU_fraction.
```

The host-to-device ledger is

```text
B_H2D_design = 4 dN (2C + 3K) bytes
T_H2D_design = B_H2D_design / H2D_rate.
```

The factor four again denotes float32. The ledger counts two covariate-sized transfers and three phenotype-sized transfers for the residualization/design path. This is a logical transfer census used by the candidate; it should be re-audited whenever the preprocessing implementation changes.

The residualized and standardized phenotype is returned once to host memory before being consumed by the subsequent design-building path:

```text
B_D2H_design = 4 dN K bytes
T_D2H_design = B_D2H_design / D2H_rate.
```

All five services occur in the setup phase and are added:

```text
T_GPU_design = tiny_setup_CPU_seconds / q
             + tiny_setup_non_CPU_seconds / GPU_fraction
             + max(
                   4 dN C K / GPU_FP32_rate,
                   4 dN (6C + 18K + 4) / HBM_rate
               ) / GPU_fraction
             + 4 dN (2C + 3K) / H2D_rate
             + 4 dN K / D2H_rate.
```

Equivalently,

```text
T_GPU_design = T_tiny_CPU + T_tiny_non_CPU
             + T_extra_GPU + T_H2D_design + T_D2H_design.
```

The `max` is used only inside `T_extra_GPU` because arithmetic and device-memory throughput limit the same GPU operations. Transfers and fixed setup services are added because the implementation performs them as separate required phases rather than overlapping them in the steady-state genotype pipeline.

#### 5.3.6 Lazy first use

```text
T_first_use = first_use_seconds.
```

This independently measured scalar represents lazy library initialization that occurs once before steady repeated GPU work. It is added once, not once per block. Because it is an empirical elapsed primitive rather than a source-derived work/rate quotient, it is another explicitly identified hybrid term.

#### 5.3.7 Per-block decoder and transfer work

For block `i`, let `b_i = min(B, M - iB)` and `f_i = b_i/M`.

`iB` is the number of variants preceding block `i`. `M-iB` is therefore the number still unprocessed. Taking the minimum with the configured capacity `B` gives a full block except at the tail. `f_i` is the fraction of all variants represented by this block; aggregate source-work counts are distributed to blocks using this fraction under the explicit uniform-work scenario.

The block’s identified decoder CPU work is

```text
U_decode,i = U_decode,total × f_i
           + I(i < W) × M × u_pgen_index,
```

where `W = min(D, decode_workers, ceil(M/B))` and `I(condition)` is 1 when the condition is true. The second term initializes one reader index for each decoder worker’s first job.

`U_decode,total × f_i` assigns block `i` its proportional share of the complete input's identified decoder CPU service. `ceil(M/B)` is the number of blocks, so `W` cannot exceed the number of blocks, available decoder workers, or pipeline slots. For the first `W` blocks, `I(i<W)=1`; each starts a distinct worker and charges one complete `M`-record index initialization. Later blocks reuse an initialized worker and the indicator is zero.

The modeled decoder memory traffic is

```text
H_decode,i = 2 f_i × record_payload_bytes
           + 2 ceil(N/4) × b_i
           + 2 N × b_i.
```

This is logical host-memory traffic, not file bytes. `f_i × record_payload_bytes` is the encoded record payload assigned to the block; the factor two represents a staging write and subsequent read. Hard calls require `ceil(N/4)` packed bytes per variant because four two-bit calls fit in one byte; multiplying by `b_i` gives the block's packed staging size, and the factor two again represents write plus read. The decoder produces an `N × b_i` int8 matrix, one byte per cell; its write and subsequent consumer read contribute `2Nb_i`. Cache reuse may reduce physical DRAM traffic, so this is a conservative logical census.

Decoder wall service is the larger CPU or memory requirement:

```text
T_decode,i = max(
    U_decode,i / q,
    H_decode,i / shared_DRAM_bytes_per_second
).
```

`U_decode,i/q` is the wall service implied by the serial decoder CPU work at scheduling fraction `q`. `H_decode,i/shared_DRAM_bytes_per_second` is the wall service implied by moving the logical bytes through host memory. They constrain the same decoder stage, so a roofline maximum is used rather than adding them.

Input-read service is

```text
T_read,i = f_i × record_payload_bytes / read_bytes_per_second.
```

Under the uniform-record scenario, block `i` owns fraction `f_i` of the encoded payload. Dividing those file bytes by sustained input bandwidth gives its storage service. Metadata bytes are not included in this term.

Host-to-device service is

```text
T_H2D,i = N × b_i / H2D_bytes_per_second.
```

The decoded hard-call block contains `Nb_i` signed 8-bit values, hence `Nb_i` transferred bytes. Division by sustained host-to-device bandwidth gives transfer service. A float32 dosage representation would require four times as many bytes and must use a different census.

Returned result bytes and device-to-host service are

```text
result_bytes_i = (16K + 5)b_i
T_D2H,i = result_bytes_i / D2H_bytes_per_second.
```

For each of the `b_i K` association cells, device-to-host transfer contains a 4-byte coefficient, a 4-byte t statistic, and an 8-byte float64 `-log10(P)`, giving `16K b_i` bytes. Each variant additionally carries a 1-byte status code and a 4-byte residual degrees-of-freedom value, giving the remaining `5b_i` bytes.

The CPU result-finishing service is

```text
U_finish,i = u_finish_fixed
           + max(0, result_bytes_i - 21×32) × u_numpy_copy_byte

T_finish,i = U_finish,i / q.
```

The fixed finishing primitive is defined at 32 variants and one phenotype. Under the new dense-output contract, each reference variant carries `16×1+5=21` returned bytes, giving the baseline `21×32`. Bytes above that baseline are charged at the independently measured NumPy copy service. The `max` prevents a negative incremental charge for a smaller tail block. `T_finish,i` divides CPU service by the serial task's scheduling fraction.

The frozen `u_finish_fixed` measurement predates the float64 `-log10(P)` result and must be remeasured with the five-array finish path before this updated expression is used numerically. Merely changing `13×32` to `21×32` corrects the byte algebra but does not update that empirical fixed primitive.

#### 5.3.8 Per-block GPU kernel work

The eager tensor graph is obtained from `eager_statistics_work()`. For every active tensor step, the candidate records:

- input and output tensor shapes and storage aliases;
- logical bytes;
- modeled L2 and HBM bytes under the declared ideal-LRU cache;
- operation count;
- compiled kernel grid and block geometry;
- host API service.

For the matrix multiplication, useful arithmetic is

```text
F_useful = 2 N b_i (K + C + 1).
```

The matrix product multiplies a block containing `b_i` variants by a resident design of width `K + C + 1`: `K` phenotype columns, `C` covariate-basis columns, and one intercept column. Each output element reduces across `N` samples. A length-`N` dot product performs `N` multiplications and approximately `N` additions, counted as `2N` FLOPs. There are `b_i(K+C+1)` output elements, yielding `2Nb_i(K+C+1)` useful FLOPs.

The runtime model uses issued tile work rather than only useful work:

```text
issued_K = ceil[N / (split_K × tile_K)] × split_K × tile_K
padded_M = grid_M × tile_M
padded_P = grid_P × tile_P

F_issued = 2 × issued_K × padded_M × padded_P.
```

The names in this GEMM equation follow matrix-multiplication terminology and should not be confused with the manuscript's phenotype count `K`. `issued_K` is the padded **reduction dimension**, whose useful length is `N`; `padded_M` is the padded variant dimension, whose useful length is `b_i`; and `padded_P` is the padded design-width dimension, whose useful length is `K+C+1`. `tile_K`, `tile_M`, and `tile_P` are compiled tile dimensions, `grid_M` and `grid_P` are grid extents, and `split_K` is the number of partial reductions along the sample dimension. Ceiling and padding matter because edge tiles issue instructions for a full tile even when some lanes correspond to no useful matrix element. The factor two again counts a multiply and an add.

If split-K is used, workspace writes, reads, and reduction additions are included. The kernel service equation is the `T_kernel` equation in Section 5.2. The block’s GPU operations remain individually represented in the execution graph; their services are not collapsed into one measured per-block coefficient.

The current dense-output source also evaluates one exact two-sided Student-t tail for every association cell:

```text
L_logp,i = b_i K tail evaluations.
```

For each cell, `upper_tail_log10_from_t_torch()` converts the statistic to float64, forms `a=df/2`, `b=1/2`, and `x=df/(df+t^2)`, evaluates three `lgamma` terms and logarithmic front factors, and then evaluates the regularized incomplete beta through a fixed 40-iteration Lentz continued fraction. Every iteration contains an even and an odd rational update. If any cell enters the reflected branch near `t=0`, a second continued fraction and its exponential/logarithmic front are evaluated before `torch.where` selects the appropriate result.

The complete per-block device service must therefore be

```text
T_GPU_block,i = T_OLS_graph,i + T_logp,i.
```

`T_logp,i` cannot defensibly be set to zero. The current source/resource candidate predates dense `-log10(P)` output and does not yet contain instruction, special-function, launch, or logical-memory counts for this tail graph. Consequently, any curve generated from that candidate is provisional for the new output contract until the fixed 40-iteration tensor graph is added to the source-work census and priced with the same device-resource model. This is a known missing term, not an empirical residual to be hidden elsewhere.

#### 5.3.9 Binary output work

The default dense writer stores three float32 arrays: coefficient, t statistic, and `-log10(P)`. Define `I_beta = 1` when coefficients are retained and `0` for the screening-only `t` field set. The number of stored arrays is

```text
arrays = 2 + I_beta
payload_per_array = 4 M K
total_binary_payload = 4(2 + I_beta) M K.
```

The two mandatory arrays are `tstat.f32` and `neglog10p.f32`; `beta.f32` is present by default. Therefore, the default `beta+t` field set has `I_beta=1`, three arrays, and `12MK` payload bytes. The optional `t` field set has `I_beta=0`, two arrays, and `8MK` payload bytes. Headers, the manifest, variant identifiers, and filesystem allocation overhead are separate.

The auto-selected writer block size per array is

```text
first_payload = 4 K × min(M, B)

writer_block_bytes = min(
    16 MiB,
    max(1 MiB, first_payload)
).
```

`min(M,B)` is the number of variants in the first chunk, including the case where the entire scan is smaller than one configured chunk. One stored float32 array contributes `4K` bytes per variant, so `first_payload` is that chunk's payload for one array. The inner `max` prevents blocks smaller than 1 MiB; the outer `min` caps them at 16 MiB. The same block size is used independently for each stored array.

With writer queue depth `Q_w`, allocated and zero-initialized staging memory is

```text
zero_initialization_bytes
    = arrays × (Q_w + 1) × writer_block_bytes.
```

For each array, the writer owns `Q_w` queued blocks plus one block currently being filled. All are allocated and zero-initialized, hence `Q_w+1` blocks per array. Multiplication by `arrays` accounts for the mandatory t-statistic and `-log10(P)` streams and, by default, the coefficient stream.

If a chunk arrives when no partial staging block exists and its payload is at least one writer block, the writer borrows that chunk and avoids a staging copy. Otherwise bytes are copied into a staging block. A copied byte contributes one read and one write, so

```text
staging_logical_memory_bytes = 2 × staging_copy_bytes.
```

Every staged payload byte is read once from the incoming array and written once into a staging block. The logical host-memory traffic is therefore twice the payload copied. Borrowed chunks bypass this copy and contribute zero staging-copy bytes.

Each queued output block has

```text
T_writer_copy = copied_bytes × u_numpy_copy_byte / q
T_writer_write = queued_bytes / write_bytes_per_second.
```

The first equation converts logical copied bytes to CPU service using the independently measured copy cost and then to wall time using `q`. The second divides bytes submitted to storage by sustained output bandwidth. Copying and writing are separate nodes in the finite writer schedule; their overlap is determined there rather than by adding all writer times globally.

The array writers share one aggregate output capacity. A bounded queue can make the scan consumer wait. Partial blocks are queued during close, and the coefficient, t-statistic, and `-log10(P)` streams are each drained and fsynced before the store is complete.

#### 5.3.10 Exact finite scan graph

`torch_scan_schedule()` builds a directed acyclic graph with nodes for:

```text
submit decode -> input read -> decode -> publish
              -> host transfer submission -> H2D
              -> ordered host API dispatch and GPU kernels
              -> result submission -> D2H -> finish
              -> writer append/copy/queue/write -> close/fsync.
```

Dependencies enforce:

- cyclic assignment to `W` decoder workers;
- at most `D` in-flight ring slots;
- serial H2D, GPU-stream, D2H, and ordered result consumption where required;
- bounded output staging and borrowed-block credit;
- ordered close and fsync.

When several active nodes share a resource, the graph uses proportional fluid sharing. For active node `j`, let `d_jr` be its nominal demand rate for resource `r`, and let `C_r` be the resource capacity. Total active demand is

```text
D_r = sum over active j of d_jr.
```

`d_jr` is node `j`'s nominal demand for resource `r`, such as CPU core-equivalents or bytes/s. Summing only currently active nodes gives simultaneous demand `D_r`. Waiting or dependency-blocked nodes contribute nothing.

The progress rate of node `j` relative to its nominal service is

```text
lambda_j = min(
    1,
    C_r / D_r for every resource r used by node j
).
```

`C_r/D_r` is the fraction of requested service that resource `r` can supply when demand exceeds capacity. The leading `1` prevents a node from progressing faster than its own nominal service when capacity is abundant. A node using several resources must accept the most restrictive fraction, hence the minimum across its resources.

The next completion event occurs after

```text
delta = min over active j of remaining_service_j / lambda_j.
```

At progress fraction `lambda_j`, node `j` needs `remaining_service_j/lambda_j` wall seconds to finish. The smallest such time identifies the next completion event. Advancing by this `delta` cannot skip an earlier dependency release or resource-demand change.

After `delta`, finished nodes release their capacity and dependent nodes may start. The process repeats until every node completes. The final graph makespan is

```text
T_scan_binary_close = final graph time.
```

The final graph time is the timestamp of the last required node, including final partial-block queuing, draining, and fsync. It is a makespan, not the sum of node durations; operations on independent resources may overlap when dependencies and capacities allow.

This is where overlap is calculated. The decoder, GPU, transfers, and writer are not added as if they were sequential.

#### 5.3.11 Sidecars and exit

```text
sidecar_bytes = pvar_field_characters + 4096

T_sidecars_exit = M × u_variant_ID_row / q
                  + sidecar_bytes / write_bytes_per_second.
```

`pvar_field_characters` counts the identifier-related characters copied from the variant table; 4,096 bytes provide the manifest/header allowance used by the candidate. `M × u_variant_ID_row` is CPU service for formatting or publishing one identifier row per variant and is divided by `q`. `sidecar_bytes/write_bytes_per_second` is storage service for those small files. CPU construction precedes or accompanies the final write in this outer phase, so the candidate adds them. This term is separate from the dense-array payload and its fsyncs already included in the writer schedule.

### 5.4 Exact PLINK 2.0 decomposition

The current source/resource candidate is `plink_runtime()` in `src/torchgwas/mechanistic_plink.py`. Its validated branch is one analyzed phenotype with matching sample order.

The final equation is

```text
T_PLINK = T_setup + T_scan.
```

PLINK completes table parsing, matching, and fixed regression preparation before entering the variant scan, so these two outer phases are sequential and are added. `T_scan` already includes internal overlap among reading, worker computation, formatting, and writing.

#### 5.4.1 Setup

Let `P = C + 2`, representing intercept, tested genotype, and covariates in the small regression system.

```text
metadata_bytes = sum of all input table bytes

token_characters = sum over tables of
                   (field_characters + number_of_fields)

U_parse = token_characters × u_token16 / 16
        + N(C+1) × u_scan_double
        + M × u_scan_uint

retained_bytes = 64 × total_number_of_fields
T_metadata_memory = retained_bytes / DRAM_rate

matching_characters = 3 × psam_field_characters
U_matching = matching_characters × u_token16 / 16

F_setup = 2N(C+1)(C+2) + 2N(C+1) + 2(C+1)^3.
```

`metadata_bytes` is the sum of PVAR, PSAM, phenotype, and covariate table sizes; its cold-read time is omitted under the selected calculator boundary. `token_characters` counts field characters plus one delimiter allowance per field. `u_token16/16` prices token handling in 16-character units. `N(C+1)` numeric phenotype/covariate entries are charged at `u_scan_double`, and `M` variant-position values at `u_scan_uint`. `retained_bytes=64×total_number_of_fields` is a fixed per-field retained-object approximation; division by DRAM bandwidth converts it to memory service. The matching ledger traverses PSAM identifiers three times, hence `3×psam_field_characters`. In `F_setup`, `C+1` includes the intercept with the covariates and `C+2` additionally includes the tested-genotype column. The first two terms count dense cross-product and centering/scaling work over `N` samples; the cubic term represents the small fixed-system factorization/update. These are source-work approximations, not exact retired-instruction counts.

Therefore

```text
T_setup = [U_parse + T_metadata_memory + U_matching
           + F_setup / CPU_FP64_rate] / q.
```

`U_parse` and `U_matching` are CPU-seconds. `F_setup/CPU_FP64_rate` converts setup FLOPs to seconds. `T_metadata_memory` is already seconds and is treated as serial setup service in this candidate. All four contributions occur before scanning and are added, then divided by the serial CPU scheduling fraction `q`.

Here `T_metadata_memory` is expressed in seconds before division by `q`, matching the code’s serial-CPU setup treatment.
Process/environment startup and metadata-read time are outside this calculator boundary. Parsing, matching, and setup arithmetic remain included, and genotype payload reading remains inside the finite scan schedule.

#### 5.4.2 Missingness branch counts

Let `r_miss` be the independent per-sample missing-call probability. A variant is complete with probability

```text
p_complete = (1 - r_miss)^N.
```

Each sample call is assumed independently missing with probability `r_miss`, so it is present with probability `1-r_miss`. A variant is complete only when all `N` calls are present; independence gives the `N`th power. This is a declared synthetic-data assumption, not an empirical claim that real missingness is independent.

Let `S_restart` be the number of worker/block restart segments. The expected branch counts are

```text
M_gram   = M(1 - p_complete)

M_opener = S_restart p_complete
         + (M - S_restart)p_complete(1 - p_complete)

M_sparse = (M - S_restart)p_complete^2.
```

`M_gram`, `M_opener`, and `M_sparse` are expected variant counts for PLINK's three modeled regression paths. Any incomplete variant takes the Gram path, giving `M(1-p_complete)`. Among complete variants, each restart segment begins with an opener, contributing `S_restart p_complete`; elsewhere a complete variant is an opener when its predecessor is incomplete, contributing `(M-S_restart)p_complete(1-p_complete)`. A non-restart complete variant following another complete variant takes the sparse/update path, giving `(M-S_restart)p_complete^2`. These mutually exclusive expected counts sum to `M`.

These sum to `M`. A Gram variant has missing calls. An opener is a complete variant at a restart or after a missing variant. A sparse-path variant follows another complete variant.

For Gram variants, the expected total observed-sample count is

```text
observed_sum = N(M_gram - M r_miss).
```

Across all variants, the expected number of observed calls is `NM(1-r_miss)`. Complete variants contribute `NM p_complete` observed calls. Subtracting those leaves `N[M(1-r_miss)-Mp_complete] = N[M(1-p_complete)-Mr_miss] = N(M_gram-Mr_miss)` observed calls among Gram-path variants.

The expected observed sample count conditional on entering the Gram branch is

```text
N_observed = observed_sum / M_gram,
```

Dividing the total expected observed calls on incomplete variants by their expected count gives the conditional mean sample count for one Gram-path regression. This expression is used only when `M_gram>0`; a no-missingness scenario has no Gram branch to price.

with `N_observed = N` when no Gram variants occur.

#### 5.4.3 Per-branch worker trace

For every stage in a Gram, opener, or sparse trace, the candidate calculates

```text
T_memory,stage = max(bytes_L1/rate_L1,
                     bytes_L2/rate_L2,
                     bytes_L3/rate_L3,
                     bytes_DRAM/rate_DRAM)

T_arithmetic,stage = FLOPs / arithmetic_rate

T_stage = max(T_instruction,stage,
              T_memory,stage,
              T_arithmetic,stage).
```

For each source stage, cache-level byte counts are divided by their measured bandwidths and the slowest cache/memory level defines `T_memory,stage`. FLOPs divided by the appropriate arithmetic rate define `T_arithmetic,stage`; source-operation counts priced by primitive instruction services define `T_instruction,stage`. These are simultaneous capacity limits on the same stage, so their maximum is used.

The arithmetic rate can be a separate independently measured BLAS rate for SYRK, vector GEMM, or GEMV. The branch cost is the sequential sum

```text
T_branch = sum over its stages of T_stage.
```

A single worker executes the stages of one variant branch in source order. They cannot overlap within that worker trace, so stage services are added to obtain per-variant branch service.

The Gram stages are mask, expansion, gather, triangular Gram matrix, variance-inflation check, `X'y`, solve, and post-regression work. The opener omits the full Gram/solve path and uses fast cross-products plus rank-one updates. The sparse branch uses carrier work plus the rank-one update and post-regression steps.

#### 5.4.4 Decode service

Decoder instruction service is the sum of exact census counts times matching primitive services:

```text
U_decode_instruction = sum_j count_j × u_j.
```

For decoder primitive `j`, `count_j` is its exact or scenario-derived occurrence count and `u_j` is CPU-seconds per occurrence. Summing across primitives gives decoder instruction service. The subscript `j` enumerates operations such as record dispatch, group headers, and category extraction; it does not denote a worker.

Let `H_decode` be the encoded payload, packed output, copies, fills, and inversions modeled for the decoder. Then

```text
U_decode = max(
    U_decode_instruction,
    H_decode / L2_rate
).
```

`H_decode/L2_rate` is the cache-traffic service of the same decoder work. The larger of instruction service and cache service is used as a decoder roofline rather than adding two overlapping descriptions of one execution.

The quantity is CPU service before distribution among workers.

#### 5.4.5 Compute-pool wall time

Let `T_thread` be the number of calculation threads and `D_j` the CPU work assigned to worker `j`:

```text
D_j = [sum over branches b of M_b × T_branch,b,j] / T_thread
    + U_decode / T_thread
    + [max(0, T_gram_first,j - T_gram_steady,j)
       × S_restart / T_thread].
```

`D_j` is normalized demand assigned to worker `j`. For each branch `b`, its expected variant count `M_b` is multiplied by worker `j`'s per-variant service and divided by the nominal thread capacity `T_thread`. Decoder work is distributed in the same units. The final term adds the excess first-Gram cost over steady Gram cost once per restart segment, clipped at zero so a faster first case cannot create negative work.

The last term charges first-touch/cache warmup at every worker/block restart.

Total pool CPU work is

```text
U_compute = sum_j D_j.
```

Summing normalized demand over workers gives total CPU service for the compute pool. This is distinct from the busiest-worker critical path used next.

The compute-pool wall time is

```text
T_compute_pool = max(
    max_j(D_j) / q,
    U_compute / C_cpu
).
```

`max_j(D_j)/q` is the time required by the most heavily loaded serial worker at scheduling share `q`. `U_compute/C_cpu` is the whole-pool CPU-capacity floor. Both must be satisfied, so the compute pool takes their maximum.

The first term enforces the slowest worker under its scheduling fraction. The second enforces total CPU capacity.

#### 5.4.6 Formatting and output

For the current `--glm hide-covar` text scenario,

```text
U_format = M × [4u_dtoa + u_line_format + u_t_tail + 2u_uint32].
```

Each of the `M` output rows is modeled as formatting four floating-point values, one complete line, one Student-t tail, and two unsigned integer fields. The corresponding `u_*` values are CPU-seconds per operation; the bracket is service per variant and multiplication by `M` gives total formatting CPU service.

Let `identifier_bytes_per_variant` be derived from the `.pvar` character count. The modeled output size is

```text
output_bytes = M × [identifier_bytes_per_variant
                    + 4×12 + 14 + digits(N) + 15].

T_output = output_bytes / write_rate.
```

The bracket estimates bytes per row: identifier text, four 12-character floating fields, fixed separators and labels, the decimal width of `N`, and a final fixed allowance. Multiplication by `M` gives total output bytes, which are divided by sustained write bandwidth. This models the requested PLINK text output, not TorchGWAS binary output.

These character allowances describe the evaluated output scenario; they are not a maximum for arbitrary identifiers or numeric values.

#### 5.4.7 Finite main-thread and compute schedule

PLINK processes blocks of `B_P = 65,536` variants. For each block fraction `f_i = b_i/M`:

```text
T_read,i = pgen_file_bytes × f_i / read_rate
T_compute,i = T_compute_pool × f_i
T_format_write,i = U_format × f_i / q + T_output × f_i.
```

Each block receives fraction `f_i` of total file bytes, compute-pool service, formatting CPU service, and output storage service. Formatting and its corresponding write are serial within this block node and are added. Division by `q` converts formatting CPU service to wall time.

The ordering is:

```text
read current block
join previous calculation
launch current calculation
format previous block
```

The final block is formatted only after its calculation finishes. `plink_block_schedule()` builds this dependency graph exactly; reading and formatting share the main thread, while the worker pool can overlap them.

Two additional whole-scan floors are

```text
T_CPU_floor = (U_compute + U_format) / C_cpu

T_DRAM_floor = total_modeled_DRAM_bytes / shared_DRAM_rate.
```

`U_compute+U_format` is total CPU-seconds across the scan; division by available core-equivalents prevents the finite schedule from exceeding aggregate CPU capacity. The second expression divides the scan's modeled host-memory traffic by shared DRAM bandwidth. They are independent lower bounds on the same finite scan.

The scan time is

```text
T_scan = max(
    T_finite_PLINK_schedule,
    T_CPU_floor,
    T_DRAM_floor
).
```

The explicit block schedule supplies the dependency-limited makespan. The CPU and DRAM terms supply whole-run capacity floors that the finite abstraction might otherwise violate. All three constraints must hold, so the largest determines modeled scan time.

Finally,

```text
T_PLINK = T_setup + T_scan.
```

### 5.5 Exact fastGWA-lr decomposition

The current candidate is `fastgwa_runtime()` in `src/torchgwas/mechanistic_cpu.py`. The tested GCTA 1.95.3 build uses one sequential reader and one effectively serial analysis consumer.

The final equation is

```text
T_fastGWA = T_setup + T_scan.
```

`T_setup` is the one-time table parsing, sample matching, and covariate-design work. `T_scan` is the finite genotype reader/decoder/analysis/output schedule. The two are added because the scan requires the completed setup objects and therefore cannot overlap setup in the modeled implementation.

#### 5.5.1 Setup

Let `P = C + 1`, containing the intercept and covariates before the tested genotype is processed.

For each input table `t`, parser service is

```text
U_text += rows_t × u_split_mode(t).
```

For table `t`, `rows_t` is its number of parsed rows. `u_split_mode(t)` is the independently measured CPU-seconds per row for the parser branch selected by that table's column count and field-length pattern. Multiplying rows by CPU-seconds per row gives CPU-seconds for the table; `+=` means the service is accumulated across the sample, phenotype, covariate, and variant tables.

The mode depends on the table’s column count and whether it contains many long fields. Retained string-copy service is

```text
U_copy += short_fields_t × u_short_string
        + long_fields_t × u_string_construct.
```

`short_fields_t` and `long_fields_t` count retained text fields in table `t`. `u_short_string` is CPU-seconds for copying a short field through the small-string path, whereas `u_string_construct` is CPU-seconds for constructing a longer heap-backed string. The two classes are added because both kinds of fields occur and require separate copy work.

Marker ID/ref/alt retention adds

```text
U_copy += 3M × u_short_string.
```

Every one of the `M` variants retains three short marker strings: identifier, reference allele, and alternate allele. This creates `3M` short-string copy operations. The added service uses the same CPU-seconds-per-copy primitive as the table census.

Numeric parsing is

```text
U_numeric = N(K_file + C) × u_strtod
          + 2M × u_stoi.
```

`K_file` is the number of phenotype columns present in the input table, including columns not selected for the analyzed phenotype. Each of the `N` sample rows contains `K_file + C` floating-point phenotype/covariate fields, so `N(K_file+C)` values pass through string-to-double conversion at cost `u_strtod` CPU-seconds each. The variant table contributes two integer fields per variant, represented by `2M` string-to-integer conversions at cost `u_stoi` each. The services add because both conversion classes must run.

Sample matching and sorting use the exact saved comparison counts:

```text
U_matching = 2 × string_sort_compares × u_string_compare
           + 3 × integer_sort_compares × u_string_compare
           + [2(N-1) + 3N] × u_string_compare.
```

`string_sort_compares` and `integer_sort_compares` are comparison counts saved from the exact sorting/matching census for the workload. The factors two and three count the repeated comparison passes performed by the source path. The final bracket counts linear matching passes over ordered sample identifiers: two adjacent-row passes contribute `2(N-1)` comparisons and three complete sample passes contribute `3N`. The implementation uses the same measured comparison-service coefficient `u_string_compare` for these retained comparison proxies. Every count is multiplied by CPU-seconds per comparison and then added, yielding total matching CPU-seconds. These multipliers are source-path counts, not fitted timing coefficients.

Setup arithmetic and retained-table traffic are

```text
F_setup = 4NP^2 + (2/3)P^3 + 4NP
H_setup = 64 × total_fields + 16NP.
```

`P=C+1` is the fixed-effect design width after adding the intercept. `F_setup` is the setup arithmetic census in floating-point operations: `4NP^2` represents the leading sample-by-design cross-product/factorization work, `(2/3)P^3` the small dense factorization, and `4NP` the remaining sample-by-design transformations. `H_setup` is the logical host-memory traffic in bytes: `total_fields` is the number of parsed table fields and the 64-byte coefficient represents the retained per-field object/metadata traffic; `16NP` represents two eight-byte movements for each entry of the `N × P` double-precision design. These coefficients describe the audited source scenario and should be changed if the parser or design representation changes.

The CPU setup work uses the larger arithmetic or memory service:

```text
U_setup = U_text + U_copy + U_numeric + U_matching
        + max(F_setup / CPU_FP64_rate,
              H_setup / DRAM_rate).
```

The first four terms are serial CPU services already expressed in CPU-seconds. `F_setup/CPU_FP64_rate` converts setup FLOPs to arithmetic-limited seconds, while `H_setup/DRAM_rate` converts logical bytes to memory-limited seconds. Arithmetic and memory constrain the same dense setup operations, so their larger value is used rather than adding both. That roofline service is then added to parsing, copying, numeric conversion, and matching because those are separate required operations.

Thus

```text
T_setup = U_setup / q.
```

`q` is the fraction of one logical CPU available to this serial setup path. Dividing the total CPU service by `q` converts it to elapsed setup time; for example, `q=0.5` doubles the elapsed time.

As for PLINK 2.0, process/environment startup and metadata-read time are outside the calculator boundary. Table parsing and matching remain included, and genotype payload reading remains in the finite reader-analysis schedule.

#### 5.5.2 Decoder

Decoder instruction service is

```text
U_decode_instruction = sum_j count_j × u_j,
```

Decoder operation class `j` might be a record dispatch, group-header parse, category extraction, or another audited primitive. `count_j` is the number of times the input census invokes that class, and `u_j` is its independently measured CPU-seconds per invocation. The products are summed because all invoked decoder operations must be performed.

including record dispatch and explicit proxies for group headers and category extraction. Let `H_decode` be modeled packed-input and decoder-copy traffic. Then

```text
U_decode = max(
    U_decode_instruction,
    H_decode / L2_rate
).
```

`H_decode` is the decoder's modeled packed-input and copy traffic in bytes. Dividing by sustained `L2_rate` gives the cache-bandwidth service for that traffic. Instruction execution and L2 traffic describe two capacity limits on the same decoding work, so the roofline takes their maximum. The result `U_decode` is serial decoder CPU service expressed in seconds at full scheduling share.

The sequential reader service for the full file is

```text
T_reader_total = pgen_file_bytes / read_rate + U_decode / q.
```

`pgen_file_bytes/read_rate` is the time needed to deliver the genotype payload from storage. `U_decode/q` is elapsed decoder time at serial CPU share `q`. They are added within the sequential reader because this candidate does not grant storage/decode overlap inside that reader service. The reader may still overlap the downstream analysis consumer in the finite queue.

#### 5.5.3 Per-variant analysis trace

The source-level stages are:

1. count packed genotypes;
2. allocate/zero the genotype vector;
3. expand genotypes;
4. make the Eigen alias-safe temporary copy;
5. allocate/zero `Hy`;
6. compute `H y`;
7. compute `X(H y)` and subtract it;
8. copy the result;
9. calculate `x'x`;
10. calculate `x'phenotype`.

For every stage,

```text
T_stage = max(T_instruction, T_arithmetic, T_cache_memory),
```

For one source-level stage, `T_instruction` is instruction-throughput service, `T_arithmetic` is floating-point service, and `T_cache_memory` is the largest applicable cache/DRAM service. They are competing bounds on the same stage, so the slowest capacity determines the stage estimate. The ten distinct stages listed above remain sequential and are therefore added after each stage's maximum is found.

and the stages are summed. The ideal-LRU trace is repeated until its cost becomes stationary. Therefore

```text
U_math = T_variant_first
       + (M - 1)T_variant_steady.
```

The first analyzed variant pays `T_variant_first`, which includes cold allocation/cache effects in the modeled trace. Each of the remaining `M-1` variants pays the stationary trace cost `T_variant_steady`. The expression is CPU service for the effectively serial analysis loop; the later division by `q` converts it to wall time.

These are CPU-seconds, despite the `T` notation, because the serial worker’s scheduling fraction is applied later.

Tail probability and formatting service are

```text
U_tail = M × u_chi_square_tail

U_format = M × [4u_ostream_double
                + u_ostream_uint
                + 2u_to_string_uint
                + 8u_short_string].
```

Every variant performs one chi-square-tail evaluation, so `M` is multiplied by the measured CPU-seconds per tail call `u_chi_square_tail`. Formatting emits four floating-point fields, one unsigned-integer field, two integer-to-string conversions, and eight retained short strings per variant. The coefficients 4, 1, 2, and 8 are source-level call counts; the corresponding `u_*` values are CPU-seconds per call. Multiplication by `M` and summation therefore produce total tail and formatting CPU-seconds.

For the current output-character scenario,

```text
output_bytes = M × [marker_characters_per_variant
                    + digits(N) + 4(6+6) + 10]

T_output = output_bytes / write_rate.
```

The bracket is the modeled number of output characters per variant. `marker_characters_per_variant` counts identifier and allele text. `digits(N)` counts the printed sample-count width. `4(6+6)` represents four numeric fields with six value characters plus six formatting/separator characters each in the declared scenario, and `10` covers the remaining fixed delimiters and line structure. Multiplying by `M` gives bytes because the output is single-byte text. Dividing by sustained `write_rate` gives storage-write seconds.

Total serial analysis service is

```text
T_analysis_total = (U_math + U_format + U_tail) / q
                 + T_output.
```

The three `U` terms are serial CPU-seconds for regression arithmetic, formatting, and tail probabilities. They are added and divided by the analysis worker's CPU share `q`. Text-output storage time is then added because this consumer writes its produced line stream as part of the same sequential analysis service.

#### 5.5.4 Three-slot finite queue

Variants are divided into blocks of 1,024. For block `i`, the reader and analysis services are proportional to its actual variant count. The first block additionally receives the difference between first-touch and steady variant cost.

For queue depth `D_q = 3`:

```text
reader_start_i = max(
    reader_end_(i-1),
    analysis_end_(i-D_q) if i >= D_q else 0
)

reader_end_i = reader_start_i + reader_service_i

analysis_start_i = max(
    analysis_end_(i-1),
    reader_end_i
)

analysis_end_i = analysis_start_i + analysis_service_i.
```

`reader_start_i` waits for both the preceding reader block and, once the three-slot queue is full, analysis of block `i-D_q`; the maximum enforces both constraints. Reader completion adds the block's reader service. `analysis_start_i` waits for both the preceding analysis block and publication of reader block `i`. Analysis completion then adds that block's analysis service. Thus blocks can overlap across the two stages, while neither stage reorders its own blocks and the reader cannot overwrite an occupied queue slot.

The finite queue makespan is the last `analysis_end_i`.

The whole-scan resource floors are

```text
T_CPU_floor = (U_decode + U_math + U_format + U_tail) / C_cpu

T_DRAM_floor = total_modeled_DRAM_bytes / shared_DRAM_rate.
```

The CPU numerator is all modeled CPU work, in CPU-seconds, across decoding and analysis. Even with ideal overlap it cannot finish faster than division by the process-wide capacity `C_cpu` in core-equivalents. The DRAM floor similarly divides all modeled host-memory bytes by shared sustained bandwidth. These are whole-run lower bounds and do not replace the dependency schedule.

Therefore

```text
T_scan = max(
    final_analysis_end,
    T_CPU_floor,
    T_DRAM_floor
)

T_fastGWA = T_setup + T_scan.
```

`final_analysis_end` is the two-stage, three-slot finite-queue makespan. The scan must also respect the process-wide CPU and DRAM floors, so `T_scan` is their maximum. Finally, one-time setup is added because scan execution begins after setup is complete.

### 5.6 Exact arithmetic example from a frozen H100 profile

The following example applies the paper-facing calculator boundary to the retained H100 candidate for `N = 2,048`, `M = 237,824`, `K = 1`, and eight covariates. It is included only to demonstrate the arithmetic and should not be quoted as a measured crossover.

#### TorchGWAS

The saved top-level terms are:

| Term | Seconds |
| --- | ---: |
| CUDA context | 0.4961648276 |
| Metadata QC and covariate basis | 0.2430765700 |
| Pinned allocation | 0.0048833829 |
| GPU design setup | 0.0002084360 |
| Lazy first use | 0.0975402421 |
| Finite scan plus binary close | 0.0714627641 |
| Sidecars and exit | 0.0331473489 |

The top-level sum is

```text
T_TorchGWAS
  = 0.4961648276 + 0.2430765700
  + 0.0048833829 + 0.0002084360
  + 0.0975402421 + 0.0714627641 + 0.0331473489
  = 0.9464835716 s.
```

The 0.0714627641-s scan term already contains overlap among 117 blocks, four decoder workers, depth four, GPU work, result transfer, 1,902,592 output bytes, two writer calls, and final close. None of those internal stage times should be added again to the top-level sum.

#### PLINK 2.0

The setup terms sum to

```text
T_setup
  = metadata parse
  + metadata memory
  + sample match
  + setup arithmetic

  = 0.0300856166 + 0.0038631370
  + 0.0001099716 + 0.0000102470
  = 0.0340689722 s.
```

The scan alternatives are

```text
finite four-block schedule = 2.6555582649 s
whole-process CPU floor    = 1.3239929700 s
whole-process DRAM floor   = 0.0028968246 s.
```

Therefore

```text
T_scan = max(2.6555582649,
             1.3239929700,
             0.0028968246)
       = 2.6555582649 s

T_PLINK = T_setup + T_scan
        = 0.0340689722 + 2.6555582649
        = 2.6896272371 s.
```

The displayed `decoder_pool`, `regression_pool`, `format_main`, and `buffered_output` fields diagnose components. They are not all added directly because compute, reading, and formatting overlap in the finite schedule.

#### fastGWA-lr

The setup terms sum to

```text
T_setup
  = 0.1087673120
  + 0.0057057684 + 0.0128955691 + 0.0004232980
  + 0.0038769862
  = 0.1316689337 s.
```

The scan alternatives are

```text
finite 233-block queue = 2.1720096985 s
whole-process CPU floor = 0.2804398034 s
whole-process DRAM floor = 0.0045897634 s.
```

Therefore

```text
T_scan = max(2.1720096985,
             0.2804398034,
             0.0045897634)
       = 2.1720096985 s

T_fastGWA = 0.1316689337 + 2.1720096985
          = 2.3036786322 s.
```

Again, `reader_seconds = 0.1184167283` and `analysis_seconds = 2.1714998310` are not added. Their overlap and queue backpressure produce the 2.1720096985-s finite makespan.

### 5.7 Exact current Figure 2 predictor and boundary search

The regenerated crossover panel is produced by `_scratch/plot_figure2_current_source_resource.py`. It calls the same current source/resource functions described above:

```text
T_Torch(N,M)   = torch_runtime(work(N,M), H100_Torch_profile)
T_PLINK(N,M)   = plink_runtime(work(N,M), H100_PLINK_profile)
T_fastGWA(N,M) = fastgwa_runtime(work(N,M), H100_fastGWA_profile).
```

Each function receives the same workload dimensions `N` and `M`, constructs the method-specific source-work census, and combines it with the corresponding frozen H100 resource profile. The three returned values are predicted elapsed seconds. This notation is a function definition, not an assertion that the tools perform identical operations internally.

No end-to-end association runtime is used as a multiplicative or additive fitting coefficient. For each `N` and `M`, `scenario_data()` constructs the source-work census for one analyzed phenotype, eight covariates, MAF 0.2, iid genotype missingness 0.001, 32 phenotype columns in the input file, and matching sample order. The exact TorchGWAS, PLINK 2.0 and fastGWA-lr decompositions are the equations in Sections 5.3, 5.4 and 5.5.

#### 5.7.1 Figure-specific H100 resource scenario

The figure uses the following explicit conditions:

- fresh-run analytical timing for TorchGWAS, with general environment bootstrap omitted but CUDA context and lazy first use included;
- genotype-input bandwidth from the retained H100 storage profile, with no separate metadata-read term;
- full modeled GPU availability;
- 48 requested PLINK threads and 47 calculation workers;
- `40.125` effective PLINK core-equivalents, independently measured as the median of `41.75` and `38.50` in the quiet-host `K=2` runs;
- two 24-core sockets with 60 MiB L3 per socket;
- a shared 60-MiB PLINK reference-stream cache scenario, so common phenotype/covariate residency is not divided 47 times among otherwise identical worker traces.

The final item is an explicit cache-sharing approximation, not a fitted timing coefficient. The equal-partition alternative is retained in `cache_policy_probe.json`; it is much too pessimistic for the measured 47-worker PLINK run because it duplicates common-source eviction in every worker trace.

#### 5.7.2 Discrete equality search

For a fixed sample count and competitor, define

```text
Delta(M) = T_Torch(N,M) - T_competitor(N,M).
```

At fixed `N`, `Delta(M)` is the predicted TorchGWAS runtime minus the predicted comparator runtime. A positive value means TorchGWAS is slower; zero is predicted equality; a negative value means TorchGWAS is faster.

`discrete_crossing()` starts at one TorchGWAS marker block, doubles `M` until the sign changes, and then performs integer bisection on the marker-block grid. If

```text
Delta(M_low)  > 0
Delta(M_high) <= 0,
```

These inequalities bracket the first sign change on the discrete chunk grid. The lower marker count is still comparator-faster, while the upper marker count is TorchGWAS-faster or tied.

then the model predicts the transition from comparator-faster to TorchGWAS-faster inside `[M_low, M_high]`. The plotted dashed line uses the geometric midpoint

```text
M_curve = sqrt(M_low × M_high).
```

The geometric midpoint is used because marker counts are searched and plotted multiplicatively. It lies halfway between the two bracket endpoints on a logarithmic `M` axis; it is a plotting representative, not an additional modeled or measured equality point.

The search also evaluates larger marker counts to check for a later sign reversal. The bracket is a numerical equality interval at the chosen block resolution; it is not a statistical confidence interval and is not proof of universal dominance.

#### 5.7.3 Retired Figure 2 draft

The earlier draft generated by `_scratch/plot_panel_b_plink48_projected.py` used the retained component predictor in `benchmarks/direct_marker_ties.py` and `benchmarks/direct_nm_model.py`. Its equations are still useful audit history, but that curve is superseded and must not be cited as the current Figure 2 predictor.

### 5.8 Exact TorchGWAS pipeline-planner equation

`src/torchgwas/pipeline_model.py` contains a separate planning model. It does not generate the cross-tool Figure 2 curve.

For a candidate plan, define:

```text
R = ceil(M / read_tile)              number of reads
D_tiles = ceil(M / decode_tile)      number of decode tiles
Q = ceil(M / B)                      number of compute chunks

transfer_bytes = M × transfer_bytes_per_variant
decoded_bytes = M × decoded_row_bytes
output_bytes = M K × output_bytes_per_test
result_bytes = M(16K + 5).
```

`read_tile`, `decode_tile`, and `B` are numbers of variants handled by one storage read, one decode operation, and one compute chunk, respectively. Applying `ceil(M/tile)` counts the final partial unit as an additional read, decode tile, or chunk. `transfer_bytes_per_variant` is the actual host-to-device representation width, so multiplying it by `M` gives total transfer bytes. `decoded_row_bytes` is the in-memory decoded width per variant. `output_bytes_per_test` is 12 for the default three-float32 binary output (coefficient, t statistic, and `-log10(P)`) and 8 for the t-statistic plus `-log10(P)` field set, so `MK` tests determine the durable payload. During the default dense scan, each phenotype-variant cell returned from the GPU carries 4-byte beta, 4-byte t, and 8-byte float64 `-log10(P)`, totaling `16MK`; every variant also returns a 1-byte status and 4-byte residual degrees of freedom, totaling `5M`. Therefore the default in-memory result transfer is `M(16K+5)`. A no-output scan that does not calculate log significance should instead use `M(8K+5)`.

The resource-service dictionary is calculated term by term.

#### Storage

```text
T_input = genotype_record_read_bytes / disk_rate
T_output = output_bytes / output_rate
T_read_latency = R × read_latency_seconds

T_storage = T_input + T_output + T_read_latency.
```

`genotype_record_read_bytes` is the number of genotype payload bytes actually read, and division by sustained `disk_rate` gives transfer time. `output_bytes/output_rate` gives durable-output service. Every one of the `R` reads also pays the fixed `read_latency_seconds`, so multiplication gives total read-call latency. The three terms add because the coarse planner assigns input, output, and read-call overhead to the same sequential storage resource. More detailed finite scheduling may overlap them only when it represents separate resources explicitly.

#### CPU decode

```text
U_decode = supplied_total_decode_CPU_seconds
```

When a source-specific census has already summed decoder primitive work, that measured total is used directly. It is CPU service, not elapsed time and not a complete-program timing.

or, when only a per-variant service is supplied,

```text
U_decode = M × decode_CPU_seconds_per_variant.
```

Otherwise, `decode_CPU_seconds_per_variant` is the independently measured average decoder CPU service for one variant under the declared encoding. Multiplication by `M` produces total decoder CPU-seconds. This fallback is less detailed than a record-form census and must not be reused for a materially different encoding.

CPU decode concurrency is

```text
C_decode = min(workers, depth),
```

At most `workers` decoder tasks can execute simultaneously, but at most `depth` pipeline slots can hold their in-flight products. The smaller count is therefore the maximum useful decoder concurrency. Explicit operating conditions may reduce this ideal count further.

possibly reduced by explicit operating conditions. Then

```text
T_CPU_decode = [U_decode
                + I(CPU_decode) × D_tiles × decode_launch_seconds]
               / C_decode.
```

`I(CPU_decode)` is one when decoding occurs on the CPU and zero otherwise. Each of the `D_tiles` CPU decode tiles pays one launch/dispatch service, so their product is added to the useful decode CPU work. Division by effective concurrency `C_decode` converts total parallelizable CPU service to an optimistic elapsed decode-resource demand. The equation assumes reasonably balanced tiles; the finite scheduler handles tail effects separately.

#### Host memory

```text
host_traffic = genotype_record_read_bytes
             + (1 + 2×host_staging_copies) × transfer_bytes
             + 3 × result_bytes
             + M × extra_host_bytes_per_variant
             + I(CPU_decode) × decoded_bytes

T_host_memory = host_traffic / host_memory_rate.
```

`genotype_record_read_bytes` represents delivery of encoded records into host memory. If `s=host_staging_copies`, `(1+2s)transfer_bytes` counts one DMA read of the final transfer buffer plus one host read and one host write for each staging copy. `3×result_bytes` counts the device-to-host DMA write into pinned memory followed by a host read and write when results are copied into owned arrays. `M×extra_host_bytes_per_variant` prices any explicitly declared auxiliary per-variant traffic. `I(CPU_decode)×decoded_bytes` includes decoded-matrix traffic only when a CPU decoder materializes it. These logical byte streams share host memory and are added, then divided by sustained `host_memory_rate`. Cache reuse can make physical DRAM traffic smaller, so the interpretation is conditional on this census.

#### PCIe transfers

```text
T_H2D = transfer_bytes / H2D_rate
T_D2H = result_bytes / D2H_rate.
```

The first equation divides the bytes in the selected genotype transport representation by sustained host-to-device bandwidth. The second divides returned result bytes—including the float64 significance matrix in default dense-output mode—by sustained device-to-host bandwidth. These are separate directional PCIe/NVLink services; the resource dictionary keeps them separate so the schedule can represent overlap when supported.

#### GPU statistics

Let `C_design = C + 1`, because the planner adds the intercept. The operation count is

```text
F_GPU = M × [2N(K + C_design) + 4N + 2C_design + 12K].
```

`C_design=C+1` includes the intercept. For each variant, `2N(K+C_design)` counts multiply-add work for genotype products with all phenotype and design columns. `4N` represents genotype centering/sum-of-squares reductions, `2C_design` the small covariate-projection products, and `12K` the retained per-phenotype statistic formation. Multiplication by `M` gives total modeled FLOPs. This is a source-level operation census for the fixed-effect path, not a claim about exact issued GPU instructions.

For the eager Torch kernel, logical genotype traffic is 32 bytes per genotype. The small-array term is 16 bytes per design/result value:

```text
H_GPU = M × [32N + 16(K + C_design)].
```

For each variant in the eager implementation, the logical genotype intermediates contribute 32 bytes per sample, hence `32N`. Small design/result arrays contribute 16 bytes for each phenotype and design column, hence `16(K+C_design)`. Multiplication by `M` gives total logical device-memory bytes. A native fused kernel substitutes `(2w+8)N` for the genotype term because its raw row width `w` and materialization behavior differ.

For a native fused kernel with raw input width `w`, the genotype term becomes `(2w+8)N`. Packed PGEN uses its physical padded row width explicitly.

The roofline reference is

```text
T_GPU_roofline = max(
    F_GPU / GPU_FLOP_rate,
    H_GPU / GPU_memory_rate
).
```

The first quotient is the arithmetic-capacity time and the second is the device-memory-capacity time. They limit the same GPU statistic work, so the roofline uses the larger rather than adding both. The result is a reference service under the declared effective rates; it does not include launch, decode, conversion, or significance-tail work unless those counts are explicitly present.

The statistics service is either an independently supplied resident-component value or the roofline reference:

```text
T_statistics = M × supplied_statistics_seconds_per_variant
```

If a separately measured resident statistics primitive is supplied, its seconds per variant are multiplied by `M`. The primitive must use the same `N`, `K`, design width, kernel path, and precision or be scaled by an explicit source model; otherwise this shortcut is not portable.

or

```text
T_statistics = T_GPU_roofline.
```

When no compatible resident primitive is supplied, the operation/byte roofline becomes the statistics service. Exactly one of these two alternatives is selected, so they are never added together.

Compute-launch service is then added:

```text
T_compute = T_statistics + Q × compute_launch_seconds.
```

There are `Q` compute chunks, and each pays one fixed launch/dispatch service. Those launches are added to the useful statistics service because both are required. This term must also include a separately modeled Student-t-tail service when dense `-log10(P)` is calculated; the older retained coefficient does not make that special-function work free.

GPU decode and input conversion share the same GPU and are therefore sequential resource demands:

```text
T_GPU_decode = M × GPU_decode_seconds_per_variant
             + I(GPU_decode) × D_tiles × decode_launch_seconds

T_GPU_conversion = M × conversion_seconds_per_variant

T_GPU = T_compute + T_GPU_decode + T_GPU_conversion.
```

GPU decoding pays useful per-variant service plus one launch for every decode tile when `I(GPU_decode)=1`. Conversion pays its per-variant service for all `M` variants. Decode, conversion, and association kernels consume the same GPU timeline in this coarse resource model, so their total demands are added. A more detailed multi-stream execution graph may overlap kernels only when it has explicit evidence and dependencies for doing so.

#### Resource maximum and fill/drain approximation

The planner forms

```text
resource = {
    storage:    T_storage,
    cpu_decode: T_CPU_decode,
    host_memory:T_host_memory,
    h2d:        T_H2D,
    d2h:        T_D2H,
    gpu:        T_GPU
}.
```

Every dictionary value is a whole-run demand in seconds on a distinct coarse resource: storage, CPU decoding, host memory, host-to-device transfer, device-to-host transfer, or GPU execution. Keeping the entries separate prevents a byte-rate bound from being mistaken for serialized elapsed time.

The steady resource service is

```text
T_bottleneck = max(resource values).
```

In an infinitely long, perfectly overlapped pipeline, distinct resources can work concurrently. Its throughput cannot be better than the slowest whole-run demand, so the largest entry is the steady-state lower bound.

Let

```text
f_batch = min[M, max(read_tile, decode_tile, B)] / M.
```

`max(read_tile,decode_tile,B)` is the largest variant batch used by any pipeline stage. It is capped at `M` so a one-batch run has fraction one. Dividing by `M` estimates what fraction of total work belongs to one fill/drain batch: the fraction approaches zero for many batches and one for a one-batch workload.

The coarse fill/drain term is

```text
T_fill_drain = f_batch × [sum(resource values) - T_bottleneck].
```

`sum(resource values)-T_bottleneck` is all non-bottleneck service that ideal steady-state overlap would hide. Multiplying by one batch fraction restores an approximate startup/drain share of that hidden service. This is an explicit structural approximation, not a fitted overlap factor.

The planning score is

```text
T_planning = T_bottleneck + T_fill_drain.
```

The plan-ranking score adds the steady-state bottleneck demand and the approximate fill/drain penalty. It is suitable for comparing candidate plans under the same assumptions, not for claiming exact end-to-end runtime.

This is a plan-ranking quantity, not the newer method-specific end-to-end candidate and not a validated runtime promise.

#### Exact finite flow-shop helper

When the planner has explicit sequential resource totals `S_j` and uniform chunks, `finite_pipeline_seconds()` uses an exact flow-shop recurrence.

Let

```text
q = floor(M/B)
tail = M mod B
s_j = S_j × B/M               service of one full chunk at stage j.
```

`q` is the number of complete `B`-variant chunks and `tail` is the remaining variant count. `S_j` is the total whole-run service assigned to stage `j`; under the uniform-work assumption, one full chunk owns fraction `B/M`, so its stage service is `s_j=S_jB/M`.

For `q` full chunks, completion at stage `j` is

```text
C_full,j = sum_(h=1..j) s_h
         + (q-1) × max_(h=1..j) s_h.
```

The first full chunk must traverse stages 1 through `j`, giving the prefix sum. Every additional full chunk advances at the slowest stage encountered in that prefix, giving `(q-1)` times the prefix bottleneck. This closed form applies to a uniform serial flow shop; it is not valid when chunks have different work or share more complicated resource constraints.

If a partial tail exists, append it in stage order:

```text
C_tail,j = max(C_tail,j-1, C_full,j)
         + s_j × tail/B.
```

At stage `j`, the partial tail must wait for both its own completion at the preceding stage and the last full chunk's completion at the current stage. The maximum enforces those dependencies. Its service scales the full-chunk service `s_j` by the actual fraction `tail/B`.

The makespan is the last stage’s completion. This equation handles one-chunk workloads and partial final chunks without a fitted overlap multiplier.

#### Plan selection

For every candidate combination of read tile, decode tile, compute chunk, workers, and depth:

1. compute resource service;
2. compute host and device memory;
3. reject infeasible candidates;
4. rank feasible candidates lexicographically by

```text
(
  T_planning,
  host_buffer_bytes + device_buffer_bytes,
  workers
).
```

Candidate plans are ordered lexicographically. The smallest predicted planning time wins. Exact time ties prefer the smaller combined host/device buffer footprint, and remaining ties prefer fewer workers. Memory and worker counts therefore break ties; they do not override a faster modeled plan.

Thus the fastest modeled plan wins first; ties prefer less memory, then fewer workers.

## 6. Step 1: state the question before calculating

A performance estimate is not meaningful until the scenario is defined. Record at least:

- `N`, `M`, `K`, and `C`;
- genotype format and representation;
- hard calls or dosages;
- missing-call rate and allele-frequency assumptions if exact records are unavailable;
- number of phenotype columns physically present in text input;
- requested output fields and output format;
- software version;
- CPU thread or worker configuration;
- GPU model;
- storage/cache condition;
- timing boundary.

Two timings with different boundaries must not be compared as though they measured the same work.

### 6.1 Timing boundaries used in this project

The full-genome manuscript benchmark uses process-start end-to-end timing. It includes Python/PyTorch startup, CUDA-context creation, input opening and decoding, association testing, result transfer, and output completion.

The frozen three-server calculator validation uses an **environment-ready** TorchGWAS boundary. Imports, thread setup, and CUDA-context readiness occur before the Torch timer; input loading, data-dependent setup, lazy numerical-library first use, scan, and binary output remain inside it. The native PLINK 2.0 and GCTA executables are timed from launch through exit.

This asymmetric boundary was chosen to isolate unstable network-mounted Python startup in that experiment. It must be stated whenever those predictions or measurements are reported.

“Cold input pages” means that private Linux page residency was checked before launch. It does not prove that SSD-controller caches, filesystem metadata caches, or every lower storage layer were cold.

## 7. Step 2: count the source work

The `work` command produces an inspectable ledger. It does not produce a complete runtime prediction.

From the repository root:

```bash
export PYTHONPATH=src:.deps

$PY runtime_calculator.py work \
  --method torchGWAS \
  --samples 8192 \
  --markers 65536 \
  --covariates 8 \
  --traits 1 \
  --missing-rate 0.001 \
  --chunk-markers 2048 \
  --tensor-work \
  --output work_torchgwas_n8192_m65536.json
```

Here `$PY` should name the project Python environment. The `--tensor-work` option imports PyTorch, so a minimal system Python without PyTorch is insufficient.

For a real PGEN file, include the file census:

```bash
$PY runtime_calculator.py work \
  --method torchGWAS \
  --samples 8192 \
  --markers 65536 \
  --covariates 8 \
  --traits 1 \
  --pgen data/example.pgen \
  --genotype-bytes 123456789 \
  --metadata-bytes 1234567 \
  --native-output-bytes 524288 \
  --chunk-markers 2048 \
  --tensor-work \
  --output work_real_pgen.json
```

Replace the byte counts with actual file sizes. Do not estimate PGEN storage from a generic compression ratio when the file exists.

### 7.1 TorchGWAS arithmetic and traffic

For the current eager fixed-effect path, the main genotype-design matrix multiplication is counted as

```text
F_GEMM = 2 N M (K + C + 1) floating-point operations.
```

The `+1` is the intercept. Omitting it would undercount the design width.

For the validated hard-call PGEN path, genotype codes cross the host-to-device link as one-byte integer values, giving

```text
H2D bytes = N M.
```

The result transfer is

```text
D2H bytes = M (16K + 5).
```

The `16K` term contains a 4-byte float32 coefficient, a 4-byte float32 t statistic, and an 8-byte float64 `-log10(P)` value for every phenotype. The additional five bytes per variant contain a 1-byte status code and a 4-byte residual degrees-of-freedom value.

The default dense binary payload contains float32 coefficient, t-statistic, and `-log10(P)` arrays:

```text
output bytes = 12 M K for the default beta+t field set.
```

The optional `t` field set omits coefficients and therefore writes `8MK` bytes. Headers, manifests, and metadata sidecars add smaller terms.

The eager PyTorch implementation also creates intermediate tensors for missingness masks, conversion to float32, centering, reductions, residual sums of squares, and result formation. The ledger records their logical traffic separately. Logical tensor traffic is not automatically equal to physical high-bandwidth-memory traffic because caches and kernel fusion can change how often bytes reach device memory.

### 7.2 PGEN record census

PGEN is not one uniform record type. A useful census includes:

- total stored bytes;
- record-form counts;
- difflist entries;
- variable-length integer counts by encoded length;
- LD-reference records;
- records that must be replayed at chunk boundaries;
- phase, dosage, or multiallelic tracks.

The current hard-call census rejects unsupported tracks instead of assigning them an invented cost. Exact record statistics should be obtained from the actual file when possible. The synthetic crossover scenario instead states an explicit expected one-bit-plus-difflist model; it does not claim to reproduce every real cohort’s encoding.

### 7.3 PLINK 2.0 missing-call branches

For a simple independent missing-call scenario with missing rate `r`, the probability that a variant is complete in all `N` samples is

```text
p_complete = (1 - r)^N.
```

The expected number of variants requiring the missing-call Gram-matrix branch is

```text
M (1 - p_complete).
```

Other PLINK branches depend on whether neighboring variants are complete and whether a worker or block restarts. Therefore, a missing-variant percentage alone is insufficient. The calculator carries restart and adjacency terms rather than assigning one average cost to every variant.

The validated PLINK candidate is currently restricted to `K = 1`. It should return unsupported outside its modeled scope instead of extrapolating silently.

### 7.4 fastGWA-lr work

The fastGWA-lr candidate counts:

- setup and text parsing;
- sample matching;
- PGEN decoding;
- per-variant covariate projection and dot products;
- test-statistic tail calculation;
- formatting and output.

For the tested official GCTA 1.95.3 executable, the association-analysis loop is effectively serial. A reader can work ahead, but requesting many threads does not divide the regression loop across all requested cores. The finite schedule therefore represents a reader and a serial analysis consumer rather than assuming ideal multithreaded scaling.

## 8. Step 3: measure resource capacities independently

Each work term needs a compatible capacity. The dimensional rule is simple:

```text
service time = work / capacity.
```

Examples are:

```text
storage seconds = stored bytes / sustained storage bytes per second
H2D seconds = transfer bytes / H2D bytes per second
D2H seconds = returned bytes / D2H bytes per second
GPU compute seconds = floating-point operations / effective operations per second
CPU service seconds = source-operation count × seconds per operation
```

Use independently measured rates under the intended environment. A complete GWAS runtime must not be divided by a work count and then reused as a supposedly independent hardware rate; that would make validation circular.

### 8.1 CPU availability is not storage bandwidth

The resource profiles separate:

- scheduling fraction for a serial worker;
- total cores available to the process;
- shared host-memory bandwidth;
- effective cache capacity;
- input and output storage bandwidth.

A busy CPU can slow parser or analysis work without reducing an already specified storage rate by the same factor. Applying one global “contention multiplier” to all resources is therefore incorrect.

### 8.2 Cache capacity and cache traffic

Data that fit in a cache may be serviced differently from data that repeatedly reach main memory. The calculator uses explicit cache sizes and source-derived working sets to select a regime. It does not place a kink in the runtime curve merely because a measured timing curve happened to bend there.

Cache models remain approximations. Real behavior also depends on associativity, data placement, simultaneous multithreading, NUMA placement, and competing jobs.

### 8.3 GPU rates

Peak specification-sheet throughput is not automatically the effective rate of the implemented kernels. The profile uses relevant component measurements and records launch/dispatch costs, cache assumptions, and compiled kernel geometry where required.

If a requested `(N, B, K, C)` geometry is outside the compiled or measured scope, the candidate should refuse or explicitly identify the extrapolation. Missing resource evidence is not equivalent to zero cost.

## 9. Step 4: convert work into method-specific stage service

After counting work and supplying capacities, the calculator forms stage-level service demands.

For TorchGWAS, a block may require:

1. storage reading;
2. CPU record parsing and decoding;
3. host-memory copies;
4. host-to-device transfer;
5. GPU conversion and association work;
6. device-to-host result transfer;
7. result copying and formatting;
8. binary output writing and final durability.

PLINK 2.0 and fastGWA-lr have different stages and different concurrency. The calculator does not use the TorchGWAS pipeline structure for the CPU tools.

At this point, dividing total work by total capacity gives useful resource service demands, but not yet the final elapsed time. The remaining question is which services can overlap.

## 10. Step 5: schedule the finite pipeline

Suppose a run has `q = ceil(M/B)` chunks. For chunk `i`, every stage has a start time and an end time.

A stage can start only when:

- its input dependency has completed;
- the required resource is available;
- a bounded buffer slot is free;
- an ordered consumer is ready if output order must be preserved.

A simplified recurrence for a three-stage line is:

```text
read_end[i] = read_start[i] + read_service[i]

compute_start[i] = max(read_end[i], compute_end[i-1])
compute_end[i] = compute_start[i] + compute_service[i]

write_start[i] = max(compute_end[i], write_end[i-1])
write_end[i] = write_start[i] + write_service[i]
```

The implementation is more detailed because CPU, memory, storage, GPU, and queues may be shared, and because several decoder workers can operate concurrently. The principle is the same: dependencies use `max`, not simple addition.

### 10.1 Startup and drain

The first chunk cannot overlap with work that has not begun. The last chunk must finish its remaining result and output work. These fill-and-drain costs matter most for small `M` or very large `B`, where only a few chunks exist.

### 10.2 Partial final chunks

If `M` is not a multiple of `B`, the last chunk is smaller. The schedule uses its actual size. Charging a full chunk for every block would overpredict work.

### 10.3 One-chunk runs

A one-chunk run has essentially no steady-state overlap. Any model that applies a long-run throughput formula to one chunk will be too optimistic.

### 10.4 Bounded output queues and final close

The binary writer has a finite staging queue. When it fills, result production waits. A partial block may remain buffered until close, so those bytes cannot overlap earlier computation. The default dense writer also performs two payload durability calls for the coefficient and t-statistic files. Sidecar handling has a separate boundary.

## 11. Step 6: produce a runtime candidate

The `candidate` command combines one prepared input-work census with one method-specific resource profile.

Examples for the retained method profiles are:

```bash
$PY runtime_calculator.py candidate \
  --method torchGWAS \
  --inputs scenario_n8192_m65536.json \
  --profile paper/calculator_torch_process_20260917/predictions_v2.json \
  --server H100 \
  --output torch_candidate.json

$PY runtime_calculator.py candidate \
  --method PLINK2 \
  --inputs scenario_n8192_m65536.json \
  --profile paper/calculator_plink_mechanistic_20260917/predictions_v3.json \
  --server H100 \
  --output plink_candidate.json

$PY runtime_calculator.py candidate \
  --method fastGWA \
  --inputs scenario_n8192_m65536.json \
  --profile paper/calculator_fastgwa_mechanistic_20260917/predictions_v3.json \
  --server H100 \
  --output fastgwa_candidate.json
```

The input-work census used by `candidate` is richer than the public `work` ledger. It includes text-table statistics, sample-sort counts, and encoded PGEN work expected by the method-specific candidates. The frozen crossover bundle generates this schema internally. Do not pass a `work` output directly to `candidate` and assume that the schemas are interchangeable.

Optional overrides include:

- `--cpu-fraction` for serial-worker scheduling availability;
- `--cpu-cores` for total CPU capacity;
- `--read-gbps` for the declared input-cache condition;
- `--available-l3-mib` for available last-level cache;
- `--missing-rate` and `--carrier-fraction` for the PLINK scenario.

Every override changes the scenario. It should be recorded with the result.

### 11.1 How to read the output

Important output fields include:

- `estimated_seconds`: conditional elapsed-time estimate;
- stage or resource service fields: the modeled contributors;
- topology and scheduling fields: workers, chunks, and queues;
- assumptions or scope: conditions under which the result is interpretable;
- unpriced or unsupported mechanisms: known omissions;
- validation status: reminder that the candidate is not a guarantee.

If a capacity is zero, missing, or outside the supported scope, a reliable calculator should fail or return unavailable. It should not invent a finite runtime.

## 12. Step 7: search for a predicted crossover

The `crossing` command compares TorchGWAS with PLINK 2.0 and fastGWA-lr for an explicit scenario.

```bash
$PY runtime_calculator.py crossing \
  --bundle paper/calculator_environment_ready_20260917/inputs.json \
  --server H100 \
  --subjects 8192 \
  --maf 0.2 \
  --missing-rate 0.001 \
  --traits-in-file 32 \
  --numeric-characters 19.4 \
  --timing-boundary environment-ready \
  --max-markers 10000000 \
  --output crossing_h100_n8192.json
```

This validated scenario analyzes one phenotype with eight covariates. The PGEN record work is generated from the stated minor-allele frequency and independent missingness assumptions. Thirty-two phenotype columns exist in the physical input file, even though only one is analyzed, because parsing those columns costs time.

### 12.1 Search algorithm

The algorithm operates in whole TorchGWAS marker chunks.

1. Evaluate both methods at one chunk.
2. If TorchGWAS is already predicted faster, return a left-censored result.
3. Otherwise, double the marker count until the predicted runtime difference changes sign or the search maximum is reached.
4. When a sign change is found, use integer bisection to reduce the interval to one marker chunk.
5. Evaluate two larger points, at approximately two and four times the upper endpoint, to look for an immediate reversal.

In pseudocode:

```text
step = TorchGWAS chunk size
low = high = 1 block

if torch_time(step) <= competitor_time(step):
    return "already faster at search minimum"

while high < maximum:
    high = 2 * high
    if torch_time(high * step) <= competitor_time(high * step):
        break
    low = high

while high - low > 1:
    middle = floor((low + high) / 2)
    if torch_time(middle * step) > competitor_time(middle * step):
        low = middle
    else:
        high = middle

return [low * step, high * step]
```

The returned interval brackets equality in the model. It is not a statistical confidence interval. The extra larger-marker checks do not prove that the crossing is unique for every possible workload.

### 12.2 Why the boundary is not one universal `KMN` number

`KMN` is a useful scale summary, but programs contain work that does not scale identically with all three dimensions. Examples include:

- process startup;
- metadata parsing;
- phenotype-file parsing;
- genotype-record decoding;
- output formatting;
- cache transitions;
- finite chunk startup and drain;
- hardware-specific CPU and GPU rates.

Two workloads with the same `KMN` can therefore have different runtimes. The manuscript’s order-of-magnitude transitions describe the evaluated settings, not a theorem that all datasets above one product favor TorchGWAS.

## 13. Step 8: validate against held-out complete runs

A useful validation plan contains both small and large workloads, more than one `N`, and workloads on both sides of the predicted ordering change.

For every run, record:

- exact command and software version;
- workload dimensions;
- input and output files;
- thread/worker settings;
- timing boundary;
- cache/page-residency condition;
- CPU and GPU load preconditions;
- wall time and, when useful, process CPU service;
- correctness checks and output row counts.

Then compare predicted and measured time without changing the predictor from the same observations.

For measured time `t` and predicted time `p`, absolute percentage error is

```text
APE = 100 × |p - t| / t.
```

Mean absolute percentage error is the average APE across the retained validation runs.

The environment-ready validation produced method/server MAPE values from 16.6% to 59.4%. These values are too large to justify reporting a sharp numerical crossover as exact. They are compatible with using the model to explain resource mechanisms and general ordering transitions, provided measured brackets and limitations are shown.

## 14. How the current Figure 2 crossover panel was assembled

The current H100 crossover panel contains two kinds of visual evidence:

- **points**, which are measured outcomes;
- **dashed boundaries**, which are predicted model crossovers.

The shaded regions show the side of the dashed boundary on which the model predicts each tool to be faster. Shading is not a confidence region.

### 14.1 PLINK 2.0 points

The accepted quiet-host PLINK comparison requested 48 threads. PLINK used approximately 38–40 CPU cores across the accepted runs. Each point is the median of two fresh-process runs for each tool with page-cache-preconditioned PGEN input.

The measured brackets were:

| Samples | Lower point | Higher point | Measured transition |
| ---: | ---: | ---: | --- |
| 8,192 | 1.5 million variants: PLINK faster | 3.0 million variants: TorchGWAS faster | Between the two points |
| 32,768 | 166,400 variants: PLINK faster | 260,000 variants: TorchGWAS faster | Between the two points |

Linear interpolation within each measured pair places local empirical ties near 2.04 million and 201 thousand variants, respectively. These are descriptive interpolations between nearby measurements, not the continuous predictor curve.

### 14.2 PLINK 2.0 dashed curve

The regenerated PLINK curve uses the current source/resource model in `src/torchgwas/mechanistic_plink.py`, not the retained component predictor. It uses 47 calculation workers, the independently measured quiet-host capacity of 40.125 core-equivalents, the retained H100 genotype-input bandwidth and the explicit shared 60-MiB reference-stream cache scenario described in Section 5.7. No plotted crossover runtime is used to fit the curve.

This figure-specific quiet-host normalization is recorded in:

```text
paper/figure2_current_predictor_20260919/
    Figure2_panel_b_current_source_resource.json
```

It should be described as a conditional source/resource prediction, not as a universal hardware curve. The quiet measured points remain the direct evidence, and the provenance file retains the cache-sharing assumption.

### 14.3 fastGWA-lr curve and points

The fastGWA-lr curve is also regenerated by the current source/resource model. The tested GCTA 1.95.3 analysis loop was effectively serial, so increasing the requested GCTA thread count does not imply proportional analysis speedup. Measured outcomes remain separate point overlays.

### 14.4 Defensible interpretation

The current evidence supports describing the fastGWA-lr transition as being on the order of `MN ≈ 10^9` and the quiet 48-thread PLINK 2.0 transition as being on the order of `MN ≈ 10^10` in the evaluated single-phenotype H100 settings.

The data do not establish an exact universal threshold, a confidence band around the dashed curve, or guaranteed dominance for every workload above the line.

## 15. The TorchGWAS memory and chunk planner

The cross-tool predictor does not choose a safe or efficient TorchGWAS chunk by itself. That task belongs to `src/torchgwas/pipeline_model.py`.

### 15.1 Why variants are chunked

Keeping all `M` variants on the GPU would require memory proportional to `NM`, which is unnecessary. TorchGWAS keeps the processed phenotype/design data resident and streams variants in blocks of `B`.

Therefore, peak GPU memory is largely independent of the total number of variants `M`. It depends mainly on:

- sample count `N`;
- resident phenotype count `K`;
- covariate rank;
- numerical precision;
- transfer representation;
- chunk size `B`;
- decode tile size;
- pipeline depth `D`;
- decoder and kernel workspace.

If the resident phenotype/design matrix itself does not fit, reducing `B` cannot solve the problem. The phenotype panel must be divided into non-overlapping groups.

### 15.2 Major memory terms

The planner accounts for several terms separately.

1. **Resident design data.** Approximately proportional to `N(K + C + 1)` and the numerical byte width.
2. **Transferred-input ring.** Approximately

   ```text
   D × decode_tile × transferred bytes per variant.
   ```

   The transferred bytes per variant depend strongly on representation. Packed two-bit hard calls are much smaller than float32 dosages.
3. **Centered genotype workspace.** The eager path may hold one or two float32 genotype chunks, approximately proportional to `4NB`. A fused packed kernel can avoid materializing this array.
4. **Decoder output buffers.** Device decoding may need a depth-sized decoded ring. CPU decoding may instead use pinned host buffers.
5. **Small products and result buffers.** These scale with `B`, `K`, and the design width.
6. **Pinned host staging and result rings.** These affect host memory and transfer overlap.

The exact implementation uses `device_ring_bytes()` so the same accounting is used both to estimate memory and to choose a chunk. Maintaining two different memory formulas would allow the planner to choose a configuration that its own report says does not fit.

### 15.3 Chunk selection

When the user does not specify a chunk size, the current automatic path:

1. determines the available GPU-memory budget;
2. subtracts or reserves known resident and workspace requirements;
3. evaluates explicit candidate chunk sizes;
4. rejects candidates whose estimated rings do not fit;
5. limits automatic chunks to the current policy range;
6. scores feasible plans using the modeled finite pipeline and then memory/worker tie-breakers.

The automatic minimum is 128 variants and the current cap is 4,096 variants. The cap is not “use as much memory as possible.” A quiet H100 sweep showed a U-shaped runtime curve with a measured optimum near 1,024 variants in one tested setting. The published comparison configuration retained the 4,096 cap for consistency, and the model has separate validation evidence through chunk 16,384. Very large chunks can consume much more memory without improving runtime.

This is why the planner should eventually be presented as an autotuning aid rather than as a solved universal optimizer.

### 15.4 Plan search

`choose_plan()` searches combinations of:

- read tile;
- decode tile;
- compute chunk;
- number of workers;
- pipeline depth.

Each candidate first passes memory and implementation checks. The planner then estimates storage, decode, host-memory, transfer, GPU, result, and output service; builds a finite-pipeline estimate; and ranks feasible plans. Unsupported source controls remain advisory rather than being silently written to arbitrary object attributes.

### 15.5 Important limitation

The planner can determine that a plan is inconsistent with its memory budget, but it cannot prove that an accepted plan will never trigger an allocation failure. Private decoder workspace, allocator fragmentation, retained caller outputs, library workspaces, and other processes can consume additional memory. Leave a reserve and validate the selected plan on the target machine.

## 16. Common mistakes

### Mistake 1: adding all stage clocks

If stages overlap, their clocks are not additive. Use end-to-end wall time for the primary observation and the dependency schedule for explanation.

### Mistake 2: treating a measured total as a hardware rate

This makes the model circular and can hide two wrong terms that cancel.

### Mistake 3: calling a model bracket a confidence interval

The crossover interval is one-chunk numerical resolution in the model. It does not quantify statistical uncertainty.

### Mistake 4: using one `KMN` boundary for every format and machine

Startup, decoding, output, missingness, storage, and resource availability do not all scale only with `KMN`.

### Mistake 5: assuming requested CPU threads are used efficiently

Observed or source-verified concurrency matters. The tested fastGWA-lr loop is effectively serial, and PLINK’s realized CPU use depends on workload and host state.

### Mistake 6: equating logical tensor bytes with physical GPU-memory traffic

Caches, kernel fusion, and reuse can reduce physical traffic. Logical work is still useful, but it must be paired with an explicit cache/kernel model.

### Mistake 7: calling a warm page-cache run cold disk

State exactly which cache level was controlled.

### Mistake 8: assuming the largest fitting chunk is fastest

Memory is a constraint, not the optimization target. Chunk size affects launch overhead, queue behavior, kernel efficiency, and unexplained large-tile penalties.

### Mistake 9: hiding unsupported work as zero

Unknown work should be listed as unpriced or unsupported. Zero means the operation is known to be absent.

## 17. Recommended reporting checklist

Before publishing a prediction or crossover figure, verify that the report includes:

- [ ] software versions and source commit;
- [ ] `N`, `M`, `K`, `C`, genotype format, missingness, and output mode;
- [ ] timing boundary;
- [ ] storage and page-cache condition;
- [ ] requested and realized CPU concurrency;
- [ ] GPU model and numerical precision;
- [ ] independently measured resource inputs and units;
- [ ] predictor/profile version and hash or archived path;
- [ ] measured points distinguished from predicted curves;
- [ ] model bracket distinguished from confidence interval;
- [ ] validation errors retained rather than corrected away;
- [ ] unsupported mechanisms listed;
- [ ] no universal-threshold claim from one conditional scenario.

## 18. Suggested concise manuscript description

A technically accurate short description is:

> Resource-based runtime predictors for TorchGWAS, PLINK 2.0 and fastGWA-lr combined source-derived operation and byte counts with independently measured hardware service rates and finite pipeline schedules. The predictors were used to estimate conditional workload scales at which the runtime ordering changed. Measured workloads were retained as direct evidence, and predicted crossover boundaries were treated as approximate, hardware- and scenario-dependent transitions rather than universal thresholds.

For the present Figure 2 PLINK curve, add that the current source/resource model used 47 calculation workers, 40.125 independently measured quiet-host core-equivalents and the explicit shared-cache reference-stream scenario, while the plotted points were direct quiet-host end-to-end measurements.

## 19. Final interpretation

The calculator is most valuable as an auditable explanation of why runtime changes:

- it makes source work explicit;
- it separates work from machine capacity;
- it represents finite overlap and output backpressure;
- it exposes memory and queue constraints;
- it shows how a predicted ordering depends on the scenario.

Its output is not an oracle. The most defensible workflow is to use the calculator to formulate a conditional prediction, validate that prediction with matched end-to-end measurements on both sides of the expected transition, report the observed brackets, and preserve the model’s numerical errors and scope limitations.
