# TorchGWAS Benchmark Report

## Overview

This document summarizes the currently tracked TorchGWAS benchmark assets and explains how to interpret them.

The benchmarks included here are intended to answer three different questions:

1. Does TorchGWAS reproduce expected statistical behavior on small, fully materialized runs?
2. What are the observed CPU and GPU memory footprints for realistic in-memory linear GWAS workloads?
3. What is the expected steady-state memory and runtime behavior for a much larger streaming-scale workload?

This report is based on benchmark outputs stored under `results/benchmarks/`.

## Included Benchmark Assets

- `accuracy_plink.tsv`
  PLINK concordance summary against an external baseline workflow.
- `runtime_scaling.tsv`
  Historical runtime comparison points retained from the original prototype notebook.
- `multivariate_case.tsv`
  Toy multivariate case summary used as a reproducibility check for the exploratory multivariate workflow.
- `linear_memory_profile.tsv`
  Simulated small-to-medium linear GWAS CPU/GPU memory profile with full result materialization.
- `linear_streaming_profile_large.tsv`
  Simulated large-scale streaming memory profile at 20k samples / 1M SNPs / 2k traits.

## Benchmark Environment

The GPU-aware memory profiling benchmarks were executed in a non-sandbox environment with:

- GPU model: `NVIDIA H100 80GB HBM3`
- CUDA-visible device during the run: `CUDA_VISIBLE_DEVICES=5`
- PyTorch CUDA availability: `True`
- PyTorch CUDA version: `12.1`

The benchmark scripts used were:

- `benchmarks/profile_linear_memory.py`
- `benchmarks/profile_linear_streaming_large.py`

## Small-to-Medium Linear Memory Profile

### Purpose

`linear_memory_profile.tsv` measures the memory footprint and runtime of the current linear TorchGWAS path when the entire simulated workload is materialized in process. These runs exercise the existing `run_linear_gwas()` implementation directly and therefore reflect both:

- the core compute path
- the cost of storing intermediate and output-facing arrays in memory

### Shape Convention

For this benchmark, shape tokens are interpreted as:

- `n_samples x n_markers x n_traits`

All runs used:

- `n_covariates = 8`
- `chunk_size = 1024`

### Results

| device | n_samples | n_markers | n_traits | cpu_peak_rss_mb | gpu_peak_allocated_mb | gpu_peak_reserved_mb | runtime_seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| cpu | 1024 | 4096 | 8 | 601.398 | 0.0 | 0.0 | 0.603605 |
| cuda | 1024 | 4096 | 8 | 898.043 | 112.196 | 166.0 | 1.228984 |
| cpu | 2048 | 8192 | 8 | 799.746 | 0.0 | 0.0 | 1.262558 |
| cuda | 2048 | 8192 | 8 | 1021.559 | 192.258 | 306.0 | 1.553763 |
| cpu | 4096 | 16384 | 16 | 1664.797 | 0.0 | 0.0 | 2.186701 |
| cuda | 4096 | 16384 | 16 | 1682.273 | 352.758 | 354.0 | 2.457663 |

### Interpretation

- `cpu_peak_rss_mb` is the peak resident set size of the entire worker process.
- `gpu_peak_allocated_mb` is the peak PyTorch-allocated GPU memory.
- `gpu_peak_reserved_mb` is the peak PyTorch-reserved GPU memory.
- The CUDA runs still show substantial CPU RSS because host-side arrays remain resident even when the compute kernel executes on GPU.
- For the two smaller configurations, CUDA is slower than CPU. This is consistent with transfer and launch overhead dominating when the problem is not large enough.
- In this rerun, CPU remained faster than CUDA for all three tested fully materialized shapes, which indicates that these in-memory problem sizes are still below the point where GPU launch and transfer overhead is amortized.

## Large-Scale Streaming Profile

### Purpose

`linear_streaming_profile_large.tsv` addresses a different question:

> How much memory does TorchGWAS require in steady state for a large linear GWAS workload that is too large to benchmark by fully materializing the complete marker-by-trait result table?

To answer this, the streaming benchmark:

- simulates a workload with `20,000` samples, `1,000,000` SNPs, and `2,000` traits
- processes repeated genotype chunks of size `4096`
- uses the current large-scale `float32` compute path
- measures peak CPU and GPU memory during steady-state chunk execution
- estimates the full runtime by extrapolating from the mean chunk time

This benchmark intentionally does **not** materialize the complete output table for all `1,000,000 x 2,000` tests. It is therefore a compute-kernel and steady-state memory profile, not a full output-materialization benchmark.

### Configuration

- `n_samples = 20,000`
- `n_markers_total = 1,000,000`
- `n_traits = 2,000`
- `n_covariates = 8`
- `chunk_size = 4096`
- `compute_dtype = float32`
- `warmup_chunks = 1`
- `profile_chunks = 2`

### Results

| device | n_samples | n_markers_total | n_traits | chunk_size | cpu_peak_rss_mb | gpu_peak_allocated_mb | gpu_peak_reserved_mb | avg_chunk_seconds | estimated_full_runtime_seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cpu | 20000 | 1000000 | 2000 | 4096 | 2767.645 | 0.0 | 0.0 | 1.539444 | 377.163868 |
| cuda | 20000 | 1000000 | 2000 | 4096 | 2546.039 | 1400.104 | 1972.0 | 1.244163 | 304.820029 |

### Interpretation

- Under this streaming benchmark, the CUDA path reduced the estimated full runtime from about `377 s` to about `305 s`.
- The observed GPU memory footprint remained small relative to H100 capacity:
  - about `1.40 GB` allocated
  - about `1.97 GB` reserved
- CPU RSS remained in the multi-gigabyte range even on the CUDA path, which indicates that current scaling is not purely GPU-memory bound.
- The float32 large-scale path substantially reduced both host and device memory relative to the earlier float64-style benchmark, while also improving the estimated CPU and CUDA runtimes.

## What These Benchmarks Do And Do Not Show

### What they do show

- The current linear TorchGWAS implementation can use GPU execution without requiring very large device memory at the tested chunk sizes.
- GPU acceleration becomes more favorable as the problem size grows.
- Peak memory can be characterized separately for:
  - small-to-medium fully materialized runs
  - large steady-state streaming workloads

### What they do not show

- They are not end-to-end production benchmarks on real genotype files.
- They do not measure the full cost of writing or storing a complete result table for the `20k / 1M / 2k` workload.
- They do not characterize BGEN decoding overhead, phenotype-table alignment overhead, or relatedness filtering overhead.
- They do not yet include repeated trials, confidence intervals, or cross-GPU reproducibility analysis.

## Recommended Citation In Repository Narratives

If this benchmark report is referenced in repository documentation or a manuscript supplement, the safest summary is:

> TorchGWAS was profiled on simulated linear GWAS workloads in both fully materialized and large-scale streaming settings. On an NVIDIA H100 80GB GPU, the large-scale float32 streaming path showed modest GPU memory requirements at the tested chunk size and improved estimated runtime relative to CPU, while smaller fully materialized workloads remained CPU-competitive and host-side memory remained a substantial component of total footprint.

## Reproduction Commands

Small-to-medium profile:

```bash
CUDA_VISIBLE_DEVICES=5 PYTHONPATH=src python benchmarks/profile_linear_memory.py \
  --shapes 1024x4096x8 2048x8192x8 4096x16384x16 \
  --devices cpu cuda \
  --chunk-size 1024 \
  --output results/benchmarks/linear_memory_profile.tsv
```

Large-scale streaming profile:

```bash
CUDA_VISIBLE_DEVICES=5 PYTHONPATH=src python benchmarks/profile_linear_streaming_large.py \
  --n-samples 20000 \
  --n-markers-total 1000000 \
  --n-traits 2000 \
  --chunk-size 4096 \
  --dtype float32 \
  --devices cpu cuda \
  --output results/benchmarks/linear_streaming_profile_large.tsv
```
