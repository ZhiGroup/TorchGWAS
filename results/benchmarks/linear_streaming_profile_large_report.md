# Large-Scale Linear Streaming Profile Report

## Goal

Profile TorchGWAS linear-scan memory behavior for a large simulated setting with:

- `20,000` samples
- `1,000,000` SNPs
- `2,000` phenotypes
- `8` covariates

This benchmark uses a streaming/chunked profile rather than a full end-to-end `run_linear_gwas()` result materialization, since storing the complete marker-by-trait output table at this scale is not practical for a repository benchmark.

## Setup

- Benchmark script: `benchmarks/profile_linear_streaming_large.py`
- Device visibility during GPU run:
  - `CUDA_VISIBLE_DEVICES=5`
  - visible GPU: `NVIDIA H100 80GB HBM3`
  - PyTorch CUDA version: `12.1`
- Chunk size: `4096`
- Warmup chunks: `1`
- Timed chunks: `2`
- Output table: `results/benchmarks/linear_streaming_profile_large.tsv`

## Interpretation

- `cpu_peak_rss_mb`: peak worker-process resident set size while processing repeated genotype chunks.
- `gpu_peak_allocated_mb`: peak PyTorch allocated memory on the visible GPU.
- `gpu_peak_reserved_mb`: peak PyTorch reserved memory on the visible GPU.
- `avg_chunk_seconds`: average runtime per timed chunk.
- `estimated_full_runtime_seconds`: extrapolated runtime for `1,000,000` SNPs at the same chunk size.

## Results

| device | n_samples | n_markers_total | n_traits | chunk_size | cpu_peak_rss_mb | gpu_peak_allocated_mb | gpu_peak_reserved_mb | avg_chunk_seconds | estimated_full_runtime_seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cpu | 20000 | 1000000 | 2000 | 4096 | 4428.297 | 0.0 | 0.0 | 2.073389 | 507.980265 |
| cuda | 20000 | 1000000 | 2000 | 4096 | 3470.035 | 2771.032 | 3936.0 | 1.515059 | 371.189553 |

## Takeaways

- At this scale, the CUDA path reduced estimated runtime from about `508 s` to about `371 s` under this benchmark setup.
- The GPU run still consumed substantial host memory because phenotype, covariate, and chunk staging arrays remain resident on CPU.
- The observed GPU peak was about `2.77 GB allocated` and `3.94 GB reserved`, which is modest relative to H100 capacity for this chunk size.
- CPU peak RSS remained in the multi-gigabyte range even for the GPU path, so current TorchGWAS scaling is not purely GPU-memory bound.
- Because this benchmark streams chunks and does not materialize the full marker-by-trait result table, it is best interpreted as a compute-kernel and steady-state memory profile, not a full production run footprint.
