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
- Compute dtype: `float32`
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
| cpu | 20000 | 1000000 | 2000 | 4096 | 2767.645 | 0.0 | 0.0 | 1.539444 | 377.163868 |
| cuda | 20000 | 1000000 | 2000 | 4096 | 2546.039 | 1400.104 | 1972.0 | 1.244163 | 304.820029 |

## Takeaways

- At this scale, the CUDA path reduced estimated runtime from about `377 s` to about `305 s` under this benchmark setup.
- The GPU run still consumed substantial host memory because phenotype, covariate, and chunk staging arrays remain resident on CPU.
- The observed GPU peak was about `1.40 GB allocated` and `1.97 GB reserved`, which is modest relative to H100 capacity for this chunk size.
- CPU peak RSS remained in the multi-gigabyte range even for the GPU path, so current TorchGWAS scaling is still not purely GPU-memory bound, although the float32 path reduced host memory substantially relative to the earlier float64-style benchmark.
- Because this benchmark streams chunks and does not materialize the full marker-by-trait result table, it is best interpreted as a compute-kernel and steady-state memory profile, not a full production run footprint.
