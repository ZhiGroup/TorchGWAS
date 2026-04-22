# Linear Memory Profile Report

## Setup

- Benchmark script: `benchmarks/profile_linear_memory.py`
- Analysis: linear TorchGWAS on simulated dense quantitative phenotypes
- Runtime environment:
  - PyTorch CUDA availability: `True`
  - CUDA version: `12.1`
  - Visible GPU during run: `NVIDIA H100 80GB HBM3`
  - GPU selection: `CUDA_VISIBLE_DEVICES=5`
- Output table: `results/benchmarks/linear_memory_profile.tsv`

## Configurations

The benchmark was run with:

- `n_covariates = 8`
- `chunk_size = 1024`
- shapes:
  - `1024 x 4096 x 8`
  - `2048 x 8192 x 8`
  - `4096 x 16384 x 16`

The shape convention is:

- `n_samples x n_markers x n_traits`

## Results

| device | n_samples | n_markers | n_traits | cpu_peak_rss_mb | gpu_peak_allocated_mb | gpu_peak_reserved_mb | runtime_seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| cpu | 1024 | 4096 | 8 | 595.836 | 0.0 | 0.0 | 0.337398 |
| cuda | 1024 | 4096 | 8 | 899.848 | 112.196 | 166.0 | 0.980079 |
| cpu | 2048 | 8192 | 8 | 800.891 | 0.0 | 0.0 | 0.726920 |
| cuda | 2048 | 8192 | 8 | 1022.129 | 192.258 | 306.0 | 1.089599 |
| cpu | 4096 | 16384 | 16 | 1665.492 | 0.0 | 0.0 | 2.157916 |
| cuda | 4096 | 16384 | 16 | 1682.008 | 352.758 | 354.0 | 1.877989 |

## Notes

- `cpu_peak_rss_mb` is the subprocess peak resident set size reported by `resource.getrusage()`. It reflects the whole Python worker process, not only genotype/phenotype arrays.
- `gpu_peak_allocated_mb` and `gpu_peak_reserved_mb` come from PyTorch peak memory counters after `torch.cuda.reset_peak_memory_stats()`.
- The GPU runs in this benchmark still allocate large host-side arrays, so GPU execution should not be interpreted as lowering total CPU memory footprint.
- In these tested shapes, the CUDA path showed lower runtime only for the largest configuration. For the smaller configurations, dispatch and transfer overhead offset the benefit.
- These results are simulated-memory reference points, not full real-cohort production benchmarks.
