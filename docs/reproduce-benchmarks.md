# Benchmark Reproduction

## Included summaries

Generated benchmark tables are stored in `results/benchmarks/`.

## PLINK concordance from an external baseline directory

The command below expects an existing local baseline output directory and is not required for the toy workflow.

```bash
python benchmarks/summarize_plink_accuracy.py \
  --fast-gwas-root /path/to/fast_GWAS \
  --trait-indices 0,17,63,127 \
  --output results/benchmarks/accuracy_plink.tsv
```

## Historical runtime points from the original notebook

```bash
python benchmarks/materialize_runtime_table.py \
  --output results/benchmarks/runtime_scaling.tsv
```

## Native BED reader scaling

This benchmark creates a temporary variant-major BED file, reads every variant, and reports wall time, process CPU time, throughput, and peak RSS for each worker count:

```bash
PYTHONPATH=src python benchmarks/benchmark_bed_reader.py \
  --n-samples 35365 \
  --n-variants 12000 \
  --chunk-size 1000 \
  --workers 1,2,4,8
```

The A100 lab-host measurements and the corresponding public-main `pandas-plink` baseline are recorded in `results/benchmarks/bed_reader_parallel_remote_20260910.tsv`.

## Calibrated stage runtime predictor

The predictor separates compressed disk read, CPU Zstandard decode, host-to-device transfer, and GPU compute. The hardware model uses measured FP32 throughput, CPU frequency, decode thread count, and H2D bandwidth; `--pipeline-contention-factor` comes from a short joint-stage calibration and captures NUMA/PCIe contention.

The following values reproduce the measured torchGWAS1.1 A100 host ceiling (about 6.7 seconds for the scan itself):

```bash
PYTHONPATH=src python benchmarks/predict_runtime.py \
  --model hardware \
  --n-variants 8086101 \
  --n-samples 22250 \
  --n-traits 128 \
  --covariate-rank 27 \
  --gpu-fp32-tflops 19.1 \
  --gpu-count 2 \
  --cpu-frequency-ghz 2.3 \
  --decode-threads 40 \
  --disk-gbps 100 \
  --h2d-gbps 42.5 \
  --pipeline-contention-factor 1.584
```

Use measured sustained values for a new system. Peak vendor specifications alone do not capture filesystem cache state, PCIe topology, NUMA placement, or decode/H2D contention.

## Toy multivariate case summary

This benchmark is kept as an exploratory extension rather than the main reported workflow.

```bash
python benchmarks/run_multivariate_case.py \
  --output results/benchmarks/multivariate_case.tsv
```

## Simulated CPU/GPU memory profile for linear GWAS

This benchmark runs linear TorchGWAS on simulated matrices and records:

- peak CPU resident set size for each isolated subprocess
- peak GPU allocated and reserved memory reported by PyTorch
- runtime for each configuration

```bash
python benchmarks/profile_linear_memory.py \
  --shapes 1024x4096x8 2048x8192x8 \
  --devices cpu cuda \
  --chunk-size 1024 \
  --output results/benchmarks/linear_memory_profile.tsv
```

## Large-scale streaming profile at 20k samples / 1M SNPs / 2k traits

For a more realistic large-scale memory benchmark, use the streaming profile script.
This benchmark measures peak steady-state memory on repeated chunks and extrapolates runtime to the full marker count without materializing the complete result table.

```bash
python benchmarks/profile_linear_streaming_large.py \
  --n-samples 20000 \
  --n-markers-total 1000000 \
  --n-traits 2000 \
  --chunk-size 4096 \
  --devices cpu cuda \
  --output results/benchmarks/linear_streaming_profile_large.tsv
```
