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
