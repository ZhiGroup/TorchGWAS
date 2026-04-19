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
