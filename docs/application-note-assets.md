# Application Note Assets

- Figure 1: workflow schematic assembled from package I/O and command structure
- Figure 2: `results/benchmarks/accuracy_plink.tsv` plotted by `scripts/make_figures.py`
- Figure 3: `results/benchmarks/runtime_scaling.tsv` plus `results/benchmarks/multivariate_case.tsv`

Each figure is generated from a tracked TSV to keep manuscript assets auditable.

Real-case execution helpers are stored under `real_case/t1_128/` after running:

```bash
python scripts/prepare_t1_real_case.py --output-dir real_case/t1_128
python scripts/write_real_case_commands.py --output-dir real_case/t1_128
```
