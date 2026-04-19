# Real-Case Assets

This directory contains optional local workflow assets and lightweight validation materials.
It is not the default starting point for the public TorchGWAS interface; the tracked toy dataset under `examples/toy/` is.

## Tracked Content

- `t1_128/`
  Command templates, column definitions, and metadata for a local T1 workflow

## Not Tracked

Large intermediate tables generated for local execution are intentionally excluded from version control.
This includes smoke-test output directories such as `*_smoke/`, `*_smoke_fast/`, and `*_smoke_small/`.
If local T1 assets are needed, the full `t1_128` phenotype and covariate TSV files can be regenerated with:

```bash
python scripts/prepare_t1_real_case.py --output-dir real_case/t1_128
```

This keeps the public repository code-oriented while preserving optional local workflows.
