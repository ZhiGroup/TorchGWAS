# TorchGWAS

TorchGWAS is a GPU-accelerated framework for covariate-adjusted linear
association testing across quantitative phenotypes. It streams genotype data
in bounded chunks and batches association calculations on the GPU.

## Statistical scope

TorchGWAS implements a fixed-effect linear association model and is intended
for unrelated or appropriately relatedness-filtered individuals. When related
individuals are retained, a mixed-model association method should be used
instead.

## Features

- Linear association testing for one or more quantitative phenotypes.
- Internal sample alignment, covariate projection, and phenotype preprocessing.
- Direct support for PLINK 1 BED, PLINK 2 PGEN, and BGEN.
- Direct hard-call PGEN decoding, with `pgenlib` used for dosage records.
- Optional reusable zstd stores for BGEN, PGEN, and lossless packed BED input.
- Binary, tiled summary-statistic output for large scans.
- Optional device-side significance and joint-trait reductions.

## Benchmark

The benchmark scanned 35,365 subjects, 8,931,083 variants, and 512
phenotypes, producing 36.68 GB of binary beta and t-statistic results.

Hardware: NVIDIA H100 80GB HBM3 GPU (driver 560.35.05), two Intel Xeon
Gold 6442Y CPUs (48 physical cores total), and 1 TiB DDR5 system memory.
Measured cold-read throughput was 5.71-5.94 GB/s, and application-level
binary-output throughput was 1.41-1.63 GB/s.

| Input format | Genotype representation | Wall time | Scan time | Peak RSS | Peak VRAM allocated / reserved | Input size |
|---|---|---:|---:|---:|---:|---:|
| BED | Hard calls | 28.46 s | 19.90 s | 4.94 GB | 0.90 / 1.37 GB | 78.97 GB |
| PGEN | Hard calls | 29.13 s | 20.35 s | 9.67 GB | 2.47 / 3.16 GB | 35.93 GB |
| Zstd hard-call store | Hard calls | 44.01 s | 35.90 s | 5.20 GB | 0.90 / 1.37 GB | 10.37 GB |
| BGEN | Dosages | 51.48 s | 28.63 s | 6.49 GB | 5.07 / 8.05 GB | 86.41 GB |
| PGEN | Dosages | 58.95 s | 50.12 s | 38.66 GB | 19.84 / 20.54 GB | 73.16 GB |

BED, PGEN, and BGEN values are medians of two full-scale runs. The zstd
hard-call result is the median of three full-scale runs.

## Installation

TorchGWAS requires Python 3.10 or newer.

```bash
git clone https://github.com/ZhiGroup/TorchGWAS.git
cd TorchGWAS
python -m pip install -e .
```

Install the optional PGEN dependency when reading dosage PGEN records:

```bash
python -m pip install -e '.[pgen]'
```

Install the optional reference BGEN dependency when that fallback is needed:

```bash
python -m pip install -e '.[bgen]'
```

The default CUDA scan uses PyTorch. Optional native decoders and scan kernels
can be built in place:

```bash
bash build_bgen_cpu.sh
bash build_direct_bgen.sh
bash build_direct_scan.sh
bash build_pgen_decode.sh
```

The CUDA build scripts require a CUDA toolkit. The direct GPU BGEN decoder also
requires nvCOMP. Native shared libraries are build products and are not stored
in Git.

## Autotuning status

Association chunk size is not fixed when `--chunk-size` is omitted. On CUDA,
the current release searches from 128 to 4,096 variants for the largest chunk
whose modeled buffers fit within 85% of free device memory, then aligns down
to the source frame size when applicable. The user can always override this
choice with `--chunk-size`. Source-specific read and decode tile defaults
remain separate from the association chunk.

This is feasibility sizing, not validated performance autotuning. The
analytical runtime calculator and `--pipeline-profile` interface remain
development-only. A later release is intended to expose that work as an
autotuner for association chunk size, read/decode tile size, queue depth, and
related execution parameters. It should not yet be treated as a supported
wall-clock predictor.

The source/resource equations, byte accounting, overlap schedule, and current
limitations of the research predictor are documented in
[`docs/runtime-predictor.md`](docs/runtime-predictor.md).

## Quick start

Run the bundled toy workflow:

```bash
torchgwas demo --output-dir demo_run
```

Phenotype and covariate tables are aligned to genotype sample order by `IID`:

```bash
torchgwas linear \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --output-dir linear_out
```

The same interface accepts `.pgen` and `.bgen` input. For BGEN files with
external sample identifiers, add `--sample-file /path/to/study.sample`.

## Unrelated-sample workflow

For fixed-effect analyses that exclude related individuals, the repository
includes `scripts/run_linear_unrelated.py`. The workflow uses PLINK 2.0
`--king-cutoff` to construct an unrelated genotype subset, intersects the
retained sample identifiers with the phenotype and covariate tables, and then
runs TorchGWAS on the aligned cohort.

Our analysis used a KING kinship cutoff of `(0.5)^4.5`, or approximately
`0.044194`, to exclude third-degree or closer relationships. Pass this value
explicitly:

```bash
PYTHONPATH=src python scripts/run_linear_unrelated.py \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --king-cutoff 0.04419417382415922 \
  --output-dir unrelated_linear_run
```

The workflow writes the retained sample list, aligned genotype, phenotype and
covariate inputs, and the TorchGWAS results. If related individuals are retained
instead, use a mixed-model association method.

## Output

Full scans write a tiled binary store under `OUTPUT_DIR/sumstats/`. By default,
each marker-trait cell contains float32 `beta`, `t_stat`, and `neg_log10_p`
values. `neg_log10_p` is computed in float64 from the exact two-sided Student-t
tail before being stored as float32; raw P values are not stored because they
are redundant and can underflow for strong associations. The output manifest
records array shapes, data types, field names, and tile layout.
`run.json` records the analysis configuration and timing, while `qc.json`
records sample, covariate, missingness, and invariant-variant decisions.
Runs without `--output-dir` return `-log10_p` in each result row and omit the
redundant raw `p_value` field.

Use `--sumstats-fields t` for screening-only output containing `t_stat` and
`neg_log10_p`. Use `--sumstats-format none` to run the scan without writing
association statistics. See [`docs/sumstats-format.md`](docs/sumstats-format.md).

## Genotype conventions

- BED reports dosage for BIM allele A2.
- PGEN reports ALT1 dosage.
- BGEN reports expected dosage for allele 2.
- Sample order is taken from the genotype source unless an explicit sample
  vector is supplied.
- Missing and invariant variants are handled during the scan and summarized in
  `qc.json`.

## Testing

```bash
python -m pytest -q
```

Native-library tests are skipped when the corresponding optional library is
not built.

## Repository layout

- `src/torchgwas/`: package and command-line implementation
- `tests/`: package regression tests
- `examples/`: tracked toy inputs
- `scripts/`: supported workflow helpers
- `native/`: source for the optional direct PGEN decoder

Benchmark code and outputs, calculator experiments, experiment logs,
manuscript files, local server settings, and generated native libraries are
intentionally not part of this repository.

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff).
