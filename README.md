# TorchGWAS

TorchGWAS is a GPU-accelerated framework for covariate-adjusted linear
association testing across quantitative phenotypes. It streams genotype data
in bounded chunks and batches association calculations on the GPU.

## Features

- Linear association testing for one or more quantitative phenotypes.
- Internal sample alignment, covariate projection, and phenotype preprocessing.
- Direct support for NumPy arrays, PLINK 1 BED, PLINK 2 PGEN, and BGEN.
- Direct hard-call PGEN decoding, with `pgenlib` used for dosage records.
- Optional reusable zstd stores for BGEN, PGEN, and lossless packed BED input.
- Binary, tiled summary-statistic output for large scans.
- Optional device-side significance and joint-trait reductions.

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

## Quick start

Run the bundled toy workflow:

```bash
torchgwas demo --output-dir demo_run
```

Run a linear scan from aligned NumPy arrays:

```bash
torchgwas linear \
  --genotype examples/toy/genotype.npy \
  --phenotype examples/toy/pheno.npy \
  --covariates examples/toy/covar.npy \
  --marker-ids examples/toy/markers.tsv \
  --sample-ids examples/toy/samples.tsv \
  --output-dir linear_out
```

For cohort data, phenotype and covariate tables can be aligned to genotype
sample order by `IID`:

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

## Output

Full scans write a tiled binary store under `OUTPUT_DIR/sumstats/`. By default,
each marker-trait cell contains float32 `beta` and `t_stat` values. The output
manifest records array shapes, data types, field names, and tile layout.
`run.json` records the analysis configuration and timing, while `qc.json`
records sample, covariate, missingness, and invariant-variant decisions.

Use `--sumstats-fields t` for screening-only output containing t statistics,
or `--sumstats-format none` to run the scan without writing association
statistics. See [`docs/sumstats-format.md`](docs/sumstats-format.md).

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
