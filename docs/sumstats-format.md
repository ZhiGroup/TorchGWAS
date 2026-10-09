# Binary summary-statistic format (jagwas-dev)

TorchGWAS writes association scans as binary arrays under
`OUTPUT_DIR/sumstats/`, described by `manifest.json`. The stored fields follow
the public release. This line adds a per-variant df sidecar for missing
genotype calls, and streaming stores record where their variants come from
instead of rewriting the variant IDs.

## Dense stores

With `--sumstats-fields beta+t` (the default) a store holds float32
`beta.f32`, `tstat.f32` and `neglog10p.f32`, 12 bytes per variant-trait cell.
With `t` it holds `tstat.f32` and `neglog10p.f32`, 8 bytes per cell. All three
arrays are `(n_variants, n_traits)`, little-endian, C order, and excluded
variants carry NaN in each.

`neg_log10_p` is -log10 of the exact two-sided Student-t tail. It is computed
in FP64 on the scan device and stored as float32. Raw P values are not stored,
because they underflow for strong associations. The manifest (version 2)
records the df the tail was taken at. That is a scalar, one value per trait
for missing phenotypes, or a per-variant `df.f32` column for missing genotype
calls.

**Missing phenotypes** (`missing_phenotype`, run metadata):
- The default is `'exact'`, and `'drop_subject'` for JAGWAS; run metadata
  records the policy used.
- `'drop_subject'`: a sample with any missing or outlier-masked value is left
  out of every trait, so the store is that of a run without those samples;
  `dropped_subjects` counts them.
- `'impute'` (the release): the trait's missing values take its mean; t is
  the imputed panel's t times sqrt(trait_df / df) and `neg_log10_p` is at
  pair df variant_df x trait_df / df. The manifest keeps each trait's df.
- `'exact'`: each trait is tested on its own observed samples, as
  plink2's `--glm` does per phenotype. The intercept and covariates are
  refitted on those samples, t, beta and `neg_log10_p` are that
  complete-case OLS, and the df is the pair's own count less the covariate
  rank there and the genotype. The manifest's per-trait df is the trait's
  observed count less rank and genotype, and a missing call lowers a pair's
  df by one more.
- Calls missing inside a trait's samples keep the genotype convention:
  centred at the variant's observed mean.
- JAGWAS takes `'drop_subject'` or `'impute'`: its joint test needs one
  sample set.

```python
from torchgwas.sumstats import open_binary_df, open_binary_sumstats

beta, t_stat, neg_log10_p, manifest = open_binary_sumstats("run/sumstats")
df = open_binary_df("run/sumstats")        # broadcasts against t_stat
```

`neg_log10_p` is None for a store written before it was stored (version 1).

A multi-GPU run writes one child store per variant shard (`layout:
variant_shards`) or per trait tile (`layout: trait_tiles`). The openers above
read them as one array.

## Indexed stores

Reductions write NumPy parts listed in the manifest
(`format: torchgwas-indexed-sumstats`):

- `reduce='significant'` and `p_value_threshold`: one row per passing pair,
  with `variant_index`, `trait_index`, `t_stat`, `neg_log10_p` (float32),
  `beta` (unless `--sumstats-fields t`) and, for significant pairs, the
  pair's `df`.
- `reduce='min-p'`: one row per valid variant, holding the trait with the
  smallest p: the same columns as significant pairs, `df` included, and
  `manifest["reduction"] == "min-p"`. `neg_log10_p` is the exact tail at
  that pair's df, computed on the device. With a complete phenotype panel
  every trait of a variant shares its df, so the winner is the largest |t|.
  With missing phenotypes df differs by trait (full output's convention:
  variant df times trait_df / df, t rescaled by sqrt(trait_df / df)), so
  the winner is chosen by the tail itself.
- `reduce='jagwas'`: one row per variant, with `variant_index` and `chi2`, a
  column per JAGWAS group. The manifest `df` is the joint test's df, a list
  with one entry per group when groups are used, and `groups` names them.

```python
from torchgwas.sumstats_indexed import open_indexed_sumstats

manifest, parts = open_indexed_sumstats("run/sumstats")
```

## Variant identity

Dense row `i` is input variant `variant_offset + i`, and indexed rows carry
`variant_index` in input order. IDs, positions and alleles are the genotype's
own metadata. By default a streaming store does not copy them. Its manifest
records the input:

```json
"variant_source": {"format": "pgen", "genotype": "/abs/path/input.pgen",
                   "load": {"pvar": "..."}, "variant_offset": 0,
                   "n_variants": 8090000, "id_digest": "...",
                   "id_digest_scheme": "sha256 over the variant count and 4097 evenly spaced IDs, newline-separated"}
```

`store_variants` resolves a store's IDs and metadata through that genotype
(with its metadata cache). It refuses a genotype whose variant list no longer
matches the digest:

```python
from torchgwas.variant_source import store_variant_ids, store_variants

ids = store_variant_ids("run/sumstats")                       # the recorded genotype
ids, metadata = store_variants("run/sumstats", "moved/input.pgen")
```

`--sumstats-variant-ids` (or `sumstats_variant_ids=True`) also embeds the IDs
(`variant_ids.txt` beside a dense store, `variant_ids.npy` and
`variant_metadata.npz` in an indexed one) for a store meant to travel without
its genotype. Embedding costs 356 MB per indexed store at 8.09M variants. A
store written from an in-memory genotype always embeds them, and so does the
in-memory result path, whose rows follow the QC-kept markers.

Use `--sumstats-format none` to discard statistics after the scan (for timing
the computation only).

Results returned without `output_dir` report `-log10_p` for each row, from the
same exact tail in FP64, and no raw `p_value`.
