# Binary summary-statistic format (jagwas-dev)

TorchGWAS writes association scans as binary arrays under
`OUTPUT_DIR/sumstats/`, described by `manifest.json`. This line differs from
the public release: stores hold `beta` and `t_stat` (no `neg_log10_p`), and
streaming stores record where their variants come from instead of rewriting
the variant IDs.

## Dense stores

`--sumstats-fields beta+t` (the default) writes float32 `beta.f32` and
`tstat.f32`, `t` writes `tstat.f32` only. Both are `(n_variants, n_traits)`,
little-endian, C order. Excluded variants carry NaN. P-values are not stored:
they are the two-sided Student t tail of `t_stat` at the store's df, which is a
scalar, one value per trait (missing phenotypes), or a per-variant `df.f32`
column (missing genotype calls).

```python
from torchgwas.sumstats import open_binary_df, open_binary_sumstats

beta, t_stat, manifest = open_binary_sumstats("run/sumstats")
df = open_binary_df("run/sumstats")        # broadcasts against t_stat
```

A multi-GPU run writes one child store per variant shard (`layout:
variant_shards`) or per trait tile (`layout: trait_tiles`). The openers above
read them as one array.

## Indexed stores

Reductions write NumPy parts listed in the manifest
(`format: torchgwas-indexed-sumstats`):

- `reduce='significant'`: one row per passing pair, with `variant_index`,
  `trait_index`, `t_stat`, `df` and `beta` (unless `--sumstats-fields t`).
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
