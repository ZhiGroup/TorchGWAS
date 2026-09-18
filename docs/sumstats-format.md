# Binary summary-statistic format

TorchGWAS writes large association scans as tiled binary arrays rather than a
marker-by-trait text table. The store lives under `OUTPUT_DIR/sumstats/` and is
described by `manifest.json`.

With the default `--sumstats-fields beta+t`, the store contains float32 `beta`,
`t_stat`, and `neg_log10_p` values for each marker-trait cell (12 bytes per
cell). `neg_log10_p` is calculated in float64 with the exact two-sided
Student-t tail and then stored as float32. Raw P values are not stored because
they are redundant with `neg_log10_p` and underflow for sufficiently strong
associations. The `t` setting writes `t_stat` and `neg_log10_p` (8 bytes per
cell) and is intended for screening when effect sizes are not required.
In-memory result rows use the display name `-log10_p` for the same quantity and
do not include a redundant raw `p_value` field.

The dense files are `beta.f32` (unless `--sumstats-fields t` is used),
`tstat.f32`, and `neglog10p.f32`. All have the marker-by-trait shape recorded
in the version-2 manifest and can be memory-mapped with:

```python
from torchgwas.sumstats import open_binary_sumstats

beta, t_stat, neg_log10_p, manifest = open_binary_sumstats("run/sumstats")
```

The manifest records the logical array shape, dtype, field names, tile files,
and write settings. Variant identifiers are written beside the arrays unless
`--no-sumstats-variant-ids` is supplied.

Use `--sumstats-format none` to discard association statistics after the scan.
This is useful for isolating computation but does not produce an analysis
result.

Indexed linear output uses NumPy binary parts with a JSON manifest. Each part
records the selected variant and trait indices, `t_stat`, `neg_log10_p`, and
`beta` when requested.
