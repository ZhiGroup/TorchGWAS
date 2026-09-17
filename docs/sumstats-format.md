# Binary summary-statistic format

TorchGWAS writes large association scans as tiled binary arrays rather than a
marker-by-trait text table. The store lives under `OUTPUT_DIR/sumstats/` and is
described by `manifest.json`.

With the default `--sumstats-fields beta+t`, the store contains float32 `beta`
and `t_stat` values for each marker-trait cell. The `t` setting writes only the
float32 t statistic and is intended for screening when effect sizes are not
required.

The manifest records the logical array shape, dtype, field names, tile files,
and write settings. Variant identifiers are written beside the arrays unless
`--no-sumstats-variant-ids` is supplied.

Use `--sumstats-format none` to discard association statistics after the scan.
This is useful for isolating computation but does not produce an analysis
result.

Indexed reduction output uses NumPy binary parts with a JSON manifest. Each
part records the selected variant and trait indices together with the requested
association fields.
