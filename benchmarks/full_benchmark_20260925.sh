#!/usr/bin/env bash
# The 09-15 full benchmark (torchGWAS1.1 _paper_ours.sh) on the current tree:
# 35,365 samples x 8,931,083 variants, full binary output (beta + t), cold
# page cache, K in {512, 2048}, 3 rounds round-robin over K. The fixed arm
# repeats the 09-15 flags exactly (BED via the hard-call store, 27 covariates,
# chunk 4096, 16 workers, prefetch 32, one free GPU); the autotune arm leaves
# GPUs, layout, chunk and readers to autotune. Outputs are deleted after each
# run (146 GB at K=2048).
#   benchmarks/full_benchmark_20260925.sh <results-dir> [bench-root]
set -uo pipefail
R=${1:?results dir}
B=${2:-/data/zxie3/torchgwas_bench}
W=$B/full_benchmark_20260925
LADDER="${LADDER:-512 2048}"
ROUNDS="${ROUNDS:-3}"
mkdir -p "$R/runs" "$W"
echo "host $(hostname -s)  $(date)  load $(cut -d' ' -f1-3 /proc/loadavg)"
for round in $(seq 1 "$ROUNDS"); do
  for K in $LADDER; do
    arms="fixed autotune"; [ $((round % 2)) = 0 ] && arms="autotune fixed"
    for arm in $arms; do
      tag="r${round}_${arm}_k${K}"
      cache="$W/cache_$tag"; rm -rf "$cache"; mkdir -p "$cache"
      extra=""; [ "$arm" = autotune ] && extra="--autotune"
      "$PY" -u benchmarks/direct_cold_arm.py \
          --path "$B/full_hardcall.bed" --label "$tag" --traits "$K" \
          --covariates 27 --sumstats-format binary \
          --output-dir "$W/out_$tag" --genotype-cache-dir "$cache" \
          --hardcall-store "$B/hcstore/full" \
          --chunk-size 4096 --workers 16 --prefetch 32 $extra 2>"$R/runs/$tag.stderr" | tail -1 > "$R/runs/$tag.json"
      rm -rf "$W/out_$tag" "$cache"
      "$PY" - "$R/runs/$tag.json" <<'PYEOF'
import json, sys
try:
    r = json.load(open(sys.argv[1]))
    print("%-22s open %6.1f scan %7.1f wall %7.1f  out %6.1f GB  load %s->%s  devices %s  chunk %s"
          % (r["label"], r["open"], r["scan"], r["wall"], r["output_bytes"] / 1e9, r["load_before"],
             r["load_after"], r.get("devices_used"), (r.get("autotune_chunk") or {}).get("choice", r.get("resolved_chunk"))),
          flush=True)
except Exception as error:
    print("NO JSON", sys.argv[1], error, flush=True)
PYEOF
    done
  done
done
echo FULL_BENCHMARK_DONE
