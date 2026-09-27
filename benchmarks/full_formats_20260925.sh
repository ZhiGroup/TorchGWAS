#!/usr/bin/env bash
# Full output (binary beta + t) at full scale for every input format, cold
# page cache: fixed arm = the 09-15 benchmark flags (chunk 4096, 16 workers,
# prefetch 32, one free GPU), autotune arm = autotune=True. K in {512, 2048},
# arms alternate order by round. Outputs are deleted after each run.
#   benchmarks/full_formats_20260925.sh <results-dir> <formats-file> [work-dir]
# formats-file: one "name|direct_cold_arm arguments" per line, e.g.
#   pgen_hardcall|--path /data/.../full_hardcall.pgen --pgen-mode hardcall
# ARMS (default "fixed autotune") and FIXED_CHUNK (default 4096) select arms;
# a fixed arm at another chunk is tagged fixed<chunk>.
set -uo pipefail
R=${1:?results dir}
F=${2:?formats file}
W=${3:-/data/zxie3/torchgwas_bench/full_formats_20260925}
LADDER="${LADDER:-512 2048}"
ROUNDS="${ROUNDS:-1}"
ARMS="${ARMS:-fixed autotune}"
FIXED_CHUNK="${FIXED_CHUNK:-4096}"
PREFETCH="${PREFETCH:-32}"  # fixed arm ring depth; a non-default one is tagged p<depth>
mkdir -p "$R/runs" "$W"
echo "host $(hostname -s)  $(date)  load $(cut -d' ' -f1-3 /proc/loadavg)"
for round in $(seq 1 "$ROUNDS"); do
  while IFS='|' read -r name source; do
    [ -z "$name" ] && continue
    for K in $LADDER; do
      arms="$ARMS"; [ $((round % 2)) = 0 ] && arms=$(echo $ARMS | tr ' ' '\n' | tac | tr '\n' ' ')
      for arm in $arms; do
        label=$arm; [ "$arm" = fixed ] && [ "$FIXED_CHUNK" != 4096 ] && label="fixed$FIXED_CHUNK"
        [ "$arm" = fixed ] && [ "$PREFETCH" != 32 ] && label="${label}p$PREFETCH"
        tag="${name}_r${round}_${label}_k${K}"
        cache="$W/cache_$tag"; rm -rf "$cache"; mkdir -p "$cache"
        extra=""; [ "$arm" = autotune ] && extra="--autotune"
        # shellcheck disable=SC2086
        "$PY" -u benchmarks/direct_cold_arm.py $source --label "$tag" --traits "$K" \
            --covariates 27 --sumstats-format binary \
            --output-dir "$W/out_$tag" --genotype-cache-dir "$cache" \
            --chunk-size "$FIXED_CHUNK" --workers 16 --prefetch "$PREFETCH" $extra 2>"$R/runs/$tag.stderr" | tail -1 > "$R/runs/$tag.json"
        rm -rf "$W/out_$tag" "$cache"
        "$PY" - "$R/runs/$tag.json" <<'PYEOF'
import json, sys
try:
    r = json.load(open(sys.argv[1]))
    print("%-34s open %6.1f scan %7.1f wall %7.1f  out %6.1f GB  load %s->%s  devices %s  chunk %s  backend %s/%s"
          % (r["label"], r["open"], r["scan"], r["wall"], r["output_bytes"] / 1e9, r["load_before"],
             r["load_after"], r.get("devices_used"), (r.get("autotune_chunk") or {}).get("choice", r.get("resolved_chunk")),
             r.get("statistics_backend"), r.get("decode_backend")), flush=True)
except Exception as error:
    print("NO JSON", sys.argv[1], error, flush=True)
PYEOF
      done
    done
  done < "$F"
done
echo FULL_FORMATS_DONE
