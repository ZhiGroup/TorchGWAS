#!/usr/bin/env bash
set -euo pipefail
NVCC=${NVCC:-/usr/local/cuda/bin/nvcc}
output=src/torchgwas/native/libtorchgwas_scan.so
temporary="${output}.tmp.$$"
trap 'rm -f -- "$temporary"' EXIT
# Turing, Ampere, Hopper, plus PTX so anything newer than Hopper JITs rather
# than failing. Built for sm_80/sm_90 only, this library loads fine on any card
# and then dies at launch with "no kernel image is available for execution on
# the device" -- which is what an RTX 2080 Ti (sm_75) does, and what nineteen
# tests on one showed. `scan_gpu.available()` now probes a real launch so the
# fallback is chosen correctly even against a library built without these.
"$NVCC" -O3 -std=c++17 -shared -Xcompiler -fPIC --fmad=false \
  -gencode arch=compute_75,code=sm_75 \
  -gencode arch=compute_80,code=sm_80 \
  -gencode arch=compute_90,code=sm_90 \
  -gencode arch=compute_90,code=compute_90 \
  src/torchgwas/native/scan_statistics.cu -o "$temporary"
mv -- "$temporary" "$output"
