#!/usr/bin/env bash
set -euo pipefail
# Build the optional CUDA BGEN decoder without depending on the torch C++ ABI.
PY=${PY:-python}
NVROOT=${NVROOT:-$("$PY" -c 'import pathlib,nvidia.libnvcomp; print(pathlib.Path(nvidia.libnvcomp.__file__).parent)')}
NVCC=${NVCC:-/usr/local/cuda/bin/nvcc}
output=src/torchgwas/native/libtorchgwas_bgen.so
temporary="${output}.tmp.$$"
trap 'rm -f -- "$temporary"' EXIT
# Turing, Ampere, Hopper, plus PTX for anything newer -- see the note in
# build_direct_scan.sh: sm_80/sm_90 alone loads on any card and then fails at
# launch, which is a much worse failure than not loading.
"$NVCC" -O3 -std=c++17 -shared -Xcompiler -fPIC \
  -gencode arch=compute_75,code=sm_75 \
  -gencode arch=compute_80,code=sm_80 \
  -gencode arch=compute_90,code=sm_90 \
  -gencode arch=compute_90,code=compute_90 \
  -I"$NVROOT/include" src/torchgwas/native/bgen_decode.cu \
  -L"$NVROOT/lib64" -l:libnvcomp.so.5 -Xlinker -rpath -Xlinker "$NVROOT/lib64" \
  -o "$temporary"
mv -- "$temporary" "$output"
