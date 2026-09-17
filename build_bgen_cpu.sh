#!/usr/bin/env bash
set -euo pipefail
CC=${CC:-cc}
output=src/torchgwas/native/libtorchgwas_bgen_cpu.so
temporary="${output}.tmp.$$"
trap 'rm -f -- "$temporary"' EXIT
"$CC" -O3 -std=c11 -shared -fPIC -Wall -Wextra \
  src/torchgwas/native/bgen_decode_cpu.c -o "$temporary" -lz -lm
mv -- "$temporary" "$output"
