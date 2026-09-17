#!/usr/bin/env bash
# Build the direct PGEN record decoder.
#
# Plain C99 with no dependencies, built as a shared library and called through
# ctypes, matching how native/scan_statistics.cu is consumed. No Torch
# extension ABI is involved, so the result is independent of the interpreter
# and Torch build in use.
set -euo pipefail
cd "$(dirname "$0")"

OUT=${OUT:-.build-libs}
mkdir -p "$OUT"

CC=${CC:-gcc}
# -O2 with strict aliasing off: the decoder reads the same bytes as uint8 and
# as multi-byte little-endian ids. -fno-strict-aliasing is cheaper than
# fighting the optimiser with memcpy everywhere.
FLAGS=(-std=c99 -O2 -fPIC -shared -fno-strict-aliasing -Wall -Wextra -Wpedantic)
if [ "${NATIVE_ARCH:-1}" = "1" ]; then
  # Not portable across machines; the lab hosts are homogeneous enough and
  # the library is rebuilt per host. Set NATIVE_ARCH=0 for a portable build.
  FLAGS+=(-march=native)
fi

"$CC" "${FLAGS[@]}" -o "$OUT/libtorchgwas_pgen.so" native/pgen_decode.c

echo "built $OUT/libtorchgwas_pgen.so"
"$CC" --version | head -1
ls -l "$OUT/libtorchgwas_pgen.so"
