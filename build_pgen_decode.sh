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
# O3 enables vectorization of the one-bit expansion loop. The complete native
# reader tests and matched cold-input scan comparison preserve exact genotype
# decoding and improve CPU service. OPT_LEVEL=2 retains the earlier build for
# comparisons. Neither setting enables fast-math or strict aliasing.
OPT_LEVEL=${OPT_LEVEL:-3}
case "$OPT_LEVEL" in
  2|3) ;;
  *) printf 'OPT_LEVEL must be 2 or 3\n' >&2; exit 2 ;;
esac
# The decoder reads bytes as uint8 and multi-byte little-endian ids.
FLAGS=(-std=c99 "-O$OPT_LEVEL" -fPIC -shared -fno-strict-aliasing -Wall -Wextra -Wpedantic)
if [ "${NATIVE_ARCH:-1}" = "1" ]; then
  # Not portable across machines; the lab hosts are homogeneous enough and
  # the library is rebuilt per host. Set NATIVE_ARCH=0 for a portable build.
  FLAGS+=(-march=native)
fi

# Replace the library atomically: a process already using the old inode must
# not see its mapped machine code truncated during a rebuild.
temporary_library=$(mktemp "$OUT/.libtorchgwas_pgen.XXXXXX.so")
trap 'rm -f -- "$temporary_library"' EXIT
"$CC" "${FLAGS[@]}" -o "$temporary_library" native/pgen_decode.c
mv -f -- "$temporary_library" "$OUT/libtorchgwas_pgen.so"

echo "built $OUT/libtorchgwas_pgen.so (-O$OPT_LEVEL)"
"$CC" --version | head -1
ls -l "$OUT/libtorchgwas_pgen.so"
