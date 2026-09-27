#!/usr/bin/env bash
# Build the scan's device -log10 P stages (tails.py) with AOTInductor for the
# GPU architecture of $DEVICE (default cuda:0). Without a build the scan
# compiles the stages in every process (torch.compile), and without a compiler
# it runs them eagerly. The output lands in .build-libs/device_tail_<key>/,
# keyed by torch, CUDA, architecture and the stage source, so a changed stage
# is rebuilt rather than reused.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python}
# The generated C++ wrapper includes cuda.h. Where the toolkit is not on the
# default include path (the A100 host: /usr/local/cuda-12.9/targets/x86_64-linux/include),
# pass it in CPATH; the pip CUDA runtime headers are incomplete for this.
DEVICE=${DEVICE:-cuda:0}
PYTHONPATH="src${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -c "import sys; from torchgwas.tails import build_device_tail; print(build_device_tail(sys.argv[1]))" "$DEVICE"
