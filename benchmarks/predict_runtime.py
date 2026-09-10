from __future__ import annotations

import argparse
import json

from torchgwas.runtime import predict_runtime, predict_runtime_from_hardware


def main() -> int:
    parser = argparse.ArgumentParser(description="Predict TorchGWAS runtime from calibrated pipeline rates")
    parser.add_argument("--model", choices=["stage", "hardware"], default="hardware")
    parser.add_argument("--n-variants", type=int, required=True)
    parser.add_argument("--n-samples", type=int, required=True)
    parser.add_argument("--n-traits", type=int, required=True)
    parser.add_argument("--covariate-rank", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=2500)
    parser.add_argument("--decoded-bytes-per-value", type=float, default=1.0)
    parser.add_argument("--compression-ratio", type=float, default=12.46)
    parser.add_argument("--disk-gbps", type=float, required=True)
    parser.add_argument("--decode-gbps", type=float, default=None)
    parser.add_argument("--h2d-gbps", type=float, required=True)
    parser.add_argument("--measured-gemm-tflops", type=float, default=None)
    parser.add_argument("--gpu-fp32-tflops", type=float, default=None, help="Measured per-GPU FP32 peak")
    parser.add_argument("--gpu-count", type=int, default=1)
    parser.add_argument("--gemm-efficiency", type=float, default=12.9 / 19.1)
    parser.add_argument("--cpu-frequency-ghz", type=float, default=None)
    parser.add_argument("--decode-threads", type=int, default=None)
    parser.add_argument("--decode-gbps-per-core-at-reference", type=float, default=1.52)
    parser.add_argument("--reference-cpu-frequency-ghz", type=float, default=2.3)
    parser.add_argument("--decode-thread-efficiency", type=float, default=55.3 / (1.52 * 40))
    parser.add_argument("--non-gemm-ms-per-chunk-per-gpu", type=float, default=0.818)
    parser.add_argument("--pipeline-contention-factor", type=float, default=1.0)
    parser.add_argument("--non-gemm-ms-per-chunk", type=float, default=0.0)
    parser.add_argument("--contention-factor", type=float, default=1.0)
    parser.add_argument("--setup-seconds", type=float, default=0.0)
    parser.add_argument("--postprocess-seconds", type=float, default=0.0)
    parser.add_argument("--joint-decode-chunks-per-second", type=float, default=None)
    parser.add_argument("--joint-h2d-chunks-per-second", type=float, default=None)
    parser.add_argument("--joint-compute-chunks-per-second", type=float, default=None)
    parser.add_argument("--output", default=None, help="Optional JSON output path")
    args = parser.parse_args()
    common = {
        "n_variants": args.n_variants,
        "n_samples": args.n_samples,
        "n_traits": args.n_traits,
        "covariate_rank": args.covariate_rank,
        "chunk_size": args.chunk_size,
        "compression_ratio": args.compression_ratio,
        "disk_gbps": args.disk_gbps,
        "h2d_gbps": args.h2d_gbps,
        "setup_seconds": args.setup_seconds,
        "postprocess_seconds": args.postprocess_seconds,
    }
    if args.model == "hardware":
        missing = [
            name
            for name, value in {
                "--gpu-fp32-tflops": args.gpu_fp32_tflops,
                "--cpu-frequency-ghz": args.cpu_frequency_ghz,
                "--decode-threads": args.decode_threads,
            }.items()
            if value is None
        ]
        if missing:
            parser.error("hardware model requires " + ", ".join(missing))
        result = predict_runtime_from_hardware(
            **common,
            gpu_fp32_tflops=args.gpu_fp32_tflops,
            gpu_count=args.gpu_count,
            gemm_efficiency=args.gemm_efficiency,
            cpu_frequency_ghz=args.cpu_frequency_ghz,
            decode_threads=args.decode_threads,
            decode_gbps_per_core_at_reference=args.decode_gbps_per_core_at_reference,
            reference_cpu_frequency_ghz=args.reference_cpu_frequency_ghz,
            decode_thread_efficiency=args.decode_thread_efficiency,
            non_gemm_ms_per_chunk_per_gpu=args.non_gemm_ms_per_chunk_per_gpu,
            pipeline_contention_factor=args.pipeline_contention_factor,
        )
    else:
        if args.decode_gbps is None or args.measured_gemm_tflops is None:
            parser.error("stage model requires --decode-gbps and --measured-gemm-tflops")
        result = predict_runtime(
            **common,
            decoded_bytes_per_value=args.decoded_bytes_per_value,
            decode_gbps=args.decode_gbps,
            measured_gemm_tflops=args.measured_gemm_tflops,
            non_gemm_ms_per_chunk=args.non_gemm_ms_per_chunk,
            contention_factor=args.contention_factor,
            joint_decode_chunks_per_second=args.joint_decode_chunks_per_second,
            joint_h2d_chunks_per_second=args.joint_h2d_chunks_per_second,
            joint_compute_chunks_per_second=args.joint_compute_chunks_per_second,
        )
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path

        Path(args.output).write_text(payload)
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
