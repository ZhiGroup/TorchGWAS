from __future__ import annotations

import math


def predict_runtime(
    *,
    n_variants: int,
    n_samples: int,
    n_traits: int,
    covariate_rank: int,
    chunk_size: int,
    decoded_bytes_per_value: float,
    compression_ratio: float,
    disk_gbps: float,
    decode_gbps: float,
    h2d_gbps: float,
    measured_gemm_tflops: float,
    non_gemm_ms_per_chunk: float = 0.0,
    contention_factor: float = 1.0,
    setup_seconds: float = 0.0,
    postprocess_seconds: float = 0.0,
    joint_decode_chunks_per_second: float | None = None,
    joint_h2d_chunks_per_second: float | None = None,
    joint_compute_chunks_per_second: float | None = None,
) -> dict[str, float | int | str]:
    """Predict a scan from measured stage rates.

    Bandwidths and GEMM throughput must be sustained measurements at the actual
    workload shape. If all three joint chunk rates are provided, they take
    precedence over the isolated-stage roofline and directly model contention.
    """

    positive_ints = {
        "n_variants": n_variants,
        "n_samples": n_samples,
        "n_traits": n_traits,
        "chunk_size": chunk_size,
    }
    for name, value in positive_ints.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    if covariate_rank < 0:
        raise ValueError("covariate_rank must be nonnegative")
    for name, value in {
        "decoded_bytes_per_value": decoded_bytes_per_value,
        "compression_ratio": compression_ratio,
        "disk_gbps": disk_gbps,
        "decode_gbps": decode_gbps,
        "h2d_gbps": h2d_gbps,
        "measured_gemm_tflops": measured_gemm_tflops,
        "contention_factor": contention_factor,
    }.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")

    n_chunks = math.ceil(n_variants / chunk_size)
    decoded_bytes = float(n_variants * n_samples) * decoded_bytes_per_value
    compressed_bytes = decoded_bytes / compression_ratio
    disk_seconds = compressed_bytes / (disk_gbps * 1e9)
    decode_seconds = decoded_bytes / (decode_gbps * 1e9)
    h2d_seconds = decoded_bytes / (h2d_gbps * 1e9)
    gemm_flops = 2.0 * n_variants * n_samples * (n_traits + covariate_rank)
    compute_seconds = gemm_flops / (measured_gemm_tflops * 1e12)
    compute_seconds += n_chunks * non_gemm_ms_per_chunk / 1000.0

    joint_rates = (
        joint_decode_chunks_per_second,
        joint_h2d_chunks_per_second,
        joint_compute_chunks_per_second,
    )
    if all(rate is not None for rate in joint_rates):
        rates = [float(rate) for rate in joint_rates]
        if any(rate <= 0 for rate in rates):
            raise ValueError("joint chunk rates must be positive")
        scan_seconds = n_chunks / min(rates)
        model = "measured_joint_pipeline"
    elif any(rate is not None for rate in joint_rates):
        raise ValueError("provide all three joint chunk rates or none of them")
    else:
        scan_seconds = max(disk_seconds, decode_seconds, h2d_seconds, compute_seconds) * contention_factor
        model = "isolated_stage_roofline_with_contention_factor"

    total_seconds = setup_seconds + scan_seconds + postprocess_seconds
    return {
        "model": model,
        "n_chunks": n_chunks,
        "decoded_gb": decoded_bytes / 1e9,
        "compressed_gb": compressed_bytes / 1e9,
        "disk_seconds_isolated": disk_seconds,
        "decode_seconds_isolated": decode_seconds,
        "h2d_seconds_isolated": h2d_seconds,
        "compute_seconds_isolated": compute_seconds,
        "scan_seconds": scan_seconds,
        "setup_seconds": setup_seconds,
        "postprocess_seconds": postprocess_seconds,
        "total_seconds": total_seconds,
    }


def predict_runtime_from_hardware(
    *,
    n_variants: int,
    n_samples: int,
    n_traits: int,
    covariate_rank: int,
    chunk_size: int = 2_500,
    gpu_fp32_tflops: float,
    gpu_count: int,
    gemm_efficiency: float,
    cpu_frequency_ghz: float,
    decode_threads: int,
    decode_gbps_per_core_at_reference: float,
    reference_cpu_frequency_ghz: float,
    decode_thread_efficiency: float,
    disk_gbps: float,
    h2d_gbps: float,
    compression_ratio: float = 12.46,
    non_gemm_ms_per_chunk_per_gpu: float = 0.818,
    pipeline_contention_factor: float = 1.0,
    setup_seconds: float = 0.0,
    postprocess_seconds: float = 0.0,
) -> dict[str, float | int | str]:
    """Estimate runtime from portable hardware measures plus calibrated efficiencies.

    The structure comes from the torchGWAS1.1 stage model. Peak FLOP/s alone
    is insufficient: zstd decode scales with CPU frequency and producer count,
    while compressed reads, H2D, and GEMM overlap and share the memory fabric.
    ``pipeline_contention_factor`` should therefore come from one short joint
    calibration run when NUMA and PCIe topology differ from the reference host.
    """

    for name, value in {
        "gpu_fp32_tflops": gpu_fp32_tflops,
        "gemm_efficiency": gemm_efficiency,
        "cpu_frequency_ghz": cpu_frequency_ghz,
        "decode_gbps_per_core_at_reference": decode_gbps_per_core_at_reference,
        "reference_cpu_frequency_ghz": reference_cpu_frequency_ghz,
        "decode_thread_efficiency": decode_thread_efficiency,
    }.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    if gpu_count <= 0 or decode_threads <= 0:
        raise ValueError("gpu_count and decode_threads must be positive")
    if gemm_efficiency > 1.0 or decode_thread_efficiency > 1.0:
        raise ValueError("efficiency values must be in (0, 1]")

    decode_gbps = (
        decode_gbps_per_core_at_reference
        * decode_threads
        * (cpu_frequency_ghz / reference_cpu_frequency_ghz)
        * decode_thread_efficiency
    )
    sustained_gemm_tflops = gpu_fp32_tflops * gpu_count * gemm_efficiency
    result = predict_runtime(
        n_variants=n_variants,
        n_samples=n_samples,
        n_traits=n_traits,
        covariate_rank=covariate_rank,
        chunk_size=chunk_size,
        decoded_bytes_per_value=1.0,
        compression_ratio=compression_ratio,
        disk_gbps=disk_gbps,
        decode_gbps=decode_gbps,
        h2d_gbps=h2d_gbps,
        measured_gemm_tflops=sustained_gemm_tflops,
        non_gemm_ms_per_chunk=non_gemm_ms_per_chunk_per_gpu / gpu_count,
        contention_factor=pipeline_contention_factor,
        setup_seconds=setup_seconds,
        postprocess_seconds=postprocess_seconds,
    )
    result.update(
        {
            "model": "hardware_calibrated_stage_roofline",
            "estimated_decode_gbps": decode_gbps,
            "estimated_sustained_gemm_tflops": sustained_gemm_tflops,
            "gpu_count": int(gpu_count),
            "decode_threads": int(decode_threads),
            "cpu_frequency_ghz": float(cpu_frequency_ghz),
            "pipeline_contention_factor": float(pipeline_contention_factor),
        }
    )
    return result
