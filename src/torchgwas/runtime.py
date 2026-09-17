"""Compatibility names for the canonical :mod:`torchgwas.pipeline_model`.

New callers should supply full Hardware/InputProfile records to that module.
Legacy signatures omit memory bandwidth and GPU decode placement; their results
are partial lower bounds, not validated runtime predictions. Fitted joint-rate
and global contention overrides are retired rather than silently applied.
"""
from __future__ import annotations
import math
import warnings
from .pipeline_model import Workload, InputProfile, Hardware, PipelinePlan, estimate


def predict_runtime(*, n_variants, n_samples, n_traits, covariate_rank,
                    chunk_size, decoded_bytes_per_value, compression_ratio,
                    disk_gbps, decode_gbps, h2d_gbps, measured_gemm_tflops,
                    h2d_bytes_per_value=None, non_gemm_ms_per_chunk=0.0,
                    contention_factor=1.0, setup_seconds=0.0, postprocess_seconds=0.0,
                    joint_decode_chunks_per_second=None,
                    joint_h2d_chunks_per_second=None,
                    joint_compute_chunks_per_second=None):
    if any(value is not None for value in (joint_decode_chunks_per_second,
            joint_h2d_chunks_per_second, joint_compute_chunks_per_second)):
        raise ValueError('joint-rate calibration overrides are retired; supply independent component costs')
    if contention_factor != 1.0:
        raise ValueError('contention_factor is retired; shared resources are modeled explicitly')
    for name, value in dict(compression_ratio=compression_ratio, decode_gbps=decode_gbps,
                            decoded_bytes_per_value=decoded_bytes_per_value).items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be positive and finite')
    for value in (setup_seconds, postprocess_seconds, non_gemm_ms_per_chunk):
        if not math.isfinite(value) or value < 0:
            raise ValueError('timing overheads must be finite and nonnegative')
    warnings.warn('predict_runtime is a partial compatibility lower bound; use pipeline_model with full hardware and decode profiles',
                  DeprecationWarning, stacklevel=2)
    transfer_size = decoded_bytes_per_value if h2d_bytes_per_value is None else h2d_bytes_per_value
    decoded_bytes = n_variants * n_samples * decoded_bytes_per_value
    # Rates absent from the old interface are explicitly nonbinding; no claim
    # of a host/HBM or decoder-placement model can be made from these inputs.
    unmodeled_capacity = 1e300
    result = estimate(
        Workload(n_variants, n_samples, n_traits, covariate_rank),
        InputProfile('legacy-cpu-decode', decoded_bytes / compression_ratio,
                     n_samples * transfer_size, decoded_bytes_per_value,
                     cpu_decode_core_seconds_per_variant=n_samples * decoded_bytes_per_value / (decode_gbps * 1e9)),
        Hardware(disk_gbps * 1e9, h2d_gbps * 1e9, unmodeled_capacity,
                 measured_gemm_tflops * 1e12, unmodeled_capacity,
                 2**62, 2**62, 1, d2h_bytes_per_second=h2d_gbps * 1e9, compute_launch_seconds=non_gemm_ms_per_chunk / 1000),
        PipelinePlan(chunk_size, chunk_size, chunk_size, 1, 2))
    resource = result['resource_seconds']
    scan = result['resource_lower_bound_seconds']
    result.update(model='canonical_partial_resource_lower_bound', n_chunks=math.ceil(n_variants / chunk_size),
                  decoded_gb=decoded_bytes / 1e9,
                  h2d_gb=n_variants * n_samples * transfer_size / 1e9,
                  compressed_gb=decoded_bytes / compression_ratio / 1e9,
                  disk_seconds_isolated=result['io_lower_bound_seconds'],
                  decode_seconds_isolated=resource['cpu_decode'],
                  h2d_seconds_isolated=resource['h2d'],
                  compute_seconds_isolated=result['gpu_compute_seconds'],
                  scan_seconds=scan, setup_seconds=setup_seconds,
                  postprocess_seconds=postprocess_seconds,
                  total_seconds=setup_seconds + scan + postprocess_seconds)
    result['assumptions'].append('Legacy signature omits host/HBM capacities and assumes CPU decode; bound is incomplete.')
    return result


def predict_runtime_from_hardware(*, n_variants, n_samples, n_traits, covariate_rank,
        gpu_fp32_tflops, gpu_count, gemm_efficiency, cpu_frequency_ghz,
        decode_threads, decode_gbps_per_core_at_reference, reference_cpu_frequency_ghz,
        decode_thread_efficiency, disk_gbps, h2d_gbps, chunk_size=2500,
        h2d_bytes_per_value=1.0, compression_ratio=12.46,
        non_gemm_ms_per_chunk_per_gpu=0.818, pipeline_contention_factor=1.0,
        setup_seconds=0.0, postprocess_seconds=0.0):
    if pipeline_contention_factor != 1.0:
        raise ValueError('pipeline_contention_factor is retired; shared resources are modeled explicitly')
    for value in (gpu_count, decode_threads, cpu_frequency_ghz, reference_cpu_frequency_ghz,
                  decode_gbps_per_core_at_reference, gpu_fp32_tflops):
        if not math.isfinite(value) or value <= 0:
            raise ValueError('hardware rates and counts must be finite and positive')
    if not 0 < gemm_efficiency <= 1 or not 0 < decode_thread_efficiency <= 1:
        raise ValueError('efficiencies must be in (0, 1]')
    decode_rate = decode_gbps_per_core_at_reference * decode_threads * (cpu_frequency_ghz / reference_cpu_frequency_ghz) * decode_thread_efficiency
    gemm_rate = gpu_fp32_tflops * gpu_count * gemm_efficiency
    result = predict_runtime(n_variants=n_variants, n_samples=n_samples,
        n_traits=n_traits, covariate_rank=covariate_rank, chunk_size=chunk_size,
        decoded_bytes_per_value=1.0, h2d_bytes_per_value=h2d_bytes_per_value,
        compression_ratio=compression_ratio, disk_gbps=disk_gbps,
        decode_gbps=decode_rate, h2d_gbps=h2d_gbps, measured_gemm_tflops=gemm_rate,
        non_gemm_ms_per_chunk=non_gemm_ms_per_chunk_per_gpu / gpu_count,
        setup_seconds=setup_seconds, postprocess_seconds=postprocess_seconds)
    result.update(estimated_decode_gbps=decode_rate, estimated_sustained_gemm_tflops=gemm_rate,
                  gpu_count=gpu_count, decode_threads=decode_threads)
    result['assumptions'].append('Legacy frequency/efficiency scaling is caller-supplied and unverified; prefer independently sustained rates.')
    return result
