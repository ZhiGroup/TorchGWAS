"""Compatibility checks for independent component calibration.

Identity is an assertion supplied by the caller, not automatic hardware detection.
Rates are portable only after recalibration in the requested execution context.
"""

CONTEXT_KEYS = ("machine", "processor", "accelerator", "software", "placement")


def check_calibration_context(measured, target, *, strict=False):
    """Reject mismatches; legacy unidentified inputs remain explicitly unresolved.

    Values are nonempty strings. software identifies tool/compiler/BLAS/CUDA
    builds; placement identifies affinity, NUMA policy and concurrency. Use
    'none' for accelerator on CPU-only paths, never an omitted field.
    """
    issues = []
    for label, context in (("measured", measured), ("target", target)):
        if context is None:
            issues.append(label + " calibration context missing")
            continue
        if not isinstance(context, dict):
            raise ValueError(label + " calibration context must be a dictionary")
        for key in CONTEXT_KEYS:
            if not isinstance(context.get(key), str) or not context[key].strip():
                raise ValueError(label + " calibration context requires " + key)
    if measured is not None and target is not None:
        mismatch = [key for key in CONTEXT_KEYS if measured[key] != target[key]]
        if mismatch:
            raise ValueError("recalibrate components for changed " + ", ".join(mismatch))
    if strict and issues:
        raise ValueError("; ".join(issues))
    return issues
