from __future__ import annotations

import numpy as np

from .utils import check_aligned_rows, chunk_bounds, column_std_mask, ensure_2d, validate_no_missing


def _phenotype_column_mask(phenotype: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nonconstant phenotype columns and their observed sample counts.

    Phenotype NaNs follow the genotype missing-call contract: they are allowed,
    excluded from the column mean, and later represented by a zero centred
    contribution (mean imputation).  Infinities are malformed input, not
    missing values, and are refused.
    """
    if np.isinf(phenotype).any():
        raise ValueError("phenotype contains infinite values")
    observed = ~np.isnan(phenotype)
    counts = observed.sum(axis=0).astype(np.int64)
    safe_counts = np.maximum(counts, 1)
    means = np.where(observed, phenotype, 0).sum(axis=0) / safe_counts
    centered = phenotype - means[None, :]
    np.copyto(centered, 0, where=~observed)
    np.square(centered, out=centered)
    sum_squares = np.sum(centered, axis=0)
    keep = (counts > 1) & np.isfinite(sum_squares) & (sum_squares > 0)
    return keep, counts


def _phenotype_column_mask_device(block, device):
    """_phenotype_column_mask plus the missing-cell count, on a CUDA device.

    Same decisions as the NumPy version: sums in float64, two-pass centred
    sum of squares, so a constant column is exactly zero either way.
    """
    import torch
    import warnings
    with warnings.catch_warnings():
        # A read-only memory map is fine: the block is only read, on its way to the GPU.
        warnings.filterwarnings('ignore', message='The given NumPy array is not writable')
        values = torch.as_tensor(np.ascontiguousarray(block), device=device)
    if torch.isinf(values).any():
        raise ValueError("phenotype contains infinite values")
    observed = ~torch.isnan(values)
    counts = observed.sum(dim=0)
    filled = torch.where(observed, values, torch.zeros((), dtype=values.dtype, device=values.device)).double()
    means = filled.sum(dim=0) / counts.clamp(min=1)
    filled.sub_(means).mul_(observed)
    sum_squares = (filled * filled).sum(dim=0)
    keep = (counts > 1) & torch.isfinite(sum_squares) & (sum_squares > 0)
    missing = int(observed.numel() - int(counts.sum()))
    return keep.cpu().numpy(), counts.cpu().numpy().astype(np.int64), missing


def _phenotype_qc(phenotype, *, dtype, device, block_bytes=512 << 20):
    """Keep mask, observed counts and missing cells over column blocks on `device`.

    Blocks are sized by bytes, not columns, so a million samples still fit.
    """
    n, k = phenotype.shape
    width = max(1, int(block_bytes // max(1, n * 4)))
    keep = np.empty(k, dtype=bool)
    counts = np.empty(k, dtype=np.int64)
    missing = 0
    for start in range(0, k, width):
        end = min(start + width, k)
        block = np.asarray(phenotype[:, start:end], dtype=dtype)
        keep[start:end], counts[start:end], cells = _phenotype_column_mask_device(block, device)
        missing += cells
    return keep, counts, missing


def prepare_inputs(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, dict]:
    genotype = ensure_2d(np.asarray(genotype, dtype=np.float64), "genotype")
    phenotype = ensure_2d(np.asarray(phenotype, dtype=np.float64), "phenotype")
    covariates = None if covariates is None else ensure_2d(np.asarray(covariates, dtype=np.float64), "covariates")

    validate_no_missing(genotype, "genotype")
    pheno_mask, phenotype_observed_counts = _phenotype_column_mask(phenotype)
    if covariates is not None:
        validate_no_missing(covariates, "covariates")
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype), ("covariates", covariates))
    else:
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype))

    geno_mask = column_std_mask(genotype)
    if covariates is not None:
        covar_mask = column_std_mask(covariates)
        covariates = covariates[:, covar_mask]
    else:
        covar_mask = np.array([], dtype=bool)

    qc = {
        "genotype_columns_input": int(genotype.shape[1]),
        "genotype_columns_kept": int(geno_mask.sum()),
        "phenotype_columns_input": int(phenotype.shape[1]),
        "phenotype_columns_kept": int(pheno_mask.sum()),
        "covariate_columns_input": int(0 if covariates is None else len(covar_mask)),
        "covariate_columns_kept": int(0 if covariates is None else covariates.shape[1]),
        "dropped_genotype_columns": int((~geno_mask).sum()),
        "dropped_phenotype_columns": int((~pheno_mask).sum()),
        "dropped_covariate_columns": int(0 if covariates is None else (~covar_mask).sum()),
        "n_samples": int(genotype.shape[0]),
        "phenotype_missing_cells": int(np.isnan(phenotype).sum()),
        "phenotype_observed_counts": phenotype_observed_counts[pheno_mask].tolist(),
    }
    if not geno_mask.all():
        # Retained input column indices, so callers keep marker and trait
        # labels aligned with the surviving columns rather than assuming a
        # prefix of the original order.
        qc["genotype_kept_column_indices"] = np.flatnonzero(geno_mask).tolist()
        genotype = genotype[:, geno_mask]
    if not pheno_mask.all():
        qc["phenotype_kept_column_indices"] = np.flatnonzero(pheno_mask).tolist()
        phenotype = phenotype[:, pheno_mask]
    if phenotype.shape[1] == 0:
        raise ValueError("all phenotype columns were dropped due to zero variance")
    if genotype.shape[1] == 0:
        raise ValueError("all genotype columns were dropped due to zero variance")
    return genotype, phenotype, covariates, qc


class PhenotypeColumnView:
    """Column selection without materializing a filtered whole phenotype panel."""
    def __init__(self,values,columns):
        self.values=values
        self.columns=np.asarray(columns,dtype=np.int64)
        self.shape=(values.shape[0],len(self.columns))
        self.ndim=2
        self.dtype=values.dtype

    def __getitem__(self,key):
        rows,columns=key
        selected=self.columns[columns]
        return self.values[rows][...,selected]

    def __array__(self,dtype=None,copy=None):
        if copy is False:
            raise ValueError('Filtered phenotype columns require a copy')
        return np.asarray(self.values[:,self.columns],dtype=dtype)


class PhenotypeRowView:
    """Row selection without materializing the kept rows of a whole phenotype panel.

    missing_phenotype='drop_subject' on a memory-mapped panel: each block a
    consumer asks for (`view[rows, columns]`) reads only the kept rows of
    those columns, so a voxel-scale panel is never copied whole.
    """
    def __init__(self,values,rows):
        self.values=values
        self.rows=np.asarray(rows,dtype=np.int64)
        self.shape=(len(self.rows),values.shape[1])
        self.ndim=2
        self.dtype=values.dtype

    def __getitem__(self,key):
        if not isinstance(key,tuple):
            # A row selection stays lazy (PhenotypeColumnView indexes rows first).
            return PhenotypeRowView(self.values,self.rows[key])
        rows,columns=key
        if rows is Ellipsis:
            rows=slice(None)
        # Columns first: a slice of a memmap is a view, so only the kept rows
        # of these columns are read.
        return np.asarray(self.values[:,columns])[self.rows[rows]]

    def __array__(self,dtype=None,copy=None):
        if copy is False:
            raise ValueError('Selected phenotype rows require a copy')
        return np.asarray(self.values[self.rows],dtype=dtype)


def prepare_inputs_for_prep(
    genotype,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
    genotype_chunk_size: int | None = None,
    validate_genotype: bool = True,
    dtype=np.float64,
    phenotype_block_size: int | None = None,
    qc_device=None,
) -> tuple[np.ndarray, np.ndarray | None, dict]:
    """Validate and QC the inputs the out-of-core scan will stream against.

    `dtype` is the working precision for the phenotype and covariates, and
    it defaults to float64 because every existing caller and test was
    written against that promotion. **A float32 CUDA scan should pass
    float32.** The streaming path ends at `torch.as_tensor(pheno_proc,
    dtype=torch.float32)`, so upcasting here builds a float64 copy that is
    downcast again and thrown away. At 33,417 subjects and 600,000 voxels
    that copy is 160 GB; at the full 2,085,000 voxels it is 557 GB, which
    is most of a 1 TB host spent on precision nothing downstream consumes.
    """
    if not hasattr(genotype, "shape") or len(genotype.shape) != 2:
        raise ValueError(f"genotype must be 2D, got {getattr(genotype, 'shape', None)}")
    # `asarray` is a no-op when the caller already holds this dtype, which
    # is the point: at voxel scale that copy is the whole problem.
    if phenotype_block_size is not None and (isinstance(phenotype_block_size,bool) or not isinstance(phenotype_block_size,int) or phenotype_block_size<1):
        raise ValueError('phenotype_block_size must be a positive integer')
    # A lazy view (drop_subject's kept rows of a memory-mapped panel) stays
    # lazy: `asarray` would run its __array__, which copies every kept row of
    # every trait -- 800 GB at 100,000 subjects and 2,000,000 traits. QC and
    # the scan read it a column block at a time, converting each block.
    lazy = isinstance(phenotype, (PhenotypeRowView, PhenotypeColumnView)) or (
        not isinstance(phenotype, np.ndarray) and hasattr(phenotype, '__getitem__')
        and hasattr(phenotype, 'dtype') and len(getattr(phenotype, 'shape', ())) == 2)
    if not lazy:
        phenotype = ensure_2d(np.asarray(phenotype, dtype=dtype if phenotype_block_size is None else None),
                              "phenotype")
    covariates = None if covariates is None else ensure_2d(np.asarray(covariates, dtype=dtype), "covariates")

    if qc_device is not None:
        # Column QC on the scan GPU (bounded blocks), not one CPU core.
        pheno_mask, phenotype_observed_counts, phenotype_missing_cells = _phenotype_qc(
            phenotype, dtype=dtype, device=qc_device)
    elif phenotype_block_size is None and not lazy:
        pheno_mask, phenotype_observed_counts = _phenotype_column_mask(phenotype)
        phenotype_missing_cells=int(np.isnan(phenotype).sum())
    else:
        # Bounded column blocks: the caller's width, else ~512 MiB of float64.
        block_size = phenotype_block_size or max(1, (512 << 20) // (8 * max(phenotype.shape[0], 1)))
        pheno_mask=np.empty(phenotype.shape[1],dtype=bool)
        phenotype_observed_counts=np.empty(phenotype.shape[1],dtype=np.int64)
        phenotype_missing_cells=0
        for start in range(0,phenotype.shape[1],block_size):
            end=min(start+block_size,phenotype.shape[1])
            block=np.asarray(phenotype[:,start:end],dtype=dtype)
            pheno_mask[start:end],phenotype_observed_counts[start:end]=_phenotype_column_mask(block)
            phenotype_missing_cells+=int(np.isnan(block).sum())
    if covariates is not None:
        validate_no_missing(covariates, "covariates")
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype), ("covariates", covariates))
    else:
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype))

    effective_chunk_size = (
        genotype_chunk_size
        or getattr(genotype, "preferred_chunk_size", None)
        or min(genotype.shape[1], 4096)
        or 1
    )
    if validate_genotype:
        geno_mask = _chunked_genotype_std_mask(genotype, chunk_size=effective_chunk_size)
        if not geno_mask.all():
            raise ValueError(
                f"out-of-core genotype contains {(~geno_mask).sum()} zero-variance variants; "
                "filter invariant variants before running TorchGWAS"
            )
        genotype_qc_mode = "chunked"
    else:
        # The native packed-BED CUDA iterator performs these checks while the
        # association scan is already resident on the GPU. This avoids a full
        # CPU decode and a second read of the input before the real scan.
        geno_mask = np.ones(genotype.shape[1], dtype=bool)
        genotype_qc_mode = "fused_gpu_scan"
    if covariates is not None:
        covar_mask = column_std_mask(covariates)
        covariates = covariates[:, covar_mask]
    else:
        covar_mask = np.array([], dtype=bool)

    qc = {
        "genotype_columns_input": int(genotype.shape[1]),
        "genotype_columns_kept": int(geno_mask.sum()),
        "phenotype_columns_input": int(phenotype.shape[1]),
        "phenotype_columns_kept": int(pheno_mask.sum()),
        "covariate_columns_input": int(0 if covariates is None else len(covar_mask)),
        "covariate_columns_kept": int(0 if covariates is None else covariates.shape[1]),
        "dropped_genotype_columns": int((~geno_mask).sum()),
        "dropped_phenotype_columns": int((~pheno_mask).sum()),
        "dropped_covariate_columns": int(0 if covariates is None else (~covar_mask).sum()),
        "n_samples": int(genotype.shape[0]),
        "genotype_qc_mode": genotype_qc_mode,
        "genotype_qc_chunk_size": int(effective_chunk_size),
        "phenotype_missing_cells": phenotype_missing_cells,
        "phenotype_observed_counts": phenotype_observed_counts[pheno_mask].tolist(),
    }

    if not pheno_mask.all():
        qc['phenotype_kept_column_indices']=np.flatnonzero(pheno_mask).tolist()
        phenotype = (phenotype[:, pheno_mask] if phenotype_block_size is None and not lazy else
                     PhenotypeColumnView(phenotype,np.flatnonzero(pheno_mask)))
    if phenotype.shape[1] == 0:
        raise ValueError("all phenotype columns were dropped due to zero variance")
    return phenotype, covariates, qc


def _chunked_genotype_std_mask(genotype: np.ndarray, chunk_size: int | None = None) -> np.ndarray:
    n_markers = genotype.shape[1]
    mask = np.ones(n_markers, dtype=bool)
    chunk = chunk_size or getattr(genotype, "preferred_chunk_size", None) or min(n_markers, 4096) or 1
    if hasattr(genotype, "iter_chunks"):
        iterator = genotype.iter_chunks(chunk_size=chunk, dtype=np.float64)
    else:
        iterator = (
            (start, end, np.asarray(genotype[:, start:end], dtype=np.float64))
            for start, end in chunk_bounds(n_markers, chunk)
        )
    for start, end, geno_chunk in iterator:
        if not np.isfinite(geno_chunk).all():
            raise ValueError("genotype contains missing/non-finite values; v0.1 requires complete matrices")
        mask[start:end] = np.nanstd(geno_chunk, axis=0) > 0
    return mask


def _covariate_basis(covariates: np.ndarray) -> np.ndarray | None:
    """Orthonormal basis for the covariate column space, or None if rank 0.

    Kept in numpy on the host: it is a 27-column SVD costing about 30 ms and it
    does not grow with the trait count, so there is nothing to gain by moving
    it, and the rank tolerance is the one piece of this whose exact behaviour
    the tests pin down.

    **The rank decision is made in float64 whatever the caller passes, and that
    is not a detail.** The tolerance is numpy's `matrix_rank` rule,
    `max(shape) * eps * s[0]`, so it inherits the input's `eps`: at N = 22,250
    that is a relative cut of 4.9e-12 in float64 but **2.6e-3 in float32** --
    nine orders of magnitude apart. Measured, a covariate correlated with
    another at r = 0.9999995 (distinct, and carrying real information) is kept
    under the float64 rule and *dropped* under the float32 one. Dropping it
    changes `covariate_rank`, which changes `df`, which changes every
    t-statistic and p-value the run reports -- silently. Several callers reach
    this function without passing through `prepare_inputs_for_prep`'s float64
    cast, so the promotion happens here rather than being assumed upstream. An
    (N x 27) float64 copy is a few megabytes, and the SVD is more accurate for
    it.

    Column *scale*, by contrast, cannot trigger the truncation: each column is
    divided by its own standard deviation before the SVD, so a covariate that is
    1e-6 of the others is restored to unit variance first. Only genuine
    near-collinearity produces a small singular value here.
    """
    original_dtype = np.asarray(covariates).dtype
    covariates = np.asarray(covariates, dtype=np.float64)
    cov_centered = covariates - covariates.mean(axis=0, keepdims=True)
    cov_std = covariates.std(axis=0, keepdims=True)
    cov_std[cov_std == 0] = 1.0
    cov_scaled = cov_centered / cov_std
    u, singular_values, _ = np.linalg.svd(cov_scaled, full_matrices=False)
    if not singular_values.size:
        return None
    tolerance = max(cov_scaled.shape) * np.finfo(cov_scaled.dtype).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    # The projection uses U, an orthonormal basis, so nothing here divides by a
    # small singular value: a near-collinear pair that survives the cut is still
    # numerically harmless. A rewrite through `pinv` or the V-side would not be.
    #
    # The basis goes back out in the dtype it came in. Only the *rank decision*
    # was wrong in float32; returning float64 unconditionally would also widen
    # the residualised phenotype downstream -- `np.result_type` promotes against
    # this matrix -- and that is a behaviour change, not a bug fix, so it is not
    # smuggled in here.
    if not rank:
        return None
    return np.asarray(u[:, :rank], dtype=original_dtype)


# Residualisation is a one-off pass of about 2 GFLOP through a 27-column
# basis -- milliseconds on an idle card -- so its block should be sized to
# stay out of the way, not to go as fast as possible. Sizing it from free
# memory instead made it the PEAK of the whole run: at 33,417 samples and
# 150,000 voxels it took half the card (28.2 GB) and, because it is still
# resident when the design is built (20.5 GB), the process peaked at 55.2 GB
# where the scan itself needs about 21 GB. A fixed working-set budget costs
# a few more blocks of a negligible GEMM and gives the memory back.
_RESIDUALIZE_WORKING_SET_BYTES = 4 << 30


def _device_trait_block(n_samples: int, itemsize: int, device,
                        budget_bytes: int = _RESIDUALIZE_WORKING_SET_BYTES,
                        ) -> int:
    """How many traits to residualise at once on the device.

    Three full-width temporaries are live at the peak -- the block itself,
    the projection `basis @ (basis.T @ values)`, and the allocator's slack
    while the subtraction lands -- so the per-trait cost is about three
    columns. `budget_bytes` bounds their total; it does not scale with the
    card, deliberately, because a bigger block buys no measurable time here
    and costs memory the scan needs.
    """
    per_trait = 3.0 * max(n_samples, 1) * max(itemsize, 1)
    return max(int(budget_bytes // per_trait), 1)

def _residualize_on_device(phenotype, q_matrix, device, keep_on_device=False):
    """The projection and standardisation, on the GPU, in trait blocks.

    Same operations in the same order as the numpy path below, including
    `correction=0` on the standard deviation -- numpy's `std` defaults to the
    population form and torch's to the sample form, and silently swapping those
    would rescale every phenotype and so every t-statistic this tool reports.

    Worth moving because the host version is both slow and *unstable*: the same
    call on the same machine has been measured at 0.37 s and at 10.26 s, since
    importing torch degrades numpy's threaded GEMM by about five times and the
    interaction depends on how the two thread pools happen to land. The work
    itself is roughly 2 GFLOP through a 27-column basis, which is milliseconds
    on a device that is otherwise idle at this point in the run.

    **Blocked over traits, because uploading the whole phenotype does not
    scale and failed silently when it stopped fitting.** This function used to
    send the entire matrix to the card in one `as_tensor`. At 33,417 subjects
    and 600,000 voxels that is 160 GB against an 80 GB card, so it raised, the
    caller caught the `RuntimeError` as "a device that is busy", and the run
    fell back to the numpy path -- which then ground three more full-size
    float64 temporaries through one core. Observed: 308 GB resident, one core
    at 100%, all eight GPUs at 0%, and no progress in over an hour.

    Blocking is exact rather than an approximation: centring, projection
    through `q_matrix`, and scaling are all per-column operations, so a block
    of columns computes the same values it would have computed in company.
    """
    import torch

    # Match numpy's promotion rather than the caller's dtype. The host path
    # subtracts a float64 projection from a float32 phenotype, so it returns
    # float64 whenever the covariates are float64 -- quietly, but everything
    # downstream has been built and tested against that, and narrowing it here
    # would be a behaviour change smuggled in as an optimisation.
    out_dtype = (phenotype.dtype if q_matrix is None
                 else np.result_type(phenotype.dtype, q_matrix.dtype))
    torch_dtype = torch.float64 if out_dtype == np.float64 else torch.float32
    n_samples, n_traits = phenotype.shape
    basis = (None if q_matrix is None else
             torch.as_tensor(np.ascontiguousarray(q_matrix), device=device,
                             dtype=torch_dtype))
    block = _device_trait_block(n_samples, np.dtype(out_dtype).itemsize, device)
    # One output array, written block by block. Accumulating blocks in a list
    # and concatenating would hold two full copies at the join, which is the
    # allocation this change exists to avoid.
    single_block = 0 < n_traits <= block
    # A completed CPU download already owns independent, writable storage.
    # Its NumPy view retains that CPU tensor; a second whole-tile allocation
    # and copy adds first-touch/page-release work without changing any value.
    if keep_on_device:
        # Callers that fan the panel out to other GPUs keep it here.
        result = torch.empty((n_samples, n_traits), dtype=torch_dtype, device=device)
        single_block = False
    else:
        result = None if single_block else np.empty((n_samples, n_traits), dtype=out_dtype)
    for begin in range(0, n_traits, block):
        stop = min(begin + block, n_traits)
        values = torch.as_tensor(
            np.ascontiguousarray(phenotype[:, begin:stop]), device=device,
            dtype=torch_dtype)
        values -= values.mean(dim=0, keepdim=True)
        if basis is not None:
            values -= basis @ (basis.T @ values)
        deviation = values.std(dim=0, keepdim=True, correction=0)
        deviation = torch.where(deviation == 0, torch.ones_like(deviation),
                                deviation)
        values /= deviation
        if keep_on_device:
            result[:, begin:stop].copy_(values)
            del values
            continue
        downloaded = values.cpu().numpy()
        if single_block:
            return downloaded
        result[:, begin:stop] = downloaded
        del downloaded
        del values
    return result


# Traits per block in the host residualisation. Chosen so one block's
# temporaries are a few hundred MB at realistic sample counts rather than a
# fraction of the whole matrix: at 35,365 samples a 4,096-trait block is
# 579 MB. Each trait is residualised independently, so the blocking is
# arithmetic-neutral -- it changes the peak memory and nothing else.
RESIDUALIZE_TRAIT_BLOCK = 4096


def residualize_and_standardize(
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    device=None,
    inplace: bool = False,
    trait_block: int = RESIDUALIZE_TRAIT_BLOCK,
    out_dtype=None,
    return_observed_counts: bool = False,
    _prevalidated_observed_counts=None,
    _prevalidated_covariate_basis=...,
    keep_on_device: bool = False,
    column_groups=None,
    missing: str = "impute",
    return_plan: bool = False,
):
    """Centre each trait, project out the covariates, and scale to unit variance.

    Missing values: `missing='impute'` centres each trait on its observed
    mean, fills the missing cells with it and projects the covariates out
    over all samples (JAGWAS's panel). `missing='complete_case'` instead
    residualizes each trait with missing values on [1, covariates] over its
    own observed rows, zero-fills the rest and scales it to unit variance
    there, for per-trait complete-case OLS (complete_case.CompleteCasePlan).
    return_plan=True appends the scan's plan for the missing traits: that
    CompleteCasePlan, or for 'impute' the release's rescaling
    (complete_case.ImputedPlan); None for a complete panel.

    `column_groups` (disjoint index arrays) residualises each group by its own
    call, exactly as a run of that group alone would; columns in no group are
    one more call. The arithmetic per column is the same either way, but
    reductions and GEMMs are chosen by block width and round differently, and
    a JAGWAS rank decision on a near-collinear group can turn on that rounding
    (jagwas_projection.JagwasGroups). The result is always NumPy.

    keep_on_device=True (CUDA device, complete phenotypes) returns the result
    as a tensor on `device` instead of downloading it; other paths still
    return NumPy.

    `device` is optional and advisory: pass a CUDA device to run the projection
    there. The numpy path stays the reference implementation and the two are
    compared in the tests, so a device that is missing, busy, or out of memory
    simply falls back rather than failing the scan.

    MEMORY, and this was a real defect rather than a tidiness point. The host
    path used to be three whole-matrix expressions:

        phenotype = phenotype - phenotype.mean(axis=0, keepdims=True)
        phenotype = phenotype - q_matrix @ (q_matrix.T @ phenotype)
        phenotype = phenotype / std

    Each allocates a full-size array, and the middle line holds THREE at once:
    the input, the `q @ (q.T @ P)` product, and the result. At the voxel stress
    configuration -- 35,365 samples x 600,000 traits, an 84.9 GB phenotype --
    that is ~255 GB of peak anonymous memory. The run was measured at 229 GB
    resident with 53,468 s of KERNEL time against 348 s of user time, a 153x
    ratio that is page-fault servicing rather than computation, while all eight
    GPUs sat idle 69% of the time waiting behind it.

    The work is now done in column blocks, so the transient is one block rather
    than one copy of the matrix. `inplace=True` additionally writes through the
    caller's own array, dropping the peak to the phenotype itself; it mutates
    the argument, so it is opt-in and the default remains the safe copy.
    """
    if missing not in ("impute", "complete_case"):
        raise ValueError("missing must be 'impute' or 'complete_case'")
    if column_groups is not None:
        if inplace:
            raise ValueError("column_groups cannot be combined with inplace=True")
        groups = [np.asarray(columns, dtype=np.int64).reshape(-1) for columns in column_groups]
        uncovered = np.setdiff1d(np.arange(phenotype.shape[1]),
                                 np.concatenate(groups) if groups else np.empty(0, np.int64))
        prevalidated = (None if _prevalidated_observed_counts is None
                        else np.asarray(_prevalidated_observed_counts))
        out = group_q = None
        counts = np.empty(phenotype.shape[1], dtype=np.int64)
        for columns in groups + ([uncovered] if len(uncovered) else []):
            block, group_q, block_counts = residualize_and_standardize(
                phenotype[:, columns], covariates, device=device, trait_block=trait_block,
                out_dtype=out_dtype, return_observed_counts=True,
                _prevalidated_observed_counts=None if prevalidated is None else prevalidated[columns],
                _prevalidated_covariate_basis=_prevalidated_covariate_basis)
            block = np.asarray(block)
            if out is None:
                out = np.empty(phenotype.shape, dtype=block.dtype)
            out[:, columns] = block
            counts[columns] = block_counts
        plan = None
        if missing == "complete_case" and np.any(counts != phenotype.shape[0]):
            # Per-group subject drops (JAGWAS groups under drop_subject): one
            # plan over the panel, each trait with missing values residualized
            # on its own rows -- per trait, so it is what a run of its group
            # alone would give.
            from .complete_case import CompleteCasePlan
            observed = np.ones(phenotype.shape, dtype=bool)
            for start in range(0, phenotype.shape[1], max(1, int(trait_block))):
                stop = min(start + max(1, int(trait_block)), phenotype.shape[1])
                observed[:, start:stop] = ~np.isnan(np.asarray(phenotype[:, start:stop], dtype=np.float64))
            plan = CompleteCasePlan(observed, group_q)
            for trait in plan.traits:
                values = np.array(phenotype[:, int(trait)], dtype=np.float64)
                out[:, int(trait)] = plan.residualize(values, int(trait)).astype(out.dtype, copy=False)
        elif return_plan and missing == "impute" and np.any(counts != phenotype.shape[0]):
            from .complete_case import ImputedPlan
            plan = ImputedPlan(counts, 0 if group_q is None else group_q.shape[1], phenotype.shape[0])
        result = (out, group_q)
        if return_observed_counts:
            result = (*result, counts)
        return (*result, plan) if return_plan else result

    q_matrix = None
    if _prevalidated_covariate_basis is not ...:
        # Internal reuse only: the API derived this basis from the same retained
        # covariates to determine rank/df before starting its tile/shard workers.
        # None deliberately means a validated zero-rank basis, not "recompute".
        q_matrix = _prevalidated_covariate_basis
        if q_matrix is not None:
            q_matrix=np.asarray(q_matrix)
            columns=0 if covariates is None else covariates.shape[1]
            if (q_matrix.ndim!=2 or q_matrix.shape[0]!=phenotype.shape[0]
                or not 0<q_matrix.shape[1]<=columns):
                raise ValueError('Prevalidated covariate basis must match the retained sample/column dimensions')
    elif covariates is not None and covariates.shape[1] > 0:
        q_matrix = _covariate_basis(covariates)

    if _prevalidated_observed_counts is None:
        phenotype_observed_counts = np.sum(~np.isnan(phenotype), axis=0).astype(np.int64)
    else:
        # Internal reuse of QC for this unchanged input, including its retained
        # column order. Callers without that provenance must use the default.
        counts = np.asarray(_prevalidated_observed_counts)
        if (phenotype.ndim != 2 or counts.shape != (phenotype.shape[1],)
                or counts.dtype.kind not in 'iu'
                or np.any(counts < 0) or np.any(counts > phenotype.shape[0])):
            raise ValueError('Prevalidated observed counts must match the phenotype columns and sample range')
        phenotype_observed_counts = counts.astype(np.int64, copy=False)
    if np.any(phenotype_observed_counts == 0):
        raise ValueError("phenotype contains an entirely missing trait")
    has_missing = bool(np.any(phenotype_observed_counts != phenotype.shape[0]))

    if (not has_missing and device is not None
            and getattr(device, "type", None) == "cuda"):
        try:
            processed = (_residualize_on_device(phenotype, q_matrix, device, keep_on_device=True) if keep_on_device
                         else _residualize_on_device(phenotype, q_matrix, device))
            result = (processed, q_matrix)
            if return_observed_counts:
                result = (*result, phenotype_observed_counts)
            return (*result, None) if return_plan else result
        except (RuntimeError, ImportError):
            # Out of memory, a driver problem, or no torch: the host path below
            # produces the same answer, so this is a slowdown and not a failure.
            pass

    if phenotype.ndim != 2:
        raise ValueError("phenotype must be 2-D (samples x traits)")
    # Default to numpy's promotion, which is what both paths have always
    # returned and what everything downstream was built against: subtracting a
    # float64 projection from a float32 phenotype gives float64. Narrowing it
    # by default would be a behaviour change smuggled in as an optimisation --
    # the warning `_residualize_on_device` already carries.
    #
    # It is worth being able to ask for otherwise, though, and the reason is
    # size rather than taste: at 35,365 samples and 600,000 traits the
    # phenotype is 84.9 GB as float32 and 170 GB promoted, and the scan casts
    # it back to float32 anyway. `out_dtype` makes that choice explicit at the
    # call site instead of silent here.
    resolved = (phenotype.dtype if q_matrix is None
                else np.result_type(phenotype.dtype, q_matrix.dtype))
    if out_dtype is not None:
        resolved = np.dtype(out_dtype)
    if inplace and resolved != phenotype.dtype:
        raise ValueError(
            f"inplace=True needs the result dtype to match the phenotype, but "
            f"{phenotype.dtype} would become {resolved}; pass "
            f"out_dtype={phenotype.dtype} to keep it, or inplace=False")

    plan = None
    if missing == "complete_case" and has_missing:
        from .complete_case import CompleteCasePlan
        # Before the loop: with inplace=True it overwrites the raw values.
        observed = np.ones(phenotype.shape, dtype=bool)
        for start in range(0, phenotype.shape[1], max(1, int(trait_block))):
            stop = min(start + max(1, int(trait_block)), phenotype.shape[1])
            observed[:, start:stop] = ~np.isnan(np.asarray(phenotype[:, start:stop], dtype=np.float64))
        plan = CompleteCasePlan(observed, q_matrix)
        raw_missing = {int(trait): np.array(phenotype[:, trait], dtype=np.float64) for trait in plan.traits}
    out = phenotype if inplace else np.empty(phenotype.shape, dtype=resolved)
    step = max(1, int(trait_block))
    for start in range(0, phenotype.shape[1], step):
        stop = min(start + step, phenotype.shape[1])
        # One working copy of the BLOCK, in the working precision; every step
        # after this stays inside it, so the transient is a block and not a
        # copy of the matrix.
        work = np.array(phenotype[:, start:stop], dtype=resolved, copy=True)
        observed = ~np.isnan(work)
        counts = observed.sum(axis=0, keepdims=True)
        means = np.where(observed, work, 0).sum(axis=0, keepdims=True) / counts
        # Same missing-value convention as genotype: after centring, an
        # unobserved cell contributes zero, i.e. the observed column mean.
        work = np.where(observed, work - means, 0)
        if q_matrix is not None:
            # `q.T @ work` is (rank x block) and tiny; the product back up is a
            # block-sized temporary, not a matrix-sized one.
            work -= (q_matrix @ (q_matrix.T @ work)).astype(resolved,
                                                            copy=False)
        std = work.std(axis=0, keepdims=True)
        std[std == 0] = 1.0
        work /= std
        out[:, start:stop] = work
    if plan is not None:
        for trait, values in raw_missing.items():
            out[:, trait] = plan.residualize(values, trait).astype(out.dtype, copy=False)
    if return_plan and plan is None and missing == "impute" and has_missing:
        from .complete_case import ImputedPlan
        plan = ImputedPlan(phenotype_observed_counts, 0 if q_matrix is None else q_matrix.shape[1],
                           phenotype.shape[0])
    result = (out, q_matrix)
    if return_observed_counts:
        result = (*result, phenotype_observed_counts)
    return (*result, plan) if return_plan else result


def standardize_genotype(genotype_chunk: np.ndarray) -> np.ndarray:
    centered = genotype_chunk - genotype_chunk.mean(axis=0, keepdims=True)
    std = centered.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    return centered / std


def mask_phenotype_outliers(phenotype, covariates, threshold, *, whole_rows):
    """Set phenotype values beyond `threshold` SD of their covariate-residualised,
    standardised trait to missing; return the masked FP64 copy and the affected rows.
    A trait with missing values is residualised and standardised on its
    observed values only.

    whole_rows (JAGWAS): a sample with any such value loses its whole panel row.
    On near-collinear imaging panels a few samples, extreme in many traits at
    once, made the low-variance directions heavy-tailed; masking or clipping
    only their extreme values broke the traits' linear relations for those
    samples and made those directions heavier-tailed still, while dropping the
    row made them Gaussian. A per-trait scan masks just the value. The missing
    values then take the ordinary phenotype-missingness path.
    """
    if not (threshold > 0 and np.isfinite(threshold)):
        raise ValueError("phenotype_outlier_sd must be a positive number")
    values = np.array(phenotype, dtype=np.float64)
    # A single trait may arrive 1-D; mask it as one column and return it 1-D.
    one_trait = values.ndim == 1
    if one_trait:
        values = values[:, None]
    # Each trait residualised and scaled on its own observed values. Mean-filled
    # missing cells would shrink the SD by sqrt(observed / samples): at half
    # missing a 4 SD cutoff masked about 60 times the values it should.
    standardized, _ = residualize_and_standardize(values, covariates, missing="complete_case")
    extreme = np.abs(standardized) > threshold
    rows = extreme.any(axis=1)
    if whole_rows:
        values[rows] = np.nan
    else:
        values[extreme] = np.nan
    return (values[:, 0] if one_trait else values), rows

