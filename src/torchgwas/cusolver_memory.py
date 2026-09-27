"""Consume exact, untimed cuSOLVER size queries without interpolation.

Two JAGWAS factors: the rounding cutoff's Cholesky (xpotrf) and the default
eigen cutoff's eigh and QR (xsyevd, then xgeqrf on the kept rows).
"""
import copy


def _rounded(size):
    return ((size+511)//512)*512


def _require_context(census, traits, device_profile):
    """The census must come from the same installed runtime, GPU architecture and library file."""
    if isinstance(traits, bool) or not isinstance(traits, int) or traits < 1:
        raise ValueError('Positive integer trait count required')
    matches = {'torch_version': 'torch_version', 'cuda_version': 'cuda_version',
        'preferred_linalg': 'preferred_linalg', 'compute_capability': 'compute_capability',
        'gpu': 'gpu_name', 'library': 'cusolver_library', 'host': 'host'}
    for artifact_key, profile_key in matches.items():
        if artifact_key not in census or profile_key not in device_profile or census[artifact_key] != device_profile[profile_key]:
            raise ValueError('Workspace query context mismatch: ' + profile_key)
    if census['torch_version'] != '2.5.1+cu124' or census['cuda_version'] != '12.4':
        raise ValueError('Workspace source contract requires PyTorch 2.5.1 CUDA 12.4')
    if census['preferred_linalg'] not in ('_LinalgBackend.Default', '_LinalgBackend.Cusolver'):
        raise ValueError('Default or explicit cuSOLVER factorization required')
    if census.get('durations_recorded') is not False:
        raise ValueError('Untimed workspace census required')
    library = census['library']
    if set(library) != {'path', 'version', 'bytes', 'mtime_ns'} or not isinstance(library['path'], str) or not library['path']:
        raise ValueError('Installed cuSOLVER file identity required')
    if not isinstance(library['version'], list) or len(library['version']) != 3 or any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in library['version']) or library['version'][0] < 11:
        raise ValueError('cuSOLVER 64-bit API version required')
    if any(isinstance(library[key], bool) or not isinstance(library[key], int) or library[key] <= 0 for key in ('bytes', 'mtime_ns')):
        raise ValueError('Installed cuSOLVER file metadata required')
    return library


def jagwas_factor_workspace(census, traits, device_profile):
    """Return requests only for a queried FP64 lower/default factor geometry.

    The profile must identify the same installed runtime, GPU architecture and
    library file. These requests do not include allocator reservation or an
    upper bound for all factorization allocations.
    """
    library = _require_context(census, traits, device_profile)
    if census.get('method', 'rounding') != 'rounding':
        raise ValueError('Cholesky (xpotrf) workspace census required')
    rows = {}
    for row in census['rows']:
        k = row['traits']
        if isinstance(k, bool) or not isinstance(k, int) or k < 1 or k in rows:
            raise ValueError('Unique positive queried trait counts required')
        for key in ('device_workspace_bytes', 'host_workspace_bytes'):
            if isinstance(row[key], bool) or not isinstance(row[key], int) or row[key] < 0:
                raise ValueError('Nonnegative integer workspace requests required')
        if row.get('matrix_allocated') is not False:
            raise ValueError('Size-only workspace query required')
        rows[k] = row
    if traits not in rows:
        raise ValueError('Exact trait count was not queried; interpolation is forbidden')
    row = rows[traits]
    device = row['device_workspace_bytes']
    return dict(traits=traits, device_requested_bytes=device,
        device_rounded_bytes=_rounded(device),
        host_requested_bytes=row['host_workspace_bytes'],
        library=copy.deepcopy(library),
        scope='Exact installed FP64 lower/default xpotrf workspace query. File metadata identifies the installation, not a binary hash. Device rounding uses the existing 512-byte tensor-accounting convention. Excludes info/error tensors, driver state, allocator reservation and triangular-solve workspace.')


EIGEN_KEYS = ('syevd_device_workspace_bytes', 'syevd_host_workspace_bytes',
              'geqrf_max_device_workspace_bytes', 'geqrf_max_host_workspace_bytes')


def jagwas_eigen_factor_workspace(census, traits, device_profile):
    """Requests of the default eigen factor at a queried K.

    torch.linalg.eigh runs one xsyevd (vectors, lower) on the K x K correlation
    and torch.linalg.qr(mode='r') one xgeqrf on the kept k x K rows
    (benchmarks/direct_jagwas_eigen_workspace_20260927.py). k is data, so the
    census holds the largest geqrf request over every k in [1, K]. eigh returns
    before QR starts, so the two workspaces are never live together and the
    factor needs the larger one. On the A100 xsyevd requests about 3 K^2 FP64
    (6.46 GB at K = 16,384), three times the correlation itself.
    """
    library = _require_context(census, traits, device_profile)
    if census.get('method') != 'eigen':
        raise ValueError('Eigen (xsyevd/xgeqrf) workspace census required')
    rows = {}
    for row in census['rows']:
        k = row['traits']
        if isinstance(k, bool) or not isinstance(k, int) or k < 1 or k in rows:
            raise ValueError('Unique positive queried trait counts required')
        for key in EIGEN_KEYS:
            if isinstance(row[key], bool) or not isinstance(row[key], int) or row[key] < 0:
                raise ValueError('Nonnegative integer workspace requests required')
        if row.get('geqrf_rows_queried') != k:
            raise ValueError('geqrf must be queried for every kept row count')
        if row.get('matrix_allocated') is not False:
            raise ValueError('Size-only workspace query required')
        rows[k] = row
    if traits not in rows:
        raise ValueError('Exact trait count was not queried; interpolation is forbidden')
    row = rows[traits]
    syevd, geqrf = row['syevd_device_workspace_bytes'], row['geqrf_max_device_workspace_bytes']
    return dict(traits=traits, method='eigen',
        syevd_device_requested_bytes=syevd, geqrf_device_requested_bytes=geqrf,
        device_requested_bytes=max(syevd, geqrf),
        device_rounded_bytes=max(_rounded(syevd), _rounded(geqrf)),
        host_requested_bytes=max(row['syevd_host_workspace_bytes'], row['geqrf_max_host_workspace_bytes']),
        library=copy.deepcopy(library),
        scope='Exact installed FP64 xsyevd (vectors, lower) and xgeqrf (largest over kept rows) workspace queries; '
              'the larger is live at once. File metadata identifies the installation, not a binary hash. Device '
              'rounding uses the existing 512-byte tensor-accounting convention. Excludes info tensors (512 B '
              'each observed), driver state and allocator reservation.')
