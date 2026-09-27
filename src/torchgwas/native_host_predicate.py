"""Explicitly built, optional CPU predicate for contiguous FP32 statistics.

Build with ``python -m torchgwas.native_host_predicate --build`` and enable
``TORCHGWAS_HOST_PREDICATE=native``. Compilation never happens in a scan.
The immutable build directory binds the source/flags recipe; execution context
also binds the actual compiler output bytes. NumPy handles unsupported layouts.
"""
from __future__ import annotations
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import threading

FLAGS = ('-std=c++17', '-O3', '-fPIC', '-shared', '-fvisibility=hidden',
         '-fno-fast-math', '-ffp-contract=off')
_LOADED = None
_LOCK = threading.Lock()


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _identity():
    source = Path(__file__).with_name('_host_predicate.cpp')
    identity = dict(cpp_sha256=_sha(source), binding_sha256=_sha(__file__),
                    flags=list(FLAGS), architecture=platform.machine())
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return source, identity, key


def _directory(key):
    return Path(__file__).resolve().parents[2]/'.build-libs'/('host_predicate_'+key)


def _open(directory, key):
    record = json.loads((directory/'build.json').read_text())
    path = directory/'predicate.so'
    if record['build_key'] != key or record['binary_sha256'] != _sha(path):
        raise RuntimeError('Native host predicate build identity mismatch')
    library = ctypes.CDLL(str(path))  # CDLL releases the Python GIL during calls.
    library.torchgwas_predicate_build_key.restype = ctypes.c_char_p
    if library.torchgwas_predicate_build_key().decode() != key:
        raise RuntimeError('Native host predicate source identity mismatch')
    library.torchgwas_predicate.argtypes = [ctypes.c_void_p]*3 + [ctypes.c_size_t]*2 + [ctypes.c_ssize_t]
    library.torchgwas_predicate.restype = None
    return library, record


def build():
    """Publish a complete build atomically; never replace an existing record."""
    if sys.platform != 'linux':
        raise RuntimeError('Native host predicate currently requires Linux')
    source, identity, key = _identity()
    target = _directory(key)
    target.parent.mkdir(exist_ok=True)
    if target.exists():
        return _open(target, key)[1]
    compiler = os.environ.get('CXX', 'c++')
    temporary = Path(tempfile.mkdtemp(prefix='.host-predicate-', dir=target.parent))
    try:
        binary = temporary/'predicate.so'
        command = [compiler, *FLAGS, '-DTORCHGWAS_PREDICATE_BUILD_KEY="'+key+'"',
                   str(source), '-o', str(binary)]
        subprocess.run(command, check=True)
        record = dict(build_key=key, identity=identity, binary_sha256=_sha(binary),
            compiler=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0])
        (temporary/'build.json').write_text(json.dumps(record, indent=2)+'\n')
        _open(temporary, key)
        if _identity()[2] != key:
            raise RuntimeError('Native predicate source changed during build')
        try:
            os.rename(temporary, target)
        except OSError:
            if not target.is_dir(): raise
            # A concurrent builder won; validate and reuse its immutable bytes.
        return _open(target, key)[1]
    finally:
        if temporary.exists(): shutil.rmtree(temporary)


def library():
    global _LOADED
    if _LOADED is None:
        with _LOCK:
            if _LOADED is None:
                _, _, key = _identity()
                try:
                    _LOADED = _open(_directory(key), key)
                except (OSError, ValueError, KeyError) as e:
                    raise RuntimeError('Build the requested native host predicate with '
                        'python -m torchgwas.native_host_predicate --build') from e
    return _LOADED[0]


def context(*, digest_cache=None):
    """Bind actual loaded bytes and reject a replaced binary or source recipe."""
    library()
    loaded, record = _LOADED
    if _identity()[2] != record['build_key']:
        raise ValueError('Native host predicate source changed after loading')
    digest = (_sha(loaded._name) if digest_cache is None
              else digest_cache.digests([loaded._name])[0])
    if digest != record['binary_sha256']:
        raise ValueError('Native host predicate binary changed after loading')
    return dict(build_key=record['build_key'], binary_sha256=digest)


def fill_mask(values, limits, out):
    """Fill an eligible matrix, returning False for the NumPy fallback.

    The caller has rounded FP32 thresholds upward and broadcast them. All
    pointer/shape checks occur before entering C++; buffers stay alive across
    the GIL-released call. Read-only input and reversed/scalar row limits work.
    """
    import numpy as np
    if (values.ndim != 2 or values.dtype != np.dtype('float32')
            or not values.flags.c_contiguous or not values.flags.aligned
            or limits.shape != values.shape or limits.dtype != np.dtype('float32')
            or not limits.flags.aligned or limits.strides[0] % 4
            or (values.shape[1] > 1 and limits.strides[1] != 0)):
        return False
    if (out.shape != values.shape or out.dtype != np.dtype('bool')
            or not out.flags.c_contiguous or not out.flags.writeable
            or np.shares_memory(out, values) or np.shares_memory(out, limits)):
        raise ValueError('Native predicate requires a separate writable C boolean mask')
    native = library()
    if values.size or values.shape[0] == 0:
        native.torchgwas_predicate(values.ctypes.data, limits.ctypes.data, out.ctypes.data,
            *values.shape, limits.strides[0]//4)
    return True


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(build(), indent=2))
