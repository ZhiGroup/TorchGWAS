"""Optional exact CPU event solver; the Python scheduler remains the fallback.

Build explicitly with ``python -m torchgwas.native_event_solver --build``.
The binary cache is keyed by both sources, build flags and the Python ABI.
No CUDA code, association timings or new resource coefficients are involved.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import sysconfig
import tempfile
import os

FLAGS=('-std=c++17','-O3','-fPIC','-shared','-fvisibility=hidden','-fno-fast-math','-ffp-contract=off')
_MODULE=None
_KEY=None


def _identity():
    source=Path(__file__).with_name('_event_solver.cpp')
    identity=dict(cpp_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        binding_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        flags=FLAGS,soabi=sysconfig.get_config_var('SOABI'),python=sys.version,
        platform=sys.platform)
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    return source,identity,key


def _path(key):
    return Path(__file__).resolve().parents[2]/'.build-libs'/('_event_solver_'+key+sysconfig.get_config_var('EXT_SUFFIX'))


def build():
    """Compile a CPU-only extension and atomically publish a matching binary."""
    if sys.platform!='linux':raise RuntimeError('Native event solver build currently requires Linux')
    source,identity,key=_identity();target=_path(key);target.parent.mkdir(exist_ok=True)
    compiler=os.environ.get('CXX','c++')
    fd,name=tempfile.mkstemp(prefix='.event-solver-',suffix='.so',dir=target.parent);os.close(fd)
    temporary=Path(name)
    try:
        command=[compiler,*FLAGS,'-I'+sysconfig.get_paths()['include'],
            '-DTORCHGWAS_SOLVER_BUILD_KEY="'+key+'"',str(source),'-o',str(temporary)]
        subprocess.run(command,check=True)
        module=_load(temporary,key)
        if module.BUILD_KEY!=key:raise RuntimeError('Native event solver identity mismatch')
        os.replace(temporary,target)
    finally:
        temporary.unlink(missing_ok=True)
    compiler_version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]
    record=dict(path=str(target),build_key=key,identity=identity,compiler=compiler_version,
        binary_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    target.with_suffix(target.suffix+'.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def _load(path,key):
    spec=importlib.util.spec_from_file_location('torchgwas._event_solver',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    if module.BUILD_KEY!=key:raise ImportError('Stale native event solver; rebuild')
    return module


def library():
    global _MODULE,_KEY
    if sys.platform!='linux':return None
    try:_,_,key=_identity()
    except OSError:return None
    if _MODULE is not None and _KEY==key:return _MODULE
    path=_path(key)
    if not path.is_file():return None
    try:module=_load(path,key)
    except (OSError,ImportError):return None
    _MODULE,_KEY=module,key
    return module


def solve(graph,chains,shared_tokens,trace,prefix):
    if chr(0) in prefix:return None
    module=library()
    if module is None:return None
    try:return module.solve(graph,chains,shared_tokens,trace,prefix)
    except NotImplementedError:
        # Unusual Python key/value types retain the reference implementation's
        # behavior. Invalid supported graphs still raise their validation error.
        return None


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--build',action='store_true',required=True)
    parser.parse_args();print(json.dumps(build(),indent=2))
