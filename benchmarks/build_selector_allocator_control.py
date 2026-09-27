"""Build the experiment-only NumPy allocator wrapper in a named output folder."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sysconfig
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('selector_allocator_control.c')
    binary=out/('selector_allocator_control'+sysconfig.get_config_var('EXT_SUFFIX'))
    command=['cc','-shared','-fPIC','-O2','-std=c11','-Wall','-Wextra',
        '-I'+sysconfig.get_path('include'),'-I'+np.get_include(),str(source),'-o',str(binary)]
    result=subprocess.run(command,text=True,capture_output=True)
    (out/'build.txt').write_text(result.stdout+result.stderr)
    result.check_returncode()
    (out/'build.json').write_text(json.dumps(dict(command=command,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),numpy=np.__version__),indent=2)+'\n')


if __name__=='__main__':main()
