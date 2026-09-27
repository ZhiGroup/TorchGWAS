"""PyTorch 2.5.1 CUDA cuBLAS handle/stream workspace requests."""
import re
SOURCE='https://raw.githubusercontent.com/pytorch/pytorch/v2.5.1/aten/src/ATen/cuda/CublasHandlePool.cpp'

def cublas_workspace(compute_capability, config=None, *, handle_stream_pairs=1):
    if len(compute_capability)!=2 or any(isinstance(v,bool) or not isinstance(v,int) or v<0 for v in compute_capability):
        raise ValueError('Explicit CUDA compute capability required')
    if isinstance(handle_stream_pairs,bool) or not isinstance(handle_stream_pairs,int) or handle_stream_pairs<1:
        raise ValueError('Positive handle/stream pair count required')
    default=32*1024**2 if tuple(compute_capability)==(9,0) else 4096*1024*2+16*1024*8
    if config is not None and not isinstance(config,str):raise ValueError('Config must be text or None')
    matches=[] if config is None else re.findall(r':([0-9]+):([0-9]+)',config)
    if any(int(v)>2147483647 for pair in matches for v in pair):
        raise ValueError('Config exceeds C++ stoi integer range')
    per_pair=sum(int(size)*1024*int(count) for size,count in matches) if matches else default
    return dict(bytes_per_handle_stream=per_pair,total_bytes=per_pair*handle_stream_pairs,
        handle_stream_pairs=handle_stream_pairs,
        config_status='parsed' if matches else ('default' if config is None else 'invalid_uses_default'),
        source=SOURCE,scope='PyTorch 2.5.1 CUDA workspace requests; one allocation per distinct cuBLAS handle/stream pair. Excludes allocator rounding and driver memory.')