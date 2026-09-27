"""PyTorch 2.5 Reduce.cuh workspace for contiguous FP32 column reductions."""
def _column_workspace(samples, traits, *, sm_count, max_threads_per_sm, accumulator_bytes):
    for v in (samples,traits,sm_count,max_threads_per_sm):
        if isinstance(v,bool) or not isinstance(v,int) or v<1:
            raise ValueError('Positive integer dimensions and device properties required')
    if samples*traits>=2**31:
        raise ValueError('TensorIterator 32-bit splitting requires a separate census')
    up=lambda a,b:(a+b-1)//b
    vec=4
    while traits%vec:vec//=2
    limit=512//vec
    pow2=lambda n:1<<(n.bit_length()-1)
    d0=min(pow2(traits//vec),limit);d1=min(pow2(samples),limit)
    width=min(d0,32);height=min(d1,limit//width)
    width=min(d0,limit//height)
    threads=width*height
    step_input=1;step_output=width
    reduce_y=samples>=height*16 or samples>=256
    if reduce_y:step_input*=height
    else:step_output*=height
    values=up(samples,step_input)
    grid=up(traits//vec,step_output)
    target=sm_count*(max_threads_per_sm//threads)
    ctas=1
    if reduce_y and values>=256 and grid<=target:
        ctas=max(min(up(target,grid),up(values,16)),up(values,256))
    workspace=accumulator_bytes*traits*ctas*width*vec if ctas>1 else 0
    semaphore=4*grid if ctas>1 else 0
    # staging_memory_offset indexes grid.x groups of CTA partial vectors.
    # The allocation formula is larger; unreferenced capacity is not traffic.
    partial=accumulator_bytes*grid*ctas*width*vec if ctas>1 else 0
    return dict(workspace_bytes=workspace,semaphore_bytes=semaphore,accumulator_bytes=accumulator_bytes,
        partial_sum_bytes=partial,partial_sum_read_write_bytes=2*partial,
        output_vector_width=vec,block=[width,height],grid=[grid,ctas],
        scope='PyTorch 2.5 Reduce.cuh FP32 aligned contiguous matrix reduction(dim=0); no 32-bit-index splitting. Exact source allocation requests, not reserved memory.')

def column_sum_workspace(samples, traits, *, sm_count, max_threads_per_sm):
    """FP32 sum/mean accumulator: one float (4 bytes)."""
    return _column_workspace(samples, traits, sm_count=sm_count,
        max_threads_per_sm=max_threads_per_sm, accumulator_bytes=4)


def column_std_workspace(samples, traits, *, sm_count, max_threads_per_sm):
    """FP32 Welford accumulator: mean, m2, int32 count, float count (16 bytes).

    PyTorch v2.5.1 SharedReduceOps.h WelfordData and ReduceMomentKernel.cu.
    Unrolling is two rather than four, but ReduceConfig launch dimensions and
    global allocation are unchanged. This counts capacity, not HBM traffic.
    """
    return _column_workspace(samples, traits, sm_count=sm_count,
        max_threads_per_sm=max_threads_per_sm, accumulator_bytes=16)