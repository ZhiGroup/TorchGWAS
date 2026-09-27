import pytest
from torchgwas.tensor_service import gemm_work


def kernel(name,grid):
    return {'name':name,'geometry':{'grid':grid}}


def test_swizzle_padding_is_not_arithmetic_work():
    main=kernel('cutlass_80_simt_sgemm_256x128_8x4_nn_align1',[256,2,5])
    work=gemm_work(8192,4096,2057,[main])
    assert work['swizzle']==8
    assert work['grid_ctas']==32*9*5
    assert work['launched_ctas']==256*2*5
    assert work['noop_ctas']==1120
    assert work['issued_flops']==2*8200*4096*2304
    assert work['kernel_count']==1
    assert work['split_k_mode']=='in_kernel_reduction_unresolved'
    assert work['workspace_write_bytes']==0
    assert work['accumulation_logical_bytes']==8*4096*2057*4
    assert work['unpriced_terms']


def test_split_count_does_not_invent_reduction_kernel():
    main=kernel('sm80_tilesize128x64x8_execute_kernel',[32,9,6])
    one=gemm_work(8192,4096,521,[main])
    two=gemm_work(8192,4096,521,[main,kernel('sm80_execute_split_k_kernel',[32,9,2])])
    assert one['kernel_count']==1 and two['kernel_count']==2
    assert one['issued_flops']==two['issued_flops']
    assert two['workspace_write_bytes']>0
    assert not two['unpriced_terms']


def test_swizzled_separate_reduction_and_invalid_grid():
    main=kernel('cutlass_80_simt_sgemm_256x128_8x4_nn_align1',[16,2,7])
    reduction=kernel('cublasLt::splitKreduce_kernel',[65,8,1])
    work=gemm_work(8192,256,2057,[main,reduction])
    assert work['kernel_count']==2 and work['grid_ctas']==2*9*7
    with pytest.raises(ValueError,match='swizzled'):
        gemm_work(8192,256,2057,[kernel(main['name'],[17,2,7])])


def test_unsplit_and_old_rectangular_grid():
    work=gemm_work(8192,1024,10,[kernel('sm80_tilesize64x32x8_execute_kernel',[16,1,16]),
        kernel('sm80_execute_split_k_kernel',[16,1,2])])
    assert work['grid_ctas']==256 and work['noop_ctas']==0
    assert work['kernel_count']==2 and work['swizzle']==1
    plain=gemm_work(8192,128,32,[kernel('sm80_tilesize128x32x8_execute_kernel',[1,1,1])])
    assert plain['split_k_mode']=='none' and plain['accumulation_logical_bytes']==0


def gemv(width,split):
    import math
    main=kernel('internal::gemvx::kernel<int, int, float, float, float, float, false, true, false, false, 9, false, cublasGemvParamsEx<int, cublasGemvTensorStridedBatched<float const> > >',[math.ceil(width/64),1,split])
    main['geometry']['block']=[64,8,1]
    reduction=kernel('cublasLt::splitKreduce_kernel<32,16,int,float,float,float,float,true,false,false>',[1,math.ceil(width/16),1])
    reduction['geometry']['block']=[32,16,1]
    return [main,reduction] if split>1 else [main]


@pytest.mark.parametrize('width,split',[(521,9),(8201,1)])
def test_singleton_gemv_counts_logical_work_without_inventing_tiles(width,split):
    work=gemm_work(16384,1,width,gemv(width,split))
    assert work['useful_flops']==work['issued_flops']==2*16384*width
    assert work['arithmetic_accounting']=='logical_floor_not_verified_issued_work'
    assert work['tile'] is None and not work['inner_k_tile_verified']
    assert work['grid_ctas']==work['launched_ctas']
    assert work['input_l2_bytes']==4*16384*(width+1)
    assert work['workspace_write_bytes']==work['workspace_reduce_bytes']==(4*width*split if split>1 else 0)
    assert work['reduce_adds']==width*(split-1)
    assert work['kernel_count']==1+(split>1) and work['unpriced_terms']


def test_gemv_fails_closed_for_unverified_shapes_and_launches():
    with pytest.raises(ValueError,match='singleton GEMV'):
        gemm_work(16384,2,521,gemv(521,9))
    with pytest.raises(ValueError,match='explicit reduction'):
        gemm_work(16384,1,521,gemv(521,9)[:1])
    kernels=gemv(521,9);kernels[1]['geometry']['grid'][1]+=1
    with pytest.raises(ValueError,match='reduction launch'):
        gemm_work(16384,1,521,kernels)
    kernels=gemv(521,9);kernels[0]['geometry']['block']=[16,32,1]
    with pytest.raises(ValueError,match='singleton GEMV'):
        gemm_work(16384,1,521,kernels)
    with pytest.raises(ValueError,match='launch census'):
        gemm_work(16384,1,521,gemv(521,9)+[kernel('sm80_tilesize128x32x8_execute_kernel',[1,17,1])])


@pytest.mark.parametrize('width,block,specialization,split',[
    (265,[8,32,1],8,1),(521,[8,32,1],8,1),(1034,[32,16,1],9,3),(2057,[128,4,1],9,6)])
def test_captured_n4096_singleton_gemv_lane_layouts(width,block,specialization,split):
    import math
    kernels=gemv(width,split)
    kernels[0]['name']=kernels[0]['name'].replace('false, 9, false','false, '+str(specialization)+', false')
    kernels[0]['geometry'].update(block=block,grid=[math.ceil(width/block[0]),1,split])
    work=gemm_work(4096,1,width,kernels)
    assert work['useful_flops']==2*4096*width
    assert work['input_l2_bytes']==4*4096*(width+1)
    assert work['grid_ctas']==math.ceil(width/block[0])*split
    assert work['kernel_count']==1+(split>1)
    assert work['arithmetic_accounting']=='logical_floor_not_verified_issued_work'
    kernels[0]['geometry']['grid'][0]+=1
    with pytest.raises(ValueError,match='singleton GEMV'):gemm_work(4096,1,width,kernels)


def test_captured_n4096_narrow_singleton_dot_pair():
    import copy
    main=kernel('dot_kernel<float, 128, 0, cublasDotParams<cublasGemvTensorStridedBatched<float const>, cublasGemvTensorStridedBatched<float> > >',[31,1,10])
    reduction=kernel('reduce_1Block_kernel<float,128,7,cublasGemvTensorStridedBatched<float>,cublasGemvTensorStridedBatched<float const>,cublasGemvTensorStridedBatched<float> >',[1,1,10])
    for value in [main,reduction]:value['geometry']['block']=[128,1,1]
    work=gemm_work(4096,1,10,[main,reduction])
    assert work['useful_flops']==81920
    assert work['grid_ctas']==310
    assert work['kernel_count']==2
    assert work['workspace_write_bytes']==work['workspace_reduce_bytes']==1240
    assert work['reduce_adds']==300
    assert work['arithmetic_accounting']=='logical_floor_not_verified_issued_work'
    with pytest.raises(ValueError,match='dot GEMV'):gemm_work(4096,2,10,[main,reduction])
    with pytest.raises(ValueError,match='explicit reduction'):gemm_work(4096,1,10,[main])
    malformed=copy.deepcopy(main);malformed['geometry']['grid'][2]=11
    with pytest.raises(ValueError,match='dot GEMV'):gemm_work(4096,1,10,[malformed,reduction])
    malformed=copy.deepcopy(reduction);malformed['geometry']['block']=[256,1,1]
    with pytest.raises(ValueError,match='reduction launch'):gemm_work(4096,1,10,[main,malformed])
