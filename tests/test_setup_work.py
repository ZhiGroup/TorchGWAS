import pytest
from torchgwas.setup_work import setup_work,setup_service
from torchgwas.native_scan import design_column_block


def test_setup_counts_both_phenotype_uploads_and_one_download():
    work=setup_work(8192,512)
    assert work['h2d_bytes']==4*8192*(2*512+2*8)
    assert work['d2h_bytes']==4*8192*512
    assert work['residual_block_traits']==512
    assert work['design_block_traits']==512
    assert [p['phase'] for p in work['phases']]==['residual_common','residual_block','design_common','design_block']


def test_huge_traits_follow_two_different_source_caps_and_tails():
    n,k=32768,40001
    work=setup_work(n,k)
    assert work['residual_block_traits']==10922
    assert work['design_block_traits']==2048
    residual=[p['traits'] for p in work['phases'] if p['phase']=='residual_block']
    design=[p['traits'] for p in work['phases'] if p['phase']=='design_block']
    assert sum(residual)==sum(design)==k
    assert residual[-1]==k%10922 and design[-1]==k%2048
    assert work['device_live_bytes_upper']>=4*n*(k+9)
    assert work['h2d_bytes']==4*n*(2*k+16)


def test_tiny_setup_price_is_exact_reference_and_wide_adds_work():
    bank={p:dict(reference_shape=[32,1,8],cpu_seconds=.01,non_cpu_seconds=.02)
          for p in ['residual_common','residual_block','design_common','design_block']}
    profile=dict(setup_primitives=bank,cpu_fraction=1.,gpu_resources=dict(gpu_fraction=1.,hbm_bytes_per_second=1e9,fp32_flops_per_second=1e10),h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8,process_units=dict(numpy_copy_bytes=1e-10))
    tiny=setup_service(setup_work(32,1),profile)
    assert tiny['seconds']==pytest.approx(.12)
    assert setup_service(setup_work(32,512),profile)['seconds']>tiny['seconds']
    assert setup_service(setup_work(8192,512),profile)['seconds']>tiny['seconds']
    assert tiny['unpriced_terms']
    baseline=setup_service(setup_work(4096,2048),profile)
    profile['reduction_gpu_properties']=dict(sm_count=132,max_threads_per_sm=2048)
    accounted=setup_service(setup_work(4096,2048),profile)
    assert accounted['seconds']-baseline['seconds']==pytest.approx(6*1048576/1e9)
    assert accounted['phases'][-1]['reduction']['workspace_bytes']==67108864


def test_invalid_dimensions_and_reference_fail_closed():
    with pytest.raises(ValueError):setup_work(32,0)
    with pytest.raises(ValueError):setup_work(32,1,-1)
    with pytest.raises(ValueError):design_column_block(0,5)


def test_full_tiles_reuse_qc_counts_but_raw_preparation_scans_once():
    raw=setup_work(4096,4096)
    cached=setup_work(4096,4096,reuse_observed_counts=True)
    assert raw['observation_cpu_work']['phenotype_isnan_cells']==4096**2
    assert raw['observation_cpu_work']['column_count_reduction_cells']==4096**2
    assert cached['observation_cpu_work']==dict(phenotype_isnan_cells=0,boolean_not_cells=0,column_count_reduction_cells=0,count_vector_cells=4096)
    assert cached['phases']==raw['phases']

def test_setup_host_copy_layout_and_direct_download_lifetime():
    contiguous=setup_work(4096,4096,input_contiguous=True)
    strided=setup_work(4096,4096,input_contiguous=False)
    assert contiguous['host_copy_bytes']==0
    assert strided['host_copy_bytes']==4*4096*4096
    assert contiguous['residual_result_assembly_bytes']==strided['residual_result_assembly_bytes']==0
    many=setup_work(32768,40001,input_contiguous=True)
    assert many['residual_result_assembly_bytes']==4*32768*40001
    # Each input residual sub-block and each design sub-block needs a copy,
    # plus one assembly copy of every downloaded residual element.
    assert many['host_copy_bytes']==3*4*32768*40001
    with pytest.raises(ValueError):setup_work(32,1,input_contiguous=1)
    # Unspecified input layout is an explicit upper bound, never silently C.
    assert setup_work(4096,4096)['host_copy_bytes']==strided['host_copy_bytes']
