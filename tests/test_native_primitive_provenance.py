from copy import deepcopy
import pytest
from torchgwas.decoder_work import verified_native_unit_prices


def fixture():
    flags=['gcc','-std=c99','-O2','-fPIC','-shared','-fno-strict-aliasing','-march=native']
    context=dict(cpu_before=dict(affinity=[12,13]),cpu_after=dict(affinity=[12,13]))
    source={'native/pgen_decode.c':'b'*64}
    loop=dict(production_library_sha256='a'*64,source_sha256=source,compiler='gcc11',commands=[flags],libraries={'O2':context},loop_units_per_call=1024,
              rows=[dict(primitive='onebit8_native',compiler_level='O2',loop_units_per_call=1024,cpu_seconds_per_unit=value) for value in [1e-9,1.1e-9,1.2e-9]])
    expansion=dict(library_sha256='a'*64,source_sha256=source,**context,
                   rows=[dict(primitive='expand4_int8',packed_bytes_per_call=1024,cpu_seconds_per_unit=value) for value in [2e-9,2.1e-9,2.2e-9]])
    verification=dict(executable_text_identical=True,source_sha256=source,compiler='gcc11',compile_command=flags,
                      libraries={'production':dict(sha256='a'*64,text_sha256='c'*64),'rebuilt':dict(text_sha256='c'*64)})
    return deepcopy((loop,expansion,verification))


def accept(reports):
    return verified_native_unit_prices(*reports,library_sha256='a'*64,source_sha256='b'*64,cpu_affinity=[12,13])


def test_only_independent_matching_native_units_are_selected():
    reports=fixture()
    reports[0]['rows'].append(dict(primitive='onebit8_native',compiler_level='O3',loop_units_per_call=1024,cpu_seconds_per_unit=1e-12))
    result=accept(reports)
    assert result['prices']=={'onebit8_native':1.1e-9,'expand4_int8':2.1e-9}
    assert result['unpriced_terms']  # Other loops are not silently certified.


def test_verified_O3_build_requires_and_selects_O3_primitives():
    loop,expansion,verification=fixture()
    verification['compile_command']=['-O3' if value=='-O2' else value for value in verification['compile_command']]
    with pytest.raises(ValueError,match='flags'):
        accept((loop,expansion,verification))
    loop['commands'].append(list(verification['compile_command']))
    loop['libraries']['O3']=deepcopy(loop['libraries']['O2'])
    loop['rows'].extend(dict(primitive='onebit8_native',compiler_level='O3',loop_units_per_call=1024,
                             cpu_seconds_per_unit=value) for value in [2e-10,3e-10,4e-10])
    result=accept((loop,expansion,verification))
    assert result['compiler_level']=='O3'
    assert result['prices']['onebit8_native']==3e-10
    assert result['prices']['expand4_int8']==2.1e-9


@pytest.mark.parametrize('flags',[['-O2','-O3'],['-Ofast'],[],['-O1']])
def test_unsupported_or_conflicting_verified_optimization_is_refused(flags):
    loop,expansion,verification=fixture()
    verification['compile_command']=[value for value in verification['compile_command'] if value!='-O2']+flags
    with pytest.raises(ValueError,match='optimization'):
        accept((loop,expansion,verification))


@pytest.mark.parametrize('extra',['-fno-tree-vectorize','-ffast-math','-march=x86-64'])
def test_unverified_codegen_flags_cannot_reuse_the_same_primitive_price(extra):
    reports=fixture()
    reports[0]['commands'][0].append(extra)
    with pytest.raises(ValueError,match='compiler flags'):
        accept(reports)


@pytest.mark.parametrize('mismatch',['binary','source','text','compiler','flags','affinity','extent','nonfinite','few_repeats'])
def test_mismatched_or_invalid_prices_are_refused(mismatch):
    reports=fixture();loop,expansion,verification=reports
    if mismatch=='binary':expansion['library_sha256']='d'*64
    elif mismatch=='source':loop['source_sha256']['native/pgen_decode.c']='e'*64
    elif mismatch=='text':verification['libraries']['rebuilt']['text_sha256']='d'*64
    elif mismatch=='compiler':loop['compiler']='gcc12'
    elif mismatch=='flags':loop['commands'][0].append('-O3')
    elif mismatch=='affinity':expansion['cpu_after']['affinity']=[14,15]
    elif mismatch=='extent':loop['rows'][0]['loop_units_per_call']=2048
    elif mismatch=='nonfinite':expansion['rows'][0]['cpu_seconds_per_unit']=float('nan')
    elif mismatch=='few_repeats':loop['rows']=loop['rows'][:2]
    with pytest.raises(ValueError):accept(reports)
