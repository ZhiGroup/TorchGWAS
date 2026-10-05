import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'benchmarks'))
crossings=pytest.importorskip('direct_marker_ties').crossings

def test_bisection_returns_marker_and_product_brackets():
    result=crossings(lambda m:{'samples':2000,'torchGWAS':10+m*.00001,'PLINK2':m*.00002},'PLINK2',minimum=1000,maximum=2000000,tolerance=10)
    assert len(result)==1
    r=result[0]
    assert r['lower_markers']<=1000000<=r['upper_markers']
    assert r['upper_markers']-r['lower_markers']<=10
    assert r['upper_NxM']==2000*r['upper_markers']

def test_no_crossing_and_reversal_are_not_extrapolated():
    assert crossings(lambda m:{'samples':100,'torchGWAS':1,'PLINK2':2},'PLINK2',maximum=20000)==[]
    r=crossings(lambda m:{'samples':100,'torchGWAS':m,'PLINK2':10000},'PLINK2',maximum=20000)
    assert r[0]['direction']=='torch becomes slower'

@pytest.mark.parametrize('minimum,maximum,tolerance',[(0,1000,1),(100,10,1),(1,10,0)])
def test_invalid_search(minimum,maximum,tolerance):
    with pytest.raises(ValueError):crossings(None,'PLINK2',minimum,maximum,tolerance)
