import pytest
from torchgwas.cublas_memory import cublas_workspace

def test_hopper_default_and_distinct_handle_stream_pairs():
    assert cublas_workspace((9,0),handle_stream_pairs=2)['total_bytes']==64*1024**2
    assert cublas_workspace((8,0))['total_bytes']==8519680

def test_configuration_parser_matches_source():
    assert cublas_workspace((9,0),':4096:2:16:8')['total_bytes']==8519680
    assert cublas_workspace((9,0),':16:8')['total_bytes']==131072
    assert cublas_workspace((9,0),':0:0')['total_bytes']==0
    assert cublas_workspace((9,0),'bad')['config_status']=='invalid_uses_default'
    assert cublas_workspace((9,0),'prefix:16:8suffix')['total_bytes']==131072

def test_invalid_counts_and_overflow_are_rejected():
    with pytest.raises(ValueError):cublas_workspace((9,0),handle_stream_pairs=0)
    with pytest.raises(ValueError):cublas_workspace((9,0),':2147483648:1')