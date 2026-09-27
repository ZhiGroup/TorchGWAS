"""A quickly cancelled peer must not hide a slower originating writer error."""
from pathlib import Path
import threading
import pytest


@pytest.mark.parametrize('axis',['trait','variant'])
def test_original_error_survives_peer_completion_order(tmp_path,monkeypatch,axis):
    import torchgwas.sumstats_tiled as tiled
    import torchgwas.sumstats_sharded as sharded
    module=tiled if axis=='trait' else sharded
    root_aborting=threading.Event();peer_aborted=threading.Event()
    entered=threading.Barrier(2)
    class Writer:
        def __init__(self,directory,*args,**kwargs):self.root=int(Path(directory).name.split('_')[1])==0
        def write_chunk(self,*args,**kwargs):raise OSError('original writer fault')
        def close(self):pytest.fail('failed store cannot close successfully')
        def abort(self):
            if self.root:
                root_aborting.set()
                assert peer_aborted.wait(5.)
                # Make the peer future observable before the real error.
                threading.Event().wait(.1)
            else:peer_aborted.set()
    monkeypatch.setattr(module,'BinarySumstatsWriter',Writer)
    def scan(first,width,device,workers):
        entered.wait(5.)
        if first:assert root_aborting.wait(5.)
        yield (0,1,None,None,None,None)
    options=dict(n_variants=4,trait_names=['a','b'],n_samples=20,df=18,
        devices=['cuda:0','cuda:1'],reader_workers=2,scan_factory=scan)
    with pytest.raises(OSError,match='original writer fault'):
        if axis=='trait':module.write_trait_tiled_sumstats(tmp_path,trait_block=1,**options)
        else:module.write_variant_sharded_sumstats(tmp_path,chunk_size=1,**options)
    assert root_aborting.is_set() and peer_aborted.is_set()
    assert not (tmp_path/'manifest.json').exists()
    assert not any(t.name.startswith('torchgwas-') for t in threading.enumerate())
