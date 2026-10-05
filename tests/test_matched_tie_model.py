import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'benchmarks'))
_model=pytest.importorskip('direct_matched_tie_model')
integer_crossings,validation_traits=_model.integer_crossings,_model.validation_traits

def test_discrete_crossing_and_validation_neighbors():
 rows=[{'traits':k,'torchGWAS':10.,'fast':float(3*k)} for k in range(1,9)]
 hits=integer_crossings(rows,'fast')
 assert hits==[{'lower_traits':3,'upper_traits':4,'direction':'torch becomes faster'}]
 assert validation_traits(rows,hits)==[1,2,3,4,5,8]

def test_no_crossing_does_not_extrapolate():
 rows=[{'traits':k,'torchGWAS':10.,'plink':float(k)} for k in range(1,9)]
 assert integer_crossings(rows,'plink')==[]
 assert validation_traits(rows,[])==[1,8]

def test_nonmonotone_and_exact_integer_tie():
 rows=[{'traits':k,'torchGWAS':float(v),'other':5.} for k,v in enumerate([8,5,4,7],1)]
 assert [r['upper_traits'] for r in integer_crossings(rows,'other')]==[2,4]
 with pytest.raises(ValueError):integer_crossings(rows[::2],'other')


def test_binary_schedule_and_fast_trait_restarts():
 from direct_matched_tie_model import predict
 policy={'OMP_WAIT_POLICY':'PASSIVE'}
 t={'phase_wall_seconds':1.,'metadata_seconds':0.,'producer_seconds_per_chunk':.1,'consumer_seconds_per_chunk':.2,'writer_seconds_per_chunk':.3,'first_use_seconds':0.}
 c={'policy':policy,'torch_zero':[{'process_seconds':2.,'metadata_seconds':0.,'calibration_variants':4096}],
    'torch':{str(k):dict(t) for k in range(1,33)},
    'plink':{str(k):{'worker_elapsed_per_variant':.0001,'format_seconds_per_row':.00001} for k in range(1,33)},
    'fast':{'0':{'CYCLE_SECONDS_PER_BLOCK_T4_READER_MEAN':.01}},'fast_setup':[{'fixed':1.,'pvar_seconds_per_byte':0.}]}
 s={'binary_writer':{str(k):{'variants':65536,'append_seconds':0.,'storage_seconds':0.,'final_seconds':.3,'block_bytes':1<<20,'queue_depth':3} for k in range(1,33)},'variants':65536,'policy':policy,'pvar_bytes':100,'pgen_bytes':1000,'warm_read_bytes_per_second':10000.,
    'torch_preparation':{str(k):{'seconds_excluding_metadata':1.} for k in range(1,33)},'record_forms':{0:65536},'fast_decode_by_form':{'0':.000001},'plink_preloop':{str(k):{'seconds':3.} for k in range(1,33)}}
 rows=predict(c,s)
 assert rows[0]['torchGWAS']==pytest.approx(8.8)
 # Historical text-writer fields must never affect native binary prediction.
 for profile in c['torch'].values():profile['writer_seconds_per_chunk']=1e9
 assert predict(c,s)==rows
 assert rows[1]['fastGWA-joint']==pytest.approx(2*rows[0]['fastGWA-joint'])
 assert rows[0]['PLINK2']==pytest.approx(5.39376)
 with pytest.raises(ValueError):predict(c,s,16384)
 s['policy']={}
 with pytest.raises(ValueError):predict(c,s)


def test_binary_writer_releases_blocks_and_applies_backpressure():
 from torchgwas.pipeline_model import binary_output_pipeline_seconds as schedule
 # Four chunks, one block per chunk, writer slower than the scan.
 assert schedule(1.,2.,4,1,1,0.,16.,.5,4,1)==pytest.approx(19.5)
 # One partial block cannot begin until the last chunk arrives.
 assert schedule(1.,2.,4,1,1,0.,16.,.5,32,3)==pytest.approx(25.5)
 # Zero-byte service reduces to producer/consumer finite pipeline.
 assert schedule(1.,2.,4,1,1,0.,0.,0.,4,3)==pytest.approx(9.)
 with pytest.raises(ValueError):schedule(1.,2.,4,0,1,0.,0.,0.)
