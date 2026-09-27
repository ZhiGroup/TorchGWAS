import math,unittest
from dataclasses import replace
from torchgwas.resource_model import ResourceState,StageDemand,shared_resource_floor,interpolate_unit_service,bounded_gpu_pipeline
class ResourceAccounting(unittest.TestCase):
 def setUp(self):self.r=ResourceState(1.,1.,4.,100.,50.,25.,1.)
 def test_cpu_capacity(self):
  self.assertEqual(self.r.cpu_wall(8,4),2)
  self.assertEqual(replace(self.r,cpu_fraction_workers=.5).cpu_wall(8,4),4)
  self.assertEqual(replace(self.r,cpu_capacity_cores=1).cpu_wall(8,4),8)
 def test_cache_and_physical_read(self):
  self.assertEqual(self.r.input_io_wall(100),0)
  self.assertEqual(replace(self.r,input_cache_hit_fraction=0).input_io_wall(100),2)
  self.assertEqual(replace(self.r,input_cache_hit_fraction=.25).input_io_wall(100),1.5)
 def test_shared_bandwidth(self):
  stages=[StageDemand('decode',memory_bytes=100),StageDemand('copy',memory_bytes=100)]
  self.assertEqual(max(s.seconds(self.r) for s in stages),1)
  self.assertEqual(shared_resource_floor(stages,self.r),2)
 def test_no_double_memory_charge(self):
  self.assertEqual(StageDemand('work',cpu_seconds=2,memory_bytes=100).seconds(self.r),2)
 def test_monotonic_constraints(self):
  stages=[StageDemand('cpu',cpu_seconds=8,threads=4,memory_bytes=100),StageDemand('GPU',gpu_seconds=2,memory_bytes=50,write_bytes=20)]
  before=shared_resource_floor(stages,self.r)
  for state in (replace(self.r,memory_bytes_per_second=10),replace(self.r,cpu_capacity_cores=1),replace(self.r,gpu_service_fraction=.2),replace(self.r,storage_write_bytes_per_second=1)):
   self.assertGreaterEqual(shared_resource_floor(stages,state),before)
 def test_invalid_unavailable_resources(self):
  for key in ('cpu_fraction_workers','memory_bytes_per_second','storage_read_bytes_per_second','gpu_service_fraction'):
   for value in (0,-1,float('nan')):
    with self.assertRaises(ValueError):replace(self.r,**{key:value})
 def test_work_units_not_total_fit(self):
  samples=[(100,1.),(400,4.)]
  self.assertAlmostEqual(interpolate_unit_service(samples,200,lambda n:n),2.)
  with self.assertRaises(ValueError):interpolate_unit_service(samples,800,lambda n:n)
class PipelineGeometry(unittest.TestCase):
 def test_single_chunk_serial_dependencies(self):
  r=bounded_gpu_pipeline(1024,1024,4,4,2,0,1,3,4,5)
  self.assertEqual(r['seconds'],15)
 def test_work_cannot_vanish(self):
  r=bounded_gpu_pipeline(100*1024,1024,4,4,4,0,0,1,2,0)
  self.assertGreaterEqual(r['seconds'],200)
  self.assertGreaterEqual(r['seconds'],100)
 def test_buffer_capacity_changes_throughput(self):
  small=bounded_gpu_pipeline(100*1024,1024,4,2,8,0,0,2,1,0)
  large=bounded_gpu_pipeline(100*1024,1024,4,8,8,0,0,2,1,0)
  self.assertGreater(small['seconds'],large['seconds'])
 def test_partial_tail(self):
  r=bounded_gpu_pipeline(512,1024,1,2,2,0,0,2,2,2)
  self.assertEqual(r['seconds'],4)
if __name__=='__main__':unittest.main()
