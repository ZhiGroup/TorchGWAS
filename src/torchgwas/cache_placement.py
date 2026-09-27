"""Cache-domain capacity from sysfs geometry and explicit worker CPUs."""
from collections import Counter

def cache_capacity_for_workers(geometry,worker_cpus,available_fraction=1.):
 if not worker_cpus or len(set(worker_cpus))!=len(worker_cpus):raise ValueError('Require distinct explicit worker logical CPUs')
 if not 0<=available_fraction<=1:raise ValueError('Available cache fraction must be in [0,1]')
 domains={};sizes={}
 for row in geometry:
  if row['type'] not in ['Data','Unified'] or row['level'] not in ['1','2','3']:continue
  cpu=row['cpu'];level='L'+row['level'];key=(level,row['shared_cpu_list']);value=row['size']
  size=int(value[:-1])*{'K':1024,'M':1024**2}[value[-1]]
  domains[cpu,level]=key
  if key in sizes and sizes[key]!=size:raise ValueError('Inconsistent physical cache size')
  sizes[key]=size
 counts=Counter(domains[cpu,level] for cpu in worker_cpus for level in ['L1','L2','L3'])
 capacities=[{level:int(sizes[domains[cpu,level]]*available_fraction/counts[domains[cpu,level]]) for level in ['L1','L2','L3']} for cpu in worker_cpus]
 return dict(worker_cpus=list(worker_cpus),worker_cache_bytes=capacities,
  domains=[dict(level=k[0],shared_cpu_list=k[1],physical_bytes=sizes[k],workers=counts[k]) for k in sorted(counts)],
  scope='Equal capacity share within each actual cache domain; background cache occupancy is an explicit available_fraction. This is a placement scenario, not evidence of historical worker placement.')
