"""Inspect source/resource CPU candidates; validation gaps remain explicit."""
import argparse,copy,json
from pathlib import Path
from .mechanistic_cpu import fastgwa_runtime
from .mechanistic_plink import plink_runtime
from .mechanistic_torch import torch_runtime, torch_scan_runtime
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--method',choices=['fastGWA','PLINK2','torchGWAS'],required=True)
p.add_argument('--inputs',type=Path,required=True,help='Input work census JSON, including PGEN and text statistics')
p.add_argument('--profile',type=Path,required=True,help='Independent resource profile, or a frozen predictions file containing profiles')
p.add_argument('--server',help='Profile name when --profile contains multiple servers')
p.add_argument('--cpu-fraction',type=float,help='Available scheduling fraction for each serial worker, 0 to 1')
p.add_argument('--cpu-cores',type=float,help='Total CPU capacity available to this process')
p.add_argument('--read-gbps',type=float,help='Available input bandwidth in decimal GB/s under the stated cache condition')
p.add_argument('--available-l3-mib',type=float,help='L3 capacity available to the modeled reference stream; independent of scheduling fraction')
p.add_argument('--missing-rate',type=float,help='Explicit iid missingness scenario for the PLINK branch expectation')
p.add_argument('--carrier-fraction',type=float,help='Expected nonreference carrier fraction for PLINK sparse work')
p.add_argument('--output',type=Path)
p.add_argument('--scope',choices=['process','scan'],default='process',help='scan: Torch native-int8 eager scan-to-discard only; no setup or output')
a=p.parse_args();data=json.loads(a.inputs.read_text());spec=json.loads(a.profile.read_text())
if 'profiles' in spec:
 if a.server not in spec['profiles']:p.error('--server must identify a profile in the supplied file')
 profile=copy.deepcopy(spec['profiles'][a.server])
else:profile=copy.deepcopy(spec)
for key,value in [('cpu_fraction',a.cpu_fraction),('cpu_available_cores',a.cpu_cores)]:
 if value is not None:profile[key]=value
if a.read_gbps is not None:profile['read_bytes_per_second']=a.read_gbps*1e9
if a.available_l3_mib is not None:
 if a.available_l3_mib<0:p.error('Available L3 cannot be negative')
 profile['cache_bytes']['L3']=int(a.available_l3_mib*1024**2)
 for cache in profile.get('worker_cache_bytes',[]):cache['L3']=profile['cache_bytes']['L3']
for key,value in [('missing_rate',a.missing_rate),('carrier_fraction',a.carrier_fraction)]:
 if value is not None:
  if a.method!='PLINK2' or not 0<=value<=1:p.error(key+' must be in [0,1] and applies to PLINK2')
  profile[key]=value
if a.scope=='scan' and a.method!='torchGWAS':p.error('scan scope applies only to torchGWAS')
result=(torch_scan_runtime(data,profile) if a.scope=='scan' else {'fastGWA':fastgwa_runtime,'PLINK2':plink_runtime,'torchGWAS':torch_runtime}[a.method](data,profile))
result['validation_status']='Candidate only: source/resource approximations and incomplete validation remain explicit. Not approved as a crossover guarantee.'
result['dimensions']={k:data[k] for k in ['samples','markers','covariates','traits_analyzed','traits_in_file']}
text=json.dumps(result,indent=2,allow_nan=False)
if a.output:a.output.write_text(text+'\n')
else:print(text)

