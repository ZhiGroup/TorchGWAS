"""Conditional source/resource crossover candidate, with explicit scenario inputs."""
import argparse,copy,json
from pathlib import Path
from .crossover import crossover_at_subjects
from .timing_boundary import timing_profiles
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bundle',type=Path,required=True,help='Frozen three-method profiles, sort counts and untimed geometry')
p.add_argument('--server',required=True);p.add_argument('--subjects',type=int,required=True)
p.add_argument('--maf',type=float,default=.2);p.add_argument('--missing-rate',type=float,default=.001)
p.add_argument('--traits-in-file',type=int,default=32);p.add_argument('--numeric-characters',type=float,default=19.4)
p.add_argument('--cpu-fraction',type=float);p.add_argument('--cpu-cores',type=float);p.add_argument('--read-gbps',type=float)
p.add_argument('--timing-boundary',choices=['environment-ready','process'],default='environment-ready',help='Default starts after imports, thread setup and CUDA context readiness, before input loading')
p.add_argument('--startup-seconds',type=float,help='Explicit Torch import/thread-pool/process entry-exit elapsed service at the target environment load; excludes CUDA and data-dependent setup')
p.add_argument('--max-markers',type=int,default=8388608);p.add_argument('--output',type=Path)
a=p.parse_args();bundle=json.loads(a.bundle.read_text())
if a.server not in bundle['profiles']:p.error('Unknown server')
if str(a.subjects) not in bundle['sort_counts']:p.error('Supply untimed sort/compiled-work census for this subject count')
profiles=timing_profiles(bundle['profiles'][a.server],a.timing_boundary)
for profile in profiles.values():
 if a.cpu_fraction is not None:profile['cpu_fraction']=a.cpu_fraction
 if a.cpu_cores is not None:profile['cpu_available_cores']=a.cpu_cores
 if a.read_gbps is not None:profile['read_bytes_per_second']=a.read_gbps*1e9
if a.startup_seconds is not None:
 if a.startup_seconds<0:p.error('Startup service must be nonnegative')
 profiles['torchGWAS']['environment_events']=dict(events=[dict(name='user_supplied_environment',count=1,seconds_per_event=a.startup_seconds)],condition='Explicit user-supplied independent environment service')
if a.timing_boundary=='process' and 'environment_events' not in profiles['torchGWAS']:p.error('Legacy startup accounting is withdrawn. Supply inputs_v2.json or explicit --startup-seconds.')
rows=crossover_at_subjects(a.subjects,profiles,bundle['sort_counts'][str(a.subjects)],maf=a.maf,missing_rate=a.missing_rate,traits_in_file=a.traits_in_file,numeric_characters=a.numeric_characters,max_markers=a.max_markers)
result=dict(server=a.server,timing_boundary=a.timing_boundary,rows=rows,status='development_candidate_not_validated_bound',scope='Numerical runtime equality under specified inputs and resources. No runtime interpolation or hyperbola fitting. See model validation failures before scientific use.')
text=json.dumps(result,indent=2,allow_nan=False)
if a.output:a.output.write_text(text+'\n')
else:print(text)
