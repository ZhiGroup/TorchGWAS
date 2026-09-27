"""Inspectable first-principles counts; no incomplete end-to-end predictions."""
import argparse,json,math
from pathlib import Path
from .first_principles import Cohort,Availability,work_ledger,necessary_resource_service
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--method',choices=['torchGWAS','PLINK2','fastGWA'],required=True)
p.add_argument('--samples',type=int,required=True)
p.add_argument('--markers',type=int,required=True)
p.add_argument('--covariates',type=int,default=8)
p.add_argument('--traits',type=int,default=1)
p.add_argument('--missing-rate',type=float,default=.001)
p.add_argument('--genotype-bytes',type=int,default=0)
p.add_argument('--metadata-bytes',type=int,default=0)
p.add_argument('--native-output-bytes',type=int)
p.add_argument('--pgen',type=Path,help='Census actual hardcall record work; dimensions must match the entire file')
p.add_argument('--resources',type=Path,help='Explicit available capacities and capacity ceilings; produces partial service floor only')
p.add_argument('--tensor-work',action='store_true',help='Include actual eager Torch operation/dtype/alias counts, without running GPU kernels or timing a workload')
p.add_argument('--chunk-markers',type=int,default=2048)
p.add_argument('--output',type=Path)
a=p.parse_args()
c=Cohort(a.samples,a.markers,a.covariates,a.traits,a.missing_rate,a.genotype_bytes,a.metadata_bytes,a.native_output_bytes)
r=work_ledger(a.method,c)
if a.pgen:
 from dataclasses import replace
 from .pgen_work_census import census
 encoded=census(a.pgen)
 if (encoded['samples'],encoded['markers'])!=(a.samples,a.markers):
  p.error('--pgen dimensions must match --samples and --markers; subset census is not implicit')
 if a.genotype_bytes and a.genotype_bytes!=encoded['file_bytes']:
  p.error('--genotype-bytes disagrees with actual file length')
 c=replace(c,stored_genotype_bytes=encoded['file_bytes'])
 r=work_ledger(a.method,c)
 r['encoded_data_work']=encoded
 from .decoder_work import decoder_work
 choices={'torchGWAS':['torch_native_int8'],'fastGWA':['pgenlib_sse2'],'PLINK2':['pgenlib_avx2_bmi2']}
 r['reader_source_work']=[decoder_work(encoded,impl,restart_ld_bases=(a.method!='fastGWA')) for impl in choices[a.method]]
 r['reader_isa_scope']='Selected implementation assumption, not automatic verification of an arbitrary installed binary; specify matching source/build before pricing.'
 r['unpriced_mechanisms']=[x.replace('PGEN record census, replay dependencies and decoder instruction service','PGEN reader-specific replay schedule and decoder instruction service') for x in r['unpriced_mechanisms']]

if a.tensor_work:
 if a.method!='torchGWAS':p.error('--tensor-work requires --method torchGWAS')
 if a.chunk_markers<1:p.error('--chunk-markers must be positive')
 from .tensor_work import eager_statistics_work
 full,tail=divmod(a.markers,a.chunk_markers)
 r['eager_tensor_work']=dict(full_chunks=full,tail_markers=tail,chunk_markers=a.chunk_markers,
   full_chunk=(eager_statistics_work(a.samples,a.chunk_markers,a.traits,a.covariates) if full else None),
   tail_chunk=(eager_statistics_work(a.samples,tail,a.traits,a.covariates) if tail else None),
   interpretation='Symbolic tensor work only; no GPU execution, runtime anchor or extrapolated timing.')

if a.resources:
 spec=json.loads(a.resources.read_text())
 r['resource_service']=necessary_resource_service(r,Availability(**spec['availability']),
     cpu_flops_per_core_second=spec['cpu_flops_per_core_second'],
     gpu_flops_per_second=spec['gpu_flops_per_second'],cpu_workers=spec['cpu_workers'])
 r['resource_provenance']=spec.get('provenance','UNSPECIFIED: resource floor is conditional on user-supplied ceilings')
if a.resources:
 unbounded=[name for name,value in r['resource_service']['resource_seconds'].items() if math.isinf(value)]
 r['resource_service']['zero_capacity_resources']=unbounded
 if unbounded:
  r['resource_service']['partial_resource_floor_seconds']=None
  r['resource_service']['interpretation']='Positive work has zero available capacity; no finite service time'
  r['resource_service']['resource_seconds']={key:(None if math.isinf(value) else value) for key,value in r['resource_service']['resource_seconds'].items()}
s=json.dumps(r,indent=2,allow_nan=False)
if a.output:a.output.write_text(s+'\n')
else:print(s)
