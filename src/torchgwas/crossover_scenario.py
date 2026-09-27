"""Expected source work for an explicit iid hardcall/PGEN scenario.

No runtime observation enters this input model. Expected counts can be fractional;
record counts and analysis dimensions remain integer. This is a statistical
scenario, not the exact byte census of an unspecified real dataset.
"""
from __future__ import annotations
import math
from functools import lru_cache

def curve_subjects():
 return sorted({128*round(1024*2**(i/3)/128) for i in range(19)})

def digit_sum_progression(count,start=0,step=1):
 """Exact sum of decimal character counts along a nonnegative progression."""
 if count<0 or start<0 or step<1:raise ValueError('Invalid integer progression')
 total=count;limit=10
 last=start+max(0,count-1)*step
 while limit<=last:
  first=max(0,(limit-start+step-1)//step)
  total+=max(0,count-first);limit*=10
 return total

@lru_cache(maxsize=128)
def onebit_expected_per_marker(n,maf=.2,missing_rate=.001):
 import numpy as np
 from scipy.stats import binom
 if n<1 or not 0<maf<.5 or not 0<=missing_rate<1:raise ValueError('Invalid genotype statistics')
 # Explicit two-common-category form: nonmissing 0 and 1 are represented by
 # the bit plane; category2 and missing are in the difflist.
 rare=(1-missing_rate)*maf**2+missing_rate
 entries=n*rare;nonempty=-math.expm1(n*math.log1p(-rare));q=1-rare
 groups=float(binom.sf(np.arange(0,n,64),n,rare).sum())
 category_bytes=float(binom.sf(np.arange(0,n,4),n,rare).sum())
 header={}
 for length in range(1,6):
  lo=0 if length==1 else 128**(length-1);hi=128**length-1
  prob=float(binom.cdf(min(n,hi),n,rare)-binom.cdf(lo-1,n,rare)) if lo<=n else 0.
  if prob:header[str(length)]=prob
 def gap_tail(lo):
  if lo>=n:return 0.
  a=n-lo
  return math.exp((lo-1)*math.log1p(-rare))*(rare*a-q*(-math.expm1(a*math.log1p(-rare))))
 adjacency=gap_tail(1);internal=entries-groups
 # Stationary-gap approximation for skipping absolute IDs at group starts.
 # The total retained delta count is exact in expectation; gap-length mix
 # approximates group-start selection as independent of the preceding gap.
 factor=internal/adjacency if adjacency else 0.
 delta={}
 for length in range(1,6):
  lo=1 if length==1 else 128**(length-1)
  count=max(0.,factor*(gap_tail(lo)-gap_tail(128**length)))
  if count:delta[str(length)]=count
 hist={k:header.get(k,0.)+delta.get(k,0.) for k in set(header)|set(delta)}
 varbytes=sum(int(k)*v for k,v in hist.items())
 ids=max(1,math.ceil((n-1).bit_length()/8))
 payload=1+math.ceil(n/8)+varbytes+groups*ids+max(0.,groups-nonempty)+category_bytes
 return dict(entries=entries,groups=groups,category_bytes=category_bytes,header=header,delta=delta,hist=hist,varint_bytes=varbytes,payload_bytes=payload,rare_rate=rare,
  assumption='Binomial rare-count/group/category/header expectations; geometric gap-length counts with explicit stationary group-head approximation. Forced valid one-bit record form, no LD.')

def scenario_data(n,m,sort_counts,*,maf=.2,missing_rate=.001,traits_in_file=32,numeric_characters=19.4):
 if not isinstance(n,int) or not isinstance(m,int) or n<128 or n%128 or m<1:raise ValueError('Scenario requires positive integer M and N multiple of128')
 if sort_counts['subjects']!=n:raise ValueError('Sample sort census mismatch')
 s=onebit_expected_per_marker(n,maf,missing_rate);hist={k:v*m for k,v in s['hist'].items()}
 varbytes=sum(int(k)*v for k,v in hist.items());entries=m*s['entries'];groups=m*s['groups']
 # A two-byte length index is a declared format choice for this N domain,
 # not a guess that an automatic encoder always chooses that index width.
 if s['payload_bytes']>=65536:raise ValueError('Two-byte record-length scenario exceeded')
 payload=m*s['payload_bytes'];index=12+8*math.ceil(m/65536)+math.ceil(m/2)+2*m
 encoded=dict(samples=n,markers=m,file_bytes=payload+index,record_payload_bytes=payload,record_form_counts={'1':m},
  total_difflist_groups=groups,sample_delta_integer_count=entries-groups,
  total_varint_lengths=hist,difflist_header_varint_lengths={k:v*m for k,v in s['header'].items()},sample_delta_varint_lengths={k:v*m for k,v in s['delta'].items()},
  source_work=dict(difflist_entries=entries,varint_bytes=varbytes),native_onebit_tail_high_count=0,
  ld_records_at_chunk_starts=0,additional_base_record_bytes_if_every_chunk_restarts=0,chunk_markers=2048,
  scope='Statistical expected-work scenario; no physical file claimed. '+s['assumption'])
 # Match the synthetic fixture lexical ID/POS conventions exactly. Numeric
 # text length remains a visible data-format parameter, not a runtime fit.
 sid_chars=n+digit_sum_progression(n)
 def table(rows,cols,characters,long=0,header=0):return dict(rows=rows,columns=cols,fields=rows*cols,field_characters=characters,long_fields=long,bytes=characters+rows*cols+header)
 pvchars=5*m+digit_sum_progression(m,1000,10)+digit_sum_progression(m,500000000)
 tables=dict(pvar=table(m,5,pvchars,header=22),psam=table(n,3,2*n+sid_chars,header=13),
  phenotype=table(n,traits_in_file+2,n+sid_chars+n*traits_in_file*numeric_characters,n*traits_in_file),
  covariate=table(n,10,n+sid_chars+n*8*numeric_characters,n*8))
 return dict(samples=n,markers=m,covariates=8,traits_analyzed=1,traits_in_file=traits_in_file,tables=tables,sort_counts=sort_counts,matching_sample_order=True,encoded=encoded,
  scenario=dict(maf=maf,missing_rate=missing_rate,numeric_characters=numeric_characters,pgen_length_bytes=2,ld_records=0,expected_work=True),
  scope='Expected source-work scenario with explicit genotype/format parameters; no observed runtime or interpolation.')
