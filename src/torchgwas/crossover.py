"""Numerical runtime-equality brackets; never a fitted hyperbola."""
from __future__ import annotations
import copy,math
from .crossover_scenario import scenario_data
from .mechanistic_torch import torch_runtime
from .mechanistic_plink import plink_runtime
from .mechanistic_cpu import fastgwa_runtime

def discrete_crossing(evaluate,step=2048,max_markers=8388608):
 if not isinstance(step,int) or step<1 or max_markers<step:raise ValueError('Invalid marker search domain')
 trials={}
 def get(blocks):
  if blocks not in trials:
   torch,other=evaluate(blocks*step)
   if torch is None or other is None or not all(math.isfinite(x) and x>=0 for x in [torch,other]):raise ValueError('Invalid runtime prediction')
   trials[blocks]=dict(M=blocks*step,torch_seconds=torch,competitor_seconds=other,difference_seconds=torch-other)
  return trials[blocks]['difference_seconds']
 limit=max_markers//step
 if get(1)<=0:return dict(status='already_faster_at_search_minimum',lower_M=None,upper_M=step,trials=list(trials.values()))
 low=high=1
 while high<limit:
  high=min(limit,high*2)
  if get(high)<=0:break
  low=high
 if high==low:return dict(status='no_crossing_in_search_domain',lower_M=high*step,upper_M=None,trials=list(trials.values()))
 while high-low>1:
  mid=(low+high)//2
  if get(mid)>0:low=mid
  else:high=mid
 # Guard against a later reversal at two independent larger marker extents.
 checks=[]
 for blocks in sorted({min(limit,2*high),min(limit,4*high)}):
  checks.append(dict(M=blocks*step,torch_faster=get(blocks)<=0))
 return dict(status='bracketed',lower_M=low*step,upper_M=high*step,geometric_midpoint_M=step*math.sqrt(low*high),
  lower=trials[low],upper=trials[high],larger_marker_checks=checks,trials=sorted(trials.values(),key=lambda r:r['M']),
  scope='First crossing found by doubling and integer bisection, resolved to one marker block. Bracket is model equality, not a confidence interval; no claim of uniqueness or global dominance.')

def crossover_at_subjects(n,profiles,sort_counts,*,maf=.2,missing_rate=.001,traits_in_file=32,numeric_characters=19.4,max_markers=8388608):
 p=copy.deepcopy(profiles);p['PLINK2']['missing_rate']=missing_rate;p['PLINK2']['carrier_fraction']=1-(1-maf)**2
 memo={}
 def evaluate(m):
  if m not in memo:
   data=scenario_data(n,m,sort_counts,maf=maf,missing_rate=missing_rate,traits_in_file=traits_in_file,numeric_characters=numeric_characters)
   memo[m]={method:fn(data,p[method])['estimated_seconds'] for method,fn in [('torchGWAS',torch_runtime),('PLINK2',plink_runtime),('fastGWA',fastgwa_runtime)]}
  return memo[m]
 rows=[]
 for other in ['PLINK2','fastGWA']:
  root=discrete_crossing(lambda m:(evaluate(m)['torchGWAS'],evaluate(m)[other]),step=p['torchGWAS']['chunk_markers'],max_markers=max_markers)
  root.update(N=n,competitor=other,scenario=dict(maf=maf,missing_rate=missing_rate,traits_analyzed=1,traits_in_file=traits_in_file,numeric_characters=numeric_characters))
  if root['upper_M'] is not None:root['upper_NM']=n*root['upper_M']
  if root['lower_M'] is not None:root['lower_NM']=n*root['lower_M']
  rows.append(root)
 return rows
