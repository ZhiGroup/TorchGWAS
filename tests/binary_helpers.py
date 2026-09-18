import json
from pathlib import Path
import numpy as np
from scipy import stats
from torchgwas.sumstats import open_binary_sumstats
from torchgwas.sumstats_indexed import open_indexed_sumstats

def binary_rows(out):
 directory=Path(out)/'sumstats';manifest=json.loads((directory/'manifest.json').read_text())
 rows=[];meta={}
 if (directory/'variant_metadata.npz').exists():
  with np.load(directory/'variant_metadata.npz',allow_pickle=False) as f:meta={k:f[k] for k in f.files}
 if manifest['format']=='torchgwas-indexed-sumstats':
  ids=np.load(directory/'variant_ids.npy',allow_pickle=False);_,parts=open_indexed_sumstats(directory)
 else:
  beta,t,logp,manifest=open_binary_sumstats(directory)
  ids=(directory/'variant_ids.txt').read_text().splitlines()
  vi,ti=np.indices(t.shape)
  parts=[{'variant_index':vi.ravel(),'trait_index':ti.ravel(),'beta':np.asarray(beta).ravel(),'t_stat':np.asarray(t).ravel(),'neg_log10_p':np.asarray(logp).ravel()}]
 for part in parts:
  for i,v in enumerate(part['variant_index']):
   row={'marker_id':str(ids[v]),'n':manifest['n_samples']}
   if 'chi2' in part:
    value=float(part['chi2'][i]);logp=-stats.chi2.logsf(value,manifest['df'])/np.log(10)
    row.update(chi2=value,df=manifest['df'],**{'-log10_p':logp})
   else:
    trait=int(part['trait_index'][i]);b=float(part['beta'][i]);t=float(part['t_stat'][i]);logp=float(part['neg_log10_p'][i])
    row.update(trait=manifest['traits'][trait],beta=b,t_stat=t,se=abs(b/t) if t else float('nan'),**{'-log10_p':logp})
   for field,values in meta.items():row[field]=int(values[v]) if field=='position' else str(values[v])
   rows.append(row)
 return rows
