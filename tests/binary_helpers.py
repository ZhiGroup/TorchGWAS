import json
from pathlib import Path
import numpy as np
from scipy import stats
from torchgwas.tails import upper_tail_log10_from_t
from torchgwas.sumstats import open_binary_sumstats, open_binary_df
from torchgwas.sumstats_indexed import open_indexed_sumstats
from torchgwas.variant_source import store_variants

def binary_rows(out):
 directory=Path(out)/'sumstats';manifest=json.loads((directory/'manifest.json').read_text())
 rows=[]
 # Embedded IDs, or the recorded input's own (variant_source.store_variants).
 ids,meta=store_variants(directory)
 if manifest['format']=='torchgwas-indexed-sumstats':
  _,parts=open_indexed_sumstats(directory)
 else:
  beta,t,_logp,manifest=open_binary_sumstats(directory)
  stored_df=np.broadcast_to(open_binary_df(directory),t.shape)
  vi,ti=np.indices(t.shape)
  parts=[{'variant_index':vi.ravel(),'trait_index':ti.ravel(),'beta':np.asarray(beta).ravel(),'t_stat':np.asarray(t).ravel(),
          **({} if _logp is None else {'neg_log10_p':np.asarray(_logp).ravel()})}]
 for part in parts:
  for i,v in enumerate(part['variant_index']):
   row={'marker_id':str(ids[v]),'n':manifest['n_samples']}
   if 'chi2' in part:
    value=float(part['chi2'][i]);logp=-stats.chi2.logsf(value,manifest['df'])/np.log(10)
    row.update(chi2=value,df=manifest['df'],**{'-log10_p':logp})
   else:
    trait=int(part['trait_index'][i]);b=float(part['beta'][i]);t=float(part['t_stat'][i])
    if 'neg_log10_p' in part:logp=float(part['neg_log10_p'][i])
    else:
     # A store written before -log10 P was stored.
     df=(part['df'][i] if 'df' in part else manifest['df']) if manifest['format']=='torchgwas-indexed-sumstats' else stored_df[v,trait]
     logp=float(upper_tail_log10_from_t(t,df))
    row.update(trait=manifest['traits'][trait],beta=b,t_stat=t,se=abs(b/t) if t else float('nan'),**{'-log10_p':logp})
   for field,values in meta.items():row[field]=int(values[v]) if field=='position' else str(values[v])
   rows.append(row)
 return rows
