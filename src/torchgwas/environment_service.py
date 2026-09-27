"""Independent runtime-environment event services, never association residuals.

The import/thread-pool event is a critical-path elapsed service at a stated
filesystem and scheduler condition. Aggregate process CPU is not elapsed time.
It must not be divided by a CPU share a second time: its observed waiting is
already included. Changing that environment requires new explicit event inputs.
"""
import math,statistics

def environment_profile(observations, *, condition, repeats=None):
 rows=observations['rows']
 if repeats is not None:rows=[r for r in rows if r['repeat'] in repeats]
 if not rows:raise ValueError('No independent environment observations')
 names=[p['name'] for p in rows[0]['phases']]
 if len(names)!=len(set(names)):raise ValueError('Duplicate environment event')
 if any([p['name'] for p in row['phases']]!=names for row in rows):raise ValueError('Environment event graph changed')
 events=[]
 for name in names:
  values=[next(p['wall_seconds'] for p in row['phases'] if p['name']==name) for row in rows]
  events.append(dict(name=name,count=1,seconds_per_event=statistics.median(values),observed_min_seconds=min(values),observed_max_seconds=max(values)))
 values=[row['uninstrumented_process_entry_exit_seconds'] for row in rows]
 events.append(dict(name='process_entry_exit',count=1,seconds_per_event=statistics.median(values),observed_min_seconds=min(values),observed_max_seconds=max(values)))
 result=dict(events=events,condition=condition,repeat_ids=[row['repeat'] for row in rows],
  includes_scheduler_and_filesystem_wait=True,scaling='Externally supplied environment event service at the stated load, not automatically rescaled by scan CPU/GPU/storage availability.',
  origin='Independent empty-process imports and thread-pool setup; no N, M, K, genotype or association call.')
 environment_service(result)
 return result

def environment_service(profile):
 if not isinstance(profile,dict) or not profile.get('events'):raise ValueError('Explicit independent environment event services are required; legacy aggregate process CPU pricing is withdrawn')
 events=[]
 for event in profile['events']:
  count=event['count'];rate=event['seconds_per_event']
  if not all(isinstance(x,(int,float)) and math.isfinite(x) and x>=0 for x in [count,rate]):raise ValueError('Invalid environment work/service')
  events.append(dict(name=event['name'],count=count,seconds_per_event=rate,seconds=count*rate))
 return dict(seconds=sum(e['seconds'] for e in events),events=events,condition=profile.get('condition'),
  scope='Independent environment-input service, not hardware-derived instruction timing or an association-runtime fit. Waiting is included; scan load overrides do not change this service automatically.')
