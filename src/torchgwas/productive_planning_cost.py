"""Advisory cost history for bounded productive calculator callbacks.

Prior observations may only raise a later job's expected tuning cost. They
never provide association service prices or approve a configuration change.
"""
from copy import deepcopy
import math
import time

from .calibration_cache import CalibrationParameterCache,_age


NAME='productive_planning_cost.v1'


def validate_planning_cost_history(config):
    if (not isinstance(config,dict) or set(config)!=
            {'cache_dir','max_age_seconds','publication_seconds'}):
        raise ValueError('Planning cost history requires a directory, age and publication forecast')
    if not isinstance(config['cache_dir'],str) or not config['cache_dir'].strip():
        raise ValueError('Nonempty planning cost cache directory required')
    _age(config['max_age_seconds'])
    value=config['publication_seconds']
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError('Finite nonnegative planning cost publication forecast required')


def _row(row):
    if not isinstance(row,dict) or set(row)!={'cpu_seconds','wall_seconds','observed_unix_seconds'}:
        raise ValueError('Complete planning cost observation required')
    for key,value in row.items():
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError('Positive finite planning cost '+key+' required')
    return deepcopy(row)


class ProductivePlanningCostHistory:
    def __init__(self,config,dependencies,*,max_steps):
        validate_planning_cost_history(config)
        if type(max_steps) is not int or not 1<=max_steps<=32:
            raise ValueError('Bounded planning-cost step count required')
        self.config=deepcopy(config);self.dependencies=deepcopy(dependencies)
        self.max_steps=max_steps;self.cache=CalibrationParameterCache(config['cache_dir'])
        self.rows=[];self.prior=None;self.lookup_report=None;self.publication=None
        self.lookup_wall_seconds=0.;self.preparation_wall_seconds=0.
        self.publication_wall_seconds=0.;self._looked_up=False

    def lookup(self):
        if self._looked_up:return self.lookup_report
        self._looked_up=True;began=time.perf_counter()
        try:
            found=self.cache.lookup('stage_observations',NAME,
                dependencies=self.dependencies,max_age_seconds=self.config['max_age_seconds'])
            if found['hit']:
                value=found['record']['value']
                if (not isinstance(value,dict) or set(value)!={'schema','rows'} or
                        value['schema']!=NAME or not isinstance(value['rows'],list) or
                        not 1<=len(value['rows'])<=self.max_steps or
                        found['record'].get('observed_unix_seconds')!=min(
                            _row(row)['observed_unix_seconds'] for row in value['rows'])):
                    raise ValueError('Invalid bound planner-cost history')
                self.prior=dict(record_sha256=found['record_sha256'],
                    observed_unix_seconds=found['record']['observed_unix_seconds'],
                    age_seconds=found['age_seconds'],rows=deepcopy(value['rows']))
            self.lookup_report=dict(hit=self.prior is not None,reason=found.get('reason'),
                record_sha256=None if self.prior is None else self.prior['record_sha256'])
        except (OSError,ValueError,KeyError,TypeError) as error:
            self.lookup_report=dict(hit=False,reason='invalid_or_unavailable_history',
                                    error=type(error).__name__)
        finally:self.lookup_wall_seconds=time.perf_counter()-began
        return deepcopy(self.lookup_report)

    def costs(self,declared):
        self.lookup()
        result=dict(declared)
        prior=[] if self.prior is None else self.prior['rows']
        rows=prior+self.rows
        if rows:
            result['expected_cpu_seconds']=max(result['expected_cpu_seconds'],
                                               *(row['cpu_seconds'] for row in rows))
            result['expected_wall_seconds']=max(result['expected_wall_seconds'],
                                                *(row['wall_seconds'] for row in rows))
        result['publication_seconds']+=self.config['publication_seconds']
        return result

    def charge_preparation(self,seconds):
        if isinstance(seconds,bool) or not isinstance(seconds,(int,float)) or not math.isfinite(seconds) or seconds<0:
            raise ValueError('Finite nonnegative planning history preparation time required')
        self.preparation_wall_seconds+=seconds

    def observe(self,result,*,observed_unix_seconds):
        if not result.get('evaluated') or 'wall_seconds' not in result or len(self.rows)>=self.max_steps:
            return
        self.rows.append(_row(dict(cpu_seconds=result['cpu_seconds'],
            wall_seconds=result['wall_seconds'],observed_unix_seconds=observed_unix_seconds)))

    def finish(self,*,successful):
        if type(successful) is not bool:raise ValueError('Boolean planning-cost completion required')
        if successful and self.rows:
            began=time.perf_counter()
            try:
                record=self.cache.store('stage_observations',NAME,
                    dict(schema=NAME,rows=deepcopy(self.rows)),dependencies=self.dependencies,
                    provenance=dict(scope='Observed productive planner CPU/wall costs only; advisory conservative forecast, never a resource price or switch approval.'),
                    max_age_seconds=self.config['max_age_seconds'],
                    observed_unix_seconds=min(row['observed_unix_seconds'] for row in self.rows))
                self.publication=dict(status='stored',record_sha256=record['record_sha256'],path=record['path'])
            except (OSError,ValueError,TypeError) as error:
                self.publication=dict(status='cache_error',error=type(error).__name__)
            finally:self.publication_wall_seconds=time.perf_counter()-began
        else:self.publication=dict(status='not_stored')
        return self.snapshot()

    def snapshot(self):
        return dict(lookup=deepcopy(self.lookup_report),prior=deepcopy(self.prior),
            observations=deepcopy(self.rows),lookup_wall_seconds=self.lookup_wall_seconds,
            preparation_wall_seconds=self.preparation_wall_seconds,
            publication=deepcopy(self.publication),publication_wall_seconds=self.publication_wall_seconds,
            declared_publication_seconds=self.config['publication_seconds'],
            scope='Source/config/profile-bound previous planner cost may only raise expected CPU/wall cost. Preparation including lookup and forecast publication are charged; actual publication is audited. A hit never renews observation age or authorizes a switch.')
