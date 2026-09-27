"""Bind explicit calculator values to original immutable component evidence.

This verifies only the declared fields. It cannot certify an unlisted price,
measurement protocol, prediction accuracy, or current resource availability.
"""
from copy import deepcopy
from pathlib import Path
import time

from .calibration_cache import _age, _digest, _json, read_calibration_record,inspect_calibration_record


KINDS=frozenset({'cpu_capacity','gpu_capacity','transfer_capacity','storage_capacity'})
MAX_BINDINGS=128
MAX_TARGETS=512


def _path(value):
    if (not isinstance(value,list) or not value or len(value)>16 or
            any(not ((type(v) is str and v) or (type(v) is int and v>=0)) for v in value)):
        raise ValueError('Explicit bounded JSON value path required')
    return tuple(value)


def _at(root,path):
    for key in path:
        if isinstance(root,dict) and type(key) is str and key in root:root=root[key]
        elif isinstance(root,list) and type(key) is int and key<len(root):root=root[key]
        else:raise ValueError('Price binding value path does not exist')
    return root


def validate_price_bindings(profile, *, max_age_seconds=None, deferred_bindings=()):
    """Check exact profile values, artifact identities, dependencies and age.

    Targets are paths relative to profile['contexts']; value paths are relative
    to each record's value. The caller may shorten all declared lifetimes, never
    extend them. No lookup chooses a newer record behind the immutable profile.
    A final age check includes time spent reading the rest of the evidence.
    Explicit deferred indexes inspect identity without accepting stale prices.
    Their presence returns pending_refresh, never declared_targets_verified.
    """
    from .detailed_calibration import sha256_file
    requested=None if max_age_seconds is None else _age(max_age_seconds)
    bindings=profile.get('price_bindings')
    if (not isinstance(deferred_bindings,(list,tuple)) or any(type(i) is not int or i<0 for i in deferred_bindings)
            or len(set(deferred_bindings))!=len(deferred_bindings)
            or (deferred_bindings and (not isinstance(bindings,list) or max(deferred_bindings)>=len(bindings)))):
        raise ValueError('Explicit unique existing deferred binding indexes required')
    if bindings is None:
        return dict(status='undeclared',bindings=[],verified_targets=0,
            scope='No declared price-to-record links; artifact identity does not establish empirical freshness.')
    if not isinstance(bindings,list) or not 0<len(bindings)<=MAX_BINDINGS:
        raise ValueError('Nonempty bounded price bindings required')
    loaded={};seen=[];rows=[]
    for index,binding in enumerate(bindings):
        fields={'artifact','kind','name','dependencies','targets','max_age_seconds'}
        if not isinstance(binding,dict) or set(binding)!=fields:
            raise ValueError('Explicit price binding fields required')
        if binding['kind'] not in KINDS:
            raise ValueError('Prices require independent component evidence, not loaded spans or live observations')
        deps=binding['dependencies']
        if (not isinstance(deps,dict) or deps.get('source_sha256')!=profile['source_sha256'] or
                deps.get('execution_context')!=profile['execution_context'] or
                not isinstance(deps.get('measurement_protocol'),dict) or not deps['measurement_protocol']):
            raise ValueError('Price evidence source, execution context or measurement protocol differs')
        age=None if binding['max_age_seconds'] is None else _age(binding['max_age_seconds'])
        if requested is not None:age=requested if age is None else min(age,requested)
        artifact=Path(binding['artifact']).expanduser().resolve(strict=True)
        if str(artifact)!=binding['artifact']:
            raise ValueError('Price artifact must use its canonical absolute path')
        digest=profile['component_artifacts'].get(str(artifact))
        if digest is None:raise ValueError('Price record is absent from profile artifacts')
        deferred=index in deferred_bindings
        key=(str(artifact),binding['kind'],binding['name'],_json(deps),deferred)
        if key not in loaded:
            if sha256_file(artifact)!=digest:raise ValueError('Price component artifact changed')
            reader=inspect_calibration_record if deferred else read_calibration_record
            loaded[key]=reader(artifact,kind=binding['kind'],name=binding['name'],dependencies=deps)
        evidence=loaded[key];record=evidence['record']
        lifetime=record['max_age_seconds'] if age is None else min(age,record['max_age_seconds'])
        targets=binding['targets']
        if not isinstance(targets,list) or not targets or len(seen)+len(targets)>MAX_TARGETS:
            raise ValueError('Nonempty bounded price targets required')
        locations=[]
        for target in targets:
            if not isinstance(target,dict) or set(target)!={'context_path','value_path'}:
                raise ValueError('Exact context and measurement value paths required')
            path=_path(target['context_path']);value_path=_path(target['value_path'])
            if type(path[0]) is not int:raise ValueError('Context path must start at an explicit context index')
            if any(path[:len(old)]==old or old[:len(path)]==path for old in seen):
                raise ValueError('Overlapping or duplicate price targets')
            seen.append(path)
            if _json(_at(profile['contexts'],path))!=_json(_at(record['value'],value_path)):
                raise ValueError('Calculator price differs from its immutable measurement')
            locations.append(deepcopy(target))
        rows.append(dict(artifact=str(artifact),artifact_sha256=digest,record_sha256=evidence['record_sha256'],
            kind=binding['kind'],name=binding['name'],targets=locations,
            observed_unix_seconds=record['observed_unix_seconds'],max_age_seconds=lifetime,
            created_unix_seconds=record['created_unix_seconds'],
            expires_unix_seconds=record['observed_unix_seconds']+lifetime,
            provenance=deepcopy(record['provenance'])))
    now=time.time()
    for index,row in enumerate(rows):
        row['age_seconds']=now-row['observed_unix_seconds']
        if now<row['created_unix_seconds'] or row['age_seconds']<0 or (index not in deferred_bindings and now>=row['expires_unix_seconds']):
            raise ValueError('Bound calculator price expired during validation')
        if deferred_bindings:row.update(pending_refresh=index in deferred_bindings,fresh=now<row['expires_unix_seconds'])
    if deferred_bindings:
        return dict(status='pending_refresh',bindings=rows,verified_targets=sum(len(r['targets']) for r in rows if not r['pending_refresh']),
            deferred_bindings=list(deferred_bindings),
            scope='Identity checked; named coefficients must complete refresh before any timing decision.')
    return dict(status='declared_targets_verified',bindings=rows,verified_targets=len(seen),
        scope='Exact declared values and original measurement ages only. Unlisted prices, measurement validity, live load and prediction accuracy are not certified.')


def canonical_price_bindings(bindings):
    """Copy declaration paths before the containing profile is published."""
    if not isinstance(bindings,list) or not 0<len(bindings)<=MAX_BINDINGS:
        raise ValueError('Nonempty bounded price bindings required')
    result=deepcopy(bindings)
    for row in result:
        if not isinstance(row,dict) or not isinstance(row.get('artifact'),(str,Path)):
            raise ValueError('Explicit price artifact required')
        row['artifact']=str(Path(row['artifact']).expanduser().resolve(strict=True))
    return result


def validate_comparison_prices(profile,comparisons):
    """Bind productive comparisons to these profiles, then recheck price age.

    Shapes may change chunk size and untimed kernel geometry. All other profile
    values must match one of the explicitly bound device profiles. Source/live
    context validation still belongs to admission, and unlisted prices remain
    unqualified even when the declared subset passes.
    """
    excluded=('chunk_markers','kernel_geometry','joint_kernel_geometry')
    for binding in profile.get('price_bindings',[]):
        for target in binding['targets']:
            path=_path(target['context_path'])
            if len(path)<4 or path[1]!='profiles' or path[3] in excluded:
                raise ValueError('Productive evidence must bind fixed per-device profile parameters')
    fingerprints={}
    for context in profile['contexts']:
        for device,values in context['profiles'].items():
            fixed={key:value for key,value in values.items()
                if key not in excluded}
            fingerprints.setdefault(device,set()).add(_digest(fixed))
    digest=_digest(profile)
    if not isinstance(comparisons,list) or not comparisons:
        raise ValueError('Explicit price-bound comparisons required')
    for comparison in comparisons:
        contract=comparison['comparison_contract'];identity=contract.get('model_identity')
        if (not isinstance(identity,dict) or identity.get('price_profile_sha256')!=digest or
                identity.get('source_sha256')!=profile['source_sha256']):
            raise ValueError('Comparison differs from its declared price profile')
        for name in ('baseline','candidate'):
            for window in contract[name]['windows']:
                if window.get('chunk_invariant_profile_sha256') not in fingerprints.get(window['device'],set()):
                    raise ValueError('Comparison device prices differ from bound profile values')
    checked=validate_price_bindings(profile)
    if checked['status']!='declared_targets_verified':
        raise ValueError('Productive price profile has no declared measurement evidence')
    return dict(checked,price_profile_sha256=digest)
