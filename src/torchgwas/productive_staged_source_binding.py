"""Audit exact source prices used by a first-chunk staged layout screen.

The caller must first validate the active detailed profile and original
measurement ages. This compares candidate values with that profile and checks
which of those values have declared immutable price targets. It does not bind
GPU, selector, archive, output or finite-completion service.
"""
from .calibration_cache import _digest


_SOURCE_FIELDS=('decode_units','cpu_fraction','depth','decode_workers',
                'cpu_available_cores','shared_dram_bytes_per_second',
                'read_bytes_per_second')
_CAPACITIES=('cpu','dram','input')


def _leaves(value,prefix):
    if isinstance(value,dict):
        for key,item in value.items():
            yield from _leaves(item,prefix+(key,))
    elif isinstance(value,list):
        for index,item in enumerate(value):
            yield from _leaves(item,prefix+(index,))
    else:
        yield prefix


def _covered(path,targets):
    return any(path[:len(target)]==target for target in targets)


def audit_staged_source_price_binding(profile,context_name,candidates,
                                      shared_source_capacities,*,profile_sha256=None,
                                      max_reported_paths=24):
    """Match source coefficients and report unbound measurement leaves.

    A shared capacity may be below its declared context rate as an explicit
    availability scenario. No candidate may increase it. An exact candidate
    source coefficient must match its device's active profile, including the
    optional buffered-read CPU price when present.
    """
    if (not isinstance(profile,dict) or not isinstance(profile.get('contexts'),list)
            or not isinstance(context_name,str) or not context_name or
            not isinstance(candidates,list) or not 0<len(candidates)<=4 or
            type(max_reported_paths) is not int or not 1<=max_reported_paths<=128):
        raise ValueError('Bounded active staged source price audit required')
    contexts=[(index,row) for index,row in enumerate(profile['contexts'])
              if isinstance(row,dict) and row.get('name')==context_name]
    if len(contexts)!=1:
        raise ValueError('One active staged price context required')
    context_index,context=contexts[0]
    base=context.get('shared_capacities')
    if (not isinstance(base,dict) or
            not isinstance(shared_source_capacities,dict) or
            set(shared_source_capacities)!=set(_CAPACITIES)):
        raise ValueError('Explicit matched staged source capacities required')
    checked=set()
    for name in _CAPACITIES:
        proposed=shared_source_capacities[name]
        available=base.get(name)
        if (isinstance(proposed,bool) or not isinstance(proposed,(int,float)) or
                isinstance(available,bool) or not isinstance(available,(int,float)) or
                not 0<proposed<=available):
            raise ValueError('Staged source capacity exceeds active context: '+name)
        checked.add((context_index,'shared_capacities',name))
    profiles=context.get('profiles')
    if not isinstance(profiles,dict):
        raise ValueError('Active staged device profiles required')
    for candidate in candidates:
        if not isinstance(candidate,dict) or not isinstance(candidate.get('partitions'),list):
            raise ValueError('Explicit staged candidate partitions required')
        parts=candidate['partitions']
        sources=candidate.get('source_profiles')
        if (not 0<len(parts)<=16 or not isinstance(sources,dict) or
                set(sources)!={part.get('id') for part in parts if isinstance(part,dict)}):
            raise ValueError('One staged source price per bounded partition required')
        for part in parts:
            device=part.get('device');part_id=part.get('id')
            source=sources[part_id]
            active=profiles.get(device)
            if not isinstance(source,dict) or not isinstance(active,dict):
                raise ValueError('Staged source device lacks an active profile')
            if (active.get('input_read_cpu_prices') is not None and
                    source.get('input_read_cpu_prices') is None):
                raise ValueError('Staged source omits active buffered-read CPU price: '
                                 +str(device))
            fields=list(_SOURCE_FIELDS)
            if 'input_read_cpu_prices' in source:fields.append('input_read_cpu_prices')
            for field in fields:
                if field not in source or field not in active or _digest(source[field])!=_digest(active[field]):
                    raise ValueError('Staged source price differs from active profile: '
                                     +str(device)+'.'+field)
                checked.update(_leaves(source[field],
                    (context_index,'profiles',device,field)))
    targets=[]
    for binding in profile.get('price_bindings',[]):
        for target in binding.get('targets',[]):
            path=target.get('context_path')
            if isinstance(path,list):targets.append(tuple(path))
    missing=sorted((path for path in checked if not _covered(path,targets)),
                   key=lambda path:repr(path))
    return dict(kind='torchgwas.staged_source_price_binding.v1',
        status=('declared_source_prices_verified' if not missing else
                'source_prices_match_but_unbound'),
        context=context_name,profile_sha256=(_digest(profile) if
            profile_sha256 is None else profile_sha256),
        matched_price_leaves=len(checked),declared_price_leaves=len(checked)-len(missing),
        unbound_price_leaves=len(missing),
        unbound_paths=[list(path) for path in missing[:max_reported_paths]],
        truncated_unbound_paths=max(0,len(missing)-max_reported_paths),
        candidate_count=len(candidates),
        scope='Exact candidate source-price values and declared target coverage in the previously validated active profile. This does not re-read measurement records or certify source price quality, GPU/transfer/output prices, loaded availability or completion time.')
