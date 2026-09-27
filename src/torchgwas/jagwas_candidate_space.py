"""Bounded full-panel JAGWAS proposals and source-derived preparation.

This connects finite chunk/device contexts to the shared detailed calculator.
It does not fit GWAS timings or certify calibration, ranking, or live admission.
"""
from __future__ import annotations

from .jagwas_candidate import detailed_jagwas_plan
from .jagwas_preparation import build_jagwas_preparation
from .mechanistic_plan import _integer
from .trait_candidate_space import prepare_trait_candidates


def prepare_jagwas_candidates(workload, contexts, *, chunks, output,
                              max_candidates=64, max_candidate_tiles=10000,
                              max_census_chunks=1000000):
    """Retain every phenotype on each GPU; vary only chunk and device context."""
    return prepare_trait_candidates(workload, contexts, chunks=chunks,
        trait_blocks=[workload['traits']], output=output, partition_axes=('variant',),
        reduction='jagwas', max_candidates=max_candidates,
        max_candidate_tiles=max_candidate_tiles, max_census_chunks=max_census_chunks)


def bounded_jagwas_plan(workload, contexts, *, bounds, joint, output, prices,
                       preparation_services):
    """Construct and price a finite JAGWAS space, including setup per scenario.

    preparation_services[context_name][host_scenario_name] declares the fixed
    common CPU/finalization steps and per-device scalar/tensor factor arithmetic.
    Factor/design/pinned work is derived from the candidate, using that host
    scenario's serialization fraction. Services are explicit independent inputs,
    not association durations; no empty preparation or dense-writer fallback is
    supplied. Graphs are built only after candidate memory/reader admission.
    """
    if not isinstance(bounds,dict) or 'chunks' not in bounds or set(bounds)-{
            'chunks','max_candidates','max_candidate_tiles','max_census_chunks'}:
        raise ValueError('JAGWAS bounds require chunks and optional finite search budgets; phenotype partitioning is excluded')
    if not isinstance(joint,dict) or any(key in joint for key in ('preparations','preparation_factory')):
        raise ValueError('JAGWAS joint options cannot supply a preparation graph or factory')
    hosts=joint.get('host_scenarios')
    if not isinstance(hosts,dict) or not hosts:
        raise ValueError('Explicit host scenarios required before candidate construction')
    if not isinstance(contexts,(list,tuple)) or not contexts:
        raise ValueError('Nonempty measured contexts required')
    names=[context['name'] for context in contexts]
    if len(set(names))!=len(names):
        raise ValueError('Unique measured context names required')
    if not isinstance(preparation_services,dict) or set(preparation_services)!=set(names):
        raise ValueError('Exactly one preparation-service map per measured context required')
    for context in contexts:
        services=preparation_services[context['name']]
        if not isinstance(services,dict) or set(services)!=set(hosts):
            raise ValueError('Preparation services must cover every host scenario exactly')
        for spec in services.values():
            if not isinstance(spec,dict) or set(spec)!={'library_arithmetic','shared_cpu_steps','finalize'}:
                raise ValueError('Explicit factor arithmetic, common CPU and finalization services required')
            arithmetic=spec['library_arithmetic']
            if (not isinstance(arithmetic,dict) or set(arithmetic)!=set(context['devices'])
                    or any(value not in ('scalar','tensor') for value in arithmetic.values())):
                raise ValueError('Explicit scalar/tensor factor arithmetic per context device required')
            if any(not isinstance(spec[key],list) or not spec[key] for key in ('shared_cpu_steps','finalize')):
                raise ValueError('Nonempty common CPU and finalization services required')
    bounded=dict(bounds)
    search_limit=_integer('max_candidates',joint.get('max_candidates',64))
    proposal_limit=_integer('max_candidates',bounded.get('max_candidates',64))
    bounded['max_candidates']=min(search_limit,proposal_limit)
    space=prepare_jagwas_candidates(workload,contexts,output=output,**bounded)
    assignments={row['candidate_index']:row['context'] for row in space['assignments']}
    def preparation(index,candidate,host_name,host):
        spec=preparation_services[assignments[index]][host_name]
        return build_jagwas_preparation(candidate,
            host_serial_fraction=host['host_serial_fraction'],**spec)['preparation']
    result=detailed_jagwas_plan(space['candidates'],prices,
        preparation_factory=preparation,**joint)
    result['search_space']={key:value for key,value in space.items() if key!='candidates'}
    return result
