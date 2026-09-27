"""Bounded significant-pairs proposals using the shared source/resource model."""
from .mechanistic_plan import _integer
from .trait_candidate_space import prepare_trait_candidates
from .significant_host_plan import detailed_significant_host_plan


def prepare_significant_host_candidates(workload,contexts,*,chunks,trait_blocks,output,
        max_candidates=64,max_candidate_tiles=10000,max_census_chunks=1000000):
    """Tile phenotypes and preserve each complete genotype pass and tail."""
    if not isinstance(output,dict) or output.get('block_bytes') is not None:
        raise ValueError('Significant-pairs output requires indexed parts without dense block coalescing')
    return prepare_trait_candidates(workload,contexts,chunks=chunks,trait_blocks=trait_blocks,
        output=output,partition_axes=('trait',),max_candidates=max_candidates,
        max_candidate_tiles=max_candidate_tiles,max_census_chunks=max_census_chunks)


def bounded_significant_host_plan(workload,contexts,*,bounds,joint,output,prices,
                                  significance_threshold=None):
    """Minimize worst supplied resource/occupancy scenario over bounded tiles.

    Every memory decision uses dense retention. Occupancy is explicit and never
    inferred from the significance threshold. Source, profile provenance and
    live resource checks remain the execution bridge's responsibility.
    """
    if not isinstance(bounds,dict) or not {'chunks','trait_blocks'}<=set(bounds) or set(bounds)-{
            'chunks','trait_blocks','max_candidates','max_candidate_tiles','max_census_chunks'}:
        raise ValueError('Significant-pairs bounds require chunks, trait_blocks and optional finite budgets')
    if not isinstance(joint,dict) or 'significance_threshold' in joint:
        raise ValueError('Supply the significance threshold separately from joint resource options')
    limited=dict(bounds)
    limited['max_candidates']=min(_integer('max_candidates',bounds.get('max_candidates',64)),
                                  _integer('max_candidates',joint.get('max_candidates',64)))
    space=prepare_significant_host_candidates(workload,contexts,output=output,**limited)
    result=detailed_significant_host_plan(space['candidates'],prices,
        significance_threshold=significance_threshold,**joint)
    result['search_space']={key:value for key,value in space.items() if key!='candidates'}
    return result
