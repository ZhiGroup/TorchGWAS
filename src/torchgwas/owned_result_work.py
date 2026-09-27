"""Owned native-result allocations and explicit copy-memory-state scenarios."""
import math


def owned_result_work(markers,traits,*,reduction=None,return_beta=True):
    if any(isinstance(v,bool) or not isinstance(v,int) or v<1 for v in (markers,traits)):
        raise ValueError('Positive integer result dimensions required')
    if reduction not in (None,'jagwas'):
        raise ValueError('Unsupported native result reduction')
    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced result layout')
    if reduction=='jagwas':
        arrays=dict(beta=4*markers,t=4*markers,trait_index=4*markers,status=markers,df=4*markers)
        return dict(array_bytes=arrays,array_elements={name:markers for name in arrays},
            allocation_calls=5,copy_bytes=sum(arrays.values()),reduction='jagwas',
            consumer_arrays=['beta','t','trait_index'])
    arrays={'beta':4*markers*traits,'t':4*markers*traits,'status':markers,'df':4*markers}
    if not return_beta:arrays.pop('beta')
    return dict(array_bytes=arrays,array_elements={name:markers*traits if name in ('beta','t') else markers for name in arrays},allocation_calls=len(arrays),copy_bytes=sum(arrays.values()))


def owned_result_copy_service(work,scenario,*,baseline_copy_bytes=None):
    """Additional bulk service beyond the existing tiny finish primitive.

    Fresh-page rates include touching pages, but not creating/freeing allocations.
    Memory-state fractions are explicit assumptions, never inferred from scans.
    """
    if baseline_copy_bytes is None:
        if work.get('reduction') is not None or 'beta' not in work['array_bytes']:
            raise ValueError('Reduced or t-only results require an explicit matching finish baseline')
        baseline_copy_bytes=13*32
    if type(baseline_copy_bytes) is not int or baseline_copy_bytes<0:
        raise ValueError('Nonnegative integer finish baseline required')
    for key in ('resident_cpu_seconds_per_byte','fresh_cpu_seconds_per_byte','fresh_fraction'):
        value=scenario[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
            raise ValueError('Invalid owned-result scenario: '+key)
    fraction=scenario['fresh_fraction']
    if fraction>1:raise ValueError('Fresh fraction must be in [0,1]')
    rate=(1-fraction)*scenario['resident_cpu_seconds_per_byte']+fraction*scenario['fresh_cpu_seconds_per_byte']
    additional=max(0,work['copy_bytes']-baseline_copy_bytes)
    return dict(additional_copy_cpu_seconds=additional*rate,additional_copy_bytes=additional,
                fresh_fraction=fraction,prediction_complete=False,
                unpriced_terms=['owned-result allocation/free service and actual page-reuse state'])


def owned_result_allocation_policy(work, policy=None):
    """Source-derived advice exposure, not a compaction latency prediction.

    Only NumPy 2.2.6's allocator threshold has been source-verified here.
    Advice eligibility is per array, never the aggregate chunk output size.
    Kernel state and free-page topology do not give a finite stall bound.
    """
    policy = {} if policy is None else policy
    enabled = policy.get('numpy_madvise_hugepage')
    if enabled is not None and not isinstance(enabled, bool):
        raise ValueError('numpy_madvise_hugepage must be boolean or unknown')
    version = policy.get('numpy_version')
    kernel = policy.get('linux_thp_enabled')
    if kernel not in (None, 'always', 'madvise', 'never'):
        raise ValueError('Supply the selected Linux THP policy, not the whole sysfs line')
    known_threshold = version == '2.2.6'
    threshold = (1 << 22) if known_threshold else None
    eligible = ({name: size for name, size in work['array_bytes'].items()
                 if size >= threshold} if known_threshold else None)
    if enabled is False or kernel == 'never' or eligible == {}:
        exposure = False
    elif enabled is True and kernel in ('always', 'madvise') and eligible:
        exposure = True
    else:
        exposure = None
    unpriced = []
    if exposure is not False:
        unpriced.append('NumPy-owned result hugepage advice, allocation compaction and page-reuse state')
    return dict(numpy_version=version, numpy_madvise_hugepage=enabled,
                linux_thp_enabled=kernel, linux_thp_defrag=policy.get('linux_thp_defrag'),
                verified_advice_threshold_bytes=threshold,
                advice_eligible_array_bytes=eligible,
                advice_compaction_exposure=exposure,
                compaction_latency_bound_seconds=None,
                prediction_complete=False, unpriced_terms=unpriced,
                scope='NumPy advice exposure only; kernel always-mode promotion, other allocators, NUMA migration and system memory pressure remain outside this diagnostic.',
                source='https://raw.githubusercontent.com/numpy/numpy/v2.2.6/numpy/_core/src/multiarray/alloc.c' if known_threshold else None)

def owned_result_copy_threading(work, policy=None):
    """NumPy 2.2.6 plain numeric copies release the GIL above 500 elements.

    Allocation, dispatch and finalization are outside this loop. Thread support
    is an explicit environment input; unknown versions/builds remain unknown.
    """
    policy={} if policy is None else policy
    support=policy.get('numpy_allow_threads')
    if support is not None and not isinstance(support,bool):
        raise ValueError('numpy_allow_threads must be boolean or unknown')
    known=policy.get('numpy_version')=='2.2.6' and support is not None
    released={name:(bool(support and count>500) if known else None)
              for name,count in work['array_elements'].items()}
    return dict(gil_released_by_array=released,verified_element_threshold=500 if known else None,
                unpriced_terms=[] if known else ['NumPy numeric-copy GIL policy requires a verified version and threading support'],
                scope='Plain numeric copy loop only; allocation and Python control remain separate.',
                source='https://raw.githubusercontent.com/numpy/numpy/v2.2.6/numpy/_core/include/numpy/ndarraytypes.h' if known else None)
