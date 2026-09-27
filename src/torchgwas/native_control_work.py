"""Native eager scan control calls, separated from payload and statistics work.

Counts follow native_scan.py and PinnedDosageLoader.produce_direct. CPU prices
must come from independent tiny operation probes, never association durations.
Blocking waits are dependencies; their ready-state CPU overhead is priced here.
"""
from .first_principles import positive


def native_control_work(reused=False,*,reduction=None,return_beta=True):
    if reduction not in (None,'jagwas','device_significant'):
        raise ValueError('Unsupported native result reduction')
    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced control layout')
    arrays=5 if reduction=='jagwas' else 3+int(return_beta)
    calls = {
        'decode_submit': {'queue_get': 1, 'tensor_slice': 1, 'tensor_numpy': 1,
                          'executor_submit': 1},
        'publish': {'future_result': 1, 'tensor_slice': 1, 'queue_put': 1},
        'transfer_submit': {'queue_get': 1, 'future_result': int(reused),
                            'stream_wait_event': 1 + int(reused),
                            'stream_context': 1, 'tensor_slice': 2,
                            'copy_h2d': 1, 'event_record': 1},
        'result_submit': {'event_record': 2, 'stream_context': 1,
                          'stream_wait_event': 1, 'tensor_slice': arrays,
                          'copy_d2h': arrays, 'tensor_record_stream': arrays,
                          'executor_submit': 2},
        'release': {'event_synchronize_ready': 1, 'queue_put': 1},
        'finish': {'event_synchronize_ready': 1, 'tensor_slice': arrays, 'tensor_numpy': arrays},
        'resolve': {'future_result': 1},
    }
    if reduction=='device_significant':
        # The main iterator records compute completion and submits only the
        # independent input-release worker. Blocking status/selection APIs
        # have their own services; no dense result-ring future is submitted.
        calls['result_submit']={'event_record':1,'executor_submit':1}
        calls['finish']={};calls['resolve']={}
    return calls


def native_control_service(primitives, reused=False,*,reduction=None,return_beta=True):
    work = native_control_work(reused,reduction=reduction,return_beta=return_beta)
    missing = sorted({name for counts in work.values() for name, count in counts.items()
                      if count and name not in primitives})
    if missing:
        raise ValueError('Missing native control primitives: ' + ', '.join(missing))
    seconds = {phase: sum(count * positive(name, primitives[name], True)
                          for name, count in counts.items() if count)
               for phase, counts in work.items()}
    return dict(cpu_seconds=seconds, calls=work,
                scope='Native int8 eager direct-fill pipeline without scan profiling. '
                      'Payload copy time, statistics, NumPy finish bookkeeping and writer '
                      'are counted separately. Python loop, decoder wrapper, allocator '
                      'and scheduling/GIL latency are not fully priced.')
