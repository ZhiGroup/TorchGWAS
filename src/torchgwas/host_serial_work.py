"""Map independently supplied shared-interpreter CPU to actual source calls."""
import math


def attach_host_serial_work(component, primitives, cpu_fraction):
    """Conserve CPU work at existing host submission boundaries.

    Prices describe held-or-unknown CPU service, not elapsed GIL waiting.
    Each source API is charged once, including views with no device kernel.
    A multi-kernel API therefore shares one cumulative submission endpoint.
    """
    if isinstance(cpu_fraction, bool) or not math.isfinite(cpu_fraction) or not 0 < cpu_fraction <= 1:
        raise ValueError('Positive host CPU fraction no greater than one required')
    if 'host_calls' not in component:
        raise ValueError('Source host-call census required for serial CPU work')
    cpu = serial = 0.; ready = {}; calls = []
    for call in component['host_calls']:
        name = call['primitive']
        if name not in primitives:
            raise ValueError('Missing independent serial CPU primitive: ' + name)
        value = primitives[name]
        service = call['cpu_seconds']
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= service:
            raise ValueError('Serial CPU service must lie within total CPU service: ' + name)
        cpu += service; serial += value
        if not math.isclose(call['submit_finish'], cpu / cpu_fraction, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError('Host call offsets must conserve CPU work')
        ready[call['submit_finish']] = serial
        calls.append(dict(call, serial_cpu_seconds=value))
    if not math.isclose(cpu, component['host_dispatch_cpu_seconds'], rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError('Host call census must conserve total CPU work')
    operations = []
    for operation in component['operations']:
        offset = operation['host_submit_finish']
        if offset not in ready:
            raise ValueError('Kernel submission must match a source host-call endpoint')
        operations.append(dict(operation, host_serial_cpu_finish=ready[offset]))
    return dict(component, operations=operations, host_calls=calls,
                host_dispatch_serial_cpu_seconds=serial)
