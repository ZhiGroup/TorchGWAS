"""Independent resident-copy CPU samples collected after useful job output.

This coefficient is not DRAM bandwidth or a loaded GWAS stage capacity. All
other prices in the integration audit remain synthetic.
"""
from copy import deepcopy
import time

import numpy as np

from torchgwas.cpu_service_refresh import CpuServiceRefresh
from torchgwas.detailed_calibration import bind_detailed_profile, sha256_file
from torchgwas.price_binding import validate_price_bindings


class CopyPriceProbe:
    def __init__(self, directory, contexts, source, execution, *, max_age_seconds=300.):
        self.contexts, self.source, self.execution = deepcopy(contexts), deepcopy(source), deepcopy(execution)
        self.protocol = dict(operation='numpy.copyto', dtype='uint8', elements=8<<20, repeats=7, loops=4,
            source='arange', destination='pre-touched', clock='thread_time',
            aggregation='median_repeat_cpu_seconds_per_byte', implementation_sha256=sha256_file(__file__))
        self.dependencies = dict(source_sha256=source, execution_context=execution, measurement_protocol=self.protocol)
        self.refresh = CpuServiceRefresh(directory, 'resident_numpy_copy.v2', dependencies=self.dependencies,
            work_units=self.protocol['elements']*self.protocol['loops'], max_age_seconds=max_age_seconds)
        self.values = self.destination = None
        self.max_age_seconds = max_age_seconds

    def sample(self):
        if self.values is None:
            self.values = np.arange(self.protocol['elements'], dtype=np.uint8)
            self.destination = np.empty_like(self.values)
            np.copyto(self.destination, self.values)
        observed = time.time(); wall = time.perf_counter(); cpu = time.thread_time()
        for _ in range(self.protocol['loops']):
            np.copyto(self.destination, self.values)
        cpu = time.thread_time()-cpu; wall = time.perf_counter()-wall
        assert np.array_equal(self.destination, self.values)
        return dict(cpu_seconds=cpu, wall_seconds=wall, work_units=self.values.nbytes*self.protocol['loops'],
                    observed_unix_seconds=observed)

    def close(self):
        self.values = self.destination = None

    def advance(self):
        began = time.perf_counter()
        try:
            state = self.refresh.advance(self.sample)
            profile = None
            if state['state'] == 'ready':
                result = state['result']; artifact = result['path']; coefficient = result['record']['value']['cpu_seconds_per_unit']
                contexts = deepcopy(self.contexts); targets = []
                for index, context in enumerate(contexts):
                    for device, entry in context['profiles'].items():
                        entry['owned_result_copy_scenario']['resident_cpu_seconds_per_byte'] = coefficient
                        targets.append(dict(context_path=[index,'profiles',device,'owned_result_copy_scenario','resident_cpu_seconds_per_byte'],
                                            value_path=['cpu_seconds_per_unit']))
                binding = dict(artifact=artifact, kind='cpu_capacity', name='resident_numpy_copy.v2',
                    dependencies=self.dependencies, max_age_seconds=self.max_age_seconds, targets=targets)
                profile = bind_detailed_profile(contexts, self.execution, sources=self.source,
                    component_artifacts={artifact:sha256_file(artifact)}, price_bindings=[binding],
                    limitations=['Only resident NumPy copy CPU service is measured and age-checked; other audit prices are synthetic.'])
                state.update(status=result['status'], cpu_seconds_per_byte=coefficient, evidence=validate_price_bindings(profile))
            if state['state'] not in ('checking', 'measuring'):
                self.close()
            state['total_callback_wall_seconds'] = time.perf_counter()-began
            return profile, state
        except Exception:
            self.close()
            raise
