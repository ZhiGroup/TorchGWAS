"""Conditional storage schedule for explicit binary-writer range requests."""
import math


class RangeWritebackSchedule:
    def __init__(self, ledger, service, *, cpu_resources=None):
        required=('pagecache_seconds_per_byte','storage_seconds_per_byte',
                  'submit_seconds','wait_seconds','fadvise_seconds')
        for key in required+tuple(key for key in ['fadvise_eviction_seconds_per_byte'] if key in service):
            value=service[key]
            if not math.isfinite(value) or value<0 or (key=='storage_seconds_per_byte' and value==0):
                raise ValueError('Invalid writeback service: '+key)
        self.ledger=ledger;self.service=service;self.cpu_resources=cpu_resources
        self.storage_resources={'output':1/service['storage_seconds_per_byte']}
        self.pending={};self.bytes=0

    def after_write(self,g,array,index,after):
        last=after
        for step,action in enumerate(self.ledger['events'][index]['actions']):
            kind=action['kind'];offset=action['offset'];length=action['bytes']
            prefix=f'writer:{array}:range:{index}:{step}'
            if kind=='submit':
                last=g.add(prefix+':submit',self.service['submit_seconds'],[last],resources=self.cpu_resources)
                self.pending[offset]=g.add(prefix+':storage',length*self.service['storage_seconds_per_byte'],[last],resources=self.storage_resources)
                self.bytes+=length
            elif kind=='wait':
                last=g.add(prefix+':wait',self.service['wait_seconds'],[last,self.pending[offset]],resources=self.cpu_resources)
            elif kind=='drop_cache':
                last=g.add(prefix+':drop',self.service['fadvise_seconds']+length*self.service.get('fadvise_eviction_seconds_per_byte',0.),[last],resources=self.cpu_resources)
        return last

    def before_fsync(self,g,array,after):
        deps=[after]+list(self.pending.values())
        tail=self.ledger['unsubmitted_tail_bytes']
        if tail:
            deps.append(g.add(f'writer:{array}:tail:storage',tail*self.service['storage_seconds_per_byte'],[after],resources=self.storage_resources))
            self.bytes+=tail
        return g.add(f'writer:{array}:storage:complete',0.,deps)
