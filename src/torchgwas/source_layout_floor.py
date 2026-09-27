"""Whole-source read/decode floors across fixed GPU GWAS partitions.

Each partition supplies a priced, source-bound native PGEN schedule. Shared
CPU, memory and input capacity is counted once across the layout, while
reader-worker chains for tiles assigned to one GPU are serial. This is a
source-stage floor, not a GPU/output continuation or switch decision.
"""
import math

from .analytical_plan_cache import input_identity


def _span(value,name):
    if (not isinstance(value,(tuple,list)) or len(value)!=2 or
            any(type(x) is not int for x in value) or not 0<=value[0]<value[1]):
        raise ValueError('Nonempty '+name+' range required')
    return tuple(value)


def _interval(value,name):
    if (not isinstance(value,(tuple,list)) or len(value)!=2 or
            any(isinstance(x,bool) or not isinstance(x,(int,float)) or
                not math.isfinite(x) or x<0 for x in value) or value[0]>value[1]):
        raise ValueError('Finite nonnegative '+name+' interval required')
    return tuple(value)


def native_layout_source_floor(partitions, *, total_traits, reduction,
                               partition_axis, max_partitions=10000):
    """Combine mandatory source work under one shared-capacity scenario.

    The caller binds these partitions to its exact unissued executor frontier.
    Dense and significant trait tiles may reread overlapping variant ranges;
    the repeated physical source work is then added. A GPU's tiles are assumed
    serial as in the public executor. JAGWAS retains the full panel per GPU.
    """
    if (type(max_partitions) is not int or max_partitions<1 or
            not isinstance(partitions,(tuple,list)) or
            not 0<len(partitions)<=max_partitions):
        raise ValueError('Bounded nonempty source partitions required')
    if type(total_traits) is not int or total_traits<1:
        raise ValueError('Positive total phenotype count required')
    if reduction not in (None,'significant','jagwas'):
        raise ValueError('Unsupported source output mode')
    if partition_axis not in ('trait','variant') or (reduction=='jagwas' and partition_axis!='variant') or (reduction=='significant' and partition_axis!='trait'):
        raise ValueError('Output mode differs from fixed partition axis')
    ids=set();devices=set();traits=[];variants=[];source=None;capacities=None;chunk=None;samples=None
    cpu=[[],[]];dram=[];input_bytes=[];by_device={};rows=[]
    for row in partitions:
        if not isinstance(row,dict) or set(row)!={'id','device','trait_range','floor'}:
            raise ValueError('Explicit partition id, device, traits and source floor required')
        key=row['id'];device=row['device'];floor=row['floor']
        if not isinstance(key,str) or not key or key in ids:
            raise ValueError('Unique source partition ids required')
        ids.add(key)
        if not isinstance(device,str) or not device.startswith('cuda:') or not device[5:].isdigit():
            raise ValueError('Explicit CUDA source device required')
        if partition_axis=='variant' and device in devices:
            raise ValueError('Variant shards require one source partition per GPU')
        devices.add(device)
        trait=_span(row['trait_range'],'phenotype')
        if trait[1]>total_traits or (partition_axis=='variant' and trait!=(0,total_traits)):
            raise ValueError('Variant shards retain the complete phenotype panel')
        if (not isinstance(floor,dict) or
                floor.get('kind')!='torchgwas.pgen_schedule_source_floor.v1'):
            raise ValueError('Typed native PGEN source floor required')
        variant=_span(floor.get('variant_range'),'source variant')
        traits.append(trait);variants.append(variant)
        if not isinstance(floor.get('input_identity'),dict) or not floor['input_identity']:
            raise ValueError('Bound source identity required')
        if source is None:source=floor['input_identity']
        elif source!=floor['input_identity']:raise ValueError('Source identities differ across GPUs')
        count=floor.get('samples')
        if type(count) is not int or count<1:
            raise ValueError('Bound native PGEN sample count required')
        if samples is None:samples=count
        elif samples!=count:raise ValueError('PGEN sample count differs across partitions')
        if (not isinstance(floor.get('capacities'),dict) or
                set(floor['capacities'])!={'cpu','dram','input'} or
                any(isinstance(v,bool) or not isinstance(v,(int,float)) or
                    not math.isfinite(v) or v<=0 for v in floor['capacities'].values())):
            raise ValueError('Complete positive shared source capacities required')
        if capacities is None:capacities=floor['capacities']
        elif capacities!=floor['capacities']:raise ValueError('Shared source capacities differ')
        if type(floor.get('chunk_markers')) is not int or floor['chunk_markers']<1:
            raise ValueError('Bound source chunk size required')
        if chunk is None:chunk=floor['chunk_markers']
        elif chunk!=floor['chunk_markers']:raise ValueError('One fixed layout chunk size required')
        if (type(floor.get('chunk_count')) is not int or
                floor['chunk_count']!=(variant[1]-variant[0]+chunk-1)//chunk):
            raise ValueError('Source chunk geometry differs from its partition')
        work=floor.get('resource_work')
        if not isinstance(work,dict):raise ValueError('Source resource work required')
        core=_interval(work.get('cpu_seconds'),'CPU work')
        worker=_interval(floor.get('reader_worker_floor_seconds'),'reader worker')
        for name,target in [('dram_bytes',dram),('input_bytes',input_bytes)]:
            value=work.get(name)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
                raise ValueError('Finite nonnegative '+name+' required')
            target.append(value)
        for i in (0,1):cpu[i].append(core[i])
        chain=by_device.setdefault(device,[[],[]])
        for i in (0,1):chain[i].append(worker[i])
        rows.append(dict(id=key,device=device,variant_range=list(variant),
                         trait_range=list(trait),chunks=floor['chunk_count']))
    spans=sorted(traits if partition_axis=='trait' else variants)
    if any(previous[1]>following[0] for previous,following in zip(spans,spans[1:])):
        raise ValueError(('Trait' if partition_axis=='trait' else 'Variant')+' partitions overlap')
    total_cpu=[math.fsum(part) for part in cpu]
    total_dram=math.fsum(dram);total_input=math.fsum(input_bytes)
    resource=[max(value/capacities['cpu'],total_dram/capacities['dram'],
                  total_input/capacities['input']) for value in total_cpu]
    chains={device:[math.fsum(part) for part in pair]
            for device,pair in by_device.items()}
    worker=[max(pair[i] for pair in chains.values()) for i in (0,1)]
    floor=[max(resource[i],worker[i]) for i in (0,1)]
    if not all(math.isfinite(value) for value in (*total_cpu,total_dram,total_input,*floor)):
        raise ValueError('Whole-layout source work overflow')
    if 'path' not in source or input_identity(source['path'])!=source:
        raise ValueError('PGEN input changed before source-layout composition')
    return dict(kind='torchgwas.pgen_layout_source_floor.v1',
        input_identity=dict(source),samples=samples,chunk_markers=chunk,partition_axis=partition_axis,
        reduction=reduction,total_traits=total_traits,partitions=rows,
        shared_capacities=dict(capacities),
        resource_work=dict(cpu_seconds=total_cpu,dram_bytes=total_dram,
                           input_bytes=total_input),
        shared_resource_floor_seconds=resource,
        per_device_reader_floor_seconds=chains,
        source_stage_floor_seconds=floor,
        prediction_complete=False,selection_validated=False,
        scope='Conditional mandatory PGEN read/decode source-stage floor across fixed GPU partitions. Shared CPU/DRAM/input work is counted once; serial tiles on one GPU contribute separate reader-worker chains. It omits GPU, transfer, selection, output, in-flight state and final drain, and is not a hardware-time guarantee or chunk-switch decision.')
