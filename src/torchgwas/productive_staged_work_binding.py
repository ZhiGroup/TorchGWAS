"""Audit the transfer, GPU and output prices of a live staged screen.

The active detailed profile and its original measurement ages must already
have been validated. Missing active rates or immutable targets remain explicit
unbound inputs; a supplied rate that disagrees with an active rate is rejected.
"""
from .calibration_cache import _digest
from .productive_staged_source_binding import _covered, _leaves


def _nested(value,path):
    for key in path:
        if not isinstance(value,dict) or key not in value:
            return False,None
        value=value[key]
    return True,value


def audit_staged_work_price_binding(profile,context_name,candidates,
                                    *,reduction_price_evidence=None,
                                    profile_sha256=None,
                                    max_reported_paths=24):
    """Check candidate work prices against active context values and targets.

    Structural settings such as chunk width, phenotype tile, output fields and
    covariate rank are not prices. Shared H2D/D2H and reduced-output primitive
    banks are reported unbound unless the detailed context/record explicitly
    exposes them; this function never invents a shared-bus capacity.
    """
    if (not isinstance(profile,dict) or not isinstance(profile.get('contexts'),list)
            or not isinstance(context_name,str) or not context_name or
            not isinstance(candidates,list) or not 0<len(candidates)<=4 or
            type(max_reported_paths) is not int or not 1<=max_reported_paths<=128):
        raise ValueError('Bounded active staged work price audit required')
    found=[(index,row) for index,row in enumerate(profile['contexts'])
           if isinstance(row,dict) and row.get('name')==context_name]
    if len(found)!=1:
        raise ValueError('One active staged work context required')
    context_index,context=found[0]
    devices=context.get('profiles')
    if not isinstance(devices,dict):
        raise ValueError('Active staged GPU profiles required')
    checked=set();missing=set();external_bound=set()
    def digest_identity(value):
        return (isinstance(value,str) and len(value)==64 and
                all(character in '0123456789abcdef' for character in value))
    bound_external=(isinstance(reduction_price_evidence,dict) and
        digest_identity(reduction_price_evidence.get('record_sha256')) and
        digest_identity(reduction_price_evidence.get('artifact_sha256')))
    external=(reduction_price_evidence.get('record',{}).get('value')
              if bound_external else None)

    def compare(label,value,context_path,*,scenario=False):
        present,active=_nested(context,context_path)
        if not present:
            missing.add(label)
            return
        if scenario:
            if (isinstance(value,bool) or not isinstance(value,(int,float)) or
                    isinstance(active,bool) or not isinstance(active,(int,float)) or
                    not 0<value<=active):
                raise ValueError('Staged work capacity exceeds active context: '+label)
        elif _digest(value)!=_digest(active):
            raise ValueError('Staged work price differs from active profile: '+label)
        checked.update(_leaves(active,(context_index,)+tuple(context_path)))

    for candidate in candidates:
        if (not isinstance(candidate,dict) or
                not isinstance(candidate.get('partitions'),list) or
                not 0<len(candidate['partitions'])<=16 or
                not isinstance(candidate.get('compute_options'),dict) or
                not isinstance(candidate.get('output_options'),dict) or
                not isinstance(candidate.get('mode_service_options'),dict)):
            raise ValueError('Explicit bounded staged work candidate required')
        parts=candidate['partitions']
        active_devices={part.get('device') for part in parts if isinstance(part,dict)}
        if len(active_devices)==0 or not active_devices<=set(devices):
            raise ValueError('Staged work uses an inactive GPU')
        compute=candidate['compute_options'];output=candidate['output_options']
        for field,profile_path in (
                ('per_device_h2d_bytes_per_second',('h2d_bytes_per_second',)),
                ('peak_fp32_flops_per_second',('gpu_resources','fp32_flops_per_second')),
                ('peak_fp64_flops_per_second',('gpu_resources','fp64_flops_per_second'))):
            bank=compute.get(field)
            if field=='peak_fp64_flops_per_second' and bank is None:continue
            if not isinstance(bank,dict) or set(bank)!=active_devices:
                raise ValueError('One staged '+field+' per active GPU required')
            for device,value in bank.items():
                compare(field+'.'+device,value,('profiles',device)+profile_path)
        bank=output.get('per_device_d2h_bytes_per_second')
        if not isinstance(bank,dict) or set(bank)!=active_devices:
            raise ValueError('One staged D2H rate per active GPU required')
        for device,value in bank.items():
            compare('per_device_d2h_bytes_per_second.'+device,value,
                    ('profiles',device,'d2h_bytes_per_second'))
        for options,field,name in (
                (compute,'shared_h2d_bytes_per_second','h2d'),
                (output,'shared_d2h_bytes_per_second','d2h'),
                (output,'output_bytes_per_second','output')):
            if field not in options:
                raise ValueError('Explicit staged '+field+' required')
            context_path=(('shared_capacities',name) if name=='output' else
                          ('shared_transfer_capacities',name))
            compare(field,options[field],context_path,scenario=True)
        compute_links=compute.get('shared_links',())
        output_links=output.get('shared_links',())
        if compute_links!=output_links:
            raise ValueError('Staged H2D/D2H link declarations differ')
        active_links=context.get('shared_links',())
        if active_links and not compute_links:
            raise ValueError('Staged work omits active shared transfer links')
        if compute_links:
            compare('shared_links',compute_links,('shared_links',))
        elif len(active_devices)>1 and not active_links:
            missing.add('shared_link_topology')
        shapes=candidate.get('shape_profiles')
        if shapes is None:
            missing.add('gpu_shape_service')
        else:
            if not isinstance(shapes,dict) or set(shapes)!=active_devices:
                raise ValueError('One staged GPU shape profile per active GPU required')
            excluded={'chunk_markers','kernel_geometry','joint_kernel_geometry'}
            for device,shape in shapes.items():
                if not isinstance(shape,dict):
                    raise ValueError('Staged GPU shape profile required')
                active=devices[device]
                fixed={key:value for key,value in shape.items() if key not in excluded}
                original={key:value for key,value in active.items() if key not in excluded}
                if _digest(fixed)!=_digest(original):
                    raise ValueError('Staged GPU shape prices differ from active profile: '+device)
                # _shape_component reads these independently priced fields.
                # The compiled kernel census is duration-free and is checked
                # for exact identity below, not treated as a timing target.
                for field in ('cpu_fraction','gpu_resources','host_primitives')+(
                        ('joint_host_primitives',) if shape.get('reduction')=='jagwas'
                        else ()):
                    value=active.get(field)
                    if value is None or (field!='cpu_fraction' and
                                         (not isinstance(value,dict) or not value)):
                        missing.add('gpu_shape_service.'+device+'.'+field)
                    else:
                        checked.update(_leaves(value,
                            (context_index,'profiles',device,field)))
                if active.get('host_serial_primitives') is not None:
                    checked.update(_leaves(active['host_serial_primitives'],
                        (context_index,'profiles',device,'host_serial_primitives')))
                for field in ('kernel_geometry','joint_kernel_geometry'):
                    chosen=shape.get(field,[])
                    bank=active.get(field,[])
                    if not isinstance(chosen,list) or not isinstance(bank,list):
                        raise ValueError('Staged GPU geometry bank required')
                    counts={}
                    for item in bank:
                        fingerprint=_digest(item)
                        counts[fingerprint]=counts.get(fingerprint,0)+1
                    for row in chosen:
                        if counts.get(_digest(row))!=1:
                            raise ValueError('Staged GPU shape geometry differs from active profile: '+device)
        services=candidate['mode_service_options']
        reduction=('jagwas' if 'jagwas_writer_fsync' in output else
                   'significant' if output.get('significant_backend') else None)
        if reduction is None:
            writers=services.get('writer_profiles')
            if not isinstance(writers,dict) or set(writers)!=active_devices:
                raise ValueError('One staged dense writer profile per active GPU required')
            for device,writer in writers.items():
                if not isinstance(writer,dict):
                    raise ValueError('Staged dense writer profile required')
                for field in ('cpu_fraction','executor_cpu_seconds',
                              'writeback_service','fsync_seconds'):
                    if field not in writer:
                        raise ValueError('Staged dense writer price missing: '+field)
                    compare('writer.'+device+'.'+field,writer[field],
                            ('profiles',device,field))
                for field in ('bytearray_zero_bytes',):
                    present,value=_nested(writer,('process_units',field))
                    if not present:raise ValueError('Staged dense writer process price missing')
                    compare('writer.'+device+'.process_units.'+field,value,
                            ('profiles',device,'process_units',field))
                if writer.get('writer_copy_service') is not None:
                    compare('writer.'+device+'.writer_copy_service',
                            writer['writer_copy_service'],
                            ('profiles',device,'writer_copy_service'))
                else:
                    active=devices[device]
                    if active.get('writer_copy_service') is not None:
                        raise ValueError('Staged writer uses a legacy copy price despite an active call price')
                    present,value=_nested(writer,('process_units','numpy_copy_bytes'))
                    if not present:raise ValueError('Staged writer fallback copy price missing')
                    compare('writer.'+device+'.process_units.numpy_copy_bytes',value,
                            ('profiles',device,'process_units','numpy_copy_bytes'))
        else:
            def compare_bank(name,fields):
                bank=services.get(name)
                if not isinstance(bank,dict) or set(bank)!=active_devices:
                    raise ValueError('One staged '+name+' per active GPU required')
                for device,row in bank.items():
                    if not isinstance(row,dict):
                        raise ValueError('Staged '+name+' profile required')
                    for path in fields:
                        present,value=_nested(row,path)
                        if not present:
                            raise ValueError('Staged '+name+' price missing: '+'.'.join(path))
                        compare(name+'.'+device+'.'+'.'.join(path),value,
                                ('profiles',device)+path)

            archive_fields=(('cpu_fraction',),
                ('shared_dram_bytes_per_second',),('fsync_seconds',),
                ('writeback_service','pagecache_seconds_per_byte'),
                ('writeback_service','storage_seconds_per_byte'))
            compare_bank('archive_profiles',archive_fields+(
                (('process_units','numpy_copy_bytes'),)
                if reduction=='significant' else ()))
            if reduction=='jagwas':
                fractions=services.get('cpu_fraction_by_device')
                if not isinstance(fractions,dict) or set(fractions)!=active_devices:
                    raise ValueError('One JAGWAS CPU fraction per active GPU required')
                for device,value in fractions.items():
                    compare('jagwas_cpu_fraction.'+device,value,
                            ('profiles',device,'cpu_fraction'))
                writer=(external.get('writer_prices') if isinstance(external,dict)
                        else None)
                for field,source_field in (('selection_prices','prices'),
                                           ('archive_price','archive')):
                    if not isinstance(writer,dict) or source_field not in writer:
                        missing.add('reduced_service.'+field)
                    elif _digest(services.get(field))!=_digest(writer[source_field]):
                        raise ValueError('Staged JAGWAS '+field+' differs from bound reduction record')
                    else:external_bound.add('jagwas.'+field)
            elif output['significant_backend']=='host':
                compare_bank('selection_profiles',(
                    ('cpu_fraction',),('shared_dram_bytes_per_second',)))
                for field,expected in (('selection_prices',external),
                    ('archive_prices',external.get('archive') if
                     isinstance(external,dict) else None)):
                    if expected is None:
                        missing.add('reduced_service.'+field)
                    elif _digest(services.get(field))!=_digest(expected):
                        raise ValueError('Staged significant '+field+' differs from bound reduction record')
                    else:external_bound.add('significant.'+field)
            else:
                launches=services.get('launch_profiles')
                if not isinstance(launches,dict) or set(launches)!=active_devices:
                    raise ValueError('One device launch profile per active GPU required')
                for device,row in launches.items():
                    if not isinstance(row,dict):
                        raise ValueError('Staged device launch profile required')
                    for field in ('kernel_launch_seconds','gpu_fraction'):
                        if field not in row:
                            raise ValueError('Staged device launch price missing: '+field)
                        compare('device_launch.'+device+'.'+field,row[field],
                                ('profiles',device,'gpu_resources',field))
                archive=(external.get('archive') if isinstance(external,dict)
                         else None)
                if archive is None:
                    missing.add('reduced_service.archive_prices')
                elif _digest(services.get('archive_prices'))!=_digest(archive):
                    raise ValueError('Staged device significant archive_prices differs from bound reduction record')
                else:external_bound.add('significant.archive_prices')
                missing.add('reduced_service.count_transfer_prices')
    targets=[]
    for binding in profile.get('price_bindings',[]):
        for target in binding.get('targets',[]):
            path=target.get('context_path')
            if isinstance(path,list):targets.append(tuple(path))
    unbound=sorted((path for path in checked if not _covered(path,targets)),
                   key=lambda path:repr(path))
    absent=sorted(missing)
    return dict(kind='torchgwas.staged_work_price_binding.v1',
        status=('declared_work_prices_verified' if not unbound and not absent
                else 'work_prices_match_but_unbound'),
        context=context_name,profile_sha256=(_digest(profile) if
            profile_sha256 is None else profile_sha256),
        matched_price_leaves=len(checked),
        declared_price_leaves=len(checked)-len(unbound),
        unbound_price_leaves=len(unbound),
        unbound_paths=[list(path) for path in unbound[:max_reported_paths]],
        truncated_unbound_paths=max(0,len(unbound)-max_reported_paths),
        missing_active_prices=absent[:max_reported_paths],
        truncated_missing_active_prices=max(0,len(absent)-max_reported_paths),
        external_price_fields_bound=sorted(external_bound),
        external_record_sha256=(None if not bound_external else
                                reduction_price_evidence.get('record_sha256')),
        external_artifact_sha256=(None if not bound_external else
                                  reduction_price_evidence.get('artifact_sha256')),
        candidate_count=len(candidates),
        scope='Exact candidate transfer/GPU/output values where the previously validated detailed profile declares them; host significant/JAGWAS primitive banks may additionally match the previously validated immutable reduction record. Missing shared rates, shape service and device-selector prices remain explicit. This does not re-read records or certify live capacity, finite completion or selection.')
