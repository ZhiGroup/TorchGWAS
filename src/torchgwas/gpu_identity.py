"""Fresh physical GPU identity without spawning nvidia-smi on the normal path.

Read-only NVML calls, resolved by CUDA's physical UUID, never NVML device index.
No device/load/property cache and no Python NVML dependency. The subprocess
fallback preserves context binding on hosts without the required NVML symbols.
"""
import ctypes
import re
import subprocess
import uuid as uuid_module


def canonical_cuda_device(value):
    if not isinstance(value,str) or not value.startswith('cuda:') or not value[5:].isdigit():
        raise ValueError('Explicit cuda:N device names required')
    if str(int(value[5:]))!=value[5:]:
        raise ValueError('Device names must use canonical CUDA indices')
    return value


class _PciInfo(ctypes.Structure):
    # NVIDIA's nvmlPciInfo_t ABI for nvmlDeviceGetPciInfo_v3, not the legacy API.
    # https://docs.nvidia.com/deploy/nvml-api/api/structnvmlPciInfo__t.html
    _fields_=[('busIdLegacy',ctypes.c_char*16),('domain',ctypes.c_uint),
        ('bus',ctypes.c_uint),('device',ctypes.c_uint),('pciDeviceId',ctypes.c_uint),
        ('pciSubSystemId',ctypes.c_uint),('busId',ctypes.c_char*32)]


def _record(bus,driver):
    if not re.fullmatch(r'[0-9a-fA-F]{4,8}:[0-9a-fA-F]{2}:[0-9a-fA-F]{2}\.[0-7]',bus):
        raise ValueError('Invalid physical GPU PCI identity')
    if not driver or driver.strip()!=driver or any(c.isspace() for c in driver):
        raise ValueError('Invalid GPU driver identity')
    domain,bus,slot=bus.split(':')
    return dict(pci_bus_id=f'{int(domain,16):08X}:{bus.upper()}:{slot.upper()}',driver_version=driver)


def _nvml_identity(uuids):
    library=ctypes.CDLL('libnvidia-ml.so.1')
    signatures={
        'nvmlInit_v2':[], 'nvmlShutdown':[],
        'nvmlSystemGetDriverVersion':[ctypes.POINTER(ctypes.c_char),ctypes.c_uint],
        'nvmlDeviceGetHandleByUUID':[ctypes.c_char_p,ctypes.POINTER(ctypes.c_void_p)],
        'nvmlDeviceGetPciInfo_v3':[ctypes.c_void_p,ctypes.POINTER(_PciInfo)],
        'nvmlDeviceGetUUID':[ctypes.c_void_p,ctypes.POINTER(ctypes.c_char),ctypes.c_uint]}
    calls={}
    for name,arguments in signatures.items():
        function=getattr(library,name);function.argtypes=arguments;function.restype=ctypes.c_int
        calls[name]=function
    def checked(name,*args):
        result=calls[name](*args)
        if result!=0:raise OSError(f'{name} returned NVML error {result}')
    checked('nvmlInit_v2')
    try:
        driver=ctypes.create_string_buffer(80)
        checked('nvmlSystemGetDriverVersion',driver,len(driver))
        result={}
        for uuid in uuids:
            handle=ctypes.c_void_p();requested=('GPU-'+uuid).encode('ascii')
            checked('nvmlDeviceGetHandleByUUID',requested,ctypes.byref(handle))
            found=ctypes.create_string_buffer(96)
            checked('nvmlDeviceGetUUID',handle,found,len(found))
            if found.value.decode('ascii').lower()!='gpu-'+uuid:
                raise ValueError('NVML handle differs from physical CUDA UUID')
            pci=_PciInfo();checked('nvmlDeviceGetPciInfo_v3',handle,ctypes.byref(pci))
            result[uuid]=_record(pci.busId.decode('ascii'),driver.value.decode('ascii'))
        return result
    finally:
        # Balance this initialization, including on a partial-query failure.
        checked('nvmlShutdown')


def _smi_identity(uuids):
    records=subprocess.check_output(['nvidia-smi','--query-gpu=uuid,pci.bus_id,driver_version',
        '--format=csv,noheader'],text=True)
    physical={}
    for line in records.splitlines():
        uuid,bus,driver=[x.strip() for x in line.split(',')]
        key=uuid.removeprefix('GPU-').lower()
        if key in physical:raise ValueError('Duplicate physical GPU UUID')
        physical[key]=_record(bus,driver)
    if any(uuid not in physical for uuid in uuids):
        raise ValueError('CUDA UUID must identify a measured physical GPU')
    return {uuid:physical[uuid] for uuid in uuids}


def physical_gpu_identity(uuids):
    uuids=list(uuids)
    if not uuids or len(set(uuids))!=len(uuids):raise ValueError('Unique physical GPU UUIDs required')
    for value in uuids:
        if not isinstance(value,str) or str(uuid_module.UUID(value))!=value:
            raise ValueError('Canonical physical GPU UUID required')
    try:return _nvml_identity(uuids)
    except (OSError,AttributeError,UnicodeError,ValueError):
        return _smi_identity(uuids)
