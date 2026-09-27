"""Live physical UUID binding must be backend-independent and never memoized."""
import ctypes
from types import SimpleNamespace
from unittest.mock import Mock
import pytest

from torchgwas import gpu_identity as identity


A='01234567-89ab-cdef-0123-456789abcdef'
B='11234567-89ab-cdef-0123-456789abcdef'


def fake_library(monkeypatch):
    state=dict(driver=b'550.54.15',requests=[],handles={A:7,B:3})
    def driver(buffer,length):buffer.value=state['driver'];return 0
    def handle(value,out):
        uuid=value.decode().removeprefix('GPU-');state['requests'].append(uuid)
        ctypes.cast(out,ctypes.POINTER(ctypes.c_void_p))[0]=state['handles'][uuid]
        return 0
    def uuid(handle,buffer,length):
        key=next(k for k,v in state['handles'].items() if v==handle.value)
        buffer.value=('GPU-'+key).encode();return 0
    def pci(handle,out):
        value=ctypes.cast(out,ctypes.POINTER(identity._PciInfo)).contents
        value.busId=b'00000000:3B:00.0' if handle.value==7 else b'00000000:AF:00.0'
        return 0
    lib=SimpleNamespace(nvmlInit_v2=Mock(return_value=0),nvmlShutdown=Mock(return_value=0),
        nvmlSystemGetDriverVersion=Mock(side_effect=driver),nvmlDeviceGetHandleByUUID=Mock(side_effect=handle),
        nvmlDeviceGetUUID=Mock(side_effect=uuid),nvmlDeviceGetPciInfo_v3=Mock(side_effect=pci))
    monkeypatch.setattr(identity.ctypes,'CDLL',Mock(return_value=lib))
    return lib,state


def test_nvml_queries_by_uuid_preserve_order_and_fresh_driver(monkeypatch):
    lib,state=fake_library(monkeypatch)
    monkeypatch.setattr(identity,'_smi_identity',Mock(side_effect=AssertionError('unexpected fallback')))
    first=identity.physical_gpu_identity([B,A]);state['driver']=b'580.99.02'
    second=identity.physical_gpu_identity([A,B])
    assert state['requests']==[B,A,A,B]
    assert first[B]==dict(pci_bus_id='00000000:AF:00.0',driver_version='550.54.15')
    assert first[A]['pci_bus_id']==second[A]['pci_bus_id']=='00000000:3B:00.0'
    assert second[A]['driver_version']=='580.99.02'
    assert lib.nvmlInit_v2.call_count==lib.nvmlShutdown.call_count==2
    assert ctypes.sizeof(identity._PciInfo)==68


@pytest.mark.parametrize('failure',['nvmlInit_v2','nvmlSystemGetDriverVersion',
    'nvmlDeviceGetHandleByUUID','nvmlDeviceGetUUID','nvmlDeviceGetPciInfo_v3','nvmlShutdown'])
def test_query_failure_discards_partial_result_and_balances_init(monkeypatch,failure):
    lib,state=fake_library(monkeypatch)
    function=getattr(lib,failure);function.side_effect=None;function.return_value=999
    fallback=Mock(return_value={'from':'fallback'});monkeypatch.setattr(identity,'_smi_identity',fallback)
    assert identity.physical_gpu_identity([A])=={'from':'fallback'}
    fallback.assert_called_once_with([A])
    assert lib.nvmlShutdown.call_count==int(failure!='nvmlInit_v2')


def test_wrong_handle_and_missing_symbols_do_not_pass_binding(monkeypatch):
    lib,_=fake_library(monkeypatch)
    def wrong(handle,buffer,length):buffer.value=('GPU-'+B).encode();return 0
    lib.nvmlDeviceGetUUID.side_effect=wrong
    with pytest.raises(ValueError,match='differs'):
        identity._nvml_identity([A])
    assert lib.nvmlShutdown.call_count==1
    del lib.nvmlDeviceGetPciInfo_v3
    lib.nvmlInit_v2.reset_mock()
    with pytest.raises(AttributeError):identity._nvml_identity([A])
    lib.nvmlInit_v2.assert_not_called()


def test_smi_fallback_matches_uuid_not_index_and_normalizes_bus(monkeypatch):
    monkeypatch.setattr(identity.ctypes,'CDLL',Mock(side_effect=OSError('no library')))
    command=Mock(return_value=f'GPU-{A}, 0000:3b:00.0, 550.54.15\nGPU-{B}, 00000000:AF:00.0, 550.54.15\n')
    monkeypatch.setattr(identity.subprocess,'check_output',command)
    result=identity.physical_gpu_identity([B,A])
    assert list(result)==[B,A] and result[A]['pci_bus_id']=='00000000:3B:00.0'
    command.assert_called_once()
    command.return_value=f'GPU-{B}, 0000:AF:00.0, 550.54.15\n'
    with pytest.raises(ValueError,match='measured physical'):identity.physical_gpu_identity([A])
    command.return_value*=2
    with pytest.raises(ValueError,match='Duplicate'):identity.physical_gpu_identity([B])


@pytest.mark.parametrize('values',[[],[A,A],[''],['GPU-'+A],['MIG-'+A],['unknown'],[None]])
def test_invalid_physical_uuid_rejected_before_io(monkeypatch,values):
    native=Mock(side_effect=AssertionError('must not query'))
    monkeypatch.setattr(identity,'_nvml_identity',native)
    with pytest.raises(ValueError):identity.physical_gpu_identity(values)


@pytest.mark.parametrize('value',['cuda','cuda:-1','cuda:01','cuda:0.0','cpu',None,0])
def test_noncanonical_device_rejected(value):
    with pytest.raises(ValueError):identity.canonical_cuda_device(value)


def test_canonical_device():
    assert identity.canonical_cuda_device('cuda:12')=='cuda:12'
