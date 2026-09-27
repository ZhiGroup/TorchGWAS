"""Fixed operation prices must represent the actual typed source APIs."""
import importlib.util
from pathlib import Path
import pytest
import torch
from torch.utils._python_dispatch import TorchDispatchMode
from torchgwas.reduction_tensor_work import jagwas_tensor_work,jagwas_host_primitive_name

path=Path(__file__).parents[1]/'benchmarks'/'direct_jagwas_host_primitives.py'
spec=importlib.util.spec_from_file_location('joint_primitive_bank',path)
bank_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(bank_module)


@pytest.mark.parametrize('dtype',['float32','float64'])
def test_fixed_bank_preserves_typed_operation_sequence_and_gemm_strides(dtype):
    bank=bank_module.build_bank('meta')
    trace=jagwas_tensor_work(64,32,32,phase='reduce',compute_dtype=dtype)
    observed=[]
    class Record(TorchDispatchMode):
        def __torch_dispatch__(self,func,types,args=(),kwargs=None):
            result=func(*args,**(kwargs or {}))
            if str(func)!='aten.detach.default':
                tensors=[v for v in args if isinstance(v,torch.Tensor)]
                observed.append((str(func),[str(v.dtype) for v in tensors]))
                if str(func)=='aten.mm.default':
                    assert args[0].stride()==(32,1)
                    assert args[1].stride()==(1,32)
            return result
    for call in trace['host_calls']:
        name=jagwas_host_primitive_name(call)
        assert name in bank
        observed.clear()
        with Record():bank[name]()
        expected=[(row['op'],[v['dtype'] for v in row['inputs']])
            for i in call['step_indices'] if (row:=trace['steps'][i])['op']!='aten.detach.default']
        assert observed==expected,(name,observed,expected)
