"""Fixed selector dispatch and blocking-copy operations, independent of GWAS."""
import torch


def build_bank(device):
    fp = torch.ones(32, device=device, dtype=torch.float32)
    integer = torch.arange(32, device=device, dtype=torch.int64)
    status = torch.zeros(32, device=device, dtype=torch.uint8)
    matrix = torch.ones((32,32), device=device, dtype=torch.float32)
    boolean = torch.ones((32,32), device=device, dtype=torch.bool)
    empty = torch.zeros_like(boolean)
    nonempty = empty.clone(); nonempty[0,0] = True
    coordinates = torch.empty_strided((32,2),(1,32), device=device, dtype=torch.int64)
    column = fp[:,None]
    rows = torch.arange(32, device=device, dtype=torch.int64)
    valid = torch.ones(32, device=device, dtype=torch.bool)
    mask = torch.empty((32,32), device=device, dtype=torch.bool)
    packed = torch.arange(5*32, device=device, dtype=torch.int32).reshape(5, 32)
    fields = tuple(packed.unbind(0))
    bank = {
        'loop_control': lambda: None,
        'view_alias': lambda: matrix[...],
        'view_slice': lambda: fp[:16],
        'view_two_slices': lambda: matrix[:16,:16],
        'view_column': lambda: fp[:,None],
        'view_coordinate': lambda: coordinates[:,0],
        'df_cast_int64': lambda: fp.to(torch.int64),
        'clamp_int64': lambda: integer.clamp(0,40),
        'gather_vector': lambda: fp[rows],
        'gather_matrix': lambda: matrix[rows,rows],
        'status_equal_zero': lambda: status==0,
        'df_greater_zero': lambda: fp>0,
        'and_bool': lambda: boolean & boolean,
        'isfinite_fp32': lambda: torch.isfinite(matrix),
        'abs_fp32': lambda: matrix.abs(),
        'compare_cutoff': lambda: matrix>=column,
        'nonzero_empty': lambda: empty.nonzero(as_tuple=False),
        'nonzero_nonempty': lambda: nonempty.nonzero(as_tuple=False),
        'index_add_int64': lambda: integer+17,
        'copy_cpu_numpy_fp32': lambda: fp.cpu().numpy(),
        'copy_cpu_numpy_int64': lambda: integer.cpu().numpy(),
        'copy_cpu_numpy_uint8': lambda: status.cpu().numpy(),
        'copy_cpu_numpy_int32': lambda: packed.cpu().numpy(),
        'coordinate_cast_int32': lambda: integer.to(torch.int32),
        'where_limit': lambda: torch.where(valid, fp, torch.inf),
        'empty_mask': lambda: fp.new_empty((32,32), dtype=torch.bool),
        # A full-range slice before None dispatches only the unsqueeze.
        'view_unsqueeze': lambda: fp[0:32, None],
        'compare_cutoff_out': lambda: torch.ge(matrix, column, out=mask),
        'finite_upper': lambda: matrix < torch.inf,
        'and_bool_inplace': lambda: mask.__iand__(boolean),
        'view_unbind': lambda: coordinates.unbind(1),
        'view_dtype': lambda: fp.view(torch.int32),
        'stack_int32': lambda: torch.stack(fields),
        'ready_stream_sync': lambda: torch.cuda.current_stream(device).synchronize(),
    }
    return bank
