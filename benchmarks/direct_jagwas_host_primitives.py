"""Fixed tiny Torch API bank for joint-reduction dispatch, not GPU throughput."""
import torch


def build_bank(device):
    # Shape/dtype/view setup stays outside every measured call. The fixed
    # dimensions are independent of the association shape and kernel census.
    matrix=torch.ones((32,32),dtype=torch.float64,device=device)
    right=torch.ones_like(matrix)
    transposed=right.T
    single=torch.ones((32,32),dtype=torch.float32,device=device)
    mask=torch.ones((32,32),dtype=torch.bool,device=device)
    rows=torch.ones(32,dtype=torch.bool,device=device)
    column=rows[:,None]
    statistic=torch.ones((32,1),dtype=torch.float64,device=device)
    status=torch.zeros(32,dtype=torch.uint8,device=device)
    # Block-triangular projection (jagwas_projection): out= product into a
    # row slice, in-place square, column sum.
    product=torch.empty((32,32),dtype=torch.float64,device=device)
    scratch=torch.ones((32,32),dtype=torch.float64,device=device)
    sums=torch.ones(32,dtype=torch.float64,device=device)
    # Score form z = t / sqrt(1 + t^2 / df) (jagwas_projection), in the scan
    # precision: out-of-place square, then in-place divide by the df column,
    # add, rsqrt and multiply; plus the FP32 guards ahead of it. Operands keep
    # their values across repeats (divide by 1, add and multiply identities,
    # rsqrt of 1), so no call meets denormals.
    score_form={}
    for dtype,width in ((torch.float32,'fp32'),(torch.float64,'fp64')):
        ones=torch.ones((32,32),dtype=dtype,device=device)
        df_column=torch.ones((32,1),dtype=dtype,device=device)
        score_form.update({
            'joint_square_out_'+width:lambda a=ones:a.square(),
            'joint_div_rows_inplace_'+width:lambda a=ones.clone(),c=df_column:a.div_(c),
            'joint_add_scalar_inplace_'+width:lambda a=ones.clone():a.add_(0.0),
            'joint_rsqrt_inplace_'+width:lambda a=ones.clone():a.rsqrt_(),
            'joint_mul_inplace_'+width:lambda a=ones.clone(),b=ones:a.mul_(b)})
    df=torch.ones(32,dtype=torch.float32,device=device)
    return {
        **score_form,
        'joint_isfinite_fp32':lambda:torch.isfinite(single),
        'joint_nan_to_num_fp32':lambda:torch.nan_to_num(single,nan=0.,posinf=0.,neginf=0.),
        'joint_unsqueeze_view_fp32':lambda:df.unsqueeze(1),
        'joint_fp32_to_fp64':lambda:single.double(),
        'joint_fp64_noop':lambda:matrix.double(),
        'joint_isfinite_fp64':lambda:torch.isfinite(matrix),
        'joint_all_bool_rows':lambda:mask.all(dim=1),
        'joint_invert_bool':lambda:~rows,
        'joint_nan_to_num_fp64':lambda:torch.nan_to_num(matrix,nan=0.,posinf=0.,neginf=0.),
        'joint_transpose_view':lambda:right.T,
        'joint_gemm_out_fp64':lambda:torch.mm(matrix,transposed,out=product),
        'joint_empty_fp64':lambda:torch.empty((32,32),dtype=torch.float64,device=device),
        'joint_slice_view':lambda:matrix[:16],
        'joint_square_inplace_fp64':lambda:scratch.square_(),
        'joint_sum_fp64_columns':lambda:matrix.sum(dim=0),
        'joint_unsqueeze_view_fp64':lambda:sums.unsqueeze(1),
        'joint_status_compare':lambda:status!=0,
        'joint_or_bool':lambda:rows|rows,
        'joint_unsqueeze_view':lambda:rows.unsqueeze(1),
        'joint_masked_fill_fp64':lambda:statistic.masked_fill(column,float('nan')),
        'joint_fill_fp64':lambda:torch.full_like(statistic,float('nan'),dtype=torch.float64),
        'joint_fill_fp32':lambda:torch.full_like(statistic,float('nan'),dtype=torch.float32),
        'joint_fp64_to_fp32':lambda:statistic.to(torch.float32),
        'joint_fill_int32':lambda:torch.zeros_like(statistic,dtype=torch.int32),
    }
