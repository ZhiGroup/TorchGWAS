"""Host memory at voxel scale: what the trait-block planner counts, and the limits it obeys."""
import json

import numpy as np
import pytest
import torch

from torchgwas.pipeline_model import auto_trait_block, host_pinned_bytes

# 100,000 subjects and 2,000,000 traits, packed PGEN rows, four 80 GB GPUs, a 500 GB host.
SCALE = dict(n_samples=100000, n_traits=2000000, covariate_rank=11, chunk_variants=4096, depth=16,
             transfer_bytes_per_variant=25024.0, device_memory_bytes=80 * 1024**3,
             host_memory_bytes=500 * 1024**3, trait_devices=4)


def test_significant_blocks_fit_the_host_with_their_residualised_phenotypes():
    selected = auto_trait_block(**SCALE, device_selection=True)
    ringed = auto_trait_block(**SCALE)
    # Significant pairs pin no result ring, so the host allows a wider block:
    # here the device's (~163,000 traits) rather than the ring's.
    assert selected == auto_trait_block(**dict(SCALE, host_memory_bytes=None)) > ringed
    # Each GPU's block holds its residualised phenotypes (float32 samples x
    # block) on the host besides its pinned rings; four of them fit.
    held = host_pinned_bytes(chunk_variants=4096, depth=16, n_traits=selected, transfer_bytes_per_variant=25024.0,
                             result_ring=False) + 4.0 * 100000 * selected
    assert held * 4 <= 0.85 * SCALE['host_memory_bytes']
    # A smaller host narrows the block rather than overcommitting it.
    assert auto_trait_block(**dict(SCALE, host_memory_bytes=200 * 1024**3), device_selection=True) < selected


def test_cgroup_limit_bounds_host_memory(tmp_path):
    from torchgwas.api import _cgroup_available_bytes
    gib = 1024**3
    # cgroup v2: Slurm limits the job; the step below it sets no limit of its own.
    job = tmp_path / 'slurm' / 'job_1'
    step = job / 'step_batch'
    step.mkdir(parents=True)
    (job / 'memory.max').write_text(str(480 * gib))
    (job / 'memory.current').write_text(str(100 * gib))
    (job / 'memory.stat').write_text(f'anon {90 * gib}\ninactive_file {6 * gib}\n')
    (step / 'memory.max').write_text('max')
    (step / 'memory.current').write_text(str(100 * gib))
    (step / 'memory.stat').write_text('')
    membership = tmp_path / 'cgroup'
    membership.write_text('0::/slurm/job_1/step_batch\n')
    assert _cgroup_available_bytes(root=tmp_path, membership=membership) == (480 - 100 + 6) * gib
    # No limit anywhere: None, and MemAvailable decides.
    (job / 'memory.max').write_text('max')
    assert _cgroup_available_bytes(root=tmp_path, membership=membership) is None


class LazyPanel:
    """A phenotype read only by slicing (as h5py, zarr or a generated panel); never converted whole."""

    def __init__(self, values, strict=True):
        self.values, self.shape, self.dtype, self.ndim = values, values.shape, values.dtype, values.ndim
        self.reads, self.strict = 0, strict

    def __getitem__(self, key):
        self.reads += 1
        return np.array(self.values[key])

    def __array__(self, dtype=None, copy=None):
        if self.strict:
            raise AssertionError('the whole lazy panel was converted')
        return np.asarray(self.values, dtype=dtype)


@pytest.mark.skipif(not torch.cuda.is_available(), reason='CUDA required')
@pytest.mark.parametrize('case', ['tiles', 'drop_subject', 'autotune', 'autotune_tiles'])
def test_a_lazy_panel_is_read_by_blocks(tmp_path, monkeypatch, case):
    if case == 'autotune_tiles' and torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    from test_min_p import _inputs, _read
    from torchgwas.api import run_linear_gwas
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    path, y, covariates, _ = _inputs(tmp_path, missing_calls=False, missing_pheno=False, k=11)
    if case == 'drop_subject':
        y[[4, 30], [2, 7]] = np.nan
    layout = dict(
        autotune=dict(autotune=True, autotune_options=dict(devices=['cuda:0'], chunk_sizes=[16], max_tile_traits=4,
                                                           cpus_per_device=1)),
        # Phenotype tiles over two GPUs, as a panel too wide for one GPU gets.
        autotune_tiles=dict(autotune=True, autotune_options=dict(devices=['cuda:0', 'cuda:1'], chunk_sizes=[16],
                                                                 min_tile_traits=2, max_tile_traits=4,
                                                                 split='traits', cpus_per_device=1)),
    ).get(case, dict(device='cuda:0', chunk_size=16, trait_block=4))
    options = dict(genotype_format='pgen', pgen_mode='hardcall', compute_dtype='float32', reduce='significant',
                   significance_threshold=0.5, **layout)
    # Tiled runs read blocks only; autotune scans 11 traits untiled, which
    # converts the panel by design.
    panel = LazyPanel(y, strict=case != 'autotune')
    run_linear_gwas(path, panel, covariates, output_dir=tmp_path / 'lazy', **options)
    run_linear_gwas(path, y, covariates, output_dir=tmp_path / 'array', **options)
    assert panel.reads > 0
    got, want = _read(tmp_path / 'lazy' / 'sumstats')[1], _read(tmp_path / 'array' / 'sumstats')[1]
    assert got.keys() == want.keys()
    for key in want:
        np.testing.assert_array_equal(got[key], want[key])


@pytest.mark.skipif(not torch.cuda.is_available(), reason='CUDA required')
def test_the_default_device_still_blocks_traits(tmp_path, monkeypatch):
    # device='auto' (the default) resolves to a GPU before the blocking
    # decision; it used to skip it, and a panel needing blocks was
    # residualised whole on the host.
    from test_min_p import _inputs
    from torchgwas.api import run_linear_gwas
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    path, y, covariates, _ = _inputs(tmp_path, missing_calls=False, missing_pheno=False, k=11)
    monkeypatch.setattr('torchgwas.pipeline_model.auto_trait_block', lambda **kwargs: 4)
    options = dict(genotype_format='pgen', pgen_mode='hardcall', compute_dtype='float32', chunk_size=16,
                   reduce='significant', significance_threshold=0.5)
    run_linear_gwas(path, y, covariates, output_dir=tmp_path / 'auto', **options)
    run = json.loads((tmp_path / 'auto' / 'run.json').read_text())
    # The 4-trait block, narrowed to one block per visible GPU.
    cards = torch.cuda.device_count()
    assert run['trait_block'] == (min(4, -(-11 // cards)) if cards > 1 else 4)
