"""One process per GPU for a trait-tiled scan: run_linear_gwas(trait_workers='process').

Threads share one interpreter, so with many GPUs each card's per-chunk Python
-- kernel launches, the complete-case correction, pair selection -- waits for
the others' under the GIL. Seven GPUs on small tiles took 25.8 s as threads
and 19.3 s as processes, and a thread asking for the GIL every millisecond
waited 5.6 ms at the 99th percentile against 0.3 ms (2026-10-10). The 2M-voxel
scan kept its GPUs 30-45% busy as threads and 100% as one process per GPU.

Each device gets a contiguous range of the phenotype columns and a child
process running the ordinary single-device scan on it, with
CUDA_VISIBLE_DEVICES set to that card. Before scanning, the children report
their QC-kept trait counts and receive the total, so reduce='significant'
uses one Bonferroni threshold over every kept trait, as a thread run does.
The parent then merges what they wrote into the store a thread run writes:

- pair stores (reduce='significant', or a p_value_threshold filter): every
  part with trait_index rebased onto the kept traits of all partitions, then
  ordered by (variant, trait);
- reduce='min-p': per variant, the partition winner with the largest
  -log10 P (a tie goes to the lower trait index);
- full output: one trait-tiled manifest over the children's tiles.

Every per-trait step is column-separable under missing_phenotype 'exact' or
'impute' (QC, per-value outlier masking, residualization, the scan), so a
partition computes exactly what the whole panel would for its columns.
'drop_subject' is not -- a sample missing any trait leaves every trait -- and
is refused, as is JAGWAS.

The cost: each child decodes the genotypes itself, where threads share one
decoder (shared_decode). That matters little when the scan is GPU-bound --
many traits per card, as in a voxel scan -- which is where threads contend.
"""
from __future__ import annotations

import json
import os
import pickle
import selectors
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

WORK_DIRECTORY = '.trait_processes'
# Rows per part when min-p winners are written back.
MERGED_PART_ROWS = 1 << 20
# Arguments a process run cannot pass on to its children (p_value_threshold and
# autotune are untested across processes).
UNSUPPORTED = ('variant_devices', 'autotune_profile', 'autotune_config', 'initial_calibration',
               'pipeline_profile', '_internal_reduction', 'jagwas_groups', 'jagwas_rcond',
               'jagwas_min_residual', 'p_value_threshold', 'autotune')
# trait_workers='auto' takes processes from this many trait devices. Seven
# H100s, 7 x 5,000 voxels: scan 23.6 s as threads, 17.5 s as processes
# (2026-10-10); two GPUs showed no GIL contention.
AUTO_PROCESS_DEVICES = 4


def resolve_trait_workers(arguments):
    """'process' or 'thread' for run_linear_gwas's trait_workers ('auto', 'thread' or 'process')."""
    requested, devices = arguments['trait_workers'], arguments['trait_devices']
    if requested == 'thread' or devices is None or len(devices) < 2:
        return 'thread'
    if requested == 'process':
        return 'process'
    return 'process' if len(devices) >= AUTO_PROCESS_DEVICES and not process_problems(arguments) else 'thread'


def run_trait_processes(arguments):
    """run_linear_gwas's arguments (its locals at entry) -> GWASResult, one child process per trait device."""
    from .types import GWASResult
    started = time.perf_counter()
    arguments = dict(arguments)
    devices = _validate(arguments)
    output_dir = Path(arguments['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    work = output_dir / WORK_DIRECTORY
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    names, phenotype = _trait_columns(arguments, work)
    bounds = _partitions(len(names), len(devices))
    devices = devices[:len(bounds)]
    budget = 24 if arguments['reader_workers'] is None else int(arguments['reader_workers'])
    if budget < len(devices):
        raise ValueError('reader_workers must provide at least one reader per trait device')
    report_kept = arguments['reduce'] == 'significant' and arguments['significance_threshold'] is None
    children = []
    try:
        for index, (device, (first, last)) in enumerate(zip(devices, bounds)):
            child = _child_arguments(arguments, phenotype, names, first, last, work / f'part_{index:02d}',
                                     budget // len(devices) + (index < budget % len(devices)))
            children.append(_launch(index, device, child, report_kept, work))
        _wait(children, report_kept)
    except BaseException:
        for child in children:
            _stop(child['process'])
        raise
    merged = _merge(arguments, children, output_dir)
    shutil.rmtree(work)
    merged['run']['runtime_seconds'] = time.perf_counter() - started
    from .utils import write_json
    write_json(merged['run'], output_dir / 'run.json')
    write_json(merged['qc'], output_dir / 'qc.json')
    return GWASResult(table=[], run_metadata=merged['run'], qc_summary=merged['qc'])


def process_problems(arguments):
    """Why this run cannot take one process per trait device; empty when it can."""
    problems = []
    devices = arguments['trait_devices']
    names = [str(device) for device in devices] if isinstance(devices, (list, tuple)) else []
    if len(names) < 2:
        problems.append("trait_workers='process' needs two or more trait_devices")
    elif len(set(names)) != len(names) or not all(name.startswith('cuda') for name in names):
        problems.append("trait_workers='process' needs distinct CUDA trait_devices")
    if arguments['output_dir'] is None or arguments['sumstats_format'] != 'binary':
        problems.append("trait_workers='process' writes a binary store: output_dir and sumstats_format='binary'")
    reduce = arguments['reduce']
    if isinstance(reduce, str) and reduce.replace('_', '-').lower() == 'min-p':
        reduce = 'min-p'
    if reduce not in (None, 'significant', 'min-p'):
        problems.append(f"trait_workers='process' supports full output, reduce='significant' and "
                        f"reduce='min-p', not reduce={reduce!r}")
    if arguments['missing_phenotype'] not in (None, 'exact', 'impute'):
        problems.append("trait_workers='process' needs missing_phenotype 'exact' or 'impute': under "
                        "'drop_subject' a sample missing any trait leaves every partition")
    for name in UNSUPPORTED:
        if arguments.get(name) not in (None, False):
            problems.append(f"trait_workers='process' does not support {name}")
    if 'devices' in (arguments.get('autotune_options') or {}):
        problems.append("trait_workers='process' takes its devices from trait_devices, not autotune_options")
    if not isinstance(arguments['genotype'], (str, Path)):
        problems.append("trait_workers='process' needs the genotype as a file path")
    if arguments['phenotype_table'] is None and not isinstance(arguments['phenotype'], (str, Path, np.ndarray)):
        problems.append("trait_workers='process' needs the phenotype as a path, an array or phenotype_table")
    return problems


def _validate(arguments):
    problems = process_problems(arguments)
    if problems:
        raise ValueError(problems[0])
    if isinstance(arguments['reduce'], str) and arguments['reduce'].replace('_', '-').lower() == 'min-p':
        arguments['reduce'] = 'min-p'
    return [str(device) for device in arguments['trait_devices']]


def _trait_columns(arguments, work):
    """(names of every phenotype column, the phenotype as a child should open it)."""
    if arguments['phenotype_table'] is not None:
        names = arguments['trait_columns']
        if names is None:
            from .io import load_table
            table = load_table(arguments['phenotype_table'])
            identifier = arguments['sample_id_column']
            if identifier not in table.columns:
                identifier = next((c for c in ('sampleid', 'sample_id', 'ID', 'id') if c in table.columns), identifier)
            names = [c for c in table.columns if c not in {'FID', 'fid', identifier}]
        return [str(name) for name in names], None
    phenotype = arguments['phenotype']
    if isinstance(phenotype, (str, Path)) and Path(phenotype).suffix.lower() == '.npy':
        count = np.load(phenotype, mmap_mode='r', allow_pickle=False).shape[1]
    else:
        # A text panel or an array: one .npy every child memory-maps.
        from .io import load_array
        values = load_array(phenotype) if isinstance(phenotype, (str, Path)) else np.asarray(phenotype)
        phenotype = work / 'phenotype.npy'
        np.save(phenotype, values, allow_pickle=False)
        count = values.shape[1]
    names = arguments['trait_columns']
    if names is not None and len(names) != count:
        raise ValueError('trait_columns must name every phenotype column')
    return ([str(name) for name in names] if names is not None else [f'trait_{i}' for i in range(count)]), phenotype


def _partitions(n_traits, n_devices):
    """Contiguous column ranges, one per device, sizes differing by at most one."""
    used = min(n_traits, n_devices)
    if used < 1:
        raise ValueError('At least one phenotype required')
    size, extra = divmod(n_traits, used)
    bounds, first = [], 0
    for index in range(used):
        last = first + size + (index < extra)
        bounds.append((first, last))
        first = last
    return bounds


def _child_arguments(arguments, phenotype, names, first, last, output_dir, readers):
    child = {key: value for key, value in arguments.items() if not key.startswith('_')}
    child.update(device='cuda:0', trait_devices=None, trait_workers='thread', output_dir=str(output_dir),
                 reader_workers=readers, trait_columns=names[first:last])
    if phenotype is not None:
        child.update(phenotype=str(phenotype), _phenotype_columns=(first, last))
    if isinstance(child['covariates'], np.ndarray):
        path = output_dir.parent / 'covariates.npy'
        if not path.exists():
            np.save(path, child['covariates'], allow_pickle=False)
        child['covariates'] = str(path)
    return dict(arguments=child, columns=(first, last))


def _physical(device):
    """The CUDA_VISIBLE_DEVICES entry a child needs to see `device` as cuda:0."""
    index = int(device.split(':')[1]) if ':' in device else 0
    visible = os.environ.get('CUDA_VISIBLE_DEVICES')
    if visible is None:
        return str(index)
    entries = [entry.strip() for entry in visible.split(',') if entry.strip()]
    if index >= len(entries):
        raise ValueError(f'{device} is not visible (CUDA_VISIBLE_DEVICES={visible})')
    return entries[index]


def _launch(index, device, child, report_kept, work):
    spec = work / f'part_{index:02d}.spec'
    spec.write_bytes(pickle.dumps(dict(arguments=child['arguments'], report_kept=report_kept)))
    spec.chmod(0o600)
    log = work / f'part_{index:02d}.log'
    # The child imports this package from where the parent did.
    package_root = str(Path(__file__).resolve().parents[1])
    search = [package_root] + [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]
    environment = dict(os.environ, CUDA_VISIBLE_DEVICES=_physical(device), PYTHONPATH=os.pathsep.join(search))
    read, write = os.pipe()
    with log.open('wb') as handle:
        # Its own session: stopping a child stops its reader and writer threads' process group.
        process = subprocess.Popen([sys.executable, '-m', 'torchgwas.trait_processes', str(spec), str(write)],
                                   stdin=subprocess.PIPE, stdout=handle, stderr=subprocess.STDOUT,
                                   env=environment, pass_fds=(write,), start_new_session=True)
    os.close(write)
    os.set_blocking(read, False)
    return dict(index=index, device=device, process=process, channel=read, buffer=b'', log=log,
                columns=child['columns'], output=Path(child['arguments']['output_dir']), kept=None, done=False)


def _wait(children, report_kept):
    """Relay the kept-trait total and wait for every child; a failure stops the rest."""
    selector = selectors.DefaultSelector()
    for child in children:
        selector.register(child['channel'], selectors.EVENT_READ, child)
        child['closed'] = False
    sent = not report_kept
    try:
        while True:
            for key, _ in selector.select(timeout=0.5):
                child = key.data
                data = os.read(child['channel'], 4096)
                if not data:
                    selector.unregister(child['channel'])
                    child['closed'] = True
                    continue
                child['buffer'] += data
                *lines, child['buffer'] = child['buffer'].split(b'\n')
                for line in lines:
                    word, _, value = line.decode().partition(' ')
                    if word == 'kept':
                        child['kept'] = int(value)
                    elif word == 'done':
                        child['done'] = True
            if not sent and all(child['kept'] is not None for child in children):
                total = sum(child['kept'] for child in children)
                for child in children:
                    child['process'].stdin.write(f'total {total}\n'.encode())
                    child['process'].stdin.flush()
                sent = True
            for child in children:
                code = child['process'].poll()
                # A child's messages are all read once its channel is closed.
                if code == 0 and child['closed'] and child['kept'] is None and report_kept:
                    # It scanned without the shared threshold (a path that never reported its
                    # count): its pairs used its own partition's Bonferroni.
                    raise RuntimeError(f"trait process {child['index']} finished without reporting its "
                                       f"trait count (log {child['log']})")
                if code not in (None, 0) or (code == 0 and child['closed'] and not child['done']):
                    tail = child['log'].read_text(errors='replace').splitlines()[-25:]
                    raise RuntimeError(f"trait process {child['index']} on {child['device']} exited with {code} "
                                       f"(log {child['log']}):\n" + '\n'.join(tail))
            if all(child['closed'] and child['process'].poll() == 0 for child in children):
                return
    finally:
        selector.close()
        for child in children:
            os.close(child['channel'])
            if child['process'].stdin is not None:
                child['process'].stdin.close()


def _stop(process, grace=10.0):
    if process.poll() is not None:
        return
    for sent in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(process.pid, sent)
        except ProcessLookupError:
            return
        try:
            process.wait(grace)
            return
        except subprocess.TimeoutExpired:
            continue


# -- merging -------------------------------------------------------------------

def _merge(arguments, children, output_dir):
    from .sumstats import read_manifest
    runs = [json.loads((child['output'] / 'run.json').read_text()) for child in children]
    qcs = [json.loads((child['output'] / 'qc.json').read_text()) for child in children]
    kept = [len(_kept_indices(qc)) for qc in qcs]
    offsets = np.concatenate([[0], np.cumsum(kept)])[:-1].tolist()
    names = [name for run in runs for name in run['trait_columns']]
    thresholds = {run['significance_threshold'] for run in runs}
    if len(thresholds) != 1:
        raise RuntimeError(f'trait processes used different significance thresholds: {sorted(thresholds)}')
    stores = [child['output'] / 'sumstats' for child in children]
    target = output_dir / 'sumstats'
    if target.exists():
        shutil.rmtree(target)
    fsync = bool(arguments['sumstats_fsync'])
    started = time.perf_counter()
    manifest = read_manifest(stores[0]) if (stores[0] / 'manifest.json').exists() else None
    if manifest is None:
        summary = dict(store=None)
    elif manifest.get('format') == 'torchgwas-indexed-sumstats':
        merge = _merge_min_p if manifest.get('reduction') == 'min-p' else _merge_pairs
        summary = merge(stores, offsets, names, target, fsync)
    else:
        summary = _merge_tiles(stores, offsets, names, target, fsync, [child['device'] for child in children])
    summary['merge_seconds'] = time.perf_counter() - started
    run = dict(runs[0])
    run.update(
        trait_workers='process', trait_devices=[child['device'] for child in children],
        device_used=[child['device'] for child in children], gpu_name=[run.get('gpu_name') for run in runs],
        trait_columns=names, phenotype_shape=[runs[0]['phenotype_shape'][0], sum(c['columns'][1] - c['columns'][0]
                                                                                for c in children)],
        reader_workers=sum(int(run['reader_workers'] or 0) for run in runs),
        n_result_rows=summary.get('rows', sum(int(run.get('n_result_rows') or 0) for run in runs)),
        sumstats_write=summary,
        trait_processes=[dict(device=child['device'], columns=list(child['columns']), kept_traits=count,
                              runtime_seconds=run.get('runtime_seconds'), phase_seconds=run.get('phase_seconds'),
                              sumstats_write=run.get('sumstats_write'))
                         for child, run, count in zip(children, runs, kept)])
    return dict(run=run, qc=_merge_qc(qcs, [child['columns'][0] for child in children]))


def _kept_indices(qc):
    if 'phenotype_kept_column_indices' in qc:
        return list(qc['phenotype_kept_column_indices'])
    return list(range(int(qc['phenotype_columns_input'])))


def _merge_qc(qcs, firsts):
    """One QC summary: per-trait lists joined, phenotype counts summed, the rest shared."""
    merged = {}
    for key in dict.fromkeys(key for qc in qcs for key in qc):
        values = [qc.get(key) for qc in qcs]
        per_trait = all(isinstance(value, list) and len(value) in (qc.get('phenotype_columns_kept'),
                                                                 qc.get('phenotype_columns_input'))
                        for value, qc in zip(values, qcs))
        if key == 'phenotype_kept_column_indices':
            merged[key] = [int(i) + first for qc, first in zip(qcs, firsts) for i in _kept_indices(qc)]
        elif key.startswith('phenotype_') and per_trait:
            merged[key] = [item for value in values for item in value]
        elif (key.startswith('phenotype_') or key == 'dropped_phenotype_columns') and \
                all(isinstance(value, (int, float)) and not isinstance(value, bool) for value in values):
            merged[key] = sum(values)
        elif all(value == values[0] for value in values):
            merged[key] = values[0]
        else:
            merged[key] = values[0]
            merged[f'{key}_by_partition'] = values
    if 'phenotype_kept_column_indices' not in merged and merged.get('phenotype_columns_kept') != \
            merged.get('phenotype_columns_input'):
        merged['phenotype_kept_column_indices'] = [i + first for qc, first in zip(qcs, firsts)
                                                   for i in _kept_indices(qc)]
    return merged


def _copy_variant_files(source, target):
    for name in ('variant_ids.npy', 'variant_metadata.npz', 'variant_ids.txt'):
        if (source / name).exists():
            shutil.copy2(source / name, target / name)


def _merge_pairs(stores, offsets, names, target, fsync):
    """Pair stores: trait_index rebased per partition, then one (variant, trait) order."""
    from .sumstats import write_manifest
    from .sumstats_indexed import _order_parts, open_indexed_sumstats
    target.mkdir(parents=True)
    parts, template = [], None
    for store, offset in zip(stores, offsets):
        manifest, values = open_indexed_sumstats(store)
        template = template or manifest
        for part, rows in zip(manifest['parts'], values):
            rows['trait_index'] = rows['trait_index'] + offset
            name = f'part_{len(parts):06d}.npz'
            _save_part(target / name, rows, fsync)
            parts.append(dict(file=name, rows=int(part['rows']), variant_range=list(part['variant_range'])))
    parts, ordering = _order_parts(target, parts, fsync, merge=True)
    manifest = dict(template, traits=names, shape=[template['shape'][0], len(names)], parts=parts,
                    rows=sum(part['rows'] for part in parts))
    manifest.pop('row_order', None)
    if ordering['globally_ordered']:
        manifest['row_order'] = 'variant_index then trait_index, across parts in manifest order'
    _copy_variant_files(stores[0], target)
    write_manifest(target, manifest, fsync=fsync)
    return dict(store='indexed', directory=str(target), rows=manifest['rows'], ordering=ordering)


def _merge_min_p(stores, offsets, names, target, fsync):
    """min-p stores: per variant, the partition winner with the largest -log10 P."""
    from .sumstats import write_manifest
    from .sumstats_indexed import open_indexed_sumstats
    target.mkdir(parents=True)
    columns, template = {}, None
    for store, offset in zip(stores, offsets):
        manifest, values = open_indexed_sumstats(store)
        template = template or manifest
        for rows in values:
            rows['trait_index'] = rows['trait_index'] + offset
            for key, value in rows.items():
                columns.setdefault(key, []).append(value)
    rows = {key: np.concatenate(value) for key, value in columns.items()}
    variant, logp = rows['variant_index'], rows['neg_log10_p'].astype(np.float64)
    # Largest -log10 P first within a variant (NaN last), then the lower trait index.
    order = np.lexsort((rows['trait_index'], np.where(np.isnan(logp), np.inf, -logp), variant))
    first = np.ones(order.size, dtype=bool)
    first[1:] = variant[order][1:] != variant[order][:-1]
    winners = order[first]
    rows = {key: value[winners] for key, value in rows.items()}
    parts = []
    for start in range(0, len(winners), MERGED_PART_ROWS):
        piece = {key: value[start:start + MERGED_PART_ROWS] for key, value in rows.items()}
        name = f'part_{len(parts):06d}.npz'
        _save_part(target / name, piece, fsync)
        span = piece['variant_index']
        parts.append(dict(file=name, rows=int(span.size), variant_range=[int(span[0]), int(span[-1]) + 1]))
    manifest = dict(template, traits=names, shape=[template['shape'][0], len(names)], parts=parts,
                    rows=int(len(winners)))
    _copy_variant_files(stores[0], target)
    write_manifest(target, manifest, fsync=fsync)
    return dict(store='indexed', directory=str(target), rows=int(len(winners)))


def _merge_tiles(stores, offsets, names, target, fsync, devices):
    """Full output: the children's tiles moved under one trait-tiled manifest."""
    from .sumstats import read_manifest, write_manifest
    target.mkdir(parents=True)
    tiles, trait_df, passes = [], [], 0
    for store, offset, device in zip(stores, offsets, devices):
        manifest = read_manifest(store)
        if manifest.get('layout') == 'trait_tiles':
            children = [(store / tile['directory'], tile['trait_range']) for tile in manifest['tiles']]
            passes += int(manifest.get('genotype_passes', len(children)))
        else:
            children = [(store, [0, manifest['shape'][1]])]
            passes += 1
        # Per-trait df (missing phenotypes) are listed; otherwise each tile has a variant df sidecar.
        trait_df.append(manifest['df'] if isinstance(manifest.get('df'), list) else None)
        for path, (first, last) in children:
            child = read_manifest(path)
            directory = f'traits_{first + offset:08d}_{last + offset:08d}'
            shutil.move(str(path), str(target / directory))
            tiles.append(dict(trait_range=[first + offset, last + offset], directory=directory, device=device))
            arrays, shape, n_samples = list(child['arrays']), child['shape'], child['n_samples']
    manifest = dict(format='torchgwas-binary-sumstats', version=2, layout='trait_tiles',
                    shape=[int(shape[0]), len(names)], dtype='float32', byte_order='little', arrays=arrays,
                    traits=names, n_samples=int(n_samples), tiles=tiles,
                    excluded_convention='NaN marks an excluded variant in stored beta/t arrays',
                    genotype_passes=passes,
                    scope='Each tile is a sequential variant-major store. No full-matrix assembly.')
    if all(df is not None for df in trait_df):
        manifest.update(df=[int(value) for df in trait_df for value in df], p_value="neg_log10_p is exact at each pair's complete-case df; df lists "
                        "each trait's observed samples less rank and genotype, which missing calls lower further")
    else:
        manifest.update(df=dict(layout='per_tile'),
                        p_value='not stored; two-sided Student t using the matching tile df sidecar')
    write_manifest(target, manifest, fsync=fsync)
    return dict(store='trait_tiles', directory=str(target), tiles=len(tiles),
                cells=int(shape[0]) * len(names))


def _save_part(path, rows, fsync):
    with path.open('wb') as handle:
        np.savez(handle, **rows)
        handle.flush()
        if fsync:
            os.fsync(handle.fileno())


# -- the child -----------------------------------------------------------------

def _child(spec_path, channel_fd):
    spec = pickle.loads(Path(spec_path).read_bytes())
    channel = os.fdopen(int(channel_fd), 'w', buffering=1)
    arguments = spec['arguments']
    if spec['report_kept']:
        def total(kept):
            channel.write(f'kept {int(kept)}\n')
            channel.flush()
            line = sys.stdin.readline()
            if not line.startswith('total '):
                raise RuntimeError('trait process: no trait total from the parent')
            return int(line.split()[1])
        arguments['_significance_traits'] = total
    from torchgwas.api import run_linear_gwas
    run_linear_gwas(**arguments)
    channel.write('done\n')
    channel.close()


if __name__ == '__main__':
    _child(*sys.argv[1:3])
