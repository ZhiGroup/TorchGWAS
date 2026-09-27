"""Does the empirical autotune pick a good layout and chunk size?

Compares run_linear_gwas(autotune=True) with fixed layouts (1/2/4/8 GPUs,
phenotype tiles or variant shards, fixed chunk) on one synthetic dataset,
fresh process per observation, interleaved repeats. Every configuration's
output is checked against the first fixed configuration. Output-inclusive API
wall time is the score. Genotypes come from `plink2 --dummy`; phenotypes are
Gaussian with a few planted signals so significant output is nonempty.
"""
import argparse, hashlib, json, os, random, statistics, subprocess, sys, time
from pathlib import Path

import numpy as np

PLINK2 = '/data/zxie3/torchgwas_bench/accuracy_3servers_20260916/assets/plink2'
ENVIRONMENT = dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
                   TORCHGWAS_SCAN_PROFILE='0', TORCHGWAS_BLOCKING_EVENTS='1', NUMPY_MADVISE_HUGEPAGE='0',
                   OMP_NUM_THREADS='4', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_WAIT_POLICY='PASSIVE')


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=1, allow_nan=False); stream.write('\n')
        stream.flush(); os.fsync(stream.fileno())


def genotype_files(path, fmt):
    """(data files to warm, sample count, variant count) of an existing genotype source."""
    import struct
    path = Path(path)
    if fmt == 'pgen':
        count = lambda p: sum(1 for line in p.open() if not line.startswith('#'))
        return [path], count(path.with_suffix('.psam')), count(path.with_suffix('.pvar'))
    if fmt == 'plink':
        count = lambda p: sum(1 for _ in p.open())
        return [path], count(path.with_suffix('.fam')), count(path.with_suffix('.bim'))
    if fmt == 'bgen':
        with path.open('rb') as stream:
            # First 4 bytes: offset; header block: length, variants, samples.
            _, _, variants, samples = struct.unpack('<IIII', stream.read(16))
        return [path], samples, variants
    if fmt == 'zstd':
        index = np.load(f'{path}.idx.npz')
        return [Path(f'{path}.zst')], int(index['nsamp']), int(index['nsnp'])
    raise ValueError(f'unknown genotype format {fmt}')


def make(data, samples, variants, traits, covariates, seed, genotype_from=None,
         genotype_path=None, genotype_format=None, pgen_mode=None):
    data.mkdir(parents=True, exist_ok=True)
    if genotype_path is not None:
        # An existing source in any format: phenotypes for its samples, in file order.
        files, samples, variants = genotype_files(genotype_path, genotype_format)
        rng = np.random.default_rng(seed)
        np.save(data/'covariates.npy', rng.normal(size=(samples, covariates)).astype(np.float32))
        np.save(data/'phenotype.npy', rng.normal(size=(samples, traits)).astype(np.float32))
        return dict(samples=samples, variants=variants, traits=traits, covariates=covariates, seed=seed,
                    genotype=str(genotype_path), genotype_format=genotype_format,
                    pgen_mode=pgen_mode or 'hardcall', genotype_files=[str(f) for f in files])
    prefix = data/'input'
    if genotype_from is not None and not prefix.with_suffix('.pgen').exists():
        for suffix in ('.pgen', '.pvar', '.psam'):
            (data/('input'+suffix)).symlink_to(Path(genotype_from)/('input'+suffix))
    if not prefix.with_suffix('.pgen').exists():
        subprocess.run([PLINK2, '--dummy', str(samples), str(variants), '0', '0', 'acgt',
                        '--make-pgen', '--out', str(prefix), '--threads', '8', '--seed', str(seed)],
                       check=True, stdout=subprocess.DEVNULL)
    rng = np.random.default_rng(seed)
    c = rng.normal(size=(samples, covariates)).astype(np.float32)
    y = rng.normal(size=(samples, traits)).astype(np.float32)
    np.save(data/'covariates.npy', c)
    np.save(data/'phenotype.npy', y)
    return dict(samples=samples, variants=variants, traits=traits, covariates=covariates, seed=seed)


def child(data, out, config):
    import torch
    from torchgwas.api import run_linear_gwas
    kwargs = dict(config['kwargs'])
    os.environ.update(config.get('env', {}))
    started = time.perf_counter()
    manifest = json.loads((data/'manifest.json').read_text())
    source = dict(genotype=manifest.get('genotype', str(data/'input.pgen')),
                  genotype_format=manifest.get('genotype_format', 'auto'))
    if source['genotype_format'] in ('auto', 'pgen'):
        source['pgen_mode'] = config.get('pgen_mode', manifest.get('pgen_mode', 'hardcall'))
    # One metadata cache per dataset: without it every run re-parses the
    # 8.09M-line pvar (8-14 s of api_seconds at full scale).
    if 'genotype_cache_dir' not in kwargs:
        (data / 'metadata_cache').mkdir(exist_ok=True)
        source['genotype_cache_dir'] = str(data / 'metadata_cache')
    result = run_linear_gwas(**source, phenotype=data/'phenotype.npy',
                             covariates=data/'covariates.npy',
                             compute_dtype='float32',
                             output_dir=out, **kwargs)
    api = time.perf_counter()-started
    meta = result.run_metadata
    timing = meta.get('sumstats_write') or {}
    return dict(api_seconds=api, executor_seconds=timing.get('setup_scan_and_write_seconds'),
                rows=meta['n_result_rows'], trait_block=meta.get('trait_block'),
                trait_devices=meta.get('trait_devices'), variant_devices=meta.get('variant_devices'),
                reader_workers=meta.get('reader_workers'), prefetch_chunks=meta.get('prefetch_chunks'),
                autotune=meta.get('autotune'), shared_decode=meta.get('shared_decode'),
                jagwas_projection=timing.get('jagwas_projection'),
                sumstats_write={key: value for key, value in timing.items()
                                if isinstance(value, (int, float, str, bool)) or key in ('ordering', 'publication')},
                numpy_hugepage_advice=meta.get('numpy_hugepage_advice'),
                process_cpu=list(os.times()[:2]))


def digest_output(directory, mode):
    """Order-free fingerprint of the association output, for cross-checking."""
    from torchgwas.sumstats import open_binary_sumstats
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    if mode == 'full':
        beta, t, _ = open_binary_sumstats(directory/'sumstats')
        return dict(kind='dense', t=np.asarray(t, dtype=np.float64))
    _, parts = open_indexed_sumstats(directory/'sumstats')
    rows = {}
    for part in parts:
        value = part['chi2'] if mode == 'jagwas' else part['t_stat']
        keys = part['variant_index'] if mode == 'jagwas' else zip(part['variant_index'], part['trait_index'])
        for key, v in zip(keys, value):
            rows[key if mode == 'jagwas' else (int(key[0]), int(key[1]))] = float(v)
    return dict(kind='indexed', rows=rows)


def compare(a, b):
    if a['kind'] == 'dense':
        finite = np.isfinite(a['t']) & np.isfinite(b['t'])
        return dict(max_abs=float(np.max(np.abs(a['t'][finite]-b['t'][finite]), initial=0.)),
                    same_keys=a['t'].shape == b['t'].shape and bool(np.array_equal(np.isnan(a['t']), np.isnan(b['t']))))
    keys = a['rows'].keys() & b['rows'].keys()
    diff = max((abs(a['rows'][k]-b['rows'][k]) for k in keys), default=0.)
    return dict(max_abs=diff, same_keys=a['rows'].keys() == b['rows'].keys(),
                only_a=len(a['rows'].keys()-b['rows'].keys()), only_b=len(b['rows'].keys()-a['rows'].keys()))


def overhead_configs(traits, devices):
    """Same 2-GPU tiled layout and 1,024 chunks; vary only what autotune adds."""
    layout = dict(trait_block=-(-traits//2), trait_devices=devices[:2], reduce='significant',
                  significance_threshold=1e-5)
    never = dict(min_job_seconds=1e9)  # tuner attached, never trials
    return [
        dict(name='base_r8_p4', kwargs=dict(layout, chunk_size=1024, reader_workers=8, prefetch_chunks=4)),
        dict(name='readers16_p8', kwargs=dict(layout, chunk_size=1024, reader_workers=16, prefetch_chunks=8)),
        dict(name='observer_cap4096', kwargs=dict(layout, reader_workers=8, prefetch_chunks=4, autotune=True,
                                                  autotune_options=dict(never, chunk_sizes=[1024, 4096]))),
        dict(name='observer_cap1024', kwargs=dict(layout, reader_workers=8, prefetch_chunks=4, autotune=True,
                                                  autotune_options=dict(never, chunk_sizes=[512, 1024]))),
    ]


def shared_configs(traits, devices, chunk, counts):
    """Significant-pair tiles with and without one shared decode pass."""
    extra = dict(reduce='significant', significance_threshold=1e-5)
    rows = [dict(name='fixed_1gpu', kwargs=dict(extra, device=devices[0], chunk_size=chunk,
                                                 reader_workers=4, prefetch_chunks=4))]
    for n in [c for c in counts if 1 < c <= len(devices)]:
        layout = dict(extra, trait_block=-(-traits//n), trait_devices=devices[:n], chunk_size=chunk,
                      reader_workers=4*n, prefetch_chunks=4)
        rows.append(dict(name=f'tiles{n}_separate', kwargs=layout, env=dict(TORCHGWAS_SHARED_DECODE='0')))
        rows.append(dict(name=f'tiles{n}_shared', kwargs=layout, env=dict(TORCHGWAS_SHARED_DECODE='1')))
    return rows


def cache_configs(traits, devices, chunk, counts):
    """More phenotype tiles than GPUs (rounds): host genotype cache on or off."""
    rows = []
    for n in [c for c in counts if c <= len(devices)]:
        layout = dict(reduce='significant', significance_threshold=1e-5, trait_block=-(-traits//(4*n)),
                      trait_devices=devices[:n], chunk_size=chunk, reader_workers=4*n, prefetch_chunks=4)
        for cache in ('1', '0'):
            rows.append(dict(name=f'rounds4_gpus{n}_cache_{"on" if cache == "1" else "off"}', kwargs=layout,
                             env=dict(TORCHGWAS_GENOTYPE_CACHE=cache, TORCHGWAS_SIGNIFICANCE_BACKEND='device')))
    return rows


def tuner_configs(traits, devices, chunk, counts):
    """One GPU, significant pairs: fixed chunk sizes against both autotuners (layout fixed)."""
    extra = dict(reduce='significant', significance_threshold=1e-5, device=devices[0],
                 reader_workers=4, prefetch_chunks=4)
    env = dict(TORCHGWAS_SIGNIFICANCE_BACKEND='device')
    rows = [dict(name=f'fixed_{size}', env=env, kwargs=dict(extra, chunk_size=size)) for size in (512, 1024, 2048)]
    for tuner in ('model', 'segments'):
        rows.append(dict(name=f'auto_{tuner}', env=env, kwargs=dict(
            extra, autotune=True, autotune_options=dict(tuner=tuner, chunk_sizes=[512, 1024, 2048],
                                                        devices=[devices[0]], min_job_seconds=0))))
    return rows


def split_configs(traits, devices, chunk, counts):
    """Significant pairs split by phenotype (tiles, one shared decode) or by variant (shards)."""
    extra = dict(reduce='significant', significance_threshold=1e-5)
    rows = [dict(name='fixed_1gpu', kwargs=dict(extra, device=devices[0], chunk_size=chunk,
                                                 reader_workers=4, prefetch_chunks=4))]
    for n in [c for c in counts if 1 < c <= len(devices)]:
        rows.append(dict(name=f'tiles{n}_shared', env=dict(TORCHGWAS_SHARED_DECODE='1', TORCHGWAS_GPU_FANOUT='pcie'),
                         kwargs=dict(extra, trait_block=-(-traits//n), trait_devices=devices[:n], chunk_size=chunk,
                                     reader_workers=4*n, prefetch_chunks=4)))
        rows.append(dict(name=f'vshards{n}', kwargs=dict(extra, variant_devices=devices[:n], chunk_size=chunk,
                                                         reader_workers=4*n, prefetch_chunks=4)))
    # Autotune over the same GPUs: does the priced split pick the faster layout?
    used = devices[:max([1] + [c for c in counts if c <= len(devices)])]
    rows.append(dict(name='autotune', kwargs=dict(extra, autotune=True, autotune_options=dict(devices=used))))
    return rows


def selection_configs(traits, devices, chunk, counts):
    """Host (full t to CPU, NumPy predicate) versus device (CUDA nonzero) pair selection."""
    extra = dict(reduce='significant', significance_threshold=1e-5)
    layouts = [('1gpu', dict(device=devices[0], reader_workers=4))]
    for n in [c for c in counts if 1 < c <= len(devices)]:
        layouts.append((f'vshards{n}', dict(variant_devices=devices[:n], reader_workers=4*n)))
        layouts.append((f'tiles{n}', dict(trait_block=-(-traits//n), trait_devices=devices[:n], reader_workers=4*n)))
    rows = []
    for name, layout in layouts:
        for backend in ('host', 'device'):
            env = dict(TORCHGWAS_SIGNIFICANCE_BACKEND=backend)
            if name.startswith('tiles'):
                env.update(TORCHGWAS_SHARED_DECODE='1', TORCHGWAS_GPU_FANOUT='pcie')
            rows.append(dict(name=f'{name}_{backend}', env=env,
                             kwargs=dict(extra, **layout, chunk_size=chunk, prefetch_chunks=4)))
    used = devices[:max([1] + [c for c in counts if c <= len(devices)])]
    rows.append(dict(name='autotune', kwargs=dict(extra, autotune=True, autotune_options=dict(devices=used))))
    return rows


def fanout_configs(traits, devices, chunk, counts):
    """Shared-decode tiles: per-GPU PCIe copies versus one PCIe copy plus NVLink peer copies."""
    extra = dict(reduce='significant', significance_threshold=1e-5)
    rows = []
    for n in [c for c in counts if 1 < c <= len(devices)]:
        layout = dict(extra, trait_block=-(-traits//n), trait_devices=devices[:n], chunk_size=chunk,
                      reader_workers=4*n, prefetch_chunks=4)
        rows.append(dict(name=f'tiles{n}_separate', kwargs=layout, env=dict(TORCHGWAS_SHARED_DECODE='0')))
        for fanout in ('pcie', 'nvlink'):
            rows.append(dict(name=f'tiles{n}_shared_{fanout}', kwargs=layout,
                             env=dict(TORCHGWAS_SHARED_DECODE='1', TORCHGWAS_GPU_FANOUT=fanout)))
    return rows


def device_selection_configs(traits, devices, chunk, counts):
    """Significant pairs with device selection only: 1 GPU, shards, tiles and autotune (full-scale runs)."""
    extra = dict(reduce='significant', significance_threshold=1e-5, chunk_size=chunk, prefetch_chunks=4)
    env = dict(TORCHGWAS_SIGNIFICANCE_BACKEND='device')
    rows = [dict(name='1gpu_device', env=env, kwargs=dict(extra, device=devices[0], reader_workers=4))]
    for n in [c for c in counts if 1 < c <= len(devices)]:
        rows.append(dict(name=f'vshards{n}_device', env=env,
                         kwargs=dict(extra, variant_devices=devices[:n], reader_workers=4*n)))
        rows.append(dict(name=f'tiles{n}_device', env=dict(env, TORCHGWAS_SHARED_DECODE='1', TORCHGWAS_GPU_FANOUT='pcie'),
                         kwargs=dict(extra, trait_block=-(-traits//n), trait_devices=devices[:n], reader_workers=4*n)))
    used = devices[:max([1] + [c for c in counts if c <= len(devices)])]
    rows.append(dict(name='autotune', kwargs=dict(reduce='significant', significance_threshold=1e-5, autotune=True,
                                                  autotune_options=dict(devices=used))))
    return rows


def projection_configs(traits, devices, chunk, counts):
    """JAGWAS on one GPU and variant shards (the dense projection was removed 2026-09-26)."""
    extra = dict(reduce='jagwas', chunk_size=chunk, prefetch_chunks=4)
    rows = []
    for n in [c for c in counts if c <= len(devices)]:
        layout = dict(device=devices[0]) if n == 1 else dict(variant_devices=devices[:n])
        rows.append(dict(name=f'jagwas_{n}gpu', kwargs=dict(extra, reader_workers=4*n, **layout)))
    used = devices[:max([1] + [c for c in counts if c <= len(devices)])]
    rows.append(dict(name='autotune', kwargs=dict(reduce='jagwas', autotune=True, autotune_options=dict(devices=used))))
    return rows


def configs(mode, traits, devices, chunk, counts=(1, 2, 4, 8)):
    """Autotune plus fixed layouts over 1..len(devices) GPUs."""
    if mode == 'projection':
        return projection_configs(traits, devices, chunk, counts)
    if mode == 'device_selection':
        return device_selection_configs(traits, devices, chunk, counts)
    if mode == 'overhead':
        return overhead_configs(traits, devices)
    if mode == 'shared':
        return shared_configs(traits, devices, chunk, counts)
    if mode == 'fanout':
        return fanout_configs(traits, devices, chunk, counts)
    if mode == 'split':
        return split_configs(traits, devices, chunk, counts)
    if mode == 'tuner':
        return tuner_configs(traits, devices, chunk, counts)
    if mode == 'cache':
        return cache_configs(traits, devices, chunk, counts)
    if mode == 'selection':
        return selection_configs(traits, devices, chunk, counts)
    base = dict(chunk_size=chunk)
    extra = dict(reduce=mode) if mode != 'full' else {}
    if mode == 'significant':
        extra['significance_threshold'] = 1e-5
    rows = [dict(name='autotune', kwargs=dict(autotune=True, **extra))]
    counts = [n for n in counts if n <= len(devices)]
    for n in counts:
        gpus = devices[:n]
        if n == 1:
            layout = dict(device=gpus[0])
        elif mode == 'significant' or (mode == 'full' and traits > 4096):
            layout = dict(trait_block=-(-traits//n), trait_devices=gpus)
        else:
            layout = dict(variant_devices=gpus)
        rows.append(dict(name=f'fixed_{n}gpu', kwargs=dict(base, reader_workers=4*n, prefetch_chunks=4,
                                                           **layout, **extra)))
    return rows


def observe(root, data, mode, devices, chunk, repeats, seed, counts=(1, 2, 4, 8), pgen_mode='hardcall', only=None):
    manifest = json.loads((data/'manifest.json').read_text())
    rows = configs(mode, manifest['traits'], devices, chunk, counts)
    if only:
        rows = [row for row in rows if row['name'] in only]
    if pgen_mode != 'hardcall':  # float32 dosage transport: 16x the bytes of 2-bit codes
        rows = [dict(row, pgen_mode=pgen_mode) for row in rows]
    rng = random.Random(seed); schedule = []
    for r in range(repeats):
        names = [c['name'] for c in rows]; rng.shuffle(names); schedule += [(n, r) for n in names]
    plan = dict(mode=mode, devices=devices, chunk=chunk, repeats=repeats, seed=seed, configs=rows,
                schedule=schedule, manifest=manifest)
    if (root/'plan.json').exists():
        assert json.loads((root/'plan.json').read_text()) == plan, 'Existing plan differs'
    else:
        save(root/'plan.json', plan)
    # Warm the page cache once so the first scheduled run is not the only cold one.
    for name in manifest.get('genotype_files', [str(data/'input.pgen')]) + [data/'phenotype.npy', data/'covariates.npy']:
        with open(name, 'rb') as stream:
            while stream.read(1 << 24):
                pass
    for name, r in schedule:
        record = root/f'{name}_r{r}.json'
        if record.exists():
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(root), '--data', str(data),
                        '--name', name, '--repeat', str(r)], check=True, env={**os.environ, **ENVIRONMENT})
        print(json.loads(record.read_text())['summary'], flush=True)
    report(root, data)


def run_child(root, data, name, repeat):
    plan = json.loads((root/'plan.json').read_text())
    config = next(c for c in plan['configs'] if c['name'] == name)
    out = data/'outputs'/root.name/f'{name}_r{repeat}'
    assert not out.exists()
    row = child(data, out, config)
    summary = dict(name=name, repeat=repeat, api_seconds=round(row['api_seconds'], 2),
                   executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                   rows=row['rows'], trait_block=row['trait_block'], trait_devices=row['trait_devices'],
                   variant_devices=row['variant_devices'],
                   chunk=((row['autotune'] or {}).get('chunk') or {}).get('choice'))
    save(root/f'{name}_r{repeat}.json', dict(row, summary=summary, output=str(out)))


def report(root, data):
    plan = json.loads((root/'plan.json').read_text())
    records = {(n, r): json.loads((root/f'{n}_r{r}.json').read_text()) for n, r in plan['schedule']}
    reference_name = plan['configs'][0 if plan['mode'] == 'overhead' or len(plan['configs']) == 1 else 1]['name']
    kind = ('significant' if plan['mode'] in ('overhead', 'shared', 'fanout', 'split', 'selection', 'tuner', 'cache',
                                              'device_selection')
            else 'jagwas' if plan['mode'] == 'projection' else plan['mode'])
    reference = digest_output(Path(records[(reference_name, 0)]['output']), kind)
    checks = {f'{n}_r{r}': compare(reference, digest_output(Path(rec['output']), kind))
              for (n, r), rec in records.items()}
    table = {}
    for config in plan['configs']:
        group = [records[(config['name'], r)] for r in range(plan['repeats'])]
        table[config['name']] = dict(
            api_seconds=[round(g['api_seconds'], 2) for g in group],
            median_api_seconds=statistics.median(g['api_seconds'] for g in group),
            executor_seconds=[g['executor_seconds'] and round(g['executor_seconds'], 2) for g in group],
            layout=group[0]['summary'], chunk_choices=[g['summary']['chunk'] for g in group])
    best = min((n for n in table if n != 'autotune'), key=lambda n: table[n]['median_api_seconds'], default=None)
    result = dict(mode=plan['mode'], table=table, best_fixed=best,
                  autotune_vs_best_fixed=(table['autotune']['median_api_seconds']/table[best]['median_api_seconds']
                                          if 'autotune' in table and best is not None else None),
                  checks=checks, reference=reference_name, manifest=plan['manifest'])
    path = root/'report.json'
    if path.exists():
        path.unlink()
    save(path, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=1))
    print('max check difference', max(c['max_abs'] for c in checks.values()),
          'all keys equal', all(c['same_keys'] for c in checks.values()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['make', 'observe', 'child', 'report'])
    parser.add_argument('--root', type=Path); parser.add_argument('--data', type=Path)
    parser.add_argument('--mode', choices=['full', 'significant', 'jagwas', 'overhead', 'shared', 'fanout', 'split', 'selection', 'tuner', 'cache', 'projection', 'device_selection'], default='significant')
    parser.add_argument('--samples', type=int, default=20000); parser.add_argument('--variants', type=int, default=200000)
    parser.add_argument('--traits', type=int, default=8192); parser.add_argument('--covariates', type=int, default=10)
    parser.add_argument('--devices', nargs='+', default=[f'cuda:{i}' for i in range(8)])
    parser.add_argument('--chunk', type=int, default=1024); parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    parser.add_argument('--genotype-from', type=Path)
    parser.add_argument('--counts', type=int, nargs='+', default=[1, 2, 4, 8])
    parser.add_argument('--pgen-mode', default='hardcall', choices=['hardcall', 'dosage'])
    parser.add_argument('--only', nargs='+', help='run only these config names')
    parser.add_argument('--genotype-path', type=Path, help='make: an existing genotype source (any format)')
    parser.add_argument('--genotype-format', choices=['pgen', 'plink', 'bgen', 'zstd'])
    parser.add_argument('--source-pgen-mode', choices=['hardcall', 'dosage'])
    args = parser.parse_args()
    if args.stage == 'make':
        manifest = make(args.data, args.samples, args.variants, args.traits, args.covariates, args.seed,
                        args.genotype_from, args.genotype_path, args.genotype_format, args.source_pgen_mode)
        (args.data/'manifest.json').write_text(json.dumps(manifest))
        print(manifest)
    elif args.stage == 'observe':
        observe(args.root, args.data, args.mode, args.devices, args.chunk, args.repeats, args.seed, tuple(args.counts),
                args.pgen_mode, args.only)
    elif args.stage == 'child':
        run_child(args.root, args.data, args.name, args.repeat)
    else:
        report(args.root, args.data)


if __name__ == '__main__':
    main()
