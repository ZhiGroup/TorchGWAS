"""Stores record which input their variant indices refer to, not the IDs themselves.

Every indexed row carries `variant_index`, and a dense store's row i is
variant `variant_offset + i` of the input. The IDs, positions and alleles are
the genotype's own metadata, so a store records only where they come from and
a fingerprint of the variant list:

    variant_source = {format, genotype, load, variant_offset, n_variants,
                      id_digest, id_digest_scheme}

`store_variant_ids` maps a store back to IDs through the genotype (with its
metadata cache) and refuses a genotype whose variant list does not match.
Embedding the IDs (variant_ids.npy / variant_ids.txt) stays available for a
store meant to travel without its genotype. Full-scale cost of embedding:
356 MB of fixed-width IDs per indexed store (8.09M variants).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

DIGEST_SAMPLES = 4097
DIGEST_SCHEME = f'sha256 over the variant count and {DIGEST_SAMPLES} evenly spaced IDs, newline-separated'


def variant_digest(ids) -> str:
    """Cheap, order-sensitive fingerprint of a variant-ID list (ms at 8M IDs)."""
    ids = np.asarray(ids)
    n = int(ids.shape[0])
    positions = np.unique(np.linspace(0, n - 1, min(n, DIGEST_SAMPLES)).astype(np.int64)) if n else []
    digest = hashlib.sha256(f'{n}\n'.encode())
    digest.update('\n'.join(str(ids[i]) for i in positions).encode('utf-8', 'surrogatepass'))
    return digest.hexdigest()


def variant_source_record(*, genotype, genotype_format, marker_ids, variant_offset, n_variants, load=None):
    """The manifest entry for a store covering `n_variants` input variants from `variant_offset`."""
    if marker_ids is None:
        return None
    ids = np.asarray(marker_ids)[:int(n_variants)]
    return dict(format=genotype_format,
                genotype=None if genotype is None else str(Path(genotype).resolve()),
                load={key: str(value) for key, value in (load or {}).items() if value is not None},
                variant_offset=int(variant_offset), n_variants=int(n_variants),
                id_digest=variant_digest(ids), id_digest_scheme=DIGEST_SCHEME)


def store_variants(directory, genotype=None, *, genotype_cache_dir=None):
    """(variant IDs, variant metadata) for a store's rows (dense) or `variant_index` values (indexed).

    Embedded IDs (and variant_metadata.npz) are read directly. Otherwise the
    genotype is opened, from the recorded path unless `genotype` (a path or a
    loaded genotype with `marker_ids`) is given, and its variant list must
    match the recorded fingerprint. Metadata is the genotype's
    `variant_metadata` over the same range, or {} if it has none.
    """
    from .sumstats import read_manifest
    directory = Path(directory)
    manifest = read_manifest(directory)
    embedded = (np.load(directory / manifest['variant_ids'], allow_pickle=False)
                if isinstance(manifest.get('variant_ids'), str) else
                np.asarray((directory / 'variant_ids.txt').read_text().splitlines())
                if (directory / 'variant_ids.txt').exists() else None)
    if embedded is not None:
        metadata = {}
        if (directory / 'variant_metadata.npz').exists():
            with np.load(directory / 'variant_metadata.npz', allow_pickle=False) as handle:
                metadata = {key: handle[key] for key in handle.files}
        return embedded, metadata
    source = manifest.get('variant_source')
    if not source:
        raise ValueError(f'{directory} records neither variant IDs nor their source')
    if genotype is None or isinstance(genotype, (str, Path)):
        from .io import load_genotype
        path = genotype if genotype is not None else source['genotype']
        if path is None:
            raise ValueError('The store was written from an in-memory genotype; pass the genotype')
        genotype, _, marker_ids, _ = load_genotype(path, genotype_format=source['format'],
                                                   genotype_cache_dir=genotype_cache_dir,
                                                   **dict(source.get('load') or {}))
    else:
        marker_ids = getattr(genotype, 'marker_ids', None)
    if marker_ids is None:
        raise ValueError('The genotype has no variant IDs')
    first, count = int(source['variant_offset']), int(source['n_variants'])
    ids = np.asarray(marker_ids)[first:first + count]
    if ids.shape[0] != count or variant_digest(ids) != source['id_digest']:
        raise ValueError('The variant list of this genotype does not match the one the store was written from')
    metadata = getattr(genotype, 'variant_metadata', None) or {}
    return ids, {key: np.asarray(value)[first:first + count] for key, value in metadata.items()}


def store_variant_ids(directory, genotype=None, *, genotype_cache_dir=None):
    """Variant IDs of a store (see store_variants)."""
    return store_variants(directory, genotype, genotype_cache_dir=genotype_cache_dir)[0]
