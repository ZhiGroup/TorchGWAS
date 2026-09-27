import threading
import numpy as np
import pytest
from torchgwas.sumstats_indexed import (write_indexed_sumstats, open_indexed_sumstats,
    IndexedOutputPartition, PartitionedIndexedChunk)
from torchgwas.indexed_writer_progress import IndexedWriterProgress


def test_selected_indices_labels_and_t_only(tmp_path):
    # Global indices must survive selection and nonzero chunk starts.
    chunk = (2, 4, np.array([3, 2]), np.array([1, 0]),
             np.array([.5, -.25]), np.array([4., -3.]), 20)
    labels = ["plain", "tab\tlabel", "newline\nlabel", "unicode-β"]
    rows, _ = write_indexed_sumstats(tmp_path, labels, ["a", "b"], 24,
        [chunk], kind="significant", df=20, store_beta=False)
    manifest, parts = open_indexed_sumstats(tmp_path)
    part = next(parts)
    assert rows == manifest["rows"] == 2
    # Stored in (variant, trait) order whatever order selection produced.
    np.testing.assert_array_equal(part["variant_index"], [2, 3])
    np.testing.assert_array_equal(part["trait_index"], [0, 1])
    np.testing.assert_array_equal(part["t_stat"], [-3., 4.])
    assert "neg_log10_p" in part
    assert "beta" not in part
    assert np.load(tmp_path / "variant_ids.npy").tolist() == labels
    assert not list(tmp_path.glob("*.tsv*"))


def _all_rows(directory):
    manifest, parts = open_indexed_sumstats(directory)
    rows = [np.column_stack([p["variant_index"], p["trait_index"], p["t_stat"]]) for p in parts]
    return manifest, np.concatenate(rows) if rows else np.empty((0, 3))


def _pairs_chunk(start, end, pairs):
    vi, ti, t = (np.array(x) for x in zip(*pairs))
    return (start, end, vi, ti, t/10, t.astype(float), 20)


def test_variant_shards_arriving_out_of_order_publish_in_order(tmp_path):
    # Two shards interleave their chunks; parts are only re-listed, not rewritten.
    chunks = [_pairs_chunk(4, 6, [(5, 1, 3.), (4, 0, 2.)]), _pairs_chunk(0, 2, [(1, 0, 4.), (0, 1, 5.)]),
              _pairs_chunk(6, 8, [(7, 0, 6.)]), _pairs_chunk(2, 4, [(3, 1, 7.), (2, 1, 8.), (2, 0, 9.)])]
    rows, summary = write_indexed_sumstats(tmp_path, [f"v{i}" for i in range(8)], ["x", "y"], 24,
                                           chunks, kind="significant", df=20, coalesce_rows=1 << 18)
    manifest, table = _all_rows(tmp_path)
    assert rows == len(table) == 8 and summary["ordering"]["merged"] is False
    np.testing.assert_array_equal(table[:, :2], [[0, 1], [1, 0], [2, 0], [2, 1], [3, 1], [4, 0], [5, 1], [7, 0]])
    # Each shard's contiguous run is one part; the list is sorted, nothing is rewritten.
    assert [p["variant_range"] for p in manifest["parts"]] == [[0, 4], [4, 8]]


def test_chunks_coalesce_into_few_durable_parts(tmp_path, monkeypatch):
    import os
    import torchgwas.sumstats_indexed as indexed
    synced = []
    original = os.fsync
    monkeypatch.setattr(os, 'fsync', lambda fd: (synced.append(os.readlink(f'/proc/self/fd/{fd}')), original(fd))[1])
    # 40 one-variant chunks, one pair each, some empty: 30 rows in 6 parts of 5.
    chunks = [_pairs_chunk(v, v+1, [(v, 0, 3.)]) if v % 4 else (v, v+1, np.array([], int), np.array([], int),
              np.array([]), np.array([]), 20) for v in range(40)]
    rows, _ = write_indexed_sumstats(tmp_path, [f"v{i}" for i in range(40)], ["x"], 24, chunks,
                                     kind="significant", df=20, coalesce_rows=5)
    manifest, table = _all_rows(tmp_path)
    assert rows == 30 and len(manifest["parts"]) == 6
    assert sum(p.endswith('.npz') for p in synced) == 6  # one fsync per part, not per chunk
    np.testing.assert_array_equal(table[:, 0], [v for v in range(40) if v % 4])


def test_overlapping_phenotype_tiles_are_merged_into_global_order(tmp_path):
    # Two phenotype tiles (traits 0-1 and 2-3) with different chunk boundaries.
    tile_a = [_pairs_chunk(0, 3, [(0, 1, 1.), (2, 0, 2.)]), _pairs_chunk(3, 6, [(5, 1, 3.), (3, 0, 4.)])]
    tile_b = [_pairs_chunk(0, 2, [(1, 3, 5.), (0, 2, 6.)]), _pairs_chunk(2, 4, [(3, 2, 7.)]),
              _pairs_chunk(4, 6, [(5, 3, 8.), (4, 2, 9.)])]
    chunks = [tile_b[0], tile_a[0], tile_b[1], tile_a[1], tile_b[2]]
    rows, summary = write_indexed_sumstats(tmp_path, [f"v{i}" for i in range(6)], list("abcd"), 24,
                                           chunks, kind="significant", df=20, coalesce_rows=1 << 18)
    manifest, table = _all_rows(tmp_path)
    assert rows == manifest["rows"] == len(table) == 9 and summary["ordering"]["merged"] is True
    keys = [tuple(map(int, row[:2])) for row in table]
    assert keys == sorted(keys) and len(set(keys)) == 9
    assert dict(zip(keys, table[:, 2])) == {(0, 1): 1., (2, 0): 2., (5, 1): 3., (3, 0): 4., (1, 3): 5.,
                                            (0, 2): 6., (3, 2): 7., (5, 3): 8., (4, 2): 9.}
    assert not list(tmp_path.glob("part_*.npz"))  # superseded by the ordered parts
    assert sum(p["rows"] for p in manifest["parts"]) == 9


def test_threshold_across_chunks(tmp_path):
    # df=20, p <= .05 two-sided: |t| >= 2.086 passes.
    chunks = [(0, 2, np.ones((2, 2)), np.array([[2., 4.], [5., np.nan]]), None),
              (2, 3, np.ones((1, 2)), np.array([[6., 8.]]), None)]
    rows, _ = write_indexed_sumstats(tmp_path, ["a", "b", "c"], ["x", "y"],
        24, chunks, kind="filtered", df=20, p_value_threshold=.05, coalesce_rows=1 << 18)
    _, parts = open_indexed_sumstats(tmp_path)
    (part,) = list(parts)  # contiguous chunks coalesce into one part
    assert rows == 4
    np.testing.assert_array_equal(part["variant_index"], [0, 1, 2, 2])
    np.testing.assert_array_equal(part["t_stat"], [4., 5., 6., 8.])
    assert "neg_log10_p" in part


def test_empty_selection_and_joint_statistics(tmp_path):
    empty = tmp_path / "empty"
    rows, _ = write_indexed_sumstats(empty, ["a"], ["x"], 24,
        [(0, 1, np.ones((1, 1)), np.full((1, 1), np.nan), None)],
        kind="filtered", df=20)
    manifest, parts = open_indexed_sumstats(empty)
    assert rows == 0 and list(parts) == [] and manifest["shape"] == [1, 1]
    joint = tmp_path / "joint"
    write_indexed_sumstats(joint, ["a", "b"], ["x", "y"], 24,
        [(0, 2, None, np.array([7., np.nan]), None)], kind="jagwas", df=20, chi2_df=2)
    manifest, parts = open_indexed_sumstats(joint)
    assert manifest["df"] == 2 and manifest["kind"] == "jagwas"
    np.testing.assert_array_equal(next(parts)["chi2"], [7.])

@pytest.mark.skipif(not __import__('pathlib').Path('/proc/self/fd').exists(),reason='Linux fd inspection')
def test_selected_pair_df_and_durable_publication(tmp_path,monkeypatch):
    import os
    import torchgwas.sumstats_indexed as indexed
    synced=[]
    original=os.fsync
    def sync(fd):
        synced.append(os.readlink('/proc/self/fd/'+str(fd)))
        return original(fd)
    monkeypatch.setattr(os,'fsync',sync)
    chunk=(0,2,np.array([0,1]),np.array([0,0]),np.ones(2,np.float32),np.full(2,3.,np.float32),np.array([2.,98.],np.float32))
    write_indexed_sumstats(tmp_path,['a','b'],['x'],100,[chunk],kind='significant',df=98)
    manifest,parts=open_indexed_sumstats(tmp_path)
    np.testing.assert_array_equal(next(parts)['df'],[2.,98.])
    assert manifest['version']==2 and manifest['nominal_df']==98
    assert any(p.endswith('variant_ids.npy') for p in synced)
    assert str(tmp_path) in synced and any('manifest.json.' in p for p in synced)


def test_failed_validation_invalidates_previous_manifest(tmp_path):
    chunks = [(0, 1, np.ones((1, 1)), np.ones((1, 1)), None)]
    write_indexed_sumstats(tmp_path, ['a'], ['x'], 24, chunks, kind='filtered', df=20)
    assert (tmp_path/'manifest.json').exists()
    def fail():
        assert (tmp_path/'variant_ids.npy').exists()
        raise RuntimeError('tile metadata mismatch')
    with pytest.raises(RuntimeError, match='tile metadata mismatch'):
        write_indexed_sumstats(tmp_path, ['a'], ['x'], 24, chunks, kind='filtered', df=20, before_publish=fail)
    assert not (tmp_path/'manifest.json').exists()

def test_indexed_writer_failure_closes_producer(tmp_path, monkeypatch):
    import torchgwas.sumstats_indexed as indexed
    closed = []
    def chunks():
        try:
            yield (0, 1, np.ones((1, 1)), np.ones((1, 1)), None)
        finally:
            closed.append(True)
    def fail(*args, **kwargs):
        raise OSError('durable output failure')
    monkeypatch.setattr(indexed.np, 'savez', fail)
    with pytest.raises(OSError, match='durable output failure'):
        write_indexed_sumstats(tmp_path, ['a'], ['x'], 24, chunks(), kind='filtered', df=20)
    assert closed == [True] and not (tmp_path/'manifest.json').exists()

def test_live_jagwas_writer_exposes_active_range_and_publication(tmp_path, monkeypatch):
    import torchgwas.sumstats_indexed as indexed
    progress = IndexedWriterProgress('jagwas')
    partition = IndexedOutputPartition('cuda:1', (6, 8), (0, 2))
    entered = threading.Event()
    release = threading.Event()
    original = indexed.np.savez
    def slow_savez(*args, **kwargs):
        entered.set()
        assert release.wait(5.)
        return original(*args, **kwargs)
    monkeypatch.setattr(indexed.np, 'savez', slow_savez)
    observed_at_callback = []
    def written(event):
        observed_at_callback.append((event, progress.snapshot()))
    errors = []
    def write():
        try:
            write_indexed_sumstats(tmp_path, ['a', 'b'], ['x', 'y'], 129,
                [(0, 2, None, np.array([7., np.nan]), None)],
                kind='jagwas', df=100, chi2_df=2, fsync=False,
                on_chunk_written=written, partition_for_range=lambda *_: partition,
                variant_offset=6, live_progress=progress)
        except BaseException as error:
            errors.append(error)
    thread = threading.Thread(target=write)
    thread.start()
    try:
        assert entered.wait(5.)
        active = progress.snapshot()
        assert active['phase'] == 'emitting_part'
        assert active['active']['variant_range'] == [6, 8]
        assert active['active']['device'] == 'cuda:1'
        assert active['completed_chunks'] == 0
    finally:
        release.set()
        thread.join(5.)
    assert not thread.is_alive() and not errors
    assert len(observed_at_callback) == 1
    event, at_callback = observed_at_callback[0]
    assert event.rows == 1 and at_callback['phase'] == 'waiting_for_result'
    assert at_callback['active'] is None and at_callback['completed_chunks'] == 1
    assert progress.snapshot()['phase'] == 'published'
    manifest, parts = open_indexed_sumstats(tmp_path)
    assert manifest['rows'] == 1
    np.testing.assert_array_equal(next(parts)['variant_index'], [0])


def test_live_significant_writer_keeps_tile_owner_without_range_resolver(tmp_path, monkeypatch):
    import torchgwas.sumstats_indexed as indexed
    progress = IndexedWriterProgress('significant')
    partition = IndexedOutputPartition('cuda:1', (0, 2), (2, 4))
    chunk = PartitionedIndexedChunk((0, 2, np.array([0, 1]),
        np.array([2, 3]), np.array([.2, .3], np.float32),
        np.array([5., 6.], np.float32), np.array([100., 100.])), partition)
    entered = threading.Event()
    release = threading.Event()
    original = indexed.np.savez
    def slow_savez(*args, **kwargs):
        entered.set()
        assert release.wait(5.)
        return original(*args, **kwargs)
    monkeypatch.setattr(indexed.np, 'savez', slow_savez)
    events = []
    failures = []
    def write():
        try:
            write_indexed_sumstats(tmp_path, ['a', 'b'],
                ['p0', 'p1', 'p2', 'p3'], 129, [chunk],
                kind='significant', df=100, fsync=False,
                on_chunk_written=events.append, live_progress=progress)
        except BaseException as error:
            failures.append(error)
    thread = threading.Thread(target=write)
    thread.start()
    try:
        assert entered.wait(5.)
        active = progress.snapshot()
        assert active['phase'] == 'emitting_part'
        assert active['active'] == dict(device='cuda:1',
            variant_range=[0, 2], trait_range=[2, 4])
    finally:
        release.set()
        thread.join(5.)
    assert not thread.is_alive() and not failures
    assert len(events) == 1 and events[0].partition == partition
    assert progress.snapshot()['phase'] == 'published'
    manifest, parts = open_indexed_sumstats(tmp_path)
    assert manifest['rows'] == 2
    np.testing.assert_array_equal(next(parts)['trait_index'], [2, 3])


def test_default_parts_follow_bytes_not_rows(tmp_path):
    # 40 chunks of 100 pairs with a 16 kB target: parts hold whole chunks until
    # the target is crossed, so their count follows bytes.
    from torchgwas.sumstats_indexed import write_indexed_sumstats, open_indexed_sumstats
    rng = np.random.default_rng(3)
    chunks = []
    for c in range(40):
        vi = np.repeat(np.arange(c * 10, c * 10 + 10), 10)
        ti = np.tile(np.arange(10), 10)
        chunks.append((c * 10, c * 10 + 10, vi, ti, rng.normal(size=100).astype(np.float32),
                       rng.normal(size=100).astype(np.float32), np.full(100, 20.0, np.float32)))
    names = [f'v{i}' for i in range(400)]
    total, summary = write_indexed_sumstats(tmp_path, names, [f't{i}' for i in range(10)], 24, iter(chunks),
                                            kind='significant', df=20, coalesce_bytes=16_000)
    manifest, parts = open_indexed_sumstats(tmp_path)
    # 32 bytes per pair (two int64 indices; beta, t, df and -log10 P as
    # float32): 3,200 per chunk, so a part closes after 5 chunks: 8 parts.
    assert total == 4000 and len(manifest['parts']) == 8
    rows = np.concatenate([part['variant_index'] for part in parts])
    assert np.all(np.diff(rows) >= 0)
