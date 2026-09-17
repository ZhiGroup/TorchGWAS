"""Synthetic BGEN fixtures exercise exact dosage bits and physical I/O order."""
import sqlite3
import struct
import zlib
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pytest
from torchgwas.bgen import BgenGenotype


def write_bgen(path, bits=8, *, missing=False, phased=False, bad_checksum=False,
               compression=1, layout=2):
    samples=['a','b','c','d']
    denom=(1<<bits)-1
    probabilities=[[(denom,0),(0,denom),(0,0),(denom//3,denom//4)],
                   [(0,0),(denom,0),(0,denom),(denom//7,denom//5)]]
    sampleblock=struct.pack('<II',8+sum(2+len(x) for x in samples),len(samples))+b''.join(struct.pack('<H',len(x))+x.encode() for x in samples)
    header=struct.pack('<IIII4sI',20+len(sampleblock),20,len(probabilities),len(samples),b'bgen',compression|(layout<<2)|(1<<31))+sampleblock
    data=bytearray(header)
    rows=[]
    expected=[]
    for i,probs in enumerate(probabilities):
        word=0
        for j,(p0,p1) in enumerate(probs):
            word |= p0 << (2*j*bits)
            word |= p1 << ((2*j+1)*bits)
        packed=word.to_bytes((len(samples)*2*bits+7)//8,'little')
        pm=bytes([130 if missing and j==1 else 2 for j in range(4)])
        raw=struct.pack('<IHBB',4,2,2,2)+pm+bytes([int(phased),bits])+packed
        if compression==0:
            compressed=bytearray(raw)
        elif compression==1:
            compressed=bytearray(zlib.compress(raw))
        else:
            import zstandard
            compressed=bytearray(zstandard.ZstdCompressor().compress(raw))
        if bad_checksum:compressed[-1]^=1
        def text(value,width='H'):
            return struct.pack('<'+width,len(value))+value.encode()
        # Reverse positions: iterator must preserve physical file order.
        position=200-i*100
        payload_header=struct.pack('<I',len(compressed)) if compression==0 else struct.pack('<II',len(compressed)+4,len(raw))
        record=text('id'+str(i))+text('rs'+str(i))+text('1')+struct.pack('<IH',position,2)+text('A','I')+text('G','I')+payload_header+compressed
        rows.append(('1',position,'rs'+str(i),2,'A','G',len(data),len(record)))
        data+=record
        expected.append([(2*denom-2*p0-p1)/denom for p0,p1 in probs])
    path.write_bytes(data)
    with sqlite3.connect(str(path)+'.bgi') as db:
        db.execute('CREATE TABLE Variant (chromosome TEXT,position INTEGER,rsid TEXT,number_of_alleles INTEGER,allele1 TEXT,allele2 TEXT,file_start_position INTEGER,size_in_bytes INTEGER)')
        db.executemany('INSERT INTO Variant VALUES (?,?,?,?,?,?,?,?)',rows[::-1])
    return np.asarray(expected,dtype=np.float64)


@pytest.mark.parametrize('bits',[1,3,8,12,16,24,32])
def test_direct_preserves_probabilities(tmp_path,bits):
    path=tmp_path/'test.bgen'
    expected=write_bgen(path,bits)
    source=BgenGenotype(path,decode_backend='cpu')
    assert source.positions.tolist()==[200,100]
    actual=np.concatenate([x.T for _,_,x in source.iter_chunks(1,np.float64)],axis=0)
    np.testing.assert_array_equal(actual,expected)
    source.select_samples(['d','a'])
    np.testing.assert_array_equal(source.read_chunk(0,2,np.float64),expected[:,[3,0]].T)
    assert source.dosage_scale==1


@pytest.mark.parametrize('kind',['phased','bad_checksum'])
def test_rejects_unsupported_or_corrupt(tmp_path,kind):
    path=tmp_path/'bad.bgen'
    write_bgen(path,**{kind:True})
    with pytest.raises((ValueError,zlib.error)):
        BgenGenotype(path).read_chunk(0,2)


def test_layout_1_is_refused_by_name(tmp_path):
    """A refused file must say *which* problem it has.

    Layout 1 and a reserved compression flag used to share one message, so
    neither told the reader what to do -- and the answers differ: Layout 1 is
    converted, a reserved codec means the file is not standard.
    """
    path=tmp_path/'layout1.bgen'
    write_bgen(path,layout=1)
    with pytest.raises(ValueError,match='Layout 1'):
        BgenGenotype(path)


def test_reserved_compression_flag_is_refused_by_name(tmp_path):
    path=tmp_path/'codec3.bgen'
    write_bgen(path,compression=3)
    with pytest.raises(ValueError,match='reserved'):
        BgenGenotype(path)


@pytest.mark.parametrize('compression',[0,1,2])
@pytest.mark.parametrize('bits',[1,8,16,32])
def test_supported_scope_decodes_exactly(tmp_path,compression,bits):
    """The documented scope, decoded exactly across its whole matrix.

    Supported is Layout 2, biallelic, unphased diploid, bit depth 1-32, under
    any of the three standard codecs. Anything outside that is refused rather
    than approximated, and the refusals have their own tests above. Asserting
    the supported set as a grid is what turns "we support Layout 2" from a
    claim in a docstring into something the suite can keep true.
    """
    if compression==2:pytest.importorskip('zstandard')
    path=tmp_path/f'scope_{compression}_{bits}.bgen'
    expected=write_bgen(path,bits,compression=compression)
    source=BgenGenotype(path,decode_backend='cpu')
    actual=np.concatenate([x.T for _,_,x in source.iter_chunks(1,np.float64)],axis=0)
    np.testing.assert_array_equal(actual,expected)


def test_reselecting_samples_invalidates_the_cached_decoder(tmp_path):
    """A cached per-thread decoder must not outlive the selection it was built for.

    The decoder holds a copy of `_sample_indices`, so a stale one keeps emitting
    the previous cohort's columns with no error at all -- the array is the right
    shape and the values are real dosages, just the wrong samples'.
    """
    path=tmp_path/'reselect.bgen'
    expected=write_bgen(path)
    source=BgenGenotype(path,decode_backend='cpu')
    np.testing.assert_array_equal(source.read_chunk(0,2,np.float32),
                                  expected.astype(np.float32).T)
    source.select_samples(['d','a'])
    np.testing.assert_array_equal(source.read_chunk(0,2,np.float32),
                                  expected[:,[3,0]].astype(np.float32).T)
    # `select_samples` is relative to the current selection, not to the file:
    # after the reorder above, 'a' is the second of the two retained samples.
    source.select_samples(['a'])
    np.testing.assert_array_equal(source.read_chunk(0,2,np.float32),
                                  expected[:,[0]].astype(np.float32).T)


def test_concurrent_chunk_reads_agree_with_serial_ones(tmp_path):
    """Reused staging buffers are per thread; prove that holds under threads.

    Sharing one would corrupt the inflate scratch and surface as
    `zlib inflate failed`, or worse, as silently interleaved rows.
    """
    from concurrent.futures import ThreadPoolExecutor

    path=tmp_path/'threaded.bgen'
    expected=write_bgen(path)
    source=BgenGenotype(path,decode_backend='cpu')
    serial=[source.read_chunk(i,i+1,np.float32) for i in range(2)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        for _ in range(20):
            results=list(pool.map(lambda i:source.read_chunk(i,i+1,np.float32),
                                  [0,1]*8))
            for index,got in zip([0,1]*8,results):
                np.testing.assert_array_equal(got,serial[index])
    np.testing.assert_array_equal(np.concatenate(serial,axis=1),
                                  expected.astype(np.float32).T)


def test_bgi_mismatch(tmp_path):
    path=tmp_path/'test.bgen'
    write_bgen(path)
    with sqlite3.connect(str(path)+'.bgi') as db:
        db.execute("UPDATE Variant SET allele2='T'")
    with pytest.raises(ValueError,match='metadata mismatch'):
        BgenGenotype(path).read_chunk(0,1)


@pytest.mark.parametrize('bits',[1,3,8,12,16,24,32])
def test_gpu_matches_exact_reference(tmp_path,bits):
    import torch
    from torchgwas.bgen_gpu import load_decoder_library
    if not torch.cuda.is_available():pytest.skip('CUDA unavailable')
    try:load_decoder_library()
    except ImportError as exc:pytest.skip(str(exc))
    path=tmp_path/'gpu.bgen'
    expected=write_bgen(path,bits)
    source=BgenGenotype(path,decode_backend='gpu',decode_batch_size=2).select_samples(['d','a','b'])
    try:
        actual=np.concatenate([x.cpu().numpy() for _,_,x in source.iter_device_chunks(1,'cuda:0')])
    except RuntimeError as exc:
        # The library loads on any card but carries code only for the
        # architectures it was built for, so an unsupported GPU fails here and
        # not at load. Seen on an RTX 2080 Ti (sm_75) against a library built
        # for sm_80/sm_90. Skipping is right for the test; the shipped fallback
        # is a separate matter and is tracked.
        if 'no kernel image' not in str(exc):raise
        pytest.skip('GPU BGEN decoder has no kernel for this architecture')
    np.testing.assert_array_equal(actual,expected[:,[3,0,1]].astype(np.float32))

def test_missing_is_nan(tmp_path):
    path=tmp_path/'missing.bgen'
    write_bgen(path,missing=True)
    assert np.isnan(BgenGenotype(path).read_chunk(0,2)[1]).all()


@pytest.mark.parametrize('compression',[0,2])
def test_cpu_compression_layouts(tmp_path,compression):
    if compression==2:pytest.importorskip('zstandard')
    path=tmp_path/'compression.bgen'
    expected=write_bgen(path,compression=compression)
    actual=BgenGenotype(path,decode_backend='auto').read_chunk(0,2,np.float64)
    np.testing.assert_array_equal(actual,expected.T)


# At 8 and 16 bits: 16 is what imputed UK Biobank-style BGEN carries, and a
# decoder fast path for it is the obvious future change, so a checksum flip
# or a missing sample must stay caught at that width too.
@pytest.mark.parametrize('bits',[8,16])
@pytest.mark.parametrize('kind',['bad_checksum','metadata','missing'])
def test_gpu_validation(tmp_path,kind,bits):
    import torch
    from torchgwas.bgen_gpu import load_decoder_library
    if not torch.cuda.is_available():pytest.skip('CUDA unavailable')
    try:load_decoder_library()
    except ImportError as exc:pytest.skip(str(exc))
    path=tmp_path/'validate.bgen'
    write_bgen(path,bits,**({kind:True} if kind!='metadata' else {}))
    if kind=='metadata':
        with sqlite3.connect(str(path)+'.bgi') as db:
            db.execute("UPDATE Variant SET allele2='T'")
    source=BgenGenotype(path,decode_backend='gpu')
    if kind=='missing':
        try:
            actual=np.concatenate([x.cpu().numpy() for _,_,x in source.iter_device_chunks(1,'cuda:0')])
        except RuntimeError as exc:
            if 'no kernel image' not in str(exc):raise
            pytest.skip('GPU BGEN decoder has no kernel for this architecture')
        assert np.isnan(actual[:,1]).all()
    else:
        with pytest.raises((ValueError,RuntimeError)):
            list(source.iter_device_chunks(1,'cuda:0'))



class BgenIndexCacheTestCase(unittest.TestCase):
    """The parsed `.bgi` must round-trip exactly, and a stale cache must not load.

    Parsing the index is 38% of a full BGEN run -- 22.7 s of 59.6 s -- and it
    is all CPU: reading the whole 391 MB index off cold disk takes 0.2 s. A
    cache takes the warm open from 19.2 s to **1.3 s, 14.9x**. The risk it
    introduces is the one these tests exist for: a cache that no longer matches
    its index would make the scan read the wrong byte offsets and produce
    confident nonsense, so staleness must be detected rather than assumed away.
    """

    def _arrays(self, count=7):
        rng = np.random.default_rng(20260913)
        lengths = rng.integers(16, 64, size=count).astype(np.uint64)
        offsets = np.cumsum(np.concatenate(([0], lengths[:-1]))).astype(np.uint64)
        return {
            "chromosomes": np.array([str(1 + i % 22) for i in range(count)],
                                    dtype=object),
            "positions": np.arange(count, dtype=np.int64) * 100,
            "marker_ids": np.array([f"rs{i}" for i in range(count)],
                                   dtype=object),
            "other_alleles": np.array(["A"] * count, dtype=object),
            "effect_alleles": np.array(["C"] * count, dtype=object),
            "offsets": offsets,
            "lengths": lengths,
        }

    def test_round_trip_is_exact_and_restores_object_dtype(self):
        from torchgwas.bgen import (_load_bgen_index_cache,
                                    write_bgen_index_cache)

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"not a real index, only its identity matters")
            arrays = self._arrays()
            write_bgen_index_cache(tmpdir, bgi, **arrays)
            loaded, _path = _load_bgen_index_cache(tmpdir, bgi)
            self.assertIsNotNone(loaded)
            for name, original in arrays.items():
                self.assertTrue(np.array_equal(loaded[name], original), name)
            # Downstream code cannot be able to tell the cache from SQLite.
            for name in ("chromosomes", "marker_ids", "other_alleles",
                         "effect_alleles"):
                self.assertEqual(loaded[name].dtype, np.dtype(object), name)

    def test_a_changed_index_invalidates_the_cache(self):
        from torchgwas.bgen import (_load_bgen_index_cache,
                                    write_bgen_index_cache)

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"original index")
            write_bgen_index_cache(tmpdir, bgi, **self._arrays())
            self.assertIsNotNone(_load_bgen_index_cache(tmpdir, bgi)[0])
            # A different index of a different size is a different file.
            bgi.write_bytes(b"a longer replacement index entirely")
            self.assertIsNone(_load_bgen_index_cache(tmpdir, bgi)[0],
                              "a stale cache was accepted")

    def test_a_same_size_rebuild_invalidates_the_cache(self):
        # The dangerous case is a rebuilt index with an identical byte count:
        # the size check cannot see it and only the mtime can. Written as its
        # own test because the size-changing case above passes even with the
        # mtime check deleted, which a mutation run showed.
        import os
        import time

        from torchgwas.bgen import (_load_bgen_index_cache,
                                    write_bgen_index_cache)

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"index one")
            write_bgen_index_cache(tmpdir, bgi, **self._arrays())
            self.assertIsNotNone(_load_bgen_index_cache(tmpdir, bgi)[0])
            bgi.write_bytes(b"index two")          # same length, new content
            stat = bgi.stat()
            os.utime(bgi, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
            self.assertEqual(bgi.stat().st_size, len(b"index one"))
            self.assertIsNone(_load_bgen_index_cache(tmpdir, bgi)[0],
                              "a same-size rebuild reused a stale cache")

    def test_a_corrupt_manifest_is_ignored_rather_than_raised(self):
        from torchgwas.bgen import (_load_bgen_index_cache,
                                    write_bgen_index_cache)

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"original index")
            path = write_bgen_index_cache(tmpdir, bgi, **self._arrays())
            (path / "manifest.json").write_text("{ this is not json")
            loaded, _ = _load_bgen_index_cache(tmpdir, bgi)
            self.assertIsNone(loaded)

    def test_the_cache_never_stores_pickles(self):
        # Loading a cache must not be able to execute code, so every array is
        # written with allow_pickle=False and must therefore be a plain dtype.
        from torchgwas.bgen import write_bgen_index_cache

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"original index")
            path = write_bgen_index_cache(tmpdir, bgi, **self._arrays())
            for entry in path.glob("*.npy"):
                array = np.load(entry, allow_pickle=False)
                self.assertNotEqual(array.dtype, np.dtype(object), entry.name)

    def test_inconsistent_lengths_are_refused(self):
        from torchgwas.bgen import write_bgen_index_cache

        with tempfile.TemporaryDirectory() as tmpdir:
            bgi = Path(tmpdir) / "x.bgen.bgi"
            bgi.write_bytes(b"original index")
            arrays = self._arrays()
            arrays["positions"] = arrays["positions"][:-1]
            with self.assertRaises(ValueError):
                write_bgen_index_cache(tmpdir, bgi, **arrays)
