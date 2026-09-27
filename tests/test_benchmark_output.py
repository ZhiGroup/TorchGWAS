import tempfile
import unittest
from pathlib import Path

import numpy as np

from benchmarks.benchmark_bed_gpu_pipeline import AsyncNpyWriter


class TestAsyncNpyWriter(unittest.TestCase):
    def test_ordered_uneven_chunks_round_trip_exactly(self):
        expected = np.arange(33, dtype=np.float32).reshape(11, 3)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "t.npy"
            writer = AsyncNpyWriter(
                output,
                shape=expected.shape,
                chunk_size=4,
                depth=2,
            )
            for start in range(0, expected.shape[0], 4):
                end = min(start + 4, expected.shape[0])
                writer.submit(start, end, expected[start:end])
            nbytes = writer.close()

            observed = np.load(output, allow_pickle=False)
            np.testing.assert_array_equal(observed, expected)
            self.assertEqual(nbytes, output.stat().st_size)

    def test_refuses_to_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "t.npy"
            output.touch()
            with self.assertRaises(FileExistsError):
                AsyncNpyWriter(output, shape=(2, 2), chunk_size=1, depth=1)


if __name__ == "__main__":
    unittest.main()
