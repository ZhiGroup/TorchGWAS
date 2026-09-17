"""The variant-id sidecar must be byte-identical to what it replaced.

It is the fixed cost of writing sumstats -- 3.72 s of a full-genome run at
K = 1, flat in K -- so it is worth making fast. But it names every marker in
the output, and silently mangling those to save a second would be a very bad
trade. These tests compare the fast path against the text path it replaced.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from torchgwas.api import _write_variant_ids


def reference(names) -> str:
    """Exactly what the writer did before the fast path."""
    return "\n".join(map(str, names)) + "\n"


class AgreementTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "variant_ids.txt"

    def assert_same(self, names):
        written = _write_variant_ids(self.path, names)
        # Read as UTF-8 explicitly. `read_text()` would use the locale default,
        # which is ASCII on the lab host -- the very dependence that made the
        # original writer crash there on a non-ASCII id.
        got = self.path.read_text(encoding="utf-8")
        self.assertEqual(got, reference(names))
        self.assertEqual(written, len(got.encode("utf-8")))
        return got

    def test_numpy_unicode_array_the_real_case(self):
        """`marker_ids` arrives as a numpy `<U16` array."""
        names = np.array(["rs1", "rs22", "chr1:100:A:T"], dtype="<U16")
        self.assert_same(names)

    def test_a_plain_list_of_strings(self):
        self.assert_same(["rs1", "rs2", "rs3"])

    def test_ids_that_look_numeric_keep_their_exact_text(self):
        """`00123` must not become `123`, the same hazard as the .bim parse."""
        names = np.array(["00123", "12345", "007"], dtype="<U16")
        got = self.assert_same(names)
        self.assertEqual(got.splitlines(), ["00123", "12345", "007"])

    def test_a_single_id(self):
        self.assert_same(np.array(["rs1"], dtype="<U8"))

    def test_ids_at_the_full_dtype_width(self):
        """`<U16` holding exactly 16 characters must not be truncated."""
        names = np.array(["A" * 16, "B" * 16], dtype="<U16")
        got = self.assert_same(names)
        self.assertEqual([len(line) for line in got.splitlines()], [16, 16])

    def test_non_ascii_falls_back_and_stays_correct(self):
        """The byte path cannot represent these, so it must not be used.

        Unusual in a marker id, but mangling one to save a second would be far
        worse than being slow.
        """
        names = np.array(["rsA", "rsé", "rs中"], dtype="<U8")
        got = self.assert_same(names)
        self.assertIn("rsé", got)
        self.assertIn("rs中", got)

    def test_an_empty_list_writes_just_a_newline(self):
        """Matching the old behaviour exactly, odd as it is."""
        self.assert_same([])


class SizeTests(unittest.TestCase):
    def test_the_reported_byte_count_matches_the_file(self):
        """The summary reports this, and a wrong number is a silent lie."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ids.txt"
            names = np.array([f"rs{i}" for i in range(1000)], dtype="<U16")
            written = _write_variant_ids(path, names)
            self.assertEqual(written, path.stat().st_size)


if __name__ == "__main__":
    unittest.main()
