"""The fast `.bim` parser must agree with the old one exactly, or not run.

Replacing a parser is a correctness change disguised as a performance one: a
column read as the wrong type, a marker id mangled, an allele swapped, and
every downstream result is labelled wrong while looking entirely plausible.
So these tests compare the two parsers rather than checking the fast one
against expectations written by the same hand that wrote it.

The speed that motivated it, measured cold on a 253.9 MB / 8,931,083-row file:
pandas with `sep=r"\\s+"` 13.08 s, pyarrow 0.49 s, and simply reading the bytes
0.35 s.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from torchgwas.bed import _bim_delimiter, _read_bim_fast


def pandas_reference(path: Path):
    """Exactly what the constructor did before the fast path existed."""
    table = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 3, 4, 5],
                        dtype={0: str, 1: str, 3: np.int64, 4: str, 5: str},
                        memory_map=True)
    return (table.iloc[:, 0].to_numpy(dtype=object),
            table.iloc[:, 1].to_numpy(dtype=object),
            table[3].to_numpy(dtype=np.int64),
            table[4].to_numpy(dtype=object),
            table[5].to_numpy(dtype=object))


ROWS = [
    ("1", "rs1", "0", "1000", "A", "G"),
    ("2", "rs2", "0", "2000", "T", "C"),
    ("X", "chrX:12345:A:T", "0", "12345", "A", "T"),
    ("MT", "rs999", "0", "16000", "AT", "A"),      # indel alleles
    ("22", "22:1-DEL", "0", "999999999", "-", "AC"),
    ("9", "rs_name_with_underscores", "0", "1", "G", "A"),
]


class AgreementTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, name: str, separator: str, rows=ROWS) -> Path:
        path = self.root / name
        path.write_text("".join(separator.join(row) + "\n" for row in rows))
        return path

    def assert_same(self, path: Path):
        fast = _read_bim_fast(path)
        self.assertIsNotNone(fast, f"fast parser declined {path.name}")
        reference = pandas_reference(path)
        names = ("chromosomes", "marker_ids", "positions", "other_alleles",
                 "effect_alleles")
        for name, got, want in zip(names, fast, reference):
            np.testing.assert_array_equal(got, want, err_msg=name)
            self.assertEqual(np.asarray(got).dtype, np.asarray(want).dtype,
                             f"{name} dtype")

    def test_tab_separated_agrees(self):
        self.assert_same(self.write("tabs.bim", "\t"))

    def test_space_separated_agrees(self):
        self.assert_same(self.write("spaces.bim", " "))

    def test_numeric_looking_marker_ids_stay_strings(self):
        """`22:1-DEL` is fine, but a bare `12345` id must not become an int.

        A type change here would break every join against the metadata, and
        it is exactly what an inferring parser does by default.
        """
        rows = [("1", "12345", "0", "500", "A", "G"),
                ("1", "00123", "0", "600", "A", "G")]
        path = self.write("numericids.bim", "\t", rows)
        fast = _read_bim_fast(path)
        self.assertIsNotNone(fast)
        self.assertEqual(list(fast[1]), ["12345", "00123"])
        self.assert_same(path)

    def test_chromosome_stays_a_string(self):
        """`1` and `X` must live in the same column without one becoming an int."""
        fast = _read_bim_fast(self.write("chrom.bim", "\t"))
        self.assertTrue(all(isinstance(value, str) for value in fast[0]))

    def test_positions_are_int64(self):
        fast = _read_bim_fast(self.write("pos.bim", "\t"))
        self.assertEqual(fast[2].dtype, np.int64)
        self.assertEqual(fast[2][-2], 999_999_999)


class DeclineTests(unittest.TestCase):
    """When the file is not simple, declining beats guessing."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_column_aligned_files_are_declined(self):
        """Runs of spaces cannot be read with a single-character delimiter.

        These exist. Reading one as space-delimited would produce empty fields
        between the runs and silently shift every column.
        """
        path = self.root / "aligned.bim"
        path.write_text("1    rs1   0     1000   A   G\n"
                        "2    rs2   0     2000   T   C\n")
        self.assertIsNone(_bim_delimiter(path))
        self.assertIsNone(_read_bim_fast(path))

    def test_an_empty_file_is_declined(self):
        path = self.root / "empty.bim"
        path.write_bytes(b"")
        self.assertIsNone(_read_bim_fast(path))

    def test_too_few_columns_is_declined(self):
        path = self.root / "short.bim"
        path.write_text("1\trs1\t0\n2\trs2\t0\n")
        self.assertIsNone(_read_bim_fast(path))

    def test_the_declined_path_still_works_through_the_reader(self):
        """A column-aligned .bim must still open, via pandas."""
        from torchgwas.bed import PlinkBedGenotype
        (self.root / "c.bim").write_text(
            "1    rs1   0     1000   A   G\n1    rs2   0  2000   T   C\n")
        (self.root / "c.fam").write_text("F1 I1 0 0 1 -9\nF2 I2 0 0 1 -9\n")
        with (self.root / "c.bed").open("wb") as handle:
            handle.write(b"\x6c\x1b\x01")
            handle.write(bytes([0b00001100, 0b00001100]))
        genotype = PlinkBedGenotype(self.root / "c")
        self.addCleanup(genotype.__del__)
        self.assertEqual(list(genotype.marker_ids), ["rs1", "rs2"])
        self.assertEqual(list(genotype.chromosomes), ["1", "1"])
        self.assertEqual(list(genotype.positions), [1000, 2000])


if __name__ == "__main__":
    unittest.main()
