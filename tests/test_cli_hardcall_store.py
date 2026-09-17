"""The store has to be reachable from the CLI, or it is not shipped.

A format nobody can build or select is a benchmark artifact, not a feature.
These tests drive the real argument parser and the real runner on a small
fixture, so a wiring mistake -- a flag added to the parser but never read, a
runner never dispatched -- fails here rather than in a user's hands.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from torchgwas.cli import _build_parser as build_parser
from torchgwas.hardcall_store import HardcallStore


def write_triplet(root: Path, variants: int, samples: int) -> Path:
    rng = np.random.default_rng(5)
    stride = (samples + 3) // 4
    packed = rng.integers(0, 256, size=(variants, stride), dtype=np.uint8)
    (root / "c.fam").write_text("".join(
        f"F{i} I{i} 0 0 1 -9\n" for i in range(samples)))
    (root / "c.bim").write_text("".join(
        f"1 rs{i} 0 {100 + i} A G\n" for i in range(variants)))
    with (root / "c.bed").open("wb") as handle:
        handle.write(b"\x6c\x1b\x01")
        handle.write(packed.tobytes())
    return root / "c"


class BuildSubcommandTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.prefix = write_triplet(self.root, 900, 61)

    def build(self, *extra):
        from torchgwas.cli import _run_build_hardcall_store
        args = build_parser().parse_args([
            "build-hardcall-store", "--genotype", str(self.prefix),
            "--out", str(self.root / "store"), *extra])
        return _run_build_hardcall_store(args)

    def test_it_builds_a_readable_store(self):
        self.assertEqual(self.build("--frame-variants", "128"), 0)
        with HardcallStore(self.root / "store") as store:
            self.assertEqual(store.shape, (61, 900))

    def test_the_store_matches_the_bed_it_came_from(self):
        self.build("--frame-variants", "128")
        raw = (self.root / "c.bed").read_bytes()[3:]
        stride = (61 + 3) // 4
        expected = np.frombuffer(raw, dtype=np.uint8).reshape(900, stride)
        with HardcallStore(self.root / "store") as store:
            np.testing.assert_array_equal(store.read_packed(0, 900), expected)

    def test_verification_runs_by_default(self):
        """It is on by default because a silent corruption is unrecoverable.

        The store stands in for every genotype byte of every scan that uses it.
        """
        parser = build_parser()
        args = parser.parse_args([
            "build-hardcall-store", "--genotype", "x", "--out", "y"])
        self.assertGreater(args.verify, 0)

    def test_defaults_match_the_measured_optimum(self):
        """zstd level 3 and 2,048-variant frames, both chosen by measurement."""
        args = build_parser().parse_args([
            "build-hardcall-store", "--genotype", "x", "--out", "y"])
        self.assertEqual(args.level, 3)
        self.assertEqual(args.frame_variants, 2048)


class LinearFlagTests(unittest.TestCase):
    def test_linear_accepts_a_store(self):
        args = build_parser().parse_args([
            "linear", "--genotype", "cohort.bed", "--phenotype", "p.npy",
            "--output-dir", "out", "--hardcall-store", "cohort"])
        self.assertEqual(args.hardcall_store, "cohort")

    def test_linear_defaults_to_no_store(self):
        args = build_parser().parse_args([
            "linear", "--genotype", "cohort.bed", "--phenotype", "p.npy",
            "--output-dir", "out"])
        self.assertIsNone(args.hardcall_store)

    def test_prep_does_not_take_the_flag(self):
        """`prep` never reads genotype bytes, and its runner has no such arg.

        Wiring it there would raise `AttributeError` at run time -- which it
        briefly did, because a blanket edit added the keyword to every
        `load_genotype` call site in this module.
        """
        args = build_parser().parse_args([
            "prep", "--genotype", "cohort.bed", "--phenotype", "p.npy",
            "--output-dir", "out"])
        self.assertFalse(hasattr(args, "hardcall_store"))


if __name__ == "__main__":
    unittest.main()
