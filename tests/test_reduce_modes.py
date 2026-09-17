"""The reduction has two user-facing modes and refuses the rest.

Standing decision: `significant` and `jagwas`. The per-variant top-k spellings
(`max-abs-t`, `max-t2`, `min-p`, `top-k`) stay as machinery -- SignificantPairs
is built on top-k and the trait-blocked paths still construct
`VariantReduction` -- but they are not an answer a caller should reach for.

These tests exist because the decision was written down and then broken by
every stress run in a single day, all of which passed `reduce='max-abs-t'`. A
decision recorded only in prose is enforced by whoever remembers to reread the
prose; a test is enforced by the test runner.
"""
import sys
import unittest

sys.path.insert(0, "src")

import numpy as np

from torchgwas.api import run_linear_gwas
from torchgwas.reduce import VariantReduction


def _tiny():
    rng = np.random.default_rng(0)
    genotype = rng.integers(0, 3, size=(40, 12)).astype(np.float32)
    phenotype = rng.normal(size=(40, 3))
    return genotype, phenotype


class RejectedModesTestCase(unittest.TestCase):
    def test_every_top_k_spelling_is_refused(self):
        genotype, phenotype = _tiny()
        for mode in ("max-abs-t", "max-t2", "min-p", "top-k"):
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError) as caught:
                    run_linear_gwas(genotype, phenotype, reduce=mode,
                                    output_dir="/tmp/should-not-be-created")
                message = str(caught.exception)
                # The error has to name the two real modes, or it sends the
                # caller looking through source for what to use instead.
                self.assertIn("significant", message)
                self.assertIn("jagwas", message)

    def test_the_refusal_names_the_offending_value(self):
        genotype, phenotype = _tiny()
        with self.assertRaises(ValueError) as caught:
            run_linear_gwas(genotype, phenotype, reduce="max-abs-t",
                            output_dir="/tmp/should-not-be-created")
        self.assertIn("max-abs-t", str(caught.exception))

    def test_reduce_top_k_is_refused_on_its_own(self):
        genotype, phenotype = _tiny()
        with self.assertRaises(ValueError):
            run_linear_gwas(genotype, phenotype, reduce_top_k=5,
                            output_dir="/tmp/should-not-be-created")

    def test_an_unknown_mode_is_refused_too(self):
        genotype, phenotype = _tiny()
        with self.assertRaises(ValueError):
            run_linear_gwas(genotype, phenotype, reduce="best-of",
                            output_dir="/tmp/should-not-be-created")


class MachineryStillWorksTestCase(unittest.TestCase):
    """The spellings remain usable INTERNALLY; only the entry point refuses."""

    def test_variant_reduction_still_constructs(self):
        # SignificantPairs is built on this, and the trait-blocked paths
        # construct it directly. Removing it would break the two modes that
        # are meant to survive.
        for mode in ("max-abs-t", "max-t2", "min-p"):
            self.assertIsNotNone(VariantReduction(mode))
        self.assertIsNotNone(VariantReduction("top-k", 4))

    def test_top_k_still_requires_its_k(self):
        with self.assertRaises(ValueError):
            VariantReduction("top-k")


class CliSurfaceTestCase(unittest.TestCase):
    def test_cli_offers_the_two_modes_and_not_the_machinery(self):
        from torchgwas.cli import _build_parser

        parser = _build_parser()
        # Parse rather than introspect argparse internals: this checks what a
        # user's command line actually does, which is the thing that was broken
        # (the CLI offered four machinery spellings and NEITHER real mode).
        for mode in ("significant", "jagwas"):
            with self.subTest(mode=mode):
                args = parser.parse_args(
                    ["linear", "--genotype", "g", "--phenotype", "p",
                     "--output-dir", "o", "--reduce", mode])
                self.assertEqual(args.reduce, mode)

    def test_cli_refuses_the_machinery_spellings(self):
        from torchgwas.cli import _build_parser

        parser = _build_parser()
        for mode in ("max-abs-t", "max-t2", "min-p", "top-k"):
            with self.subTest(mode=mode):
                with self.assertRaises(SystemExit):
                    parser.parse_args(
                        ["linear", "--genotype", "g", "--phenotype", "p",
                         "--output-dir", "o", "--reduce", mode])


if __name__ == "__main__":
    unittest.main()
