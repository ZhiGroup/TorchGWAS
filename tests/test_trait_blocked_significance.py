"""Trait blocking must not change which (variant, trait) pairs are significant.

Blocking is chosen from whatever memory the card happens to have, so if it
altered the result set the answers would depend on the hardware. These compare
a blocked run against an unblocked one on the same inputs and require the rows
to match exactly.

This path had no coverage at all, which is how two defects reached a run:
`reduce='significant'` leaves `reduction` as None, so automatic blocking never
fired for it (a 600,000-voxel scan silently put a 79.4 GB design on one card
and left seven idle); and when blocking was finally forced it crashed with
`not enough values to unpack (expected 6, got 5)`, because the blocked driver
was written for a `VariantReduction`'s one-row-per-variant tuples and
significance emits a variable number of pairs with no merge step.
"""
from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch

from torchgwas.api import run_linear_gwas
from binary_helpers import binary_rows


def read_rows(directory: Path) -> list[tuple]:
    rows = binary_rows(directory)
    # Order is not part of the contract: a blocked run emits block by block and
    # a sharded one interleaves whatever its threads finish first. The SET of
    # pairs is the contract, so both sides are sorted before comparison.
    return sorted(
        (row["marker_id"], row["trait"], round(float(row["beta"]), 6),
         round(float(row["t_stat"]), 6))
        for row in rows)


class TraitBlockedSignificanceTestCase(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20260914)
        self.n, self.m, self.k = 200, 60, 40
        self.genotype = rng.integers(0, 3, size=(self.n, self.m)).astype(np.float32)
        self.covariates = rng.normal(size=(self.n, 3))
        phenotype = rng.normal(size=(self.n, self.k))
        # Plant real signal so the threshold actually keeps rows; a test that
        # compares two empty files proves nothing.
        for trait in range(0, self.k, 4):
            phenotype[:, trait] += 3.0 * self.genotype[:, trait % self.m]
        self.phenotype = phenotype
        self.sample_ids = [f"s{i}" for i in range(self.n)]
        self.marker_ids = [f"m{j}" for j in range(self.m)]

    def run_scan(self, tmpdir, name, **kwargs):
        out = Path(tmpdir) / name
        run_linear_gwas(
            self.genotype, self.phenotype, covariates=self.covariates,
            sample_ids=self.sample_ids, marker_ids=self.marker_ids,
            output_dir=out, sumstats_format="binary", reduce="significant",
            significance_threshold=0.05, chunk_size=8, **kwargs)
        return read_rows(out)

    def test_blocked_matches_unblocked(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            whole = self.run_scan(tmpdir, "whole")
            self.assertGreater(len(whole), 0, "planted signal produced no rows")
            for block in (1, 3, 7, self.k - 1, self.k):
                with self.subTest(trait_block=block):
                    blocked = self.run_scan(tmpdir, f"block{block}",
                                            trait_block=block)
                    self.assertEqual(blocked, whole)

    def test_threshold_is_over_the_whole_trait_axis(self):
        """A block must not be Bonferroni-corrected as if it were the run.

        Correcting within the block would use a looser threshold the smaller
        the block, so a narrow block would emit strictly more pairs. Equality
        across block widths above is the strong form of this; this states the
        intent separately so a future change that reintroduces per-block
        correction fails with a message that says why.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            whole = self.run_scan(tmpdir, "whole")
            narrow = self.run_scan(tmpdir, "narrow", trait_block=1)
            self.assertEqual(len(narrow), len(whole))

    @unittest.skipUnless(torch.cuda.is_available() and torch.cuda.device_count() > 1,
                         "needs more than one CUDA device")
    def test_sharded_across_devices_matches_unblocked(self):
        devices = [f"cuda:{index}" for index in range(torch.cuda.device_count())]
        with tempfile.TemporaryDirectory() as tmpdir:
            whole = self.run_scan(tmpdir, "whole", device="cuda:0")
            sharded = self.run_scan(tmpdir, "sharded", device="cuda:0",
                                    trait_block=7, trait_devices=devices)
            self.assertEqual(sharded, whole)


if __name__ == "__main__":
    unittest.main()
