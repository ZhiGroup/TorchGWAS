from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from torchgwas.api import run_linear_gwas, run_multivariate_gwas
from torchgwas.datasets import get_toy_dataset_paths
from torchgwas.preprocess import prepare_inputs_for_prep


class APITestCase(unittest.TestCase):
    def test_linear_and_multivariate_return_rows(self):
        toy = get_toy_dataset_paths(Path(__file__).resolve().parents[1] / "examples" / "toy")
        linear = run_linear_gwas(
            genotype=toy["genotype"],
            phenotype=toy["phenotype"],
            covariates=toy["covariates"],
            marker_ids=toy["marker_ids"],
            sample_ids=toy["sample_ids"],
            chunk_size=3,
        )
        multi = run_multivariate_gwas(
            genotype=toy["genotype"],
            phenotype=toy["phenotype"],
            covariates=toy["covariates"],
            marker_ids=toy["marker_ids"],
            sample_ids=toy["sample_ids"],
            chunk_size=3,
        )
        self.assertEqual(len(linear.table), 12 * 3)
        self.assertEqual(len(multi.table), 12)
        self.assertIn("beta", linear.table[0])
        self.assertIn("chi2", multi.table[0])

    def test_tabular_inputs_align_to_sample_ids(self):
        toy = get_toy_dataset_paths(Path(__file__).resolve().parents[1] / "examples" / "toy")
        linear = run_linear_gwas(
            genotype=toy["genotype"],
            phenotype=None,
            covariates=None,
            phenotype_table=toy["phenotype_table"],
            covariates_table=toy["covariates_table"],
            sample_ids=toy["sample_ids"],
            marker_ids=toy["marker_ids"],
            trait_columns=["trait_0", "trait_1", "trait_2"],
            covariate_columns=["covar_0", "covar_1", "covar_2", "covar_3"],
            chunk_size=4,
        )
        self.assertEqual(len(linear.table), 36)
        self.assertEqual(linear.run_metadata["trait_columns"], ["trait_0", "trait_1", "trait_2"])

    def test_output_dir_writes_expected_files(self):
        toy = get_toy_dataset_paths(Path(__file__).resolve().parents[1] / "examples" / "toy")
        with tempfile.TemporaryDirectory() as tmpdir:
            run_linear_gwas(
                genotype=toy["genotype"],
                phenotype=toy["phenotype"],
                covariates=toy["covariates"],
                marker_ids=toy["marker_ids"],
                sample_ids=toy["sample_ids"],
                output_dir=tmpdir,
                chunk_size=5,
            )
            base = Path(tmpdir)
            self.assertTrue((base / "results.tsv.gz").exists())
            self.assertTrue((base / "run.json").exists())
            self.assertTrue((base / "qc.json").exists())
            with gzip.open(base / "results.tsv.gz", "rt") as handle:
                header = handle.readline().strip()
            self.assertIn("marker_id", header)
            metadata = json.loads((base / "run.json").read_text())
            self.assertEqual(metadata["analysis"], "linear")

    def test_prep_path_supports_chunked_genotype_qc(self):
        toy = get_toy_dataset_paths(Path(__file__).resolve().parents[1] / "examples" / "toy")
        genotype = np.load(toy["genotype"], mmap_mode="r")
        phenotype = np.load(toy["phenotype"])
        covariates = np.load(toy["covariates"])
        phenotype, covariates, qc = prepare_inputs_for_prep(
            genotype,
            phenotype,
            covariates,
            genotype_chunk_size=5,
        )
        self.assertEqual(phenotype.shape[0], genotype.shape[0])
        self.assertEqual(covariates.shape[0], genotype.shape[0])
        self.assertEqual(qc["genotype_qc_mode"], "chunked")
        self.assertEqual(qc["genotype_qc_chunk_size"], 5)


if __name__ == "__main__":
    unittest.main()
