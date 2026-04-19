from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class CLITestCase(unittest.TestCase):
    def test_demo_command(self):
        repo = Path(__file__).resolve().parents[1]
        env = dict(__import__("os").environ)
        env["PYTHONPATH"] = str(repo / "src")
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [sys.executable, "-m", "torchgwas", "demo", "--output-dir", tmpdir]
            subprocess.run(cmd, cwd=repo, env=env, check=True)
            summary = json.loads((Path(tmpdir) / "run_summary.json").read_text())
            self.assertGreater(summary["linear_rows"], 0)
            self.assertGreater(summary["multi_rows"], 0)

    def test_prep_command_with_tables(self):
        repo = Path(__file__).resolve().parents[1]
        env = dict(__import__("os").environ)
        env["PYTHONPATH"] = str(repo / "src")
        toy = repo / "examples" / "toy"
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                sys.executable,
                "-m",
                "torchgwas",
                "prep",
                "--genotype",
                str(toy / "genotype.npy"),
                "--phenotype-table",
                str(toy / "pheno.tsv"),
                "--covariates-table",
                str(toy / "covar.tsv"),
                "--sample-ids",
                str(toy / "samples.tsv"),
                "--output-dir",
                tmpdir,
            ]
            subprocess.run(cmd, cwd=repo, env=env, check=True)
            self.assertTrue((Path(tmpdir) / "phenotype_processed.npy").exists())
            self.assertTrue((Path(tmpdir) / "prep.json").exists())


if __name__ == "__main__":
    unittest.main()
