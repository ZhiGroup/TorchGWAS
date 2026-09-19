from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


def _write_plink_triplet(prefix: Path, dosage_a1: np.ndarray) -> Path:
    dosage_a1 = np.asarray(dosage_a1, dtype=np.float32)
    n_samples, n_markers = dosage_a1.shape
    bed = Path(f"{prefix}.bed")
    Path(f"{prefix}.fam").write_text(
        "".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_samples)),
        encoding="utf-8",
    )
    Path(f"{prefix}.bim").write_text(
        "".join(f"1 rs{j} 0 {100 + j} A C\n" for j in range(n_markers)),
        encoding="utf-8",
    )
    dosage_to_code = {0.0: 0b00, 1.0: 0b10, 2.0: 0b11}
    payload = bytearray(b"\x6c\x1b\x01")
    for marker in range(n_markers):
        for sample_start in range(0, n_samples, 4):
            byte = 0
            for offset in range(4):
                sample = sample_start + offset
                if sample >= n_samples or np.isnan(dosage_a1[sample, marker]):
                    code = 0b01
                else:
                    code = dosage_to_code[float(dosage_a1[sample, marker])]
                byte |= code << (2 * offset)
            payload.append(byte)
    bed.write_bytes(bytes(payload))
    return bed


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
            self.assertNotIn("multi_rows", summary)

    def test_prep_command_with_tables(self):
        repo = Path(__file__).resolve().parents[1]
        env = dict(__import__("os").environ)
        env["PYTHONPATH"] = str(repo / "src")
        toy = repo / "examples" / "toy"
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = _write_plink_triplet(
                Path(tmpdir) / "genotype",
                np.load(toy / "genotype.npy", allow_pickle=False),
            )
            cmd = [
                sys.executable,
                "-m",
                "torchgwas",
                "prep",
                "--genotype",
                str(bed),
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

    def test_numpy_genotype_format_is_not_public(self):
        repo = Path(__file__).resolve().parents[1]
        env = dict(__import__("os").environ)
        env["PYTHONPATH"] = str(repo / "src")
        toy = repo / "examples" / "toy"
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                sys.executable,
                "-m",
                "torchgwas",
                "linear",
                "--genotype",
                str(toy / "genotype.npy"),
                "--genotype-format",
                "npy",
                "--phenotype-table",
                str(toy / "pheno.tsv"),
                "--output-dir",
                tmpdir,
            ]
            completed = subprocess.run(
                cmd,
                cwd=repo,
                env=env,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("invalid choice: 'npy'", completed.stderr)


if __name__ == "__main__":
    unittest.main()
