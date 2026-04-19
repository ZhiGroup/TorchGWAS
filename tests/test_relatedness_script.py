from __future__ import annotations

import tempfile
import unittest
import importlib.util
from pathlib import Path

import pandas as pd

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_linear_unrelated.py"
SPEC = importlib.util.spec_from_file_location("run_linear_unrelated", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RelatednessScriptTestCase(unittest.TestCase):
    def test_filter_table_to_ids_preserves_keep_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            table_path = Path(tmpdir) / "pheno.tsv"
            out_path = Path(tmpdir) / "pheno_filtered.tsv"
            pd.DataFrame(
                {
                    "IID": ["s3", "s1", "s2"],
                    "trait": [3.0, 1.0, 2.0],
                }
            ).to_csv(table_path, sep="\t", index=False)
            filtered = MODULE._filter_table_to_ids(table_path, ["s2", "s1"], out_path)
            self.assertEqual(filtered["IID"].tolist(), ["s2", "s1"])
            self.assertTrue(out_path.exists())

    def test_write_keep_file_uses_fid_iid_pairs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            keep_path = MODULE._write_keep_file(["1001", "1002"], Path(tmpdir) / "samples.keep")
            self.assertEqual(keep_path.read_text(), "1001\t1001\n1002\t1002\n")


if __name__ == "__main__":
    unittest.main()
