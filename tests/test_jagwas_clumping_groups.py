import argparse
from pathlib import Path
import sys

sys.path.insert(0, "scripts")

import numpy as np
import pytest

import run_jagwas_clumping as wrapper


def test_phenotype_group_needs_a_safe_name():
    assert wrapper.parse_phenotype_group("CNN_v2=/x/CNN.npy") == ("CNN_v2", Path("/x/CNN.npy"))
    for bad in ("CNN", "=/x.npy", "a/b=/x.npy", "torchgwas=/x.npy", "a="):
        with pytest.raises(argparse.ArgumentTypeError):
            wrapper.parse_phenotype_group(bad)


def test_group_panel_concatenates_and_names_traits(tmp_path):
    first = np.arange(12, dtype=np.float32).reshape(4, 3)
    second = -np.arange(8, dtype=np.float32).reshape(4, 2)
    np.save(tmp_path / "first.npy", first)
    np.save(tmp_path / "second.npy", second)
    (tmp_path / "first.traits.txt").write_text("x\ny\nz\n")
    (tmp_path / "second.traits.txt").write_text("only-one\n")  # wrong length: ignored
    panel, names, columns = wrapper.group_panel(
        [("f", tmp_path / "first.npy"), ("s", tmp_path / "second.npy")]
    )
    np.testing.assert_array_equal(panel, np.column_stack([first, second]))
    assert names == ["f:x", "f:y", "f:z", "s:trait_0", "s:trait_1"]
    assert columns == [("f", [0, 1, 2]), ("s", [3, 4])]
    np.save(tmp_path / "short.npy", np.zeros((3, 2), dtype=np.float32))
    with pytest.raises(ValueError, match="has 3 rows"):
        wrapper.group_panel([("f", tmp_path / "first.npy"), ("t", tmp_path / "short.npy")])


def test_harmonize_keeps_only_eligible_rows_in_position_order():
    marker = np.array(["rs1", "rs2", "id3", "rs4", "rs5", "rs6"])
    chromosome = np.array(["2", "1", "1", "1", "6", "1"])
    position = np.array([50, 300, 200, 100, 30_000_000, 400], dtype=np.int64)
    effect = np.array(["a", "C", "G", "T", "A", "AT"])
    other = np.array(["G", "A", "C", "C", "G", "C"])
    maf = np.full(6, 0.2, dtype=np.float32)
    chi2 = np.array([40.0, 30.0, 50.0, 0.1, 60.0, 45.0])
    table, report = wrapper.harmonize(
        (marker, chromosome, position, effect, other, maf), np.arange(6), chi2,
        degrees_of_freedom=2, n_samples=100, gwas_p=0.05, maf_min=0.01, exclude_mhc=True,
    )
    # id3 (not an rsID), rs4 (P > 0.05), rs5 (MHC) and rs6 (indel) are removed.
    assert table["SNP"].tolist() == ["rs2", "rs1"]
    assert table["uniqID"].tolist() == ["1:300:A:C", "2:50:A:G"]
    assert report["jagwas_rows"] == 6 and report["harmonized_rows"] == 2
