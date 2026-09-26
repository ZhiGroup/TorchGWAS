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


def test_each_panel_drops_only_its_own_outlier_rows(tmp_path):
    rng = np.random.default_rng(2)
    first, second = rng.standard_normal((400, 3)), rng.standard_normal((400, 2))
    first[7, 1] = 40.0
    second[123, 0] = -40.0
    np.save(tmp_path / "first.npy", first)
    np.save(tmp_path / "second.npy", second)
    covariates = rng.standard_normal((400, 2))
    groups = [("f", tmp_path / "first.npy"), ("s", tmp_path / "second.npy")]
    panels = wrapper.load_panels([path for _name, path in groups])
    args = argparse.Namespace(panels=panels, outlier_rows={
        name: wrapper.phenotype_outliers(panels[path], covariates, 5.0) for name, path in groups})
    assert [np.flatnonzero(args.outlier_rows[name]).tolist() for name in "fs"] == [[7], [123]]
    panel, _names, columns = wrapper.masked_panel(args, groups)
    # A group's whole row goes, in its own columns only; the other group keeps the sample.
    assert np.isnan(panel[7, :3]).all() and np.isfinite(panel[7, 3:]).all()
    assert np.isnan(panel[123, 3:]).all() and np.isfinite(panel[123, :3]).all()
    assert int(np.isnan(panel).sum()) == 3 + 2
    assert wrapper.parse_group_rcond("g=1e-3") == ("g", 1e-3)
    with pytest.raises(argparse.ArgumentTypeError):
        wrapper.parse_group_rcond("g=2")


def _options(**changes):
    values = dict(phenotype_outlier_sd=None, jagwas_min_residual=None, jagwas_rcond=None, reuse_scan=False)
    return argparse.Namespace(**dict(values, **changes))


def test_the_pipeline_defaults_to_outlier_rows_and_eigen_truncation():
    args = _options()
    wrapper.resolve_defaults(args)
    assert (args.phenotype_outlier_sd, args.jagwas_rcond, args.jagwas_min_residual) == (5.0, 1e-3, None)
    off = _options(phenotype_outlier_sd=0.0, jagwas_rcond=0.0)
    wrapper.resolve_defaults(off)
    assert (off.phenotype_outlier_sd, off.jagwas_rcond) == (None, 0.0)  # 0: the rounding cutoff
    traits = _options(jagwas_min_residual=1e-2)
    wrapper.resolve_defaults(traits)
    assert (traits.jagwas_rcond, traits.jagwas_min_residual) == (None, 1e-2)
    reuse = _options(reuse_scan=True)
    wrapper.resolve_defaults(reuse)
    assert (reuse.phenotype_outlier_sd, reuse.jagwas_rcond) == (None, None)
    for bad in (_options(reuse_scan=True, phenotype_outlier_sd=4.0), _options(jagwas_rcond=2.0),
                _options(jagwas_rcond=1e-3, jagwas_min_residual=1e-2)):
        with pytest.raises(ValueError):
            wrapper.resolve_defaults(bad)


def test_group_cutoffs_carry_each_groups_rule():
    args = argparse.Namespace(group_rcond={"e": 1e-3}, group_min_residual={"t": 1e-2})
    assert wrapper.group_cutoffs(args, [("e", [0, 1]), ("t", [2]), ("d", [3])]) == [
        ("e", [0, 1], {"rcond": 1e-3}), ("t", [2], {"min_residual": 1e-2}), ("d", [3], None)]


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
