import numpy as np
import pytest
from torchgwas.sumstats_indexed import write_indexed_sumstats, open_indexed_sumstats


def test_selected_indices_labels_and_t_only(tmp_path):
    # Global indices must survive selection and nonzero chunk starts.
    chunk = (2, 4, np.array([3, 2]), np.array([1, 0]),
             np.array([.5, -.25]), np.array([4., -3.]), 20)
    labels = ["plain", "tab\tlabel", "newline\nlabel", "unicode-β"]
    rows, _ = write_indexed_sumstats(tmp_path, labels, ["a", "b"], 24,
        [chunk], kind="significant", df=20, store_beta=False)
    manifest, parts = open_indexed_sumstats(tmp_path)
    part = next(parts)
    assert rows == manifest["rows"] == 2
    np.testing.assert_array_equal(part["variant_index"], [3, 2])
    np.testing.assert_array_equal(part["trait_index"], [1, 0])
    assert "beta" not in part
    assert np.load(tmp_path / "variant_ids.npy").tolist() == labels
    assert not list(tmp_path.glob("*.tsv*"))


def test_topk_across_chunks_and_threshold(tmp_path):
    chunks = [(0, 2, np.ones((2, 2)), np.array([[2., 4.], [5., np.nan]]), None),
              (2, 3, np.ones((1, 2)), np.array([[6., 8.]]), None)]
    rows, _ = write_indexed_sumstats(tmp_path, ["a", "b", "c"], ["x", "y"],
        24, chunks, kind="filtered", df=20, topk_per_trait=1,
        p_value_threshold=.05)
    _, parts = open_indexed_sumstats(tmp_path)
    part = next(parts)
    assert rows == 2
    np.testing.assert_array_equal(part["variant_index"], [2, 2])
    np.testing.assert_array_equal(part["t_stat"], [6., 8.])


def test_empty_selection_and_joint_statistics(tmp_path):
    empty = tmp_path / "empty"
    rows, _ = write_indexed_sumstats(empty, ["a"], ["x"], 24,
        [(0, 1, np.ones((1, 1)), np.full((1, 1), np.nan), None)],
        kind="filtered", df=20)
    manifest, parts = open_indexed_sumstats(empty)
    assert rows == 0 and list(parts) == [] and manifest["shape"] == [1, 1]
    joint = tmp_path / "joint"
    write_indexed_sumstats(joint, ["a", "b"], ["x", "y"], 24,
        [(0, 2, None, np.array([7., np.nan]), None)], kind="jagwas", df=20, chi2_df=2)
    manifest, parts = open_indexed_sumstats(joint)
    assert manifest["df"] == 2 and manifest["kind"] == "jagwas"
    np.testing.assert_array_equal(next(parts)["chi2"], [7.])
