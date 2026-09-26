import sys

sys.path.insert(0, "scripts")

import numpy as np
import pytest

from jagwas_streaming import ClumpingPrep, GroupedClumpers, build_clumping_cache, validate_clumping_cache


def test_compact_cache_applies_only_static_clumping_filters(tmp_path):
    cache = tmp_path / "cache"
    manifest = build_clumping_cache(
        cache,
        marker_ids=np.array(["rs1", "rs2", "id3", "rs4", "rs5", "rs6"]),
        variant_metadata={
            "chromosome": np.array(["1", "1", "1", "6", "23", "2"]),
            "position": np.array([100, 200, 300, 30_000_000, 100, 500]),
            "effect_allele": np.array(["A", "AT", "C", "G", "T", "C"]),
            "other_allele": np.array(["G", "C", "A", "A", "A", "G"]),
        },
        maf=np.array([0.1, 0.2, 0.1, 0.1, 0.1, 0.001]),
        maf_min=0.01,
        exclude_mhc=True,
        mhc_start=29_614_758,
        mhc_end=33_170_276,
    )
    assert manifest["eligible_variants"] == 1
    prep = ClumpingPrep(cache)
    rows = prep.rows(
        0,
        np.array([10.0, 9.0, 8.0, 7.0, 6.0, 5.0]),
        degrees_of_freedom=2,
        gwas_p=1.0,
        n_samples=40,
    )
    assert rows[["SNP", "CHR", "POS", "A1", "A2", "N"]].to_dict("records") == [
        {"SNP": "rs1", "CHR": "1", "POS": 100, "A1": "A", "A2": "G", "N": 40}
    ]
    assert rows.iloc[0]["uniqID"] == "1:100:A:G"


def test_cache_validation_refuses_a_different_filter(tmp_path):
    cache = tmp_path / "cache"
    build_clumping_cache(
        cache,
        marker_ids=np.array(["rs1"]),
        variant_metadata={
            "chromosome": np.array(["1"]),
            "position": np.array([100]),
            "effect_allele": np.array(["A"]),
            "other_allele": np.array(["G"]),
        },
        maf=np.array([0.1]),
        maf_min=0.01,
        exclude_mhc=True,
        mhc_start=29_614_758,
        mhc_end=33_170_276,
    )
    with pytest.raises(ValueError, match="does not match"):
        validate_clumping_cache(
            cache, n_variants=1, maf_min=0.05, exclude_mhc=True
        )


class _RecordingClumper:
    def __init__(self, log, index):
        self.log, self.index, self.rank = log, index, None

    def set_rank(self, report):
        self.rank = report

    def consume(self, chunk):
        self.log.append(("consume", self.index, chunk[0], chunk[1], np.asarray(chunk[3]).tolist()))

    def close(self):
        self.log.append(("close", self.index))

    def collect(self, timeout):
        self.log.append(("collect", self.index))
        return {"ok": True, "index": self.index}

    def abort(self):
        self.log.append(("abort", self.index))


def test_grouped_clumpers_route_columns_and_finish_every_group_together():
    log = []
    clumpers = GroupedClumpers([_RecordingClumper(log, 0), _RecordingClumper(log, 1)])
    clumpers.set_rank([{"rank": 3}, {"rank": 5}])
    assert [clumper.rank for clumper in clumpers.clumpers] == [{"rank": 3}, {"rank": 5}]
    clumpers.consume((10, 12, None, np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)))
    assert log == [("consume", 0, 10, 12, [1.0, 3.0]), ("consume", 1, 10, 12, [2.0, 4.0])]
    log.clear()
    reports = clumpers.finish()
    # Every worker is told the scan is over before any is waited on.
    assert log == [("close", 0), ("close", 1), ("collect", 0), ("collect", 1)]
    assert [report["index"] for report, _done in reports] == [0, 1]
    with pytest.raises(ValueError, match="expected 2 JAGWAS group reports"):
        clumpers.set_rank([{"rank": 3}])
    with pytest.raises(ValueError, match="expected 2 JAGWAS group columns"):
        clumpers.consume((0, 1, None, np.zeros((1, 3))))
