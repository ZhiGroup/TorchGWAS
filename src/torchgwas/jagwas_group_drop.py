"""JAGWAS groups that each run on their own subjects (missing_phenotype='drop_subject').

A subject missing, or outlier-masked in, one of a group's traits leaves that
group only (api._drop_subjects_by_group marks the group's rows missing
across its traits, so every trait of a group shares one set of samples).
The scan then tests each trait on its own rows (complete_case.CompleteCasePlan):
t and each pair's df are those of the group's subjects.

Here each group's factor is prepared on its kept rows alone -- R, the
sample count and the rounding precision are exactly a run of that group
alone -- and its score z = t / sqrt(1 + t^2 / df) uses the group's pair df.

A separate module, not a change to jagwas_projection.py: that file's bytes
identify the recorded factor calibrations (reduction_tensor_work's
source_sha256), as reduce.py's identify the selector census.
"""
from __future__ import annotations

import numpy as np
import torch

from .jagwas_projection import JagwasGroups


class GroupDropJagwas(JagwasGroups):
    """JagwasGroups whose groups each keep only the subjects their own traits observe."""

    takes_pair_df = True

    def prepare(self, phenotype, device=None, plan=None):
        """Each group's factor on its kept rows. plan: the scan's CompleteCasePlan (None: every row)."""
        matrix = torch.as_tensor(phenotype)
        if device is not None:
            matrix = matrix.to(device)
        samples, traits = matrix.shape
        for name, columns in zip(self.names, self.columns):
            if columns.min() < 0 or columns.max() >= traits:
                raise ValueError(f"jagwas group {name} refers outside the {traits} scanned traits")
        self._indices = []
        # Each group's first column as a host int: reading it from the device
        # index inside reduce() would sync the compute stream once per group
        # per chunk.
        self._first_columns = []
        for name, columns, reduction in zip(self.names, self.columns, self.reductions):
            missing = _group_missing_rows(plan, name, columns)
            kept = np.setdiff1d(np.arange(samples), missing)
            rows = torch.as_tensor(kept, dtype=torch.int64, device=matrix.device)
            index = torch.as_tensor(columns, dtype=torch.int64, device=matrix.device)
            reduction.prepare(matrix.index_select(0, rows).index_select(1, index))
            self._indices.append(index)
            self._first_columns.append(int(columns[0]))
        return self

    def reduce(self, beta, t, status, variant_df, width, pair_df=None):
        """(chunk, K) t of the whole panel in, (chunk, groups) T out; each group at its own df."""
        statistic = []
        for reduction, index, first in zip(self.reductions, self._indices, self._first_columns):
            group_df = variant_df if pair_df is None else pair_df[:, first].to(variant_df.dtype)
            statistic.append(reduction.reduce(beta, t.index_select(1, index), status, group_df, 1)[1])
        statistic = torch.cat(statistic, dim=1)
        return (torch.full_like(statistic, float("nan"), dtype=beta.dtype),
                statistic,
                torch.zeros_like(statistic, dtype=torch.int32),
                status, variant_df)


def _group_missing_rows(plan, name, columns):
    """The rows a group's traits lack, which must be the same rows for all of them."""
    if plan is None:
        return np.empty(0, dtype=np.int64)
    position = {int(trait): i for i, trait in enumerate(plan.traits)}
    patterns = {position.get(int(column)) for column in columns}
    rows = {None if index is None else int(plan.pattern_of_trait[index]) for index in patterns}
    if len(rows) != 1:
        raise ValueError(f"jagwas group {name}'s traits must share their samples")
    pattern = rows.pop()
    return np.empty(0, dtype=np.int64) if pattern is None else np.asarray(plan.rows[pattern], dtype=np.int64)
