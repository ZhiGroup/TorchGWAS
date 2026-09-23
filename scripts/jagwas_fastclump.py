"""Legacy-compatible accelerators for the local FUMA-style clumper.

These are the two measured hot-path fixes from the optimized legacy JAGWAS
workflow: NumPy-backed chromosome stage 1 and cached CSR-row lookups in stage 2.
They patch the lab-local ``fuma_clump`` module only inside the clumping worker.
"""

from __future__ import annotations

import numpy as np


class Row(dict):
    @property
    def index(self):
        return self.keys()


def _uid(chromosome, position, first, second):
    if first <= second:
        return f"{chromosome}:{position}:{first}:{second}"
    return f"{chromosome}:{position}:{second}:{first}"


def fast_process_chromosome_with_ld(chrom, chr_gwas, ld_table, params):
    if chr_gwas.empty:
        return [], {}, []

    gwas_set = set(chr_gwas.loc[chr_gwas["P"] <= params["gwasP"], "SNP"])
    deduplicated = chr_gwas.drop_duplicates("SNP")
    d_snp = deduplicated["SNP"].to_numpy()
    d_pos = deduplicated["POS"].to_numpy()
    d_a1 = deduplicated["A1"].to_numpy()
    d_a2 = deduplicated["A2"].to_numpy()
    d_p = deduplicated["P"].to_numpy()
    d_chrom = deduplicated["CHR"].to_numpy()
    d_af = deduplicated["AF1"].to_numpy() if "AF1" in deduplicated.columns else None
    row_of = {snp: index for index, snp in enumerate(d_snp)}

    significant = (
        chr_gwas[chr_gwas["P"] <= params["leadP"]]
        .drop_duplicates("SNP")
        .sort_values("P")
    )
    if significant.empty:
        return [], {}, []
    s_snp = significant["SNP"].to_numpy()
    s_pos = significant["POS"].to_numpy()
    s_a1 = significant["A1"].to_numpy()
    s_a2 = significant["A2"].to_numpy()
    s_p = significant["P"].to_numpy()
    s_chrom = significant["CHR"].to_numpy()
    s_af = significant["AF1"].to_numpy() if "AF1" in significant.columns else None
    print(
        f"  [CHR{chrom}] Using pre-computed LD table "
        f"({len(significant)} candidate index SNPs)..."
    )

    def row(c, p, a1, a2, p_value, snp, af):
        result = Row(CHR=c, POS=int(p), A1=a1, A2=a2, P=float(p_value), SNP=snp)
        if af is not None:
            result["AF1"] = af
        return result

    assigned = set()
    independent = []
    candidates = {}
    ld_pairs = []
    for index, rsid in enumerate(s_snp):
        if rsid in assigned:
            continue
        assigned.add(rsid)
        independent_uid = _uid(chrom, int(s_pos[index]), s_a1[index], s_a2[index])
        candidates[rsid] = {
            "r2": 1.0,
            "ind_sig_rsID": rsid,
            "gwas_row": row(
                s_chrom[index], s_pos[index], s_a1[index], s_a2[index],
                s_p[index], rsid, None if s_af is None else s_af[index],
            ),
        }
        ld_pairs.append((independent_uid, independent_uid, 1.0))
        n_clumped = 1
        for neighbor, r2 in ld_table.get_neighbors(rsid, params["r2"]):
            if neighbor not in gwas_set or neighbor in assigned:
                continue
            neighbor_index = row_of.get(neighbor)
            if neighbor_index is None:
                continue
            assigned.add(neighbor)
            ld_pairs.append(
                (
                    independent_uid,
                    _uid(
                        chrom, int(d_pos[neighbor_index]),
                        d_a1[neighbor_index], d_a2[neighbor_index],
                    ),
                    r2,
                )
            )
            n_clumped += 1
            current = candidates.get(neighbor)
            if current is None or current["r2"] < r2:
                candidates[neighbor] = {
                    "r2": r2,
                    "ind_sig_rsID": rsid,
                    "gwas_row": row(
                        d_chrom[neighbor_index], d_pos[neighbor_index],
                        d_a1[neighbor_index], d_a2[neighbor_index],
                        d_p[neighbor_index], neighbor,
                        None if d_af is None else d_af[neighbor_index],
                    ),
                }
        independent.append(
            {
                "rsID": rsid,
                "chr": int(chrom),
                "pos": int(s_pos[index]),
                "p": float(s_p[index]),
                "A1": s_a1[index],
                "A2": s_a2[index],
                "uniqID": independent_uid,
                "nSNPs": n_clumped,
                "nGWASSNPs": n_clumped,
            }
        )
    print(f"  [CHR{chrom}] {len(independent)} IndSigSNP(s) ({len(candidates)} candidates)")
    return independent, candidates, ld_pairs


def fast_get_r2(self, snp_a, snp_b):
    first = self._s2i.get(snp_a, -1)
    second = self._s2i.get(snp_b, -1)
    if first < 0 or second < 0:
        return 0.0
    cache = self.__dict__.setdefault("_rowcache", {})
    row = cache.get(first)
    if row is None:
        start, end = int(self._indptr[first]), int(self._indptr[first + 1])
        row = dict(
            zip(
                np.asarray(self._ib[start:end]).tolist(),
                np.asarray(self._r2[start:end]).tolist(),
            )
        )
        cache[first] = row
    return float(row.get(second, 0.0))


def install(*, fast_stage1: bool = True):
    import fuma_clump

    if fast_stage1:
        fuma_clump.process_chromosome_with_ld = fast_process_chromosome_with_ld
    fuma_clump.LDTable.get_r2 = fast_get_r2
    return fuma_clump
