"""Overlap current TorchGWAS JAGWAS chunks with locus clumping.

The scan process writes each narrow chi-square chunk to ``/dev/shm`` and sends
only its path to a forked CPU worker. The worker computes tail probabilities,
builds rows from a compact precomputed variant cache, and completes clumping
stage 1 as soon as a chromosome has passed through the scan. Final clumping
stages therefore wait for neither whole-genome harmonization nor a dense
summary-statistic table.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
from pathlib import Path
import queue
import shutil
import time
import traceback
import uuid

import numpy as np
import pandas as pd
from scipy import special


CACHE_FORMAT = "torchgwas-jagwas-clump-prep-v1"


def _single_base_codes(values) -> tuple[np.ndarray, np.ndarray]:
    encoded = np.empty(len(values), dtype="S1")
    valid = np.zeros(len(values), dtype=bool)
    for index, value in enumerate(values):
        text = str(value).upper()
        if len(text) == 1 and text in "ACGT":
            encoded[index] = text.encode("ascii")
            valid[index] = True
        else:
            encoded[index] = b""
    return encoded, valid


def build_clumping_cache(
    directory: Path,
    *,
    marker_ids,
    variant_metadata: dict,
    maf: np.ndarray,
    maf_min: float,
    exclude_mhc: bool,
    mhc_start: int,
    mhc_end: int,
) -> dict:
    """Create a compact, mmap-friendly cache of clumping-eligible variants."""
    directory = Path(directory)
    complete = directory / "manifest.json"
    if complete.is_file():
        return validate_clumping_cache(
            directory,
            n_variants=len(marker_ids),
            maf_min=maf_min,
            exclude_mhc=exclude_mhc,
        )
    if directory.exists():
        raise ValueError(f"incomplete clumping cache already exists: {directory}")

    started = time.perf_counter()
    marker_ids = np.asarray(marker_ids, dtype=object)
    maf = np.asarray(maf, dtype=np.float32)
    chromosome_raw = np.asarray(variant_metadata["chromosome"])
    position = np.asarray(variant_metadata["position"], dtype=np.int64)
    effect_raw = np.asarray(variant_metadata["effect_allele"], dtype=object)
    other_raw = np.asarray(variant_metadata["other_allele"], dtype=object)
    n_variants = len(marker_ids)
    if not all(
        len(values) == n_variants
        for values in (maf, chromosome_raw, position, effect_raw, other_raw)
    ):
        raise ValueError("variant metadata, marker IDs and MAF are not aligned")

    try:
        chromosome = chromosome_raw.astype(np.int16)
    except (TypeError, ValueError):
        chromosome = pd.to_numeric(
            pd.Series(chromosome_raw), errors="coerce"
        ).fillna(0).to_numpy(dtype=np.int16)
    effect, valid_effect = _single_base_codes(effect_raw)
    other, valid_other = _single_base_codes(other_raw)
    valid_rsid = np.fromiter(
        (str(value).startswith("rs") for value in marker_ids),
        dtype=bool,
        count=n_variants,
    )
    keep = np.isfinite(maf) & (maf >= maf_min)
    keep &= (chromosome >= 1) & (chromosome <= 22)
    keep &= valid_effect & valid_other & valid_rsid
    if exclude_mhc:
        keep &= ~(
            (chromosome == 6)
            & (position >= mhc_start)
            & (position <= mhc_end)
        )
    source_index = np.flatnonzero(keep).astype(np.int64)

    staging = directory.with_name(
        f"{directory.name}.tmp-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    )
    staging.mkdir(parents=True)
    try:
        selected_marker_ids = marker_ids[source_index].astype(str).astype("S")
        np.save(staging / "source_index.npy", source_index, allow_pickle=False)
        np.save(staging / "chromosome.npy", chromosome[source_index].astype(np.int8), allow_pickle=False)
        np.save(staging / "position.npy", position[source_index], allow_pickle=False)
        np.save(staging / "maf.npy", maf[source_index], allow_pickle=False)
        np.save(staging / "snp.npy", selected_marker_ids, allow_pickle=False)
        np.save(staging / "effect_allele.npy", effect[source_index], allow_pickle=False)
        np.save(staging / "other_allele.npy", other[source_index], allow_pickle=False)
        chromosome_ends = {
            str(chrom): int(source_index[chromosome[source_index] == chrom].max()) + 1
            for chrom in range(1, 23)
            if np.any(chromosome[source_index] == chrom)
        }
        manifest = {
            "format": CACHE_FORMAT,
            "n_variants": int(n_variants),
            "eligible_variants": int(len(source_index)),
            "maf_min": float(maf_min),
            "exclude_mhc": bool(exclude_mhc),
            "mhc_start": int(mhc_start),
            "mhc_end": int(mhc_end),
            "chromosome_ends": chromosome_ends,
            "build_seconds": time.perf_counter() - started,
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        staging.replace(directory)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return manifest


def validate_clumping_cache(
    directory: Path,
    *,
    n_variants: int,
    maf_min: float,
    exclude_mhc: bool,
) -> dict:
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest.get("format") != CACHE_FORMAT:
        raise ValueError(f"unsupported clumping cache: {directory}")
    expected = {
        "n_variants": int(n_variants),
        "maf_min": float(maf_min),
        "exclude_mhc": bool(exclude_mhc),
    }
    mismatches = {
        key: (manifest.get(key), value)
        for key, value in expected.items()
        if manifest.get(key) != value
    }
    if mismatches:
        raise ValueError(f"clumping cache does not match this run: {mismatches}")
    return manifest


class ClumpingPrep:
    def __init__(self, directory: Path):
        self.directory = Path(directory)
        self.manifest = json.loads((self.directory / "manifest.json").read_text())
        self.source_index = np.load(
            self.directory / "source_index.npy", mmap_mode="r", allow_pickle=False
        )
        self.chromosome = np.load(
            self.directory / "chromosome.npy", mmap_mode="r", allow_pickle=False
        )
        self.position = np.load(
            self.directory / "position.npy", mmap_mode="r", allow_pickle=False
        )
        self.maf = np.load(self.directory / "maf.npy", mmap_mode="r", allow_pickle=False)
        self.snp = np.load(self.directory / "snp.npy", mmap_mode="r", allow_pickle=False)
        self.effect = np.load(
            self.directory / "effect_allele.npy", mmap_mode="r", allow_pickle=False
        )
        self.other = np.load(
            self.directory / "other_allele.npy", mmap_mode="r", allow_pickle=False
        )
        self.chromosome_ends = {
            int(chrom): int(end)
            for chrom, end in self.manifest["chromosome_ends"].items()
        }

    def rows(
        self,
        start: int,
        statistic: np.ndarray,
        *,
        degrees_of_freedom: int,
        gwas_p: float,
        n_samples: int,
    ) -> pd.DataFrame:
        statistic = np.asarray(statistic, dtype=np.float64).reshape(-1)
        end = start + len(statistic)
        left = int(np.searchsorted(self.source_index, start, side="left"))
        right = int(np.searchsorted(self.source_index, end, side="left"))
        columns = ["CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "P", "uniqID"]
        if right <= left:
            return pd.DataFrame(columns=columns)
        relative = np.asarray(self.source_index[left:right] - start, dtype=np.int64)
        values = statistic[relative]
        p_value = special.gammaincc(degrees_of_freedom / 2.0, values / 2.0)
        p_value = np.maximum(p_value, np.nextafter(np.float64(0), np.float64(1)))
        selected = np.flatnonzero(np.isfinite(p_value) & (p_value > 0.0) & (p_value <= gwas_p))
        if len(selected) == 0:
            return pd.DataFrame(columns=columns)
        compact = selected + left
        chromosome = np.asarray(self.chromosome[compact], dtype=np.int8)
        chromosome_text = chromosome.astype(str)
        snp = np.char.decode(np.asarray(self.snp[compact]))
        effect = np.char.decode(np.asarray(self.effect[compact]))
        other = np.char.decode(np.asarray(self.other[compact]))
        position = np.asarray(self.position[compact], dtype=np.int64)
        first = np.where(effect <= other, effect, other)
        second = np.where(effect <= other, other, effect)
        uniq = np.char.add(
            np.char.add(
                np.char.add(np.char.add(chromosome_text, ":"), position.astype(str)),
                np.char.add(":", first),
            ),
            np.char.add(":", second),
        )
        return pd.DataFrame(
            {
                "CHR": chromosome_text,
                "SNP": snp,
                "POS": position,
                "A1": effect,
                "A2": other,
                "N": np.int32(n_samples),
                "AF1": np.asarray(self.maf[compact], dtype=np.float32),
                "P": p_value[selected],
                "uniqID": uniq,
            }
        )


def _clumping_worker(control, done, config):
    try:
        import sys

        sys.path.insert(0, str(config["clumping_dir"]))
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from jagwas_fastclump import install

        # The legacy NumPy stage-1 replacement assumes a unique rsID identifies
        # one variant row. The current BGEN contains at least one duplicated
        # rsID with different alleles, for which that shortcut selected the
        # wrong record. Keep the safe CSR-row cache for stage 2, but run the
        # lab clumper's exact stage-1 implementation while it overlaps the scan.
        fuma_clump = install(fast_stage1=False)
        prep = ClumpingPrep(Path(config["cache_dir"]))
        setup_started = time.perf_counter()
        n_preloaded = fuma_clump.preload_ld(config["ld_dir"])
        setup_seconds = time.perf_counter() - setup_started
        params = {
            "leadP": config["lead_p"],
            "gwasP": config["gwas_p"],
            "r2": fuma_clump.R2,
            "r2_2": fuma_clump.R2_2,
            "maf": fuma_clump.MAF,
            "merge_dist": fuma_clump.MERGE_DIST,
            "clump_kb": fuma_clump.CLUMP_KB,
        }
        pending: dict[int, list[pd.DataFrame]] = {}
        frames: dict[int, pd.DataFrame] = {}
        stage1_cache = {}
        completed = set()
        build_seconds = 0.0
        stage1_seconds = 0.0
        rows_before_dedup = 0

        def finalize(chromosome: int):
            nonlocal build_seconds, stage1_seconds, rows_before_dedup
            if chromosome in completed:
                return
            completed.add(chromosome)
            started = time.perf_counter()
            pieces = pending.pop(chromosome, [])
            if pieces:
                frame = pd.concat(pieces, ignore_index=True)
                rows_before_dedup += len(frame)
                frame = frame.sort_values("P").drop_duplicates("uniqID")
                # Match the established harmonizer exactly, including its
                # tie order for distinct records sharing one position/rsID.
                # Sorting only POS uses a different single-key code path in
                # pandas and changed which duplicate rsID row stage 1 retained.
                frame = frame.sort_values(
                    ["CHR", "POS"],
                    key=lambda column: (
                        column.astype(int) if column.name == "CHR" else column
                    ),
                ).reset_index(drop=True)
            else:
                frame = pd.DataFrame(
                    columns=["CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "P", "uniqID"]
                )
            frames[chromosome] = frame
            build_seconds += time.perf_counter() - started
            started = time.perf_counter()
            if not frame.empty and (frame["P"] <= config["lead_p"]).any():
                stage1_cache[chromosome] = fuma_clump.process_chromosome_with_ld(
                    chromosome,
                    frame.copy(),
                    fuma_clump.LDTable().load("", chromosome),
                    params,
                )
            else:
                stage1_cache[chromosome] = ([], {}, [])
            stage1_seconds += time.perf_counter() - started

        while True:
            item = control.get()
            if item is None:
                break
            start, end, path = item
            try:
                statistic = np.load(path, allow_pickle=False)
            finally:
                Path(path).unlink(missing_ok=True)
            started = time.perf_counter()
            rows = prep.rows(
                int(start),
                statistic,
                degrees_of_freedom=config["degrees_of_freedom"],
                gwas_p=config["gwas_p"],
                n_samples=config["n_samples"],
            )
            build_seconds += time.perf_counter() - started
            if not rows.empty:
                for chromosome, frame in rows.groupby("CHR", sort=False):
                    pending.setdefault(int(chromosome), []).append(frame)
            for chromosome, chromosome_end in prep.chromosome_ends.items():
                if chromosome_end <= end:
                    finalize(chromosome)

        for chromosome in sorted(set(pending).union(prep.chromosome_ends)):
            finalize(chromosome)
        nonempty_frames = [
            frames[chromosome]
            for chromosome in sorted(frames)
            if not frames[chromosome].empty
        ]
        combined = (
            pd.concat(nonempty_frames, ignore_index=True)
            if nonempty_frames
            else pd.DataFrame(
                columns=["CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "P", "uniqID"]
            )
        )
        combined.attrs["gwasP"] = float(config["gwas_p"])
        combined.attrs["n_input"] = int(config["n_variants"])
        finalize_started = time.perf_counter()
        original = fuma_clump.process_chromosome_with_ld
        fuma_clump.process_chromosome_with_ld = (
            lambda chrom, _gwas, _table, _params: stage1_cache.get(
                int(chrom), ([], {}, [])
            )
        )
        try:
            fuma_clump.run_df(
                combined,
                config["output_dir"],
                config["ld_dir"],
                lead_p=config["lead_p"],
                gwas_p=config["gwas_p"],
                no_mhc=False,
            )
        finally:
            fuma_clump.process_chromosome_with_ld = original
        done.put(
            {
                "ok": True,
                "setup_seconds": setup_seconds,
                "ld_tables_preloaded": int(n_preloaded),
                "harmonization_seconds": build_seconds,
                "overlapped_stage1_seconds": stage1_seconds,
                "final_clumping_seconds": time.perf_counter() - finalize_started,
                "rows_before_dedup": int(rows_before_dedup),
                "harmonized_rows": int(len(combined)),
            }
        )
    except BaseException:
        done.put({"ok": False, "traceback": traceback.format_exc()})


class PipelinedClumper:
    def __init__(
        self,
        *,
        cache_dir: Path,
        output_dir: Path,
        clumping_dir: Path,
        ld_dir: Path,
        degrees_of_freedom: int,
        n_samples: int,
        n_variants: int,
        gwas_p: float,
        lead_p: float,
        variant_offset: int = 0,
    ):
        self.context = mp.get_context("fork")
        self.control = self.context.Queue()
        self.done = self.context.Queue()
        self.scratch = Path("/dev/shm") / f"torchgwas-jagwas-{os.getpid()}-{uuid.uuid4().hex[:8]}"
        self.scratch.mkdir(parents=True)
        self.counter = 0
        self.finite_rows = 0
        self.variant_offset = int(variant_offset)
        config = {
            "cache_dir": str(cache_dir),
            "output_dir": str(output_dir),
            "clumping_dir": str(clumping_dir),
            "ld_dir": str(ld_dir),
            "degrees_of_freedom": int(degrees_of_freedom),
            "n_samples": int(n_samples),
            "n_variants": int(n_variants),
            "gwas_p": float(gwas_p),
            "lead_p": float(lead_p),
        }
        self.process = self.context.Process(
            target=_clumping_worker,
            args=(self.control, self.done, config),
            daemon=False,
        )
        self.process.start()

    def consume(self, chunk):
        if not self.process.is_alive():
            try:
                report = self.done.get_nowait()
            except queue.Empty:
                report = None
            detail = "" if not report else "\n" + report.get("traceback", "")
            raise RuntimeError("clumping worker exited during the scan" + detail)
        start, end = int(chunk[0]), int(chunk[1])
        statistic = np.asarray(chunk[3], dtype=np.float64).reshape(-1)
        if len(statistic) != end - start:
            raise ValueError("JAGWAS callback received a non-dense result chunk")
        self.finite_rows += int(np.isfinite(statistic).sum())
        path = self.scratch / f"chunk_{self.counter:08d}.npy"
        self.counter += 1
        np.save(path, statistic, allow_pickle=False)
        self.control.put(
            (start + self.variant_offset, end + self.variant_offset, str(path))
        )

    def finish(self, timeout: float = 900.0) -> dict:
        self.control.put(None)
        self.control.close()
        self.control.join_thread()
        try:
            report = self.done.get(timeout=timeout)
        except queue.Empty as error:
            self.process.terminate()
            self.process.join()
            raise RuntimeError("pipelined clumping worker timed out") from error
        self.process.join()
        shutil.rmtree(self.scratch, ignore_errors=True)
        if self.process.exitcode != 0 and report.get("ok"):
            raise RuntimeError(f"clumping worker exited with code {self.process.exitcode}")
        if not report.get("ok"):
            raise RuntimeError("clumping worker failed:\n" + report["traceback"])
        report["finite_jagwas_rows"] = int(self.finite_rows)
        return report

    def abort(self):
        if self.process.is_alive():
            self.process.terminate()
        self.process.join()
        shutil.rmtree(self.scratch, ignore_errors=True)
