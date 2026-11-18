#!/usr/bin/env python3
"""Quality-first one-to-one precursor/MS2 assignment with reassignment."""

from __future__ import annotations

import argparse
import heapq
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "mZML" / "Python"))
import select_max_pif_unique_ms2_and_plot_ratios as ratio_plot  # noqa: E402
from make_strict_one_to_one_quant_and_plot import (  # noqa: E402
    ANALYSIS_CORRECTED_CHANNELS,
    ANALYSIS_CHANNELS,
    ALL_CORRECTED_CHANNELS,
    CORRECTED_CHANNELS,
    add_analysis_fields,
    make_ratio_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--candidate-pif", type=Path, required=True)
    parser.add_argument("--corrected-reporters", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def priority_tuple(row: pd.Series | dict[str, object], precursor_tie: bool = False) -> tuple[float, float, float, float, int]:
    pif = float(row["Precursor.PIF"]) if np.isfinite(row["Precursor.PIF"]) else -np.inf
    target = float(row["PIF.TargetIsotopeIntensity"]) if np.isfinite(row["PIF.TargetIsotopeIntensity"]) else -np.inf
    matched = float(row["PIF.MatchedIsotopeCount"]) if np.isfinite(row["PIF.MatchedIsotopeCount"]) else -np.inf
    rt = float(row["RT.Difference"]) if np.isfinite(row["RT.Difference"]) else np.inf
    scan = int(row["CandidateMS2.ScanID"]) if np.isfinite(row["CandidateMS2.ScanID"]) else 2**63 - 1
    return (pif, target, matched, -abs(rt), -scan)


def resolve_assignments(candidates: pd.DataFrame, precursor_ids: set[int]) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = candidates.loc[candidates["PrecursorRecordID"].isin(precursor_ids)].copy()
    frame = frame.loc[frame["PIF.Status"].eq("ok")].copy()
    frame["Selection.AbsRTDifference"] = frame["RT.Difference"].abs()
    frame = frame.sort_values(
        [
            "PrecursorRecordID",
            "Precursor.PIF",
            "PIF.TargetIsotopeIntensity",
            "PIF.MatchedIsotopeCount",
            "Selection.AbsRTDifference",
            "CandidateMS2.ScanID",
        ],
        ascending=[True, False, False, False, True, True],
        kind="mergesort",
    )
    candidate_lists: dict[int, list[dict[str, object]]] = {}
    for precursor_id, group in frame.groupby("PrecursorRecordID", sort=False):
        records = group.to_dict("records")
        for rank, record in enumerate(records, start=1):
            record["StrictMS2.ReassignmentRank"] = rank
        candidate_lists[int(precursor_id)] = records

    # Deferred acceptance: each precursor proposes in its own quality order;
    # each MS2 keeps the highest-priority competing precursor.
    next_index = {precursor_id: 0 for precursor_id in candidate_lists}
    assigned_scan: dict[int, dict[str, object]] = {}
    assigned_precursor: dict[int, int] = {}
    queue = list(candidate_lists)
    heapq.heapify(queue)
    proposal_count = 0
    while queue:
        precursor_id = heapq.heappop(queue)
        records = candidate_lists[precursor_id]
        index = next_index[precursor_id]
        if index >= len(records):
            continue
        candidate = records[index]
        next_index[precursor_id] = index + 1
        proposal_count += 1
        scan_id = int(candidate["CandidateMS2.ScanID"])
        incumbent = assigned_scan.get(scan_id)
        if incumbent is None:
            assigned_scan[scan_id] = candidate
            assigned_precursor[scan_id] = precursor_id
            continue
        incumbent_id = assigned_precursor[scan_id]
        new_priority = priority_tuple(candidate)
        old_priority = priority_tuple(incumbent)
        if new_priority > old_priority:
            assigned_scan[scan_id] = candidate
            assigned_precursor[scan_id] = precursor_id
            if next_index[incumbent_id] < len(candidate_lists[incumbent_id]):
                heapq.heappush(queue, incumbent_id)
        elif next_index[precursor_id] < len(records):
            heapq.heappush(queue, precursor_id)

    selected = pd.DataFrame.from_records(list(assigned_scan.values()))
    # The rank is attached while building candidate_lists, avoiding an
    # expensive per-row scan through every candidate list after matching.
    selected["StrictMS2.WasReassigned"] = selected["StrictMS2.ReassignmentRank"].gt(1)
    scan_competitor_count = frame.groupby("CandidateMS2.ScanID")["PrecursorRecordID"].nunique()
    selected["StrictMS2.CompetitionCount"] = selected["CandidateMS2.ScanID"].map(scan_competitor_count).fillna(1).astype(np.int32)
    selected["StrictMS2.WasSharedInCandidateGraph"] = selected["StrictMS2.CompetitionCount"].gt(1)
    selected["StrictMS2.Assignment"] = np.where(
        selected["StrictMS2.WasSharedInCandidateGraph"],
        "One-to-one winner after quality-first reassignment",
        "Unique MS2 assignment",
    )
    selected["Selection.Method"] = "Quality-first deferred acceptance; one precursor per MS2"
    selected["Selection.Rule"] = (
        "Precursor proposals by PIF, target isotope intensity, matched isotope count, "
        "min abs RT difference, min scan ID; each MS2 retains highest pair priority"
    )
    stats = {
        "CandidateGraph.PairRows": len(frame),
        "CandidateGraph.UniquePrecursors": int(frame["PrecursorRecordID"].nunique()),
        "CandidateGraph.UniqueMS2": int(frame["CandidateMS2.ScanID"].nunique()),
        "CandidateGraph.SharedMS2": int((scan_competitor_count > 1).sum()),
        "Assignment.ProposalsEvaluated": proposal_count,
        "Assignment.StrictOneToOneRows": len(selected),
        "Assignment.UniquePrecursors": int(selected["PrecursorRecordID"].nunique()),
        "Assignment.UniqueMS2": int(selected["CandidateMS2.ScanID"].nunique()),
        "Assignment.UnassignedPrecursors": int(frame["PrecursorRecordID"].nunique() - selected["PrecursorRecordID"].nunique()),
        "Assignment.ReassignedPrecursors": int(selected["StrictMS2.WasReassigned"].sum()),
    }
    return selected, stats


def main() -> None:
    args = parse_args()
    master = pd.read_csv(args.master, usecols=["PrecursorRecordID"])
    precursor_ids = set(pd.to_numeric(master["PrecursorRecordID"], errors="coerce").dropna().astype(int))

    candidate_columns = [
        "PrecursorRecordID",
        "Run",
        "Precursor.Id",
        "Protein.Names",
        "Precursor.Mz",
        "Precursor.Charge",
        "RT",
        "RT.Start",
        "RT.Stop",
        "CandidateMS2.ScanID",
        "CandidateMS2.RT",
        "Isolation.TargetMz",
        "Isolation.LowerMz",
        "Isolation.UpperMz",
        "Match.Type",
        "RT.Difference",
        "MS1.Before.ScanID",
        "MS1.Before.RT",
        "MS1.After.ScanID",
        "MS1.After.RT",
        "MS1.InterpolationWeight",
        "PIF.M0.Intensity",
        "PIF.M1.Intensity",
        "PIF.M2.Intensity",
        "PIF.M3.Intensity",
        "PIF.M4.Intensity",
        "PIF.TargetIntensity.BeforeMS1",
        "PIF.TargetIntensity.AfterMS1",
        "PIF.WindowIntensity.BeforeMS1",
        "PIF.WindowIntensity.AfterMS1",
        "PIF.TargetIsotopeIntensity",
        "PIF.TotalIsolationWindowIntensity",
        "PIF.InterferenceIntensity",
        "Precursor.PIF",
        "PIF.MatchedIsotopeCount",
        "PIF.MatchedIsotopes",
        "PIF.Status",
        "PIF.MassTolerancePPM",
        "PIF.IsotopeDefinition",
        "PIF.CalculationMethod",
        "CandidateMS2.PrecursorMultiplicity",
    ]
    candidates = pd.read_csv(args.candidate_pif, usecols=candidate_columns)
    numeric_columns = [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        "RT.Difference",
        "Precursor.PIF",
        "PIF.TargetIsotopeIntensity",
        "PIF.MatchedIsotopeCount",
    ]
    for column in numeric_columns:
        candidates[column] = pd.to_numeric(candidates[column], errors="coerce")
    selected, stats = resolve_assignments(candidates, precursor_ids)

    corrected_columns = [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        *ALL_CORRECTED_CHANNELS,
        "Corrected.Reporter.TotalIntensity",
        "Corrected.Reporter.DetectedChannelCount",
        "Corrected.NegativeBeforeClippingCount",
        "IsotopeCorrection.Method",
    ]
    corrected = pd.read_csv(args.corrected_reporters, usecols=corrected_columns)
    for column in ["PrecursorRecordID", "CandidateMS2.ScanID", *ALL_CORRECTED_CHANNELS]:
        corrected[column] = pd.to_numeric(corrected[column], errors="coerce")
    selected = selected.merge(
        corrected,
        on=["PrecursorRecordID", "CandidateMS2.ScanID"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    stats["Assignment.MissingCorrectedReporterMatch"] = int(selected["_merge"].ne("both").sum())
    selected = selected.drop(columns="_merge")
    selected = add_analysis_fields(selected)
    selected["Reporter.Intensity.Type"] = "Official full-18plex isotope corrected"
    selected["PIF.Selection.Context"] = "Uncalibrated theoretical isotope centers; 10 ppm MS1 matching"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = "strictOneToOneReassigned_uncalibrated10ppm_officialIsotopeCorrected"
    selected.to_csv(args.output_dir / f"precursor_quantification_allSpecies_{stem}_18channel.csv", index=False)
    minimal = [
        "PrecursorRecordID",
        "Precursor.Id",
        "Protein.Names",
        "CandidateMS2.ScanID",
        "Precursor.PIF",
        "PIF.TargetIsotopeIntensity",
        "PIF.MatchedIsotopeCount",
        "RT.Difference",
        "StrictMS2.WasSharedInCandidateGraph",
        "StrictMS2.CompetitionCount",
        "StrictMS2.ReassignmentRank",
        "StrictMS2.WasReassigned",
        "StrictMS2.Assignment",
        *ANALYSIS_CORRECTED_CHANNELS,
    ]
    selected[minimal].to_csv(args.output_dir / f"precursor_quantification_allSpecies_{stem}_12channel_minimal.csv", index=False)

    ecoli = selected.loc[selected["Protein.Names"].fillna("").str.contains("ECOLI", regex=False)].copy()
    ratios = make_ratio_table(ecoli)
    ratio_long = args.output_dir / f"ecoli_reporter_ratios_{stem}_long.csv"
    ratios.to_csv(ratio_long, index=False)
    ratio_plot.make_summary(ratios).to_csv(args.output_dir / f"ecoli_reporter_ratios_{stem}_summary.csv", index=False)
    from plot_unique_ms2_ratio_transparent_smooth import main as plot_main

    old_argv = sys.argv
    sys.argv = [
        "plot_unique_ms2_ratio_transparent_smooth.py",
        "--input-long", str(ratio_long),
        "--output-pdf", str(args.output_dir / f"ecoli_reporter_ratio_{stem}_transparent80_smooth_8x5.pdf"),
        "--title", "E. coli reporter ion ratios - strict one-to-one MS2 (uncalibrated 10 ppm)",
        "--method-note", "quality-first reassignment; official isotope corrected",
    ]
    plot_main()
    sys.argv = old_argv

    stats["Output.StrictECOLI.Rows"] = int(ecoli.shape[0])
    stats["Output.StrictECOLI.RatioPrecursors"] = int(ratios["PrecursorRecordID"].nunique())
    stats["Output.StrictECOLI.All12Positive"] = int(
        (ecoli[ANALYSIS_CORRECTED_CHANNELS].gt(0).all(axis=1)).sum()
    )
    pd.DataFrame([{"Metric": key, "Value": value} for key, value in stats.items()]).to_csv(
        args.output_dir / f"precursor_quantification_{stem}_QC.csv", index=False
    )
    print(pd.DataFrame([stats]).T.to_string(header=False))


if __name__ == "__main__":
    main()
