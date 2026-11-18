#!/usr/bin/env python3
"""Resolve shared MS2 scans and make a strict one-to-one quantification table."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "mZML" / "Python"))
import select_max_pif_unique_ms2_and_plot_ratios as ratio_plot  # noqa: E402


CORRECTED_CHANNELS = [f"Corrected_{channel}" for channel in ratio_plot.CHANNELS]
ANALYSIS_CHANNELS = ratio_plot.CHANNELS
ALL_CORRECTED_CHANNELS = [
    "Corrected_TMT_126",
    "Corrected_TMT_127N",
    "Corrected_TMT_127C",
    "Corrected_TMT_128N",
    "Corrected_TMT_128C",
    "Corrected_TMT_129N",
    "Corrected_TMT_129C",
    "Corrected_TMT_130N",
    "Corrected_TMT_130C",
    "Corrected_TMT_131N",
    "Corrected_TMT_131C",
    "Corrected_TMT_132N",
    "Corrected_TMT_132C",
    "Corrected_TMT_133N",
    "Corrected_TMT_133C",
    "Corrected_TMT_134N",
    "Corrected_TMT_134C",
    "Corrected_TMT_135N",
]
ANALYSIS_CORRECTED_CHANNELS = [f"Corrected_{channel}" for channel in ANALYSIS_CHANNELS]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def select_one_precursor_per_ms2(master: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    work = master.copy()
    numeric_columns = [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        "Precursor.PIF",
        "PIF.TargetIsotopeIntensity",
        "PIF.MatchedIsotopeCount",
        "RT.Difference",
    ]
    for column in numeric_columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work["Selection.AbsRTDifference"] = work["RT.Difference"].abs()
    scan_counts = work.groupby("CandidateMS2.ScanID")["PrecursorRecordID"].transform("size")
    work["StrictMS2.CompetitionCount"] = scan_counts.astype(np.int32)
    work["StrictMS2.WasShared"] = scan_counts.gt(1)
    work = work.sort_values(
        [
            "CandidateMS2.ScanID",
            "Precursor.PIF",
            "PIF.TargetIsotopeIntensity",
            "PIF.MatchedIsotopeCount",
            "Selection.AbsRTDifference",
            "PrecursorRecordID",
        ],
        ascending=[True, False, False, False, True, True],
        kind="mergesort",
    )
    winners = work.drop_duplicates("CandidateMS2.ScanID", keep="first").copy()
    winners["StrictMS2.Assignment"] = np.where(
        winners["StrictMS2.WasShared"],
        "Winner among shared MS2 by PIF priority",
        "Unique MS2 assignment",
    )
    winners = winners.sort_values("PrecursorRecordID", kind="mergesort").reset_index(drop=True)
    stats = {
        "Input.Rows": len(master),
        "Input.UniquePrecursors": int(master["PrecursorRecordID"].nunique()),
        "Input.UniqueSelectedMS2": int(master["CandidateMS2.ScanID"].nunique()),
        "Input.RowsUsingSharedMS2": int(master["SelectedMS2.IsShared"].sum()),
        "Input.SharedMS2Count": int(master.loc[master["SelectedMS2.IsShared"], "CandidateMS2.ScanID"].nunique()),
        "Output.StrictOneToOne.Rows": len(winners),
        "Output.StrictOneToOne.UniquePrecursors": int(winners["PrecursorRecordID"].nunique()),
        "Output.StrictOneToOne.UniqueMS2": int(winners["CandidateMS2.ScanID"].nunique()),
        "Removed.SharedMS2.LosingPrecursors": len(master) - len(winners),
        "Output.StrictECOLI.Rows": int(winners["Protein.Names"].fillna("").str.contains("ECOLI", regex=False).sum()),
    }
    return winners, stats


def add_analysis_fields(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame[ANALYSIS_CORRECTED_CHANNELS].to_numpy(dtype=np.float64)
    finite = np.isfinite(values)
    frame["Corrected.All12Finite"] = finite.all(axis=1)
    frame["Corrected.All12Positive"] = (finite & (values > 0)).all(axis=1)
    frame["Corrected.Min12ChannelIntensity"] = np.where(
        finite.all(axis=1), np.min(values, axis=1), np.nan
    )
    frame["Corrected.Median12ChannelIntensity"] = np.where(
        finite.all(axis=1), np.median(values, axis=1), np.nan
    )
    names = np.asarray(ANALYSIS_CHANNELS, dtype=object)
    indices = np.argmin(np.where(finite, values, np.inf), axis=1)
    frame["Corrected.Min12ChannelName"] = np.where(
        finite.all(axis=1), names[indices], ""
    )
    frame["Corrected.Min12ChannelPercentile"] = (
        frame["Corrected.Min12ChannelIntensity"].rank(method="average", pct=True) * 100.0
    )
    return frame


def make_ratio_table(ecoli: pd.DataFrame) -> pd.DataFrame:
    frame = ecoli.copy()
    for corrected, standard in zip(ANALYSIS_CORRECTED_CHANNELS, ANALYSIS_CHANNELS, strict=True):
        frame[standard] = frame[corrected]
    frame["Reporter.TotalIntensity"] = frame["Corrected.Reporter.TotalIntensity"]
    frame["Reporter.DetectedChannelCount"] = frame["Corrected.Reporter.DetectedChannelCount"]
    frame = add_analysis_fields(frame)
    frame["All12_PositiveFinite"] = frame["Corrected.All12Positive"]
    frame["Min_12Channel_Intensity"] = frame["Corrected.Min12ChannelIntensity"]
    frame["Min_12Channel_Name"] = frame["Corrected.Min12ChannelName"]
    frame["Median_12Channel_Intensity"] = frame["Corrected.Median12ChannelIntensity"]
    frame["Min12_Intensity_Percentile_WithinFinalECOLI"] = frame["Corrected.Min12ChannelPercentile"]
    frame = frame.loc[frame["Corrected.All12Positive"]].copy()
    frame["Method"] = "Strict one-to-one MS2; uncalibrated 10 ppm PIF; official isotope corrected"
    ratios = ratio_plot.calculate_ratios(frame, frame["Method"].iloc[0])
    return ratios


def main() -> None:
    args = parse_args()
    master = pd.read_csv(args.master)
    required = [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        "Precursor.PIF",
        "PIF.TargetIsotopeIntensity",
        "PIF.MatchedIsotopeCount",
        "RT.Difference",
        "Protein.Names",
        "SelectedMS2.IsShared",
        *CORRECTED_CHANNELS,
        "Corrected.Reporter.TotalIntensity",
        "Corrected.Reporter.DetectedChannelCount",
    ]
    missing = [column for column in required if column not in master.columns]
    if missing:
        raise ValueError(f"Master table is missing columns: {missing}")
    winners, stats = select_one_precursor_per_ms2(master)
    winners = add_analysis_fields(winners)
    winners["StrictMS2.PIFRankWithinScan"] = 1
    winners["StrictMS2.PIFCompetitionResolved"] = True

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = "strictOneToOne_uncalibrated10ppm_officialIsotopeCorrected"
    full_path = args.output_dir / f"precursor_quantification_allSpecies_{stem}_18channel.csv"
    winners.to_csv(full_path, index=False)

    minimal_columns = [
        "PrecursorRecordID",
        "Precursor.Id",
        "Protein.Names",
        "CandidateMS2.ScanID",
        "Precursor.PIF",
        "PIF.TargetIsotopeIntensity",
        "PIF.MatchedIsotopeCount",
        "RT.Difference",
        "StrictMS2.WasShared",
        "StrictMS2.CompetitionCount",
        "StrictMS2.Assignment",
        *ANALYSIS_CORRECTED_CHANNELS,
    ]
    winners[minimal_columns].to_csv(
        args.output_dir / f"precursor_quantification_allSpecies_{stem}_12channel_minimal.csv",
        index=False,
    )

    strict_ecoli = winners.loc[
        winners["Protein.Names"].fillna("").str.contains("ECOLI", regex=False)
    ].copy()
    stats["Output.StrictECOLI.All12CorrectedPositive"] = int(
        (strict_ecoli[CORRECTED_CHANNELS[:12]].gt(0).all(axis=1)).sum()
    )
    ratios = make_ratio_table(strict_ecoli)
    ratio_long = args.output_dir / f"ecoli_reporter_ratios_{stem}_long.csv"
    ratios.to_csv(ratio_long, index=False)
    summary = ratio_plot.make_summary(ratios)
    summary.to_csv(args.output_dir / f"ecoli_reporter_ratios_{stem}_summary.csv", index=False)

    # Match the publication-style transparent/smoothed plot used for the prior comparison.
    from plot_unique_ms2_ratio_transparent_smooth import main as plot_main

    old_argv = sys.argv
    sys.argv = [
        "plot_unique_ms2_ratio_transparent_smooth.py",
        "--input-long",
        str(ratio_long),
        "--output-pdf",
        str(args.output_dir / f"ecoli_reporter_ratio_{stem}_transparent80_smooth_8x5.pdf"),
        "--title",
        "E. coli reporter ion ratios - strict one-to-one MS2 (uncalibrated 10 ppm)",
        "--method-note",
        "one precursor per MS2; official isotope corrected",
    ]
    plot_main()
    sys.argv = old_argv

    stats["Output.StrictECOLI.RatioPrecursors"] = int(ratios["PrecursorRecordID"].nunique())
    pd.DataFrame([{"Metric": key, "Value": value} for key, value in stats.items()]).to_csv(
        args.output_dir / f"precursor_quantification_{stem}_QC.csv", index=False
    )
    print(pd.DataFrame([stats]).T.to_string(header=False))


if __name__ == "__main__":
    main()
