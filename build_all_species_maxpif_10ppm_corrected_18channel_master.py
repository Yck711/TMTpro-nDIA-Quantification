#!/usr/bin/env python3
"""Build an all-species, max-PIF, isotope-corrected 18-channel master table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


RAW_REPORTER_CHANNELS = [
    "TMT_126",
    "TMT_127N",
    "TMT_127C",
    "TMT_128N",
    "TMT_128C",
    "TMT_129N",
    "TMT_129C",
    "TMT_130N",
    "TMT_130C",
    "TMT_131N",
    "TMT_131C",
    "TMT_132N",
    "TMT_132C",
    "TMT_133N",
    "TMT_133C",
    "TMT_134N",
    "TMT_134C",
    "TMT_135N",
]
CORRECTED_REPORTER_CHANNELS = [
    f"Corrected_{channel}" for channel in RAW_REPORTER_CHANNELS
]
PAIR_KEYS = ["PrecursorRecordID", "CandidateMS2.ScanID"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pif", type=Path, required=True)
    parser.add_argument("--corrected-reporters", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qc-output", type=Path, required=True)
    return parser.parse_args()


def numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def select_maximum_pif(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    original_precursor_count = int(frame["PrecursorRecordID"].nunique())
    all_candidate_count = frame.groupby("PrecursorRecordID")[
        "CandidateMS2.ScanID"
    ].transform("size")
    frame["Selection.AllCandidateCount"] = all_candidate_count.astype(np.int32)

    ok = frame.loc[frame["PIF.Status"].eq("ok")].copy()
    ok["Selection.AbsRTDifference"] = ok["RT.Difference"].abs()
    ok["Selection.OKCandidateCount"] = (
        ok.groupby("PrecursorRecordID")["CandidateMS2.ScanID"]
        .transform("size")
        .astype(np.int32)
    )
    maximum_pif = ok.groupby("PrecursorRecordID")["Precursor.PIF"].transform("max")
    maximum_tie = np.isclose(
        ok["Precursor.PIF"].to_numpy(dtype=np.float64),
        maximum_pif.to_numpy(dtype=np.float64),
        rtol=0,
        atol=1e-12,
    )
    ok["Selection.MaxPIFTieCount"] = (
        pd.Series(maximum_tie, index=ok.index)
        .groupby(ok["PrecursorRecordID"])
        .transform("sum")
        .astype(np.int32)
    )
    ok = ok.sort_values(
        [
            "PrecursorRecordID",
            "Precursor.PIF",
            "PIF.MatchedIsotopeCount",
            "PIF.TargetIsotopeIntensity",
            "Selection.AbsRTDifference",
            "CandidateMS2.ScanID",
        ],
        ascending=[True, False, False, False, True, True],
        kind="mergesort",
    )
    selected = ok.drop_duplicates("PrecursorRecordID", keep="first").copy()
    selected["Selection.Method"] = "Maximum PIF per precursor; shared MS2 allowed"
    selected["Selection.Rule"] = (
        "PIF.Status=ok; max Precursor.PIF; max PIF.MatchedIsotopeCount; "
        "max PIF.TargetIsotopeIntensity; min abs RT.Difference; min scan ID"
    )
    stats = {
        "Input.PrecursorMS2.PairRows": len(frame),
        "Input.UniquePrecursors": original_precursor_count,
        "PIF.Status.OK.PairRows": len(ok),
        "Precursors.WithOKCandidate": int(ok["PrecursorRecordID"].nunique()),
        "Precursors.WithNoOKCandidate": original_precursor_count
        - int(ok["PrecursorRecordID"].nunique()),
        "Precursors.WithMultipleOKCandidates": int(
            (ok.groupby("PrecursorRecordID").size() > 1).sum()
        ),
        "Precursors.WithMaxPIFTie": int(
            selected["Selection.MaxPIFTieCount"].gt(1).sum()
        ),
        "Selected.OneRowPerPrecursor": len(selected),
    }
    return selected, stats


def main() -> None:
    args = parse_args()
    pif = pd.read_csv(args.pif)
    required_pif = [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        "Precursor.PIF",
        "PIF.Status",
        "PIF.MatchedIsotopeCount",
        "PIF.TargetIsotopeIntensity",
        "RT.Difference",
    ]
    missing = [column for column in required_pif if column not in pif.columns]
    if missing:
        raise ValueError(f"Missing PIF columns: {missing}")
    numeric(
        pif,
        [
            "PrecursorRecordID",
            "CandidateMS2.ScanID",
            "Precursor.PIF",
            "PIF.MatchedIsotopeCount",
            "PIF.TargetIsotopeIntensity",
            "RT.Difference",
        ],
    )
    selected, stats = select_maximum_pif(pif)

    corrected_columns = [
        *PAIR_KEYS,
        *CORRECTED_REPORTER_CHANNELS,
        "Corrected.Reporter.TotalIntensity",
        "Corrected.Reporter.DetectedChannelCount",
        "Corrected.NegativeBeforeClippingCount",
        "IsotopeCorrection.Method",
    ]
    corrected = pd.read_csv(args.corrected_reporters, usecols=corrected_columns)
    numeric(
        corrected,
        [
            *PAIR_KEYS,
            *CORRECTED_REPORTER_CHANNELS,
            "Corrected.Reporter.TotalIntensity",
            "Corrected.Reporter.DetectedChannelCount",
            "Corrected.NegativeBeforeClippingCount",
        ],
    )

    raw_columns_to_remove = [
        *RAW_REPORTER_CHANNELS,
        "Reporter.TotalIntensity",
        "Reporter.DetectedChannelCount",
    ]
    selected = selected.drop(columns=raw_columns_to_remove, errors="ignore")
    master = selected.merge(
        corrected,
        on=PAIR_KEYS,
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    stats["Selected.MissingCorrectedReporterMatch"] = int(
        master["_merge"].ne("both").sum()
    )
    master = master.drop(columns="_merge")

    corrected_values = master[CORRECTED_REPORTER_CHANNELS].to_numpy(
        dtype=np.float64
    )
    finite = np.isfinite(corrected_values)
    positive = corrected_values > 0
    master["Corrected.All18Finite"] = finite.all(axis=1)
    master["Corrected.All18Positive"] = (finite & positive).all(axis=1)
    master["Corrected.NonFiniteChannelCount"] = np.count_nonzero(~finite, axis=1)
    master["Corrected.ZeroOrNegativeChannelCount"] = np.count_nonzero(
        finite & ~positive, axis=1
    )
    all_finite = finite.all(axis=1)
    master["Corrected.Min18ChannelIntensity"] = np.where(
        all_finite, np.min(corrected_values, axis=1), np.nan
    )
    master["Corrected.Median18ChannelIntensity"] = np.where(
        all_finite, np.median(corrected_values, axis=1), np.nan
    )
    minimum_index = np.argmin(
        np.where(finite, corrected_values, np.inf), axis=1
    )
    master["Corrected.Min18ChannelName"] = np.where(
        all_finite,
        np.asarray(CORRECTED_REPORTER_CHANNELS, dtype=object)[minimum_index],
        "",
    )

    selected_scan_count = master.groupby("CandidateMS2.ScanID")[
        "PrecursorRecordID"
    ].transform("size")
    master["SelectedMS2.PrecursorCount"] = selected_scan_count.astype(np.int32)
    master["SelectedMS2.IsShared"] = selected_scan_count.gt(1)
    master["Reporter.Intensity.Type"] = (
        "Official full-18plex isotope corrected; negative inverse solutions "
        "clipped to zero"
    )
    master["PIF.Selection.Context"] = (
        "Uncalibrated theoretical isotope centers; 10 ppm MS1 matching"
    )
    master = master.sort_values("PrecursorRecordID", kind="mergesort").reset_index(
        drop=True
    )

    stats.update(
        {
            "Output.Rows": len(master),
            "Output.UniquePrecursors": int(master["PrecursorRecordID"].nunique()),
            "Output.UniqueSelectedMS2": int(master["CandidateMS2.ScanID"].nunique()),
            "Output.Rows.WithSharedSelectedMS2": int(
                master["SelectedMS2.IsShared"].sum()
            ),
            "Output.SelectedMS2.SharedByMultiplePrecursors": int(
                master.loc[
                    master["SelectedMS2.IsShared"], "CandidateMS2.ScanID"
                ].nunique()
            ),
            "Output.All18CorrectedFinite": int(master["Corrected.All18Finite"].sum()),
            "Output.All18CorrectedPositive": int(
                master["Corrected.All18Positive"].sum()
            ),
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(args.output, index=False)
    pd.DataFrame(
        [{"Metric": metric, "Value": value} for metric, value in stats.items()]
    ).to_csv(args.qc_output, index=False)
    print(f"Wrote: {args.output}")
    print(f"Wrote: {args.qc_output}")
    print(pd.DataFrame([stats]).T.to_string(header=False))


if __name__ == "__main__":
    main()
