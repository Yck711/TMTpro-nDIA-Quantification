#!/usr/bin/env python3
"""Extract raw TMTpro 18-plex reporter intensities for all candidate DIA MS2 scans."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzml


C13_MASS_DIFFERENCE = 1.00335
DEFAULT_RUN = "TMTpro18_Ecoli_HeLa_LC60min_sample4_0dot8Th_overlap_1"
PARQUET_COLUMNS = [
    "Run",
    "Precursor.Id",
    "Precursor.Charge",
    "Precursor.Mz",
    "RT",
    "RT.Start",
    "RT.Stop",
    "Protein.Names",
]
REPORTER_COLUMNS = [
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
OUTPUT_COLUMNS = [
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
    *REPORTER_COLUMNS,
    "Reporter.TotalIntensity",
    "Reporter.DetectedChannelCount",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Match DIA-NN precursor records to all candidate MS2 scans and extract "
            "uncorrected TMTpro 18-plex reporter intensities."
        )
    )
    parser.add_argument("--parquet", type=Path, required=True)
    parser.add_argument("--mzml", type=Path, required=True)
    parser.add_argument("--reporter-mz", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run", default=DEFAULT_RUN)
    parser.add_argument("--ppm", type=float, default=5.0)
    parser.add_argument(
        "--rt-bin-width",
        type=float,
        default=0.05,
        help="RT bin width in minutes used only to accelerate candidate lookup.",
    )
    parser.add_argument("--progress-every", type=int, default=50000)
    return parser.parse_args()


def load_reporter_masses(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    observed_channels = [row["Channel"] for row in rows]
    if observed_channels != REPORTER_COLUMNS:
        raise ValueError(
            "Reporter channel order does not match the expected TMTpro 18-plex order. "
            f"Observed: {observed_channels}"
        )

    masses = np.asarray([float(row["Theoretical_mz"]) for row in rows], dtype=np.float64)
    if masses.shape != (18,) or not np.all(np.isfinite(masses)):
        raise ValueError("Reporter m/z table must contain 18 finite theoretical masses.")
    return masses


def load_precursors(path: Path, run: str) -> pd.DataFrame:
    frame = pd.read_parquet(path, columns=PARQUET_COLUMNS)
    frame = frame.loc[frame["Run"].eq(run), PARQUET_COLUMNS].copy()
    frame.insert(0, "PrecursorRecordID", np.arange(1, len(frame) + 1, dtype=np.int64))

    numeric_columns = [
        "Precursor.Charge",
        "Precursor.Mz",
        "RT",
        "RT.Start",
        "RT.Stop",
    ]
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    valid = np.ones(len(frame), dtype=bool)
    for column in numeric_columns:
        valid &= np.isfinite(frame[column].to_numpy(dtype=np.float64))
    valid &= frame["Precursor.Charge"].to_numpy(dtype=np.float64) > 0
    valid &= frame["RT.Start"].to_numpy(dtype=np.float64) <= frame["RT.Stop"].to_numpy(
        dtype=np.float64
    )
    invalid_count = int((~valid).sum())
    if invalid_count:
        print(
            f"Warning: excluding {invalid_count} target-run rows with invalid matching fields.",
            file=sys.stderr,
        )
        frame = frame.loc[valid].copy()

    if frame.empty:
        raise ValueError(f"No valid precursor records found for Run={run!r}.")
    return frame.reset_index(drop=True)


def build_rt_bins(
    rt_start: np.ndarray,
    rt_stop: np.ndarray,
    bin_width: float,
) -> tuple[list[np.ndarray], float]:
    if bin_width <= 0:
        raise ValueError("RT bin width must be positive.")

    min_rt = float(np.min(rt_start))
    max_rt = float(np.max(rt_stop))
    origin = math.floor(min_rt / bin_width) * bin_width
    bin_count = int(math.floor((max_rt - origin) / bin_width)) + 2
    members: list[list[int]] = [[] for _ in range(bin_count)]

    for index, (start, stop) in enumerate(zip(rt_start, rt_stop, strict=True)):
        first_bin = max(0, int(math.floor((float(start) - origin) / bin_width)))
        last_bin = min(
            bin_count - 1,
            int(math.floor((float(stop) - origin) / bin_width)),
        )
        for bin_index in range(first_bin, last_bin + 1):
            members[bin_index].append(index)

    return [np.asarray(indices, dtype=np.int64) for indices in members], origin


def scan_number(spectrum_id: str) -> str:
    match = re.search(r"scan=(\d+)", spectrum_id)
    if match:
        return match.group(1)
    numbers = re.findall(r"\d+", spectrum_id)
    return numbers[-1] if numbers else spectrum_id


def isolation_window(spectrum: dict) -> tuple[float, float, float] | None:
    try:
        precursor = spectrum["precursorList"]["precursor"][0]
        window = precursor["isolationWindow"]
        target = float(window["isolation window target m/z"])
        if "isolation window lower offset" in window:
            lower = target - float(window["isolation window lower offset"])
            upper = target + float(window["isolation window upper offset"])
        else:
            lower = float(window["isolation window lower limit"])
            upper = float(window["isolation window upper limit"])
    except (KeyError, IndexError, TypeError, ValueError):
        return None

    if not all(math.isfinite(value) for value in (target, lower, upper)):
        return None
    return target, lower, upper


def extract_reporters(
    spectrum: dict,
    reporter_mz: np.ndarray,
    ppm: float,
) -> np.ndarray:
    mz_values = np.asarray(spectrum.get("m/z array", []), dtype=np.float64)
    intensities = np.asarray(spectrum.get("intensity array", []), dtype=np.float64)
    result = np.zeros(len(reporter_mz), dtype=np.float64)
    if mz_values.size == 0:
        return result

    errors = reporter_mz * ppm * 1e-6
    for index, (target, error) in enumerate(zip(reporter_mz, errors, strict=True)):
        left = int(np.searchsorted(mz_values, target - error, side="left"))
        right = int(np.searchsorted(mz_values, target + error, side="right"))
        if left < right:
            result[index] = float(np.sum(intensities[left:right], dtype=np.float64))
    return result


def as_csv_value(value: object) -> object:
    if pd.isna(value):
        return ""
    return value


def main() -> None:
    args = parse_args()
    started = time.time()
    reporter_mz = load_reporter_masses(args.reporter_mz)
    precursors = load_precursors(args.parquet, args.run)

    precursor_mz = precursors["Precursor.Mz"].to_numpy(dtype=np.float64)
    charge = precursors["Precursor.Charge"].to_numpy(dtype=np.float64)
    precursor_m1 = precursor_mz + C13_MASS_DIFFERENCE / charge
    rt = precursors["RT"].to_numpy(dtype=np.float64)
    rt_start = precursors["RT.Start"].to_numpy(dtype=np.float64)
    rt_stop = precursors["RT.Stop"].to_numpy(dtype=np.float64)
    rt_bins, rt_origin = build_rt_bins(rt_start, rt_stop, args.rt_bin_width)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    spectra_seen = 0
    ms2_seen = 0
    ms2_missing_window = 0
    matched_unique_ms2 = 0
    matched_precursor = np.zeros(len(precursors), dtype=bool)

    print(
        f"Loaded {len(precursors):,} precursor records for {args.run}; "
        f"extracting reporters at +/-{args.ppm:g} ppm.",
        file=sys.stderr,
    )

    with args.output.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.writer(output_handle)
        writer.writerow(OUTPUT_COLUMNS)

        with mzml.read(str(args.mzml), use_index=False, iterative=True) as reader:
            for spectrum in reader:
                spectra_seen += 1
                if spectrum.get("ms level") != 2:
                    continue
                ms2_seen += 1

                scan = spectrum.get("scanList", {}).get("scan", [{}])[0]
                scan_rt = float(scan.get("scan start time", math.nan))
                if not math.isfinite(scan_rt):
                    continue

                window = isolation_window(spectrum)
                if window is None:
                    ms2_missing_window += 1
                    continue
                isolation_target, isolation_lower, isolation_upper = window

                rt_bin_index = int(math.floor((scan_rt - rt_origin) / args.rt_bin_width))
                if rt_bin_index < 0 or rt_bin_index >= len(rt_bins):
                    continue
                possible = rt_bins[rt_bin_index]
                if possible.size == 0:
                    continue

                exact_rt = (rt_start[possible] <= scan_rt) & (scan_rt <= rt_stop[possible])
                possible = possible[exact_rt]
                if possible.size == 0:
                    continue

                m0_match = (precursor_mz[possible] >= isolation_lower) & (
                    precursor_mz[possible] <= isolation_upper
                )
                m1_match = (precursor_m1[possible] >= isolation_lower) & (
                    precursor_m1[possible] <= isolation_upper
                )
                selected_mask = m0_match | m1_match
                selected = possible[selected_mask]
                if selected.size == 0:
                    continue

                selected_m0 = m0_match[selected_mask]
                selected_m1 = m1_match[selected_mask]
                reporter_values = extract_reporters(spectrum, reporter_mz, args.ppm)
                reporter_total = float(np.sum(reporter_values, dtype=np.float64))
                detected_channel_count = int(np.count_nonzero(reporter_values > 0))
                current_scan_id = scan_number(str(spectrum.get("id", "")))
                matched_unique_ms2 += 1

                for precursor_index, is_m0, is_m1 in zip(
                    selected,
                    selected_m0,
                    selected_m1,
                    strict=True,
                ):
                    if is_m0 and is_m1:
                        match_type = "M0+M1"
                    elif is_m0:
                        match_type = "M0"
                    else:
                        match_type = "M1"

                    precursor_row = precursors.iloc[int(precursor_index)]
                    output_row = [
                        int(precursor_row["PrecursorRecordID"]),
                        as_csv_value(precursor_row["Run"]),
                        as_csv_value(precursor_row["Precursor.Id"]),
                        as_csv_value(precursor_row["Protein.Names"]),
                        float(precursor_row["Precursor.Mz"]),
                        int(precursor_row["Precursor.Charge"]),
                        float(precursor_row["RT"]),
                        float(precursor_row["RT.Start"]),
                        float(precursor_row["RT.Stop"]),
                        current_scan_id,
                        scan_rt,
                        isolation_target,
                        isolation_lower,
                        isolation_upper,
                        match_type,
                        scan_rt - rt[int(precursor_index)],
                        *reporter_values.tolist(),
                        reporter_total,
                        detected_channel_count,
                    ]
                    writer.writerow(output_row)
                    rows_written += 1
                    matched_precursor[int(precursor_index)] = True

                if args.progress_every > 0 and ms2_seen % args.progress_every == 0:
                    elapsed = time.time() - started
                    print(
                        f"MS2={ms2_seen:,}; matched rows={rows_written:,}; "
                        f"elapsed={elapsed / 60:.1f} min",
                        file=sys.stderr,
                    )

    elapsed = time.time() - started
    print(
        "Completed: "
        f"spectra={spectra_seen:,}; MS2={ms2_seen:,}; "
        f"MS2 missing isolation window={ms2_missing_window:,}; "
        f"matched unique MS2={matched_unique_ms2:,}; "
        f"precursor-MS2 rows={rows_written:,}; "
        f"precursors with candidate MS2={int(matched_precursor.sum()):,}; "
        f"precursors without candidate MS2={int((~matched_precursor).sum()):,}; "
        f"elapsed={elapsed / 60:.2f} min; output={args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
