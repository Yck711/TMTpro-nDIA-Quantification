#!/usr/bin/env python3
"""Select one maximum-PIF MS2 per precursor and plot E. coli reporter ratios."""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
import TMTpro_MAplot_analysis as base  # noqa: E402


CHANNELS = [
    "TMT_127N",
    "TMT_127C",
    "TMT_128N",
    "TMT_128C",
    "TMT_129N",
    "TMT_129C",
    "TMT_130N",
    "TMT_131N",
    "TMT_132N",
    "TMT_133N",
    "TMT_134N",
    "TMT_135N",
]
COMPARISONS = [
    ("0.6 Da", "Mix", "2:8", "TMT_127N", "TMT_128N", 0.25),
    ("0.6 Da", "Mix", "4:8", "TMT_127C", "TMT_128N", 0.50),
    ("0.6 Da", "Pure", "2:8", "TMT_128C", "TMT_129C", 0.25),
    ("0.6 Da", "Pure", "4:8", "TMT_129N", "TMT_129C", 0.50),
    ("1 Da", "Mix", "2:8", "TMT_130N", "TMT_132N", 0.25),
    ("1 Da", "Mix", "4:8", "TMT_131N", "TMT_132N", 0.50),
    ("1 Da", "Pure", "2:8", "TMT_133N", "TMT_135N", 0.25),
    ("1 Da", "Pure", "4:8", "TMT_134N", "TMT_135N", 0.50),
]
REQUIRED_COLUMNS = [
    "PrecursorRecordID",
    "Precursor.Id",
    "Protein.Names",
    "CandidateMS2.ScanID",
    "CandidateMS2.RT",
    "RT",
    "RT.Difference",
    "Precursor.PIF",
    "PIF.Status",
    "PIF.MatchedIsotopeCount",
    "PIF.TargetIsotopeIntensity",
    "Reporter.TotalIntensity",
    "Reporter.DetectedChannelCount",
    *CHANNELS,
]

COLORS = {"2:8": "#D94B45", "4:8": "#216B91"}
PAGE_W = 576.0
PAGE_H = 360.0
SCATTER_W = 166.0
BOX_W = 62.0
PANEL_GAP = 8.0
BLOCK_GAP = 28.0
PLOT_H = 108.0
BLOCK_X = {"0.6 Da": 56.0, "1 Da": 56.0 + SCATTER_W + PANEL_GAP + BOX_W + BLOCK_GAP}
SCATTER_X = {spacing: x0 for spacing, x0 in BLOCK_X.items()}
BOX_X = {spacing: x0 + SCATTER_W + PANEL_GAP for spacing, x0 in BLOCK_X.items()}
ROW_Y = {"Mix": 66.0, "Pure": 205.0}
X_LIMITS = (6.0, 18.0)
X_TICKS = [6, 9, 12, 15, 18]
Y_LIMITS = (-6.0, 2.0)
Y_TICKS = [-6, -4, -2, 0, 2]
REFERENCE_Y = [-1.0, -2.0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uncalibrated", type=Path, required=True)
    parser.add_argument("--calibrated", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def quantile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * p
    lower = int(math.floor(index))
    upper = int(math.ceil(index))
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def box_stats(values: list[float]) -> tuple[float, float, float, float, float]:
    ordered = sorted(values)
    q1 = quantile(ordered, 0.25)
    median = quantile(ordered, 0.50)
    q3 = quantile(ordered, 0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    lower = next((value for value in ordered if value >= lower_bound), ordered[0])
    upper = next((value for value in reversed(ordered) if value <= upper_bound), ordered[-1])
    return lower, q1, median, q3, upper


def select_unique_ms2(path: Path, method: str) -> tuple[pd.DataFrame, dict[str, object]]:
    frame = pd.read_csv(path, usecols=REQUIRED_COLUMNS)
    for column in [
        "PrecursorRecordID",
        "CandidateMS2.ScanID",
        "CandidateMS2.RT",
        "RT",
        "RT.Difference",
        "Precursor.PIF",
        "PIF.MatchedIsotopeCount",
        "PIF.TargetIsotopeIntensity",
        "Reporter.TotalIntensity",
        "Reporter.DetectedChannelCount",
        *CHANNELS,
    ]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    original_rows = len(frame)
    ok = frame.loc[frame["PIF.Status"].eq("ok")].copy()
    ok["Selection.AbsRTDifference"] = ok["RT.Difference"].abs()
    candidate_count = ok.groupby("PrecursorRecordID")["CandidateMS2.ScanID"].transform("size")
    maximum_pif = ok.groupby("PrecursorRecordID")["Precursor.PIF"].transform("max")
    maximum_tie = np.isclose(
        ok["Precursor.PIF"].to_numpy(dtype=np.float64),
        maximum_pif.to_numpy(dtype=np.float64),
        rtol=0,
        atol=1e-12,
    )
    tie_counts = (
        pd.Series(maximum_tie, index=ok.index)
        .groupby(ok["PrecursorRecordID"])
        .transform("sum")
    )
    ok["Selection.OKCandidateCount"] = candidate_count.astype(np.int16)
    ok["Selection.MaxPIFTieCount"] = tie_counts.astype(np.int16)
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
    selected["Selection.Method"] = method
    selected["Selection.Rule"] = (
        "PIF.Status=ok; max PIF; max matched isotope count; max target isotope "
        "intensity; min abs RT difference; min scan ID"
    )

    channel_values = selected[CHANNELS].to_numpy(dtype=np.float64)
    finite = np.isfinite(channel_values)
    all_finite = finite.all(axis=1)
    all_positive = (channel_values > 0).all(axis=1)
    complete = all_finite & all_positive
    selected["All12_PositiveFinite"] = complete
    selected["Invalid_12Channel_Count"] = np.count_nonzero(~finite | (channel_values <= 0), axis=1)
    selected["Min_12Channel_Intensity"] = np.where(
        all_finite,
        np.min(channel_values, axis=1),
        np.nan,
    )
    selected["Median_12Channel_Intensity"] = np.where(
        all_finite,
        np.median(channel_values, axis=1),
        np.nan,
    )
    minimum_indices = np.argmin(np.where(finite, channel_values, np.inf), axis=1)
    selected["Min_12Channel_Name"] = np.where(
        all_finite,
        np.asarray(CHANNELS, dtype=object)[minimum_indices],
        "",
    )
    selected["Min12_Intensity_Percentile"] = selected["Min_12Channel_Intensity"].rank(
        method="average", pct=True
    ) * 100.0
    selected = selected.sort_values(
        ["Min_12Channel_Intensity", "PrecursorRecordID"],
        ascending=[True, True],
        na_position="last",
    ).reset_index(drop=True)

    strict_ecoli = selected["Protein.Names"].fillna("").str.contains("ECOLI", regex=False)
    final = selected.loc[strict_ecoli & selected["All12_PositiveFinite"]].copy()
    final["Min12_Intensity_Percentile_WithinFinalECOLI"] = final[
        "Min_12Channel_Intensity"
    ].rank(method="average", pct=True) * 100.0
    final = final.sort_values(
        ["Min_12Channel_Intensity", "PrecursorRecordID"], ascending=[True, True]
    ).reset_index(drop=True)

    stats = {
        "Method": method,
        "Input.File": str(path),
        "Input.PairRows": original_rows,
        "PIF.Status.OK.PairRows": len(ok),
        "Unique.Precursors.WithOKCandidate": int(ok["PrecursorRecordID"].nunique()),
        "Precursors.WithMultipleOKCandidates": int(
            (ok.groupby("PrecursorRecordID").size() > 1).sum()
        ),
        "Precursors.WithMaxPIFTie": int(
            selected["Selection.MaxPIFTieCount"].gt(1).sum()
        ),
        "Selected.UniqueMS2.Rows": len(selected),
        "Selected.StrictECOLI.Rows": int(strict_ecoli.sum()),
        "Selected.All12PositiveFinite.Rows": int(selected["All12_PositiveFinite"].sum()),
        "Final.StrictECOLI.All12PositiveFinite.Rows": len(final),
        "Excluded.StrictECOLI.SelectedMS2.Not12Complete": int(
            (strict_ecoli & ~selected["All12_PositiveFinite"]).sum()
        ),
    }
    return selected, final, stats


def calculate_ratios(final: pd.DataFrame, method: str) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    metadata_columns = [
        "PrecursorRecordID",
        "Precursor.Id",
        "Protein.Names",
        "CandidateMS2.ScanID",
        "CandidateMS2.RT",
        "Precursor.PIF",
        "PIF.MatchedIsotopeCount",
        "PIF.TargetIsotopeIntensity",
        "Reporter.TotalIntensity",
        "Min_12Channel_Intensity",
        "Min_12Channel_Name",
        "Median_12Channel_Intensity",
        "Min12_Intensity_Percentile_WithinFinalECOLI",
    ]
    for row in final.itertuples(index=False, name=None):
        data = dict(zip(final.columns, row, strict=True))
        for spacing, condition, ratio_type, target, reference, theoretical in COMPARISONS:
            target_value = float(data[target])
            reference_value = float(data[reference])
            measured = target_value / reference_value
            record = {column: data[column] for column in metadata_columns}
            record.update(
                {
                    "Method": method,
                    "Reporter_Channel_Spacing": spacing,
                    "Condition": condition,
                    "Ratio_Type": ratio_type,
                    "Target_Channel": target,
                    "Reference_Channel": reference,
                    "Target_Intensity": target_value,
                    "Reference_Intensity": reference_value,
                    "Theoretical_Ratio": theoretical,
                    "Theoretical_Log2_Ratio": math.log2(theoretical),
                    "Measured_Ratio": measured,
                    "Measured_Log2_Ratio": math.log2(measured),
                    "Log2_Mean_Intensity": math.log2((target_value + reference_value) / 2.0),
                }
            )
            records.append(record)
    return pd.DataFrame.from_records(records)


def make_summary(ratios: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for spacing, condition, ratio_type, target, reference, theoretical in COMPARISONS:
        group = ratios.loc[
            ratios["Reporter_Channel_Spacing"].eq(spacing)
            & ratios["Condition"].eq(condition)
            & ratios["Ratio_Type"].eq(ratio_type)
        ]
        raw = group["Measured_Ratio"].to_numpy(dtype=np.float64)
        log2_values = group["Measured_Log2_Ratio"].to_numpy(dtype=np.float64)
        intensities = group["Log2_Mean_Intensity"].to_numpy(dtype=np.float64)
        visible = (
            (intensities >= X_LIMITS[0])
            & (intensities <= X_LIMITS[1])
            & (log2_values >= Y_LIMITS[0])
            & (log2_values <= Y_LIMITS[1])
        )
        records.append(
            {
                "Method": group["Method"].iloc[0],
                "Reporter_Channel_Spacing": spacing,
                "Condition": condition,
                "Ratio_Type": ratio_type,
                "Target_Channel": target,
                "Reference_Channel": reference,
                "Theoretical_Ratio": theoretical,
                "Theoretical_Log2_Ratio": math.log2(theoretical),
                "N_UniquePrecursors": len(group),
                "Mean_Measured_Ratio": float(np.mean(raw)),
                "Median_Measured_Ratio": float(np.median(raw)),
                "Q1_Measured_Ratio": float(np.quantile(raw, 0.25)),
                "Q3_Measured_Ratio": float(np.quantile(raw, 0.75)),
                "SD_Measured_Ratio": float(np.std(raw, ddof=1)),
                "Mean_Measured_Log2_Ratio": float(np.mean(log2_values)),
                "Median_Measured_Log2_Ratio": float(np.median(log2_values)),
                "Median_Log2_Mean_Intensity": float(np.median(intensities)),
                "N_Visible_In_MA_Axis_Range": int(np.count_nonzero(visible)),
                "N_Outside_MA_Axis_Range_NotExcluded": int(np.count_nonzero(~visible)),
            }
        )
    return pd.DataFrame.from_records(records)


def x_map(value: float, x0: float, width: float) -> float:
    return x0 + (value - X_LIMITS[0]) / (X_LIMITS[1] - X_LIMITS[0]) * width


def y_map(value: float, y0: float, height: float) -> float:
    clamped = max(Y_LIMITS[0], min(Y_LIMITS[1], value))
    return y0 + height - (clamped - Y_LIMITS[0]) / (Y_LIMITS[1] - Y_LIMITS[0]) * height


def gaussian_smooth(values: list[float], sigma: float, radius: int) -> list[float]:
    weights = [math.exp(-(i * i) / (2.0 * sigma * sigma)) for i in range(-radius, radius + 1)]
    smoothed: list[float] = []
    for index in range(len(values)):
        numerator = 0.0
        denominator = 0.0
        for offset, weight in zip(range(-radius, radius + 1), weights):
            j = index + offset
            if 0 <= j < len(values):
                numerator += values[j] * weight
                denominator += weight
        smoothed.append(numerator / denominator)
    return smoothed


def loess_trend(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    points = sorted((x, y) for x, y in points if X_LIMITS[0] <= x <= X_LIMITS[1])
    if len(points) < 25:
        return []
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    n = len(points)
    k = max(80, min(n, int(n * 0.62)))
    lower_grid = max(X_LIMITS[0], xs[0])
    upper_grid = min(X_LIMITS[1], xs[-1])
    trend: list[tuple[float, float]] = []
    for index in range(180):
        x0 = lower_grid + (upper_grid - lower_grid) * index / 179.0
        distances = sorted(abs(x - x0) for x in xs)
        bandwidth = distances[min(k - 1, len(distances) - 1)] or max(distances)
        sw = sx = sy = sxx = sxy = 0.0
        for x, y in zip(xs, ys):
            distance = abs(x - x0) / bandwidth if bandwidth > 0 else 0.0
            if distance >= 1.0:
                continue
            weight = (1.0 - distance**3) ** 3
            sw += weight
            sx += weight * x
            sy += weight * y
            sxx += weight * x * x
            sxy += weight * x * y
        if sw <= 0:
            continue
        denominator = sw * sxx - sx * sx
        if abs(denominator) < 1e-12:
            y0 = sy / sw
        else:
            slope = (sw * sxy - sx * sy) / denominator
            intercept = (sy - slope * sx) / sw
            y0 = intercept + slope * x0
        trend.append((x0, y0))
    smoothed = gaussian_smooth([point[1] for point in trend], sigma=4.5, radius=16)
    smoothed = gaussian_smooth(smoothed, sigma=3.0, radius=10)
    return list(zip([point[0] for point in trend], smoothed))


def draw_panel_frame(
    canvas: base.PDFCanvas,
    x0: float,
    y0: float,
    width: float,
    height: float,
    show_y: bool,
    show_x: bool,
    is_box: bool,
) -> None:
    for reference in REFERENCE_Y:
        y = y_map(reference, y0, height)
        canvas.line(x0, y, x0 + width, y, color="#777777", lw=0.8, dashed=True)
    if show_y:
        canvas.line(x0, y0, x0, y0 + height, color="#000000", lw=1.1)
    canvas.line(x0, y0 + height, x0 + width, y0 + height, color="#000000", lw=1.1)
    if show_y:
        for tick in Y_TICKS:
            y = y_map(tick, y0, height)
            canvas.line(x0 - 4, y, x0, y, color="#000000", lw=0.8)
            canvas.text(x0 - 8, y + 3, f"{tick:g}", size=8.0, anchor="end")
    if is_box:
        for position, label in ((0.28, "2:8"), (0.72, "4:8")):
            x = x0 + width * position
            canvas.line(x, y0 + height, x, y0 + height + 4, color="#000000", lw=0.8)
            if show_x:
                canvas.text(x, y0 + height + 18, label, size=8.8)
    else:
        for tick in X_TICKS:
            x = x_map(tick, x0, width)
            canvas.line(x, y0 + height, x, y0 + height + 4, color="#000000", lw=0.8)
            if show_x:
                canvas.text(x, y0 + height + 18, str(tick), size=8.8)


def ratio_group(ratios: pd.DataFrame, spacing: str, condition: str, ratio_type: str) -> pd.DataFrame:
    return ratios.loc[
        ratios["Reporter_Channel_Spacing"].eq(spacing)
        & ratios["Condition"].eq(condition)
        & ratios["Ratio_Type"].eq(ratio_type)
    ]


def draw_scatter(
    canvas: base.PDFCanvas,
    ratios: pd.DataFrame,
    spacing: str,
    condition: str,
    x0: float,
    y0: float,
    show_y: bool,
    show_x: bool,
) -> None:
    draw_panel_frame(canvas, x0, y0, SCATTER_W, PLOT_H, show_y, show_x, False)
    for ratio_type in ("2:8", "4:8"):
        group = ratio_group(ratios, spacing, condition, ratio_type)
        points = list(zip(group["Log2_Mean_Intensity"], group["Measured_Log2_Ratio"], strict=True))
        visible = [
            (x_map(x, x0, SCATTER_W), y_map(y, y0, PLOT_H))
            for x, y in points
            if X_LIMITS[0] <= x <= X_LIMITS[1] and Y_LIMITS[0] <= y <= Y_LIMITS[1]
        ]
        canvas.points(visible, COLORS[ratio_type], size=0.55)
        trend = loess_trend(points)
        canvas.polyline(
            [(x_map(x, x0, SCATTER_W), y_map(y, y0, PLOT_H)) for x, y in trend],
            COLORS[ratio_type],
            lw=1.25,
        )


def draw_box(
    canvas: base.PDFCanvas,
    ratios: pd.DataFrame,
    spacing: str,
    condition: str,
    x0: float,
    y0: float,
    show_x: bool,
) -> None:
    draw_panel_frame(canvas, x0, y0, BOX_W, PLOT_H, False, show_x, True)
    for index, ratio_type in enumerate(("2:8", "4:8")):
        group = ratio_group(ratios, spacing, condition, ratio_type)
        values = group["Measured_Log2_Ratio"].to_list()
        raw_median = float(group["Measured_Ratio"].median())
        lower, q1, median, q3, upper = box_stats(values)
        center_x = x0 + BOX_W * (0.28 if index == 0 else 0.72)
        box_width = BOX_W * 0.20
        cap_width = box_width * 0.72
        color = COLORS[ratio_type]
        y_lower = y_map(lower, y0, PLOT_H)
        y_q1 = y_map(q1, y0, PLOT_H)
        y_median = y_map(median, y0, PLOT_H)
        y_q3 = y_map(q3, y0, PLOT_H)
        y_upper = y_map(upper, y0, PLOT_H)
        canvas.line(center_x, y_upper, center_x, y_lower, color=color, lw=1.0)
        canvas.line(center_x - cap_width / 2, y_upper, center_x + cap_width / 2, y_upper, color=color, lw=1.0)
        canvas.line(center_x - cap_width / 2, y_lower, center_x + cap_width / 2, y_lower, color=color, lw=1.0)
        canvas.rect(center_x - box_width / 2, y_q3, box_width, max(0.5, y_q1 - y_q3), fill="#FFFFFF", stroke=color, lw=1.1)
        canvas.line(center_x - box_width / 2, y_median, center_x + box_width / 2, y_median, color=color, lw=1.4)
        label_x = center_x - 2 if index == 0 else center_x + 2
        anchor = "end" if index == 0 else "start"
        canvas.text(label_x, y_median - 5, f"{raw_median:.3f}", size=5.7, anchor=anchor, bold=True, color=color)


def draw_pdf(
    ratios: pd.DataFrame,
    final_count: int,
    method_title: str,
    output_pdf: Path,
) -> None:
    canvas = base.PDFCanvas(PAGE_W, PAGE_H)
    canvas.text(PAGE_W / 2, 14, method_title, size=10.5, bold=True)
    canvas.text(85, 36, "Ratio", size=8.0, anchor="end", bold=True)
    legend_x = 98.0
    for ratio_type in ("2:8", "4:8"):
        canvas.rect(legend_x, 31, 12, 3.0, fill=COLORS[ratio_type], stroke=None)
        canvas.text(legend_x + 17, 36, ratio_type, size=7.5, anchor="start", bold=True)
        legend_x += 68.0
    canvas.text(
        250,
        36,
        f"PIF.Status=ok; max-PIF unique MS2; all 12 > 0; n={final_count:,}",
        size=7.2,
        anchor="start",
        bold=True,
    )
    canvas.text(518, 55, "Box labels: raw median", size=6.3, anchor="end", color="#444444")
    for spacing in ("0.6 Da", "1 Da"):
        center = BLOCK_X[spacing] + (SCATTER_W + PANEL_GAP + BOX_W) / 2.0
        canvas.text(center, 55, spacing, size=10.0, bold=True, color="#1A1A1A")
    canvas.text(30.0, ROW_Y["Mix"] + PLOT_H / 2, "Mix", size=10.5, bold=True, color="#1A1A1A", rotate=-90)
    canvas.text(30.0, ROW_Y["Pure"] + PLOT_H / 2, "Pure", size=10.5, bold=True, color="#1A1A1A", rotate=-90)
    for condition in ("Mix", "Pure"):
        for spacing_index, spacing in enumerate(("0.6 Da", "1 Da")):
            draw_scatter(
                canvas,
                ratios,
                spacing,
                condition,
                SCATTER_X[spacing],
                ROW_Y[condition],
                show_y=spacing_index == 0,
                show_x=condition == "Pure",
            )
            draw_box(
                canvas,
                ratios,
                spacing,
                condition,
                BOX_X[spacing],
                ROW_Y[condition],
                show_x=condition == "Pure",
            )
    canvas.text(238.0, 348.0, "Log2 Mean Intensity", size=10.5)
    canvas.text(14.0, 181.0, "Log2 Measured Reporter Ion Ratio", size=10.5, rotate=-90)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    canvas.write_pdf(output_pdf)


def overlap_summary(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    left = old[["PrecursorRecordID", "CandidateMS2.ScanID", "Precursor.PIF"]].rename(
        columns={"CandidateMS2.ScanID": "Uncalibrated.SelectedScanID", "Precursor.PIF": "Uncalibrated.SelectedPIF"}
    )
    right = new[["PrecursorRecordID", "CandidateMS2.ScanID", "Precursor.PIF"]].rename(
        columns={"CandidateMS2.ScanID": "Calibrated.SelectedScanID", "Precursor.PIF": "Calibrated.SelectedPIF"}
    )
    merged = left.merge(right, on="PrecursorRecordID", how="outer", indicator=True)
    merged["Same.SelectedMS2"] = (
        merged["Uncalibrated.SelectedScanID"].eq(merged["Calibrated.SelectedScanID"])
        & merged["_merge"].eq("both")
    )
    return merged


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    methods = [
        (
            "Uncalibrated_10ppm",
            args.uncalibrated,
            "E. coli reporter ion ratios - max-PIF unique MS2 (uncalibrated 10 ppm)",
            "uncalibrated10ppm",
        ),
        (
            "MassCalibrated_3ppm",
            args.calibrated,
            "E. coli reporter ion ratios - max-PIF unique MS2 (mass-calibrated 3 ppm)",
            "masscalibrated3ppm",
        ),
    ]
    selected_by_method: dict[str, pd.DataFrame] = {}
    qc_rows: list[dict[str, object]] = []
    combined_summary: list[pd.DataFrame] = []
    for method, input_path, title, stem in methods:
        selected, final, stats = select_unique_ms2(input_path, method)
        selected_by_method[method] = selected
        for metric, value in stats.items():
            if metric in {"Method", "Input.File"}:
                continue
            qc_rows.append({"Method": method, "Metric": metric, "Value": value})
        qc_rows.append({"Method": method, "Metric": "Input.File", "Value": stats["Input.File"]})
        selected.to_csv(
            args.output_dir / f"maxPIF_uniqueMS2_{stem}_all_selected_sorted_by_min12.csv",
            index=False,
        )
        final.to_csv(
            args.output_dir / f"maxPIF_uniqueMS2_{stem}_strict_ECOLI_12channel_complete_sorted_by_min12.csv",
            index=False,
        )
        ratios = calculate_ratios(final, method)
        ratios.to_csv(args.output_dir / f"ecoli_reporter_ratios_maxPIF_uniqueMS2_{stem}_long.csv", index=False)
        summary = make_summary(ratios)
        summary.to_csv(args.output_dir / f"ecoli_reporter_ratios_maxPIF_uniqueMS2_{stem}_summary.csv", index=False)
        combined_summary.append(summary)
        draw_pdf(
            ratios,
            len(final),
            title,
            args.output_dir / f"ecoli_reporter_ratio_maxPIF_uniqueMS2_{stem}_8x5.pdf",
        )

    pd.DataFrame(qc_rows).to_csv(args.output_dir / "maxPIF_uniqueMS2_selection_QC.csv", index=False)
    pd.concat(combined_summary, ignore_index=True).to_csv(
        args.output_dir / "ecoli_reporter_ratio_maxPIF_uniqueMS2_combined_summary.csv",
        index=False,
    )
    overlap = overlap_summary(
        selected_by_method["Uncalibrated_10ppm"],
        selected_by_method["MassCalibrated_3ppm"],
    )
    overlap.to_csv(args.output_dir / "maxPIF_uniqueMS2_selected_scan_overlap.csv", index=False)
    overlap_stats = pd.DataFrame(
        [
            {"Metric": "Union.UniquePrecursors", "Value": len(overlap)},
            {"Metric": "Present.InBothMethods", "Value": int(overlap["_merge"].eq("both").sum())},
            {"Metric": "Same.SelectedMS2", "Value": int(overlap["Same.SelectedMS2"].sum())},
            {
                "Metric": "Different.SelectedMS2.AmongCommon",
                "Value": int((overlap["_merge"].eq("both") & ~overlap["Same.SelectedMS2"]).sum()),
            },
            {"Metric": "Only.Uncalibrated.HasOK", "Value": int(overlap["_merge"].eq("left_only").sum())},
            {"Metric": "Only.Calibrated.HasOK", "Value": int(overlap["_merge"].eq("right_only").sum())},
        ]
    )
    overlap_stats.to_csv(args.output_dir / "maxPIF_uniqueMS2_selected_scan_overlap_summary.csv", index=False)
    print(f"Output directory: {args.output_dir}")
    print(pd.DataFrame(qc_rows).to_string(index=False))
    print(overlap_stats.to_string(index=False))


if __name__ == "__main__":
    main()
