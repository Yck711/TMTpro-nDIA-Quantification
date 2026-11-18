#!/usr/bin/env python3
"""Plot transparent MA scatters, smooth trends and ratio boxplots."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d


COLORS = {"2:8": "#C84A42", "4:8": "#28698E"}
SPACINGS = ["0.6 Da", "1 Da"]
CONDITIONS = ["Mix", "Pure"]
X_LIMITS = (6.0, 18.0)
X_TICKS = [6, 9, 12, 15, 18]
X_COLUMN = "Log2_Mean_Intensity"
Y_LIMITS = (-6.0, 2.0)
Y_TICKS = [-6, -4, -2, 0, 2]
POINT_ALPHA = 0.20  # 80% transparent.
POINT_SIZE = 1.0
TREND_LINEWIDTH = 1.7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-long", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--method-note", required=True)
    parser.add_argument("--point-alpha", type=float, default=POINT_ALPHA)
    parser.add_argument("--point-size", type=float, default=POINT_SIZE)
    parser.add_argument("--trend-linewidth", type=float, default=TREND_LINEWIDTH)
    parser.add_argument("--x-column", default=X_COLUMN)
    parser.add_argument("--x-label", default="Log2 Mean Intensity")
    parser.add_argument("--x-min", type=float, default=X_LIMITS[0])
    parser.add_argument("--x-max", type=float, default=X_LIMITS[1])
    parser.add_argument("--x-ticks", default=",".join(str(value) for value in X_TICKS))
    return parser.parse_args()


def smooth_binned_median(
    x: np.ndarray,
    y: np.ndarray,
    bins: int = 42,
    sigma: float = 2.2,
) -> tuple[np.ndarray, np.ndarray]:
    valid = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= X_LIMITS[0])
        & (x <= X_LIMITS[1])
    )
    x = x[valid]
    y = y[valid]
    if len(x) < 100:
        return np.asarray([]), np.asarray([])

    edges = np.linspace(X_LIMITS[0], X_LIMITS[1], bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    indices = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    medians = np.full(bins, np.nan, dtype=np.float64)
    counts = np.zeros(bins, dtype=np.int64)
    for index in range(bins):
        values = y[indices == index]
        counts[index] = len(values)
        if len(values) >= 15:
            medians[index] = np.median(values)

    observed = np.isfinite(medians)
    if np.count_nonzero(observed) < 5:
        return np.asarray([]), np.asarray([])
    interpolated = np.interp(centers, centers[observed], medians[observed])
    smoothed = gaussian_filter1d(interpolated, sigma=sigma, mode="nearest")
    # Do not extrapolate the trend into sparse edge bins where a handful of
    # points can produce a visually unstable curve.
    supported = counts >= 15
    first = int(np.flatnonzero(supported)[0])
    last = int(np.flatnonzero(supported)[-1]) + 1
    return centers[first:last], smoothed[first:last]


def style_axis(axis: plt.Axes, show_y: bool, show_x: bool) -> None:
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xticks(X_TICKS)
    axis.set_yticks(Y_TICKS)
    axis.tick_params(axis="both", labelsize=8.5, width=0.8, length=3.5)
    axis.tick_params(labelleft=show_y, labelbottom=show_x)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(show_y)
    axis.spines["bottom"].set_linewidth(0.9)
    axis.spines["left"].set_linewidth(0.9)
    axis.axhline(-2.0, color="#777777", linewidth=0.8, linestyle=(0, (4, 3)), zorder=0)
    axis.axhline(-1.0, color="#777777", linewidth=0.8, linestyle=(0, (4, 3)), zorder=0)


def draw_ma(
    axis: plt.Axes,
    frame: pd.DataFrame,
    spacing: str,
    condition: str,
    show_y: bool,
    show_x: bool,
) -> None:
    style_axis(axis, show_y, show_x)
    for ratio_type in ["2:8", "4:8"]:
        group = frame.loc[
            frame["Reporter_Channel_Spacing"].eq(spacing)
            & frame["Condition"].eq(condition)
            & frame["Ratio_Type"].eq(ratio_type)
        ]
        x = group[X_COLUMN].to_numpy(dtype=np.float64)
        y = group["Measured_Log2_Ratio"].to_numpy(dtype=np.float64)
        visible = (
            np.isfinite(x)
            & np.isfinite(y)
            & (x >= X_LIMITS[0])
            & (x <= X_LIMITS[1])
            & (y >= Y_LIMITS[0])
            & (y <= Y_LIMITS[1])
        )
        axis.scatter(
            x[visible],
            y[visible],
            s=POINT_SIZE,
            color=COLORS[ratio_type],
            alpha=POINT_ALPHA,
            edgecolors="none",
            rasterized=True,
            zorder=1,
        )
        trend_x, trend_y = smooth_binned_median(x, y)
        axis.plot(
            trend_x,
            trend_y,
            color=COLORS[ratio_type],
            linewidth=TREND_LINEWIDTH,
            solid_capstyle="round",
            zorder=4,
        )


def draw_box(
    axis: plt.Axes,
    frame: pd.DataFrame,
    spacing: str,
    condition: str,
    show_x: bool,
) -> None:
    values: list[np.ndarray] = []
    medians: list[float] = []
    for ratio_type in ["2:8", "4:8"]:
        group = frame.loc[
            frame["Reporter_Channel_Spacing"].eq(spacing)
            & frame["Condition"].eq(condition)
            & frame["Ratio_Type"].eq(ratio_type)
        ]
        values.append(group["Measured_Log2_Ratio"].to_numpy(dtype=np.float64))
        medians.append(float(group["Measured_Ratio"].median()))

    boxes = axis.boxplot(
        values,
        positions=[1, 2],
        widths=0.48,
        showfliers=False,
        patch_artist=True,
        whis=1.5,
        medianprops={"linewidth": 1.4},
        whiskerprops={"linewidth": 0.9},
        capprops={"linewidth": 0.9},
    )
    for index, ratio_type in enumerate(["2:8", "4:8"]):
        color = COLORS[ratio_type]
        boxes["boxes"][index].set(facecolor="white", edgecolor=color, linewidth=1.1)
        boxes["medians"][index].set(color=color)
        boxes["whiskers"][2 * index].set(color=color)
        boxes["whiskers"][2 * index + 1].set(color=color)
        boxes["caps"][2 * index].set(color=color)
        boxes["caps"][2 * index + 1].set(color=color)
        median_y = float(np.median(values[index]))
        offset = -0.13 if index == 0 else 0.13
        alignment = "right" if index == 0 else "left"
        axis.text(
            index + 1 + offset,
            median_y - 0.18,
            f"{medians[index]:.3f}",
            color=color,
            fontsize=7.0,
            fontweight="bold",
            ha=alignment,
            va="top",
            clip_on=False,
        )
    axis.set_xlim(0.45, 2.65)
    axis.set_ylim(*Y_LIMITS)
    axis.set_yticks(Y_TICKS)
    axis.tick_params(axis="y", left=False, labelleft=False)
    axis.set_xticks([1, 2], ["2:8", "4:8"] if show_x else ["", ""])
    axis.tick_params(axis="x", labelsize=8.5, width=0.8, length=3.5)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.spines["bottom"].set_linewidth(0.9)
    axis.axhline(-2.0, color="#777777", linewidth=0.8, linestyle=(0, (4, 3)), zorder=0)
    axis.axhline(-1.0, color="#777777", linewidth=0.8, linestyle=(0, (4, 3)), zorder=0)


def main() -> None:
    global POINT_ALPHA, POINT_SIZE, TREND_LINEWIDTH, X_COLUMN, X_LIMITS, X_TICKS
    args = parse_args()
    POINT_ALPHA = args.point_alpha
    POINT_SIZE = args.point_size
    TREND_LINEWIDTH = args.trend_linewidth
    X_COLUMN = args.x_column
    X_LIMITS = (args.x_min, args.x_max)
    X_TICKS = [float(value) for value in args.x_ticks.split(",") if value.strip()]
    frame = pd.read_csv(args.input_long)
    precursor_count = int(frame["PrecursorRecordID"].nunique())

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(8.0, 5.0))
    grid = figure.add_gridspec(
        2,
        4,
        width_ratios=[3.15, 1.05, 3.15, 1.05],
        left=0.095,
        right=0.985,
        bottom=0.13,
        top=0.79,
        wspace=0.18,
        hspace=0.32,
    )

    for row, condition in enumerate(CONDITIONS):
        for spacing_index, spacing in enumerate(SPACINGS):
            ma_axis = figure.add_subplot(grid[row, spacing_index * 2])
            box_axis = figure.add_subplot(grid[row, spacing_index * 2 + 1])
            draw_ma(
                ma_axis,
                frame,
                spacing,
                condition,
                show_y=spacing_index == 0,
                show_x=row == 1,
            )
            draw_box(box_axis, frame, spacing, condition, show_x=row == 1)
            if row == 0:
                left = ma_axis.get_position().x0
                right = box_axis.get_position().x1
                figure.text(
                    (left + right) / 2.0,
                    0.805,
                    spacing,
                    ha="center",
                    va="bottom",
                    fontsize=10.5,
                    fontweight="bold",
                )
        figure.text(
            0.062,
            0.625 if row == 0 else 0.27,
            condition,
            rotation=90,
            ha="center",
            va="center",
            fontsize=10.5,
            fontweight="bold",
        )

    figure.suptitle(args.title, x=0.5, y=0.975, fontsize=10.8, fontweight="bold")
    handles = [
        Line2D([0], [0], color=COLORS[ratio], lw=2.0, label=ratio)
        for ratio in ["2:8", "4:8"]
    ]
    figure.legend(
        handles=handles,
        title="Ratio",
        loc="upper left",
        bbox_to_anchor=(0.095, 0.915),
        frameon=False,
        ncol=2,
        handlelength=1.7,
        columnspacing=1.5,
        fontsize=8.5,
        title_fontsize=8.5,
    )
    figure.text(
        0.64,
        0.895,
        f"{args.method_note}; 12 channels > 0; n={precursor_count:,}",
        ha="center",
        va="center",
        fontsize=7.2,
        fontweight="bold",
    )
    figure.text(
        0.54,
        0.035,
        args.x_label,
        ha="center",
        va="center",
        fontsize=10.0,
    )
    figure.text(
        0.018,
        0.45,
        "Log2 Measured Reporter Ion Ratio",
        rotation=90,
        ha="center",
        va="center",
        fontsize=10.0,
    )
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_pdf, dpi=600, bbox_inches=None)
    plt.close(figure)
    print(f"Wrote: {args.output_pdf}")


if __name__ == "__main__":
    main()
