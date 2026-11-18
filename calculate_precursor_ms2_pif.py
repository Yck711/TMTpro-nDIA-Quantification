#!/usr/bin/env python3
"""Calculate target-specific precursor ion fraction (PIF) for precursor-MS2 pairs."""

from __future__ import annotations

import argparse
import math
import re
import sys
import time
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from pyteomics import mzml


C13_MASS_DIFFERENCE = 1.00335
ISOTOPE_COUNT = 5
ESSENTIAL_COLUMNS = [
    "PrecursorRecordID",
    "Precursor.Mz",
    "Precursor.Charge",
    "CandidateMS2.ScanID",
    "CandidateMS2.RT",
    "Isolation.LowerMz",
    "Isolation.UpperMz",
]

STATUS_LABELS = {
    0: "unprocessed",
    1: "ms2_not_found",
    2: "missing_ms1_before",
    3: "missing_ms1_after",
    4: "invalid_ms1_rt_bracket",
    5: "invalid_isolation_window",
    6: "no_peaks_in_window",
    7: "no_target_isotope",
    8: "ok",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate an interpolated MS1 precursor ion fraction for every "
            "precursor-candidate MS2 pair without applying a PIF filter."
        )
    )
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--mzml", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-docx", type=Path, required=True)
    parser.add_argument("--ppm", type=float, default=10.0)
    parser.add_argument("--chunksize", type=int, default=50000)
    parser.add_argument("--progress-every", type=int, default=50000)
    return parser.parse_args()


def scan_number(spectrum_id: str) -> int | None:
    match = re.search(r"scan=(\d+)", spectrum_id)
    if match:
        return int(match.group(1))
    numbers = re.findall(r"\d+", spectrum_id)
    return int(numbers[-1]) if numbers else None


def spectrum_rt(spectrum: dict) -> float:
    scan = spectrum.get("scanList", {}).get("scan", [{}])[0]
    return float(scan.get("scan start time", math.nan))


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
    if lower >= upper:
        return None
    return target, lower, upper


def make_ms1_record(spectrum: dict) -> dict:
    mz_values = np.asarray(spectrum.get("m/z array", []), dtype=np.float64)
    intensities = np.asarray(spectrum.get("intensity array", []), dtype=np.float64)
    if mz_values.shape != intensities.shape:
        raise ValueError("MS1 m/z and intensity arrays have different lengths.")
    return {
        "scan_id": scan_number(str(spectrum.get("id", ""))),
        "rt": spectrum_rt(spectrum),
        "mz": mz_values,
        "intensity": intensities,
    }


def summed_window_intensity(ms1: dict, lower: float, upper: float) -> float:
    mz_values = ms1["mz"]
    intensities = ms1["intensity"]
    left = int(np.searchsorted(mz_values, lower, side="left"))
    right = int(np.searchsorted(mz_values, upper, side="right"))
    if left >= right:
        return 0.0
    values = intensities[left:right]
    valid = np.isfinite(values) & (values > 0)
    return float(np.sum(values[valid], dtype=np.float64)) if np.any(valid) else 0.0


def target_isotope_intensities(
    ms1: dict,
    precursor_mz: float,
    charge: int,
    lower: float,
    upper: float,
    ppm: float,
) -> np.ndarray:
    result = np.zeros(ISOTOPE_COUNT, dtype=np.float64)
    mz_values = ms1["mz"]
    intensities = ms1["intensity"]

    for isotope_index in range(ISOTOPE_COUNT):
        target = precursor_mz + isotope_index * C13_MASS_DIFFERENCE / charge
        if target < lower or target > upper:
            continue
        tolerance = target * ppm * 1e-6
        search_lower = max(lower, target - tolerance)
        search_upper = min(upper, target + tolerance)
        left = int(np.searchsorted(mz_values, search_lower, side="left"))
        right = int(np.searchsorted(mz_values, search_upper, side="right"))
        if left >= right:
            continue

        local_mz = mz_values[left:right]
        local_intensity = intensities[left:right]
        valid = np.isfinite(local_mz) & np.isfinite(local_intensity) & (local_intensity > 0)
        if not np.any(valid):
            continue
        local_mz = local_mz[valid]
        local_intensity = local_intensity[valid]
        errors = np.abs(local_mz - target)
        minimum_error = float(np.min(errors))
        nearest = np.isclose(errors, minimum_error, rtol=0, atol=1e-12)
        result[isotope_index] = float(np.max(local_intensity[nearest]))
    return result


class PIFResults:
    def __init__(self, row_count: int) -> None:
        self.ms1_before_scan = np.full(row_count, -1, dtype=np.int64)
        self.ms1_before_rt = np.full(row_count, np.nan, dtype=np.float64)
        self.ms1_after_scan = np.full(row_count, -1, dtype=np.int64)
        self.ms1_after_rt = np.full(row_count, np.nan, dtype=np.float64)
        self.weight = np.full(row_count, np.nan, dtype=np.float64)
        self.isotope_intensity = np.full((row_count, ISOTOPE_COUNT), np.nan, dtype=np.float64)
        self.target_before = np.full(row_count, np.nan, dtype=np.float64)
        self.target_after = np.full(row_count, np.nan, dtype=np.float64)
        self.window_before = np.full(row_count, np.nan, dtype=np.float64)
        self.window_after = np.full(row_count, np.nan, dtype=np.float64)
        self.target_interpolated = np.full(row_count, np.nan, dtype=np.float64)
        self.window_interpolated = np.full(row_count, np.nan, dtype=np.float64)
        self.interference = np.full(row_count, np.nan, dtype=np.float64)
        self.pif = np.full(row_count, np.nan, dtype=np.float64)
        self.matched_count = np.zeros(row_count, dtype=np.int8)
        self.matched_mask = np.zeros(row_count, dtype=np.uint8)
        self.status = np.full(row_count, 1, dtype=np.int8)


def mark_status(results: PIFResults, indices: np.ndarray, status: int) -> None:
    results.status[indices] = status


def calculate_buffered_ms2(
    buffered: list[dict],
    previous_ms1: dict | None,
    next_ms1: dict | None,
    candidates: pd.DataFrame,
    results: PIFResults,
    ppm: float,
) -> None:
    for ms2 in buffered:
        indices = ms2["indices"]
        if previous_ms1 is None:
            mark_status(results, indices, 2)
            continue
        if next_ms1 is None:
            mark_status(results, indices, 3)
            continue
        if ms2["window"] is None:
            mark_status(results, indices, 5)
            continue

        rt_before = float(previous_ms1["rt"])
        rt_after = float(next_ms1["rt"])
        ms2_rt = float(ms2["rt"])
        if (
            not all(math.isfinite(value) for value in (rt_before, ms2_rt, rt_after))
            or rt_after <= rt_before
            or ms2_rt < rt_before - 1e-9
            or ms2_rt > rt_after + 1e-9
        ):
            mark_status(results, indices, 4)
            continue

        weight = (ms2_rt - rt_before) / (rt_after - rt_before)
        weight = float(np.clip(weight, 0.0, 1.0))
        _, lower, upper = ms2["window"]
        denominator_before = summed_window_intensity(previous_ms1, lower, upper)
        denominator_after = summed_window_intensity(next_ms1, lower, upper)
        denominator_interpolated = (
            (1.0 - weight) * denominator_before + weight * denominator_after
        )

        results.ms1_before_scan[indices] = int(previous_ms1["scan_id"])
        results.ms1_before_rt[indices] = rt_before
        results.ms1_after_scan[indices] = int(next_ms1["scan_id"])
        results.ms1_after_rt[indices] = rt_after
        results.weight[indices] = weight
        results.window_before[indices] = denominator_before
        results.window_after[indices] = denominator_after
        results.window_interpolated[indices] = denominator_interpolated

        if denominator_interpolated <= 0:
            mark_status(results, indices, 6)
            continue

        for row_index in indices:
            precursor_mz = float(candidates.at[row_index, "Precursor.Mz"])
            charge = int(candidates.at[row_index, "Precursor.Charge"])
            before_isotopes = target_isotope_intensities(
                previous_ms1,
                precursor_mz,
                charge,
                lower,
                upper,
                ppm,
            )
            after_isotopes = target_isotope_intensities(
                next_ms1,
                precursor_mz,
                charge,
                lower,
                upper,
                ppm,
            )
            interpolated = (1.0 - weight) * before_isotopes + weight * after_isotopes
            target_before = float(np.sum(before_isotopes, dtype=np.float64))
            target_after = float(np.sum(after_isotopes, dtype=np.float64))
            target_interpolated = float(np.sum(interpolated, dtype=np.float64))

            results.isotope_intensity[row_index] = interpolated
            results.target_before[row_index] = target_before
            results.target_after[row_index] = target_after
            results.target_interpolated[row_index] = target_interpolated
            results.interference[row_index] = max(
                denominator_interpolated - target_interpolated,
                0.0,
            )

            present = interpolated > 0
            bitmask = sum((1 << i) for i, is_present in enumerate(present) if is_present)
            results.matched_mask[row_index] = bitmask
            results.matched_count[row_index] = int(np.count_nonzero(present))

            if target_interpolated <= 0:
                results.pif[row_index] = 0.0
                results.status[row_index] = 7
                continue

            pif = target_interpolated / denominator_interpolated
            if pif < -1e-12 or pif > 1.0 + 1e-9:
                raise ValueError(
                    f"PIF outside [0,1] for row {row_index}: {pif}; "
                    f"target={target_interpolated}, window={denominator_interpolated}"
                )
            results.pif[row_index] = float(np.clip(pif, 0.0, 1.0))
            results.status[row_index] = 8


def calculate_pif(
    candidate_path: Path,
    mzml_path: Path,
    ppm: float,
    progress_every: int,
) -> tuple[pd.DataFrame, PIFResults, np.ndarray, dict[str, int]]:
    candidates = pd.read_csv(candidate_path, usecols=ESSENTIAL_COLUMNS)
    for column in ESSENTIAL_COLUMNS:
        if column not in candidates.columns:
            raise ValueError(f"Missing required column: {column}")
    candidates = candidates.reset_index(drop=True)

    numeric_columns = [
        "Precursor.Mz",
        "Precursor.Charge",
        "CandidateMS2.ScanID",
        "CandidateMS2.RT",
        "Isolation.LowerMz",
        "Isolation.UpperMz",
    ]
    for column in numeric_columns:
        candidates[column] = pd.to_numeric(candidates[column], errors="raise")
    if candidates["PrecursorRecordID"].isna().any():
        raise ValueError("PrecursorRecordID contains missing values.")

    scan_groups = candidates.groupby("CandidateMS2.ScanID", sort=False).indices
    scan_groups = {int(scan_id): indices for scan_id, indices in scan_groups.items()}
    multiplicity = candidates.groupby("CandidateMS2.ScanID")["PrecursorRecordID"].transform(
        "nunique"
    ).to_numpy(dtype=np.int16)
    results = PIFResults(len(candidates))

    previous_ms1 = None
    buffered: list[dict] = []
    spectra_seen = 0
    ms1_seen = 0
    ms2_seen = 0
    candidate_ms2_seen = 0
    started = time.time()

    with mzml.read(str(mzml_path), use_index=False, iterative=True) as reader:
        for spectrum in reader:
            spectra_seen += 1
            ms_level = int(spectrum.get("ms level", 0))
            if ms_level == 1:
                current_ms1 = make_ms1_record(spectrum)
                ms1_seen += 1
                if buffered:
                    calculate_buffered_ms2(
                        buffered,
                        previous_ms1,
                        current_ms1,
                        candidates,
                        results,
                        ppm,
                    )
                    buffered = []
                previous_ms1 = current_ms1
            elif ms_level == 2:
                ms2_seen += 1
                current_scan_id = scan_number(str(spectrum.get("id", "")))
                if current_scan_id is None or current_scan_id not in scan_groups:
                    continue
                candidate_ms2_seen += 1
                buffered.append(
                    {
                        "scan_id": current_scan_id,
                        "rt": spectrum_rt(spectrum),
                        "window": isolation_window(spectrum),
                        "indices": scan_groups[current_scan_id],
                    }
                )
                if progress_every > 0 and candidate_ms2_seen % progress_every == 0:
                    elapsed = (time.time() - started) / 60
                    completed = int(np.count_nonzero(results.status != 1))
                    print(
                        f"Candidate MS2 seen={candidate_ms2_seen:,}; "
                        f"precursor-MS2 rows processed={completed:,}; "
                        f"elapsed={elapsed:.1f} min",
                        file=sys.stderr,
                        flush=True,
                    )

    if buffered:
        calculate_buffered_ms2(
            buffered,
            previous_ms1,
            None,
            candidates,
            results,
            ppm,
        )

    scan_stats = {
        "spectra_seen": spectra_seen,
        "ms1_seen": ms1_seen,
        "ms2_seen": ms2_seen,
        "candidate_ms2_seen": candidate_ms2_seen,
        "candidate_ms2_expected": len(scan_groups),
    }
    return candidates, results, multiplicity, scan_stats


def bitmask_to_isotopes(mask: int) -> str:
    return ";".join(f"M{i}" for i in range(ISOTOPE_COUNT) if mask & (1 << i))


def write_output_csv(
    input_path: Path,
    output_path: Path,
    results: PIFResults,
    multiplicity: np.ndarray,
    ppm: float,
    chunksize: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    first = True
    offset = 0

    for input_chunk in pd.read_csv(input_path, chunksize=chunksize):
        end = offset + len(input_chunk)
        sl = slice(offset, end)
        pif_frame = pd.DataFrame(
            {
                "MS1.Before.ScanID": pd.array(
                    np.where(results.ms1_before_scan[sl] >= 0, results.ms1_before_scan[sl], None),
                    dtype="Int64",
                ),
                "MS1.Before.RT": results.ms1_before_rt[sl],
                "MS1.After.ScanID": pd.array(
                    np.where(results.ms1_after_scan[sl] >= 0, results.ms1_after_scan[sl], None),
                    dtype="Int64",
                ),
                "MS1.After.RT": results.ms1_after_rt[sl],
                "MS1.InterpolationWeight": results.weight[sl],
                "PIF.M0.Intensity": results.isotope_intensity[sl, 0],
                "PIF.M1.Intensity": results.isotope_intensity[sl, 1],
                "PIF.M2.Intensity": results.isotope_intensity[sl, 2],
                "PIF.M3.Intensity": results.isotope_intensity[sl, 3],
                "PIF.M4.Intensity": results.isotope_intensity[sl, 4],
                "PIF.TargetIntensity.BeforeMS1": results.target_before[sl],
                "PIF.TargetIntensity.AfterMS1": results.target_after[sl],
                "PIF.WindowIntensity.BeforeMS1": results.window_before[sl],
                "PIF.WindowIntensity.AfterMS1": results.window_after[sl],
                "PIF.TargetIsotopeIntensity": results.target_interpolated[sl],
                "PIF.TotalIsolationWindowIntensity": results.window_interpolated[sl],
                "PIF.InterferenceIntensity": results.interference[sl],
                "Precursor.PIF": results.pif[sl],
                "PIF.MatchedIsotopeCount": results.matched_count[sl],
                "PIF.MatchedIsotopes": [
                    bitmask_to_isotopes(int(mask)) for mask in results.matched_mask[sl]
                ],
                "PIF.Status": [STATUS_LABELS[int(code)] for code in results.status[sl]],
                "PIF.MassTolerancePPM": ppm,
                "PIF.IsotopeDefinition": "M0-M4_within_actual_isolation_window",
                "PIF.CalculationMethod": "BracketedMS1_linear_interpolation_of_intensities",
                "CandidateMS2.PrecursorMultiplicity": multiplicity[sl],
            }
        )
        output_chunk = pd.concat(
            [input_chunk.reset_index(drop=True), pif_frame.reset_index(drop=True)],
            axis=1,
        )
        output_chunk.to_csv(
            output_path,
            mode="w" if first else "a",
            header=first,
            index=False,
            na_rep="NA",
        )
        first = False
        offset = end

    if offset != len(results.pif):
        raise ValueError(f"Output row count mismatch: wrote {offset}, expected {len(results.pif)}")


def set_run_font(run, name: str = "宋体", size: float = 10.5, bold: bool = False) -> None:
    run.font.name = name
    run._element.rPr.rFonts.set(qn("w:eastAsia"), name)
    run.font.size = Pt(size)
    run.bold = bold


def add_text(document: Document, text: str, bold_prefix: str | None = None) -> None:
    paragraph = document.add_paragraph()
    if bold_prefix and text.startswith(bold_prefix):
        first = paragraph.add_run(bold_prefix)
        set_run_font(first, bold=True)
        rest = paragraph.add_run(text[len(bold_prefix) :])
        set_run_font(rest)
    else:
        run = paragraph.add_run(text)
        set_run_font(run)


def add_bullets(document: Document, items: list[str]) -> None:
    for item in items:
        paragraph = document.add_paragraph(style="List Bullet")
        run = paragraph.add_run(item)
        set_run_font(run)


def add_numbered(document: Document, items: list[str]) -> None:
    for item in items:
        paragraph = document.add_paragraph(style="List Number")
        run = paragraph.add_run(item)
        set_run_font(run)


def add_table(document: Document, headers: list[str], rows: list[list[object]]) -> None:
    table = document.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    for index, header in enumerate(headers):
        run = table.rows[0].cells[index].paragraphs[0].add_run(str(header))
        set_run_font(run, bold=True)
    for row in rows:
        cells = table.add_row().cells
        for index, value in enumerate(row):
            run = cells[index].paragraphs[0].add_run(str(value))
            set_run_font(run, size=9)


def write_method_docx(
    output_path: Path,
    input_path: Path,
    mzml_path: Path,
    output_csv: Path,
    candidates: pd.DataFrame,
    results: PIFResults,
    multiplicity: np.ndarray,
    scan_stats: dict[str, int],
    ppm: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    document = Document()
    section = document.sections[0]
    section.top_margin = Inches(0.75)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.8)
    section.right_margin = Inches(0.8)

    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("TMTpro-nDIA precursor-MS2 PIF计算方法与分析记录")
    set_run_font(run, name="黑体", size=16, bold=True)

    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = subtitle.add_run(f"生成日期：{date.today().isoformat()}")
    set_run_font(run, size=10)

    for heading in [
        "一、分析目的",
        "二、输入与输出文件",
        "三、PIF的分析单位与定义",
        "四、计算参数",
        "五、逐步计算逻辑",
        "六、输出字段说明",
        "七、异常与缺失值处理",
        "八、数据结构与运行结果",
        "九、质量控制",
        "十、方法局限性",
        "十一、复现命令",
    ]:
        paragraph = document.add_heading(heading, level=1)
        for heading_run in paragraph.runs:
            set_run_font(heading_run, name="黑体", size=13, bold=True)

        if heading == "一、分析目的":
            add_text(
                document,
                "本分析仅为每一个 precursor-candidate MS2 配对计算采集MS2时的目标前体离子纯度（precursor ion fraction, PIF）。本阶段不进行PIF阈值过滤，不删除任何precursor-MS2记录，不选择最终定量MS2，也不修改18通道reporter ion intensity。",
            )
        elif heading == "二、输入与输出文件":
            add_table(
                document,
                ["类型", "文件"],
                [
                    ["候选定量表", str(input_path)],
                    ["原始谱图", str(mzml_path)],
                    ["PIF结果表", str(output_csv)],
                    ["计算脚本", str(Path(__file__).resolve())],
                    ["方法文档", str(output_path)],
                ],
            )
        elif heading == "三、PIF的分析单位与定义":
            add_text(
                document,
                "PIF不是一个MS2扫描固定不变的单一属性，而是一个目标precursor在一个具体MS2 isolation window中的纯度。因此分析单位是 PrecursorRecordID × CandidateMS2.ScanID。同一MS2对应多个precursor时，需要分别计算每个precursor的PIF。",
            )
            add_text(document, "PIF = N(tMS2) / D(tMS2)")
            add_bullets(
                document,
                [
                    "N(tMS2)：MS2采集时刻，isolation window内目标precursor的M0-M4同位素峰强度之和。",
                    "D(tMS2)：MS2采集时刻，isolation window内全部MS1 centroid峰的强度之和。",
                    "PIF接近1表示目标同位素峰占窗口总信号的比例较高；PIF接近0表示共隔离信号较多或未检测到目标同位素峰。",
                ],
            )
        elif heading == "四、计算参数":
            add_table(
                document,
                ["参数", "设置"],
                [
                    ["前体质量容差", f"±{ppm:g} ppm"],
                    ["同位素质量差", "1.00335 Da"],
                    ["目标同位素", "M0、M1、M2、M3、M4"],
                    ["隔离窗口", "mzML中每个MS2的实际Isolation.LowerMz至Isolation.UpperMz"],
                    ["时间估计", "MS2前后两个MS1之间进行线性插值"],
                    ["PIF filtering", "不进行"],
                    ["reporter intensity", "不参与PIF计算且不修改"],
                ],
            )
        elif heading == "五、逐步计算逻辑":
            add_numbered(
                document,
                [
                    "读取候选表，保留每一条 precursor-MS2 配对及其 Precursor.Mz、Precursor.Charge、CandidateMS2.ScanID。",
                    "按 CandidateMS2.ScanID 统计 CandidateMS2.PrecursorMultiplicity，用于描述同一MS2匹配多少个已鉴定precursor；该数量不是PIF。",
                    "流式读取mzML，定位每个候选MS2之前和之后相邻的两个MS1扫描。",
                    "读取该MS2的实际isolation window上下界以及MS2采集RT。",
                    "根据 Mi = Precursor.Mz + i × 1.00335 / charge 计算M0-M4理论m/z，只考虑理论中心位于当前isolation window内的同位素。",
                    f"分别在前后MS1中以每个理论同位素m/z为中心、按±{ppm:g} ppm搜索centroid峰；存在多个峰时先选择质量误差最小的峰，误差相同时选择强度最高的峰。",
                    "分别计算前后MS1中目标同位素强度和N_before、N_after，以及窗口全部峰强度D_before、D_after。",
                    "计算插值权重 w = (RT_MS2 - RT_before) / (RT_after - RT_before)。",
                    "逐同位素插值：I_i(tMS2) = (1-w) × I_i,before + w × I_i,after。",
                    "窗口总强度插值：D(tMS2) = (1-w) × D_before + w × D_after。",
                    "目标强度 N(tMS2) 为插值后M0-M4强度之和，干扰强度为 max[D(tMS2)-N(tMS2), 0]。",
                    "若D(tMS2)>0，则PIF=N(tMS2)/D(tMS2)；计算完成后将PIF字段追加到原始表，不改变原行顺序。",
                ],
            )
        elif heading == "六、输出字段说明":
            add_table(
                document,
                ["字段", "含义"],
                [
                    ["MS1.Before/After.ScanID、RT", "用于插值的前后MS1扫描及保留时间"],
                    ["MS1.InterpolationWeight", "MS2在两个MS1之间的线性插值权重w"],
                    ["PIF.M0-M4.Intensity", "MS2时刻插值得到的各目标同位素强度"],
                    ["PIF.TargetIntensity.Before/AfterMS1", "前后MS1的目标同位素强度和"],
                    ["PIF.WindowIntensity.Before/AfterMS1", "前后MS1的隔离窗口总强度"],
                    ["PIF.TargetIsotopeIntensity", "MS2时刻目标M0-M4插值强度和，即PIF numerator"],
                    ["PIF.TotalIsolationWindowIntensity", "MS2时刻窗口总插值强度，即PIF denominator"],
                    ["PIF.InterferenceIntensity", "窗口总强度减去目标同位素强度"],
                    ["Precursor.PIF", "目标precursor在该候选MS2中的PIF"],
                    ["PIF.MatchedIsotopes", "检测到的目标同位素，例如M0;M1"],
                    ["PIF.Status", "ok、no_target_isotope或无法计算的具体原因"],
                    ["CandidateMS2.PrecursorMultiplicity", "该MS2匹配的唯一precursor记录数"],
                ],
            )
        elif heading == "七、异常与缺失值处理":
            add_table(
                document,
                ["PIF.Status", "处理"],
                [
                    ["ok", "窗口有信号且检测到至少一个目标同位素，正常计算PIF"],
                    ["no_target_isotope", "窗口有信号但没有目标M0-M4峰，PIF=0"],
                    ["no_peaks_in_window", "插值后窗口总强度为0，PIF=NA"],
                    ["missing_ms1_before/after", "缺少用于插值的MS1，PIF=NA"],
                    ["invalid_ms1_rt_bracket", "MS1不能有效包围MS2 RT，PIF=NA"],
                    ["invalid_isolation_window", "MS2缺少有效隔离窗口，PIF=NA"],
                    ["ms2_not_found", "候选表中的ScanID未在mzML中找到，PIF=NA"],
                ],
            )
        elif heading == "八、数据结构与运行结果":
            precursor_counts = candidates.groupby("PrecursorRecordID")["CandidateMS2.ScanID"].nunique()
            ms2_counts = candidates.groupby("CandidateMS2.ScanID")["PrecursorRecordID"].nunique()
            status_counts = Counter(STATUS_LABELS[int(code)] for code in results.status)
            valid_pif = results.pif[np.isfinite(results.pif)]
            add_table(
                document,
                ["统计项目", "结果"],
                [
                    ["precursor-MS2配对数", f"{len(candidates):,}"],
                    ["唯一PrecursorRecordID", f"{candidates['PrecursorRecordID'].nunique():,}"],
                    ["唯一候选MS2", f"{candidates['CandidateMS2.ScanID'].nunique():,}"],
                    ["每个precursor候选MS2中位数", f"{precursor_counts.median():g}"],
                    ["每个MS2匹配precursor中位数", f"{ms2_counts.median():g}"],
                    ["mzML MS1扫描数", f"{scan_stats['ms1_seen']:,}"],
                    ["mzML MS2扫描数", f"{scan_stats['ms2_seen']:,}"],
                    ["成功在mzML中定位的候选MS2", f"{scan_stats['candidate_ms2_seen']:,}"],
                    ["PIF可计算记录数", f"{len(valid_pif):,}"],
                    ["PIF=ok记录数", f"{status_counts.get('ok', 0):,}"],
                    ["PIF=no_target_isotope记录数", f"{status_counts.get('no_target_isotope', 0):,}"],
                ],
            )
            if len(valid_pif):
                quantiles = np.quantile(valid_pif, [0, 0.25, 0.5, 0.75, 0.95, 1])
                add_table(
                    document,
                    ["PIF统计", "值"],
                    [
                        ["Minimum", f"{quantiles[0]:.6f}"],
                        ["25th percentile", f"{quantiles[1]:.6f}"],
                        ["Median", f"{quantiles[2]:.6f}"],
                        ["75th percentile", f"{quantiles[3]:.6f}"],
                        ["95th percentile", f"{quantiles[4]:.6f}"],
                        ["Maximum", f"{quantiles[5]:.6f}"],
                    ],
                )
            add_table(
                document,
                ["PIF.Status", "记录数"],
                [[name, f"{count:,}"] for name, count in sorted(status_counts.items())],
            )
        elif heading == "九、质量控制":
            add_bullets(
                document,
                [
                    "正常PIF必须位于0至1之间；若numerator大于denominator则脚本停止并报告错误。",
                    "同一个MS2的窗口总强度和前后MS1信息应相同，但不同precursor因理论同位素位置不同而具有不同target intensity和PIF。",
                    "输出行数必须与输入候选表完全一致，原始行顺序保持不变。",
                    "PIF计算不使用TMT reporter ion intensity，因此同位素污染校正不会改变PIF。",
                ],
            )
        elif heading == "十、方法局限性":
            add_bullets(
                document,
                [
                    "PIF基于MS1 centroid峰和理论同位素位置估计隔离纯度，不能证明窗口内信号在分子层面唯一属于目标precursor。",
                    "若其他共洗脱precursor的峰与目标同位素峰在质量容差内重叠，该峰可能被计入目标numerator，使PIF偏高。",
                    "低强度MS1峰可能因未被centroid化或低于检测阈值而缺失，使PIF不稳定。",
                    "PIF只描述前体隔离纯度，不能单独评估reporter ion区的化学噪声、同位素污染或碎片共干扰。",
                    "本阶段不根据PIF选择最终MS2；最终选择规则需要另行预先定义并与导师讨论。",
                ],
            )
        elif heading == "十一、复现命令":
            command = (
                "PYTHONPATH=.analysis_python python3 "
                "mZML/Python/calculate_precursor_ms2_pif.py "
                f"--input-csv {input_path} --mzml {mzml_path} "
                f"--output-csv {output_csv} --output-docx {output_path} "
                f"--ppm {ppm:g}"
            )
            paragraph = document.add_paragraph()
            run = paragraph.add_run(command)
            set_run_font(run, name="Courier New", size=8.5)

    document.save(output_path)


def validate_results(results: PIFResults, input_rows: int) -> dict[str, object]:
    if len(results.pif) != input_rows:
        raise ValueError("PIF result length does not equal input row count.")
    valid = np.isfinite(results.pif)
    if np.any(results.pif[valid] < 0) or np.any(results.pif[valid] > 1):
        raise ValueError("Finite PIF values outside [0,1].")
    ok = results.status == 8
    no_target = results.status == 7
    if np.any(~np.isfinite(results.pif[ok])):
        raise ValueError("Status=ok contains missing PIF values.")
    if np.any(results.pif[no_target] != 0):
        raise ValueError("Status=no_target_isotope does not have PIF=0.")
    return {
        "input_rows": input_rows,
        "finite_pif_rows": int(np.count_nonzero(valid)),
        "missing_pif_rows": int(np.count_nonzero(~valid)),
        "status_counts": Counter(STATUS_LABELS[int(code)] for code in results.status),
    }


def main() -> None:
    args = parse_args()
    if args.ppm <= 0:
        raise ValueError("Mass tolerance must be positive.")
    started = time.time()
    candidates, results, multiplicity, scan_stats = calculate_pif(
        args.input_csv,
        args.mzml,
        args.ppm,
        args.progress_every,
    )
    validation = validate_results(results, len(candidates))
    write_output_csv(
        args.input_csv,
        args.output_csv,
        results,
        multiplicity,
        args.ppm,
        args.chunksize,
    )
    write_method_docx(
        args.output_docx,
        args.input_csv,
        args.mzml,
        args.output_csv,
        candidates,
        results,
        multiplicity,
        scan_stats,
        args.ppm,
    )
    elapsed = (time.time() - started) / 60
    print(f"Output CSV: {args.output_csv}")
    print(f"Method DOCX: {args.output_docx}")
    print(f"Rows: {validation['input_rows']:,}")
    print(f"Finite PIF: {validation['finite_pif_rows']:,}")
    print(f"Missing PIF: {validation['missing_pif_rows']:,}")
    print(f"Status counts: {dict(validation['status_counts'])}")
    print(f"Elapsed: {elapsed:.2f} min")


if __name__ == "__main__":
    main()
