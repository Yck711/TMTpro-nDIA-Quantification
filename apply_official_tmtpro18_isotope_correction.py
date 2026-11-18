#!/usr/bin/env python3
"""Build the official TMTpro 18-plex isotope matrix and correct reporter intensities."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


CHANNELS = [
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

THEORETICAL_MZ = np.array(
    [
        126.127726,
        127.124761,
        127.131081,
        128.128116,
        128.134436,
        129.131471,
        129.137790,
        130.134825,
        130.141145,
        131.138180,
        131.144500,
        132.141535,
        132.147855,
        133.144890,
        133.151210,
        134.148245,
        134.154565,
        135.151600,
    ],
    dtype=np.float64,
)

ISOTOPE_NAMES = [
    "-2x13C",
    "-13C-15N",
    "-13C",
    "-15N",
    "M+",
    "+15N",
    "+13C",
    "+15N+13C",
    "+2x13C",
]
CHANNEL_INDEX_OFFSETS = np.arange(-4, 5, dtype=np.int64)

# Official reporter-ion isotopic distributions. M+ is defined as 100 and all
# other entries are relative abundances to M+; N/A entries are represented as 0.
OFFICIAL_RELATIVE_DISTRIBUTION = np.array(
    [
        [0, 0, 0, 0, 100, 0.31, 9.09, 0.02, 0.32],
        [0, 0, 0, 0.78, 100, 0, 9.41, 0, 0.33],
        [0, 0, 0.84, 0, 100, 0.23, 8.40, 0.02, 0.27],
        [0, 0.00, 0.82, 0.65, 100, 0, 8.13, 0, 0.26],
        [0.00, 0, 1.44, 0, 100, 0.34, 6.26, 0.00, 0.17],
        [0.00, 0.14, 1.30, 0.89, 100, 0, 7.52, 0, 0.12],
        [0.13, 0, 2.59, 0, 100, 0.32, 6.07, 0.01, 0.09],
        [0.13, 0.00, 2.41, 0.27, 100, 0, 5.58, 0, 0.10],
        [0.25, 0, 3.22, 0, 100, 0.28, 5.06, 0.00, 0.06],
        [0.03, 0.00, 2.78, 0.63, 100, 0, 4.57, 0, 0.12],
        [0.08, 0, 3.94, 0, 100, 0.45, 3.33, 0.00, 0.02],
        [0.07, 0.01, 3.14, 0.73, 100, 0, 3.40, 0, 0.03],
        [0.11, 0, 4.55, 0, 100, 0.43, 1.86, 0.00, 0.00],
        [0.15, 0.01, 3.58, 0.72, 100, 0, 1.80, 0, 0.00],
        [0.22, 0, 4.96, 0, 100, 0.34, 1.03, 0.00, 0],
        [0.30, 0.03, 5.49, 0.62, 100, 0, 1.14, 0, 0],
        [0.16, 0, 5.81, 0, 100, 0.39, 0, 0, 0],
        [0.21, 0.04, 5.90, 0.74, 100, 0, 0, 0, 0],
    ],
    dtype=np.float64,
)

# Display values preserve the distinction between an official N/A entry and
# an explicitly reported 0.00 value. Both contribute zero to the matrix.
OFFICIAL_DISPLAY_DISTRIBUTION = [
    [None, None, None, None, 100, 0.31, 9.09, 0.02, 0.32],
    [None, None, None, 0.78, 100, None, 9.41, None, 0.33],
    [None, None, 0.84, None, 100, 0.23, 8.40, 0.02, 0.27],
    [None, 0.00, 0.82, 0.65, 100, None, 8.13, None, 0.26],
    [0.00, None, 1.44, None, 100, 0.34, 6.26, 0.00, 0.17],
    [0.00, 0.14, 1.30, 0.89, 100, None, 7.52, None, 0.12],
    [0.13, None, 2.59, None, 100, 0.32, 6.07, 0.01, 0.09],
    [0.13, 0.00, 2.41, 0.27, 100, None, 5.58, None, 0.10],
    [0.25, None, 3.22, None, 100, 0.28, 5.06, 0.00, 0.06],
    [0.03, 0.00, 2.78, 0.63, 100, None, 4.57, None, 0.12],
    [0.08, None, 3.94, None, 100, 0.45, 3.33, 0.00, 0.02],
    [0.07, 0.01, 3.14, 0.73, 100, None, 3.40, None, 0.03],
    [0.11, None, 4.55, None, 100, 0.43, 1.86, 0.00, 0.00],
    [0.15, 0.01, 3.58, 0.72, 100, None, 1.80, None, 0.00],
    [0.22, None, 4.96, None, 100, 0.34, 1.03, 0.00, None],
    [0.30, 0.03, 5.49, 0.62, 100, None, 1.14, None, None],
    [0.16, None, 5.81, None, 100, 0.39, None, None, None],
    [0.21, 0.04, 5.90, 0.74, 100, None, None, None, None],
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-xlsx", type=Path, required=True)
    parser.add_argument("--chunksize", type=int, default=50000)
    return parser.parse_args()


def build_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    row_sums = OFFICIAL_RELATIVE_DISTRIBUTION.sum(axis=1)
    normalized = OFFICIAL_RELATIVE_DISTRIBUTION / row_sums[:, None]
    matrix = np.zeros((len(CHANNELS), len(CHANNELS)), dtype=np.float64)

    # Columns are true channels and rows are observed channels: Y = A @ X.
    for true_index, distribution in enumerate(normalized):
        for fraction, offset in zip(
            distribution,
            CHANNEL_INDEX_OFFSETS,
            strict=True,
        ):
            observed_index = true_index + int(offset)
            if 0 <= observed_index < len(CHANNELS):
                matrix[observed_index, true_index] += fraction

    inverse = np.linalg.inv(matrix)
    if not np.allclose(matrix.sum(axis=0), 1.0, rtol=0, atol=1e-12):
        raise ValueError("Correction-matrix columns do not sum to one.")
    return normalized, matrix, inverse


def style_sheet(sheet) -> None:
    fill = PatternFill("solid", fgColor="1F4E78")
    for cell in sheet[1]:
        cell.fill = fill
        cell.font = Font(color="FFFFFF", bold=True)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    for column in sheet.columns:
        values = [str(cell.value) if cell.value is not None else "" for cell in column]
        width = min(max(max(map(len, values)) + 2, 11), 32)
        sheet.column_dimensions[get_column_letter(column[0].column)].width = width


def append_dataframe(sheet, frame: pd.DataFrame, number_format: str = "0.0000000000") -> None:
    sheet.append(list(frame.columns))
    for row in frame.itertuples(index=False, name=None):
        sheet.append(list(row))
    for row in sheet.iter_rows(min_row=2):
        for cell in row:
            if isinstance(cell.value, float):
                cell.number_format = number_format
    style_sheet(sheet)


def write_matrix_workbook(
    output_path: Path,
    normalized: np.ndarray,
    matrix: np.ndarray,
    inverse: np.ndarray,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook = Workbook()
    workbook.remove(workbook.active)

    relative = pd.DataFrame(OFFICIAL_DISPLAY_DISTRIBUTION, columns=ISOTOPE_NAMES)
    relative = relative.map(lambda value: "N/A" if pd.isna(value) else value)
    relative.insert(0, "Reporter_Ion_mz", THEORETICAL_MZ)
    relative.insert(0, "Mass_Tag", CHANNELS)
    relative["Relative_Total"] = OFFICIAL_RELATIVE_DISTRIBUTION.sum(axis=1)
    relative["Normalized_Main_Fraction"] = (
        OFFICIAL_RELATIVE_DISTRIBUTION[:, 4]
        / OFFICIAL_RELATIVE_DISTRIBUTION.sum(axis=1)
    )
    append_dataframe(workbook.create_sheet("Official_Relative"), relative)

    normalized_frame = pd.DataFrame(
        normalized,
        columns=[f"Fraction_{name}" for name in ISOTOPE_NAMES],
    )
    normalized_frame.insert(0, "Reporter_Ion_mz", THEORETICAL_MZ)
    normalized_frame.insert(0, "Mass_Tag", CHANNELS)
    normalized_frame["Fraction_Sum"] = normalized.sum(axis=1)
    append_dataframe(workbook.create_sheet("Normalized_Distribution"), normalized_frame)

    mapping = pd.DataFrame(
        {
            "Isotope_Component": ISOTOPE_NAMES,
            "Channel_Index_Offset": CHANNEL_INDEX_OFFSETS,
            "Meaning": [
                "two 13C losses",
                "one 13C and one 15N loss",
                "one 13C loss",
                "one 15N loss",
                "main reporter isotopologue",
                "one 15N gain",
                "one 13C gain",
                "one 15N and one 13C gain",
                "two 13C gains",
            ],
        }
    )
    append_dataframe(workbook.create_sheet("Isotope_Mapping"), mapping)

    matrix_frame = pd.DataFrame(matrix, columns=[f"True_{name}" for name in CHANNELS])
    matrix_frame.insert(0, "Observed_Channel", CHANNELS)
    append_dataframe(workbook.create_sheet("Correction_Matrix_A"), matrix_frame)

    inverse_frame = pd.DataFrame(inverse, columns=[f"Observed_{name}" for name in CHANNELS])
    inverse_frame.insert(0, "Corrected_Channel", CHANNELS)
    append_dataframe(workbook.create_sheet("Inverse_Matrix"), inverse_frame)

    notes = workbook.create_sheet("Method_Notes")
    notes_rows = [
        ("Item", "Description"),
        ("Matrix model", "Observed vector Y = A x true vector X"),
        ("Correction", "Unclipped corrected vector X = inverse(A) x Y"),
        ("Normalization", "For each true tag, divide all official relative abundances by their row total; M+ is listed as 100 in the official distribution."),
        ("Matrix orientation", "Correction_Matrix_A rows are observed reporter channels; columns are true reporter channels."),
        ("N/A handling", "Official N/A entries are represented numerically as zero."),
        ("Negative correction", "CSV retains the number of negative elements in the unconstrained inverse solution, then reports corrected channel intensities after clipping negative values to zero."),
        ("Input source", "Official Reporter Ion Isotopic Distributions images supplied by the user, including the supplementary TMTpro-134C and TMTpro-135N rows."),
        ("Channels", "TMTpro 18-plex: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N."),
        ("Generated", date.today().isoformat()),
        ("Condition number", f"{np.linalg.cond(matrix):.12g}"),
    ]
    for row in notes_rows:
        notes.append(row)
    style_sheet(notes)
    notes.column_dimensions["A"].width = 24
    notes.column_dimensions["B"].width = 120
    for row in notes.iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical="top")

    workbook.save(output_path)


def correct_csv(
    input_path: Path,
    output_path: Path,
    inverse: np.ndarray,
    chunksize: int,
) -> dict[str, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    corrected_columns = [f"Corrected_{channel}" for channel in CHANNELS]
    first_chunk = True
    row_count = 0
    rows_with_negative = 0
    negative_cell_count = 0
    all_zero_corrected_rows = 0

    for chunk in pd.read_csv(input_path, chunksize=chunksize):
        missing = [channel for channel in CHANNELS if channel not in chunk.columns]
        if missing:
            raise ValueError(f"Input CSV is missing reporter columns: {missing}")

        raw = chunk[CHANNELS].to_numpy(dtype=np.float64)
        if not np.all(np.isfinite(raw)):
            raise ValueError("Input reporter intensities contain non-finite values.")
        if np.any(raw < 0):
            raise ValueError("Input raw reporter intensities contain negative values.")

        # Each input row is one observed vector. X = inv(A) @ Y becomes
        # X_rows = Y_rows @ inv(A).T for a two-dimensional row-major array.
        unconstrained = raw @ inverse.T
        negative_mask = unconstrained < 0
        corrected = np.clip(unconstrained, 0.0, None)

        corrected_frame = pd.DataFrame(corrected, columns=corrected_columns, index=chunk.index)
        corrected_frame["Corrected.Reporter.TotalIntensity"] = corrected.sum(axis=1)
        corrected_frame["Corrected.Reporter.DetectedChannelCount"] = np.count_nonzero(
            corrected > 0,
            axis=1,
        )
        corrected_frame["Corrected.NegativeBeforeClippingCount"] = negative_mask.sum(axis=1)
        corrected_frame["IsotopeCorrection.Method"] = (
            "OfficialFull18plex_NormalizedMatrix_Inverse_ClipNegativeToZero"
        )

        raw_output_columns = [
            *CHANNELS,
            "Reporter.TotalIntensity",
            "Reporter.DetectedChannelCount",
        ]
        metadata_frame = chunk.drop(
            columns=[column for column in raw_output_columns if column in chunk.columns]
        )
        output_chunk = pd.concat(
            [
                metadata_frame.reset_index(drop=True),
                corrected_frame.reset_index(drop=True),
            ],
            axis=1,
        )
        output_chunk.to_csv(
            output_path,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            index=False,
        )
        first_chunk = False

        row_count += len(chunk)
        rows_with_negative += int(np.count_nonzero(negative_mask.any(axis=1)))
        negative_cell_count += int(np.count_nonzero(negative_mask))
        all_zero_corrected_rows += int(np.count_nonzero(np.all(corrected == 0, axis=1)))
        print(f"Corrected {row_count:,} rows", flush=True)

    return {
        "rows": row_count,
        "rows_with_negative_before_clipping": rows_with_negative,
        "negative_cells_before_clipping": negative_cell_count,
        "all_zero_corrected_rows": all_zero_corrected_rows,
    }


def validate_workbook(path: Path) -> None:
    workbook = load_workbook(path, read_only=True, data_only=True)
    expected = [
        "Official_Relative",
        "Normalized_Distribution",
        "Isotope_Mapping",
        "Correction_Matrix_A",
        "Inverse_Matrix",
        "Method_Notes",
    ]
    if workbook.sheetnames != expected:
        raise ValueError(f"Unexpected workbook sheets: {workbook.sheetnames}")
    workbook.close()


def main() -> None:
    args = parse_args()
    normalized, matrix, inverse = build_matrices()
    write_matrix_workbook(args.output_xlsx, normalized, matrix, inverse)
    validate_workbook(args.output_xlsx)
    stats = correct_csv(args.input, args.output_csv, inverse, args.chunksize)
    print(f"Matrix condition number: {np.linalg.cond(matrix):.12g}")
    print(f"Workbook: {args.output_xlsx}")
    print(f"Corrected CSV: {args.output_csv}")
    for key, value in stats.items():
        print(f"{key}: {value:,}")


if __name__ == "__main__":
    main()
