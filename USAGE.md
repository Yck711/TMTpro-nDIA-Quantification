# Usage

This document explains how to run the five scripts and obtain a TMTpro 18-plex quantification table with one MS2 scan per precursor. The method summary is in [README.md](README.md).

Run every command from the directory that contains the scripts. Process one DIA-NN `Run` and its matching mzML at a time. Do not concatenate runs before a later step. Step 1 assigns `PrecursorRecordID` as 1…N in the order of that run's rows in the parquet. Later steps use that identifier.

## 1. Prepare the directory

Place these five quantification scripts in one directory:

```text
extract_candidate_ms2_raw_reporters.py
apply_official_tmtpro18_isotope_correction.py
calculate_precursor_ms2_pif.py
build_all_species_maxpif_10ppm_corrected_18channel_master.py
make_strict_one_to_one_reassigned_quant_and_plot.py
```

Step 5 imports three more files. Put them in the same directory:

```text
make_strict_one_to_one_quant_and_plot.py
select_max_pif_unique_ms2_and_plot_ratios.py
plot_unique_ms2_ratio_transparent_smooth.py
```

Also create `TMTpro_MAplot_analysis.py` in that directory:

```python
class PDFCanvas:
    pass
```

This repository does not contain the original `TMTpro_MAplot_analysis` module. The ratio summary used by step 5 does not call `PDFCanvas`, but the import fails if the module is missing.

Do not create `mZML/Python/` two directory levels above the scripts. Step 5 searches that path first. An old copy there overrides the scripts in the working directory.

Install the dependencies:

```bash
python -m pip install "numpy>=1.24" "pandas>=2.0" pyarrow pyteomics openpyxl python-docx matplotlib scipy
```

Python 3.10 or newer is required. `matplotlib` and `scipy` are used only by the E. coli ratio PDF at the end of step 5. Install them anyway if you run step 5: after writing the quantification tables, `main()` always attempts to draw that figure.

## 2. Inputs

### DIA-NN `report.parquet`

Use the DIA-NN main report, not a protein matrix. These columns are required:

```text
Run
Precursor.Id
Precursor.Charge
Precursor.Mz
RT
RT.Start
RT.Stop
Protein.Names
```

`--run` must match the `Run` column exactly, including case. The default run name inside step 1 is one file from this study. Pass `--run` for any other dataset.

Step 1 drops, and counts on standard error, rows with charge ≤ 0, non-finite m/z or RT, or `RT.Start` > `RT.Stop`.

### mzML

- The file must be the acquisition that corresponds to that run.
- Spectra must be centroided, with m/z arrays in ascending order. Profile spectra are not supported.
- The spectrum id should contain `scan=<integer>`. If it does not, the last integer in the id is used.
- Every MS2 spectrum must carry an isolation window: `isolation window target m/z`, plus either lower/upper offsets or lower/upper limits. An MS2 spectrum without a window is not matched.
- DIA-NN RT and the mzML `scan start time` must use the same unit. For a standard mzML, pyteomics returns minutes.

### Theoretical reporter m/z

Step 1 reads reporter masses from a CSV rather than hard-coding them. The header must be `Channel,Theoretical_mz`, and the 18 rows must be in exactly this order. Any other order or spelling raises an error. Save the table as `TMTpro18_reporter_ion_theoretical_mz.csv`:

```csv
Channel,Theoretical_mz
TMT_126,126.127726
TMT_127N,127.124761
TMT_127C,127.131081
TMT_128N,128.128116
TMT_128C,128.134436
TMT_129N,129.131471
TMT_129C,129.137790
TMT_130N,130.134825
TMT_130C,130.141145
TMT_131N,131.138180
TMT_131C,131.144500
TMT_132N,132.141535
TMT_132C,132.147855
TMT_133N,133.144890
TMT_133C,133.151210
TMT_134N,134.148245
TMT_134C,134.154565
TMT_135N,135.151600
```

### Isotope-distribution workbook

The correction values are embedded in `apply_official_tmtpro18_isotope_correction.py`. Step 2 does **not** read an Excel file. Publish this workbook next to the scripts so the matrix can be inspected:

`TMTpro18_official_isotopic_distribution_and_correction_matrix.xlsx`

SHA-256: `9af880a5a7701f44881c4bc5afeff5e825a46bc259adb117c8f893cd283eef4c`

| Sheet | Contents |
|---|---|
| `Official_Relative` | Official relative abundances. M+ = 100. N/A stays N/A. |
| `Normalized_Distribution` | Each channel divided by its row total. |
| `Isotope_Mapping` | Channel-index offset from −2×13C through +2×13C. |
| `Correction_Matrix_A` | Correction matrix. Rows are observed channels; columns are true channels. |
| `Inverse_Matrix` | Inverse of A. Quantification uses this sheet: X = inv(A) × Y. |
| `Method_Notes` | Matrix orientation, N/A treated as zero, and clipping of negative solutions. |

Do not publish `official_tmtpro18_correction_matrix_used.xlsx` from an individual run as the reference table. That file is a copy written by step 2. Step 2 still writes a workbook for the matrix it used, so the run can be checked.

## 3. Run the five steps

Replace the paths with your files. `YOUR_RUN_NAME` must match the parquet `Run` column.

```bash
python extract_candidate_ms2_raw_reporters.py \
  --parquet /data/report.parquet \
  --mzml /data/your_run.mzML \
  --reporter-mz TMTpro18_reporter_ion_theoretical_mz.csv \
  --output 01_raw_candidates.csv \
  --run "YOUR_RUN_NAME" \
  --ppm 5

python apply_official_tmtpro18_isotope_correction.py \
  --input 01_raw_candidates.csv \
  --output-csv 02_isotope_corrected.csv \
  --output-xlsx 02_correction_matrix.xlsx

python calculate_precursor_ms2_pif.py \
  --input-csv 01_raw_candidates.csv \
  --mzml /data/your_run.mzML \
  --output-csv 03_with_pif.csv \
  --output-docx 03_pif_method.docx \
  --ppm 10

python build_all_species_maxpif_10ppm_corrected_18channel_master.py \
  --pif 03_with_pif.csv \
  --corrected-reporters 02_isotope_corrected.csv \
  --output 04_maxpif_master.csv \
  --qc-output 04_maxpif_master_QC.csv

python make_strict_one_to_one_reassigned_quant_and_plot.py \
  --master 04_maxpif_master.csv \
  --candidate-pif 03_with_pif.csv \
  --corrected-reporters 02_isotope_corrected.csv \
  --output-dir 05_final
```

Steps 2 and 3 both depend only on step 1, so either may run first. PIF must read the raw candidate table from step 1. Do not pass the corrected table to step 3. Isotope correction does not enter the PIF calculation.

To reproduce this method, keep `--ppm 5` for reporter ions and `--ppm 10` for PIF. Do not switch PIF to a mass-calibrated 3 ppm or 5 ppm tolerance.

Step 1 scans the whole mzML. The number of candidate rows can be much larger than the number of precursors. Progress is written to standard error. Step 3 loads the candidate table into memory and scans the mzML again. It keeps only the MS1 scans on either side of the current MS2, but the candidate table itself can be large.

## 4. What each step does

### Step 1: match MS2 scans and extract raw reporter ions

An MS2 scan is a candidate for a precursor only when both conditions are true:

- The MS2 retention time lies between that precursor's `RT.Start` and `RT.Stop`, inclusive.
- Precursor M0, or M1 = M0 + 1.00335 / charge, lies inside the isolation window recorded for that MS2 in the mzML, inclusive.

Each reporter ion is the sum of centroid intensity within ± 5 ppm of its theoretical m/z. No peak in the window gives 0.

The output is `01_raw_candidates.csv`. One row is one precursor × one candidate MS2, not one row per precursor. `Match.Type` is `M0`, `M1`, or `M0+M1`. `RT.Difference` is MS2 RT minus precursor RT. The last 18 columns are uncorrected intensities, named `TMT_126` through `TMT_135N`.

`precursors without candidate MS2` on standard error counts precursors that have no MS2 under the RT and isolation-window rules. Those precursors do not appear in later tables.

### Step 2: isotope-impurity correction

For the 18 raw intensities Y in each row, the script solves X = inv(A) × Y. A comes from the official Thermo TMTpro 18-plex reporter isotopic distribution: M+ is 100, N/A is 0, and each row is divided by its own total before it is entered into the matrix. Rows of A are observed channels and columns are true channels. Negative values in the unclipped solution are counted in `Corrected.NegativeBeforeClippingCount`. Intensities written to the CSV have already been clipped to 0.

Corrected columns are named `Corrected_TMT_126` through `Corrected_TMT_135N`. This step does not calculate or modify PIF. `02_correction_matrix.xlsx` records the relative abundances, the normalized distribution, A, and inv(A). It is for inspection, not an input to the next step.

The step stops if any raw intensity is negative or non-finite.

### Step 3: precursor ion fraction

For every row from step 1, intensities are linearly interpolated to the MS2 retention time from the MS1 scan immediately before and the MS1 scan immediately after that MS2.

- Numerator: M0–M4. An isotope whose theoretical m/z is outside the actual isolation window is excluded. Inside the window, the peak is the positive centroid within ± 10 ppm that has the smallest mass error. If two peaks have the same mass error, the more intense one is kept.
- Denominator: the sum of all positive centroids in that isolation window.
- `Precursor.PIF` = interpolated numerator / interpolated denominator.

This step does **not** apply a PIF cutoff. `PIF.Status` means:

| Status | Meaning |
|---|---|
| `ok` | Eligible for steps 4 and 5. |
| `ms2_not_found` | The candidate scan is not in this mzML. |
| `missing_ms1_before` / `missing_ms1_after` | There is no MS1 before or after the MS2. |
| `invalid_ms1_rt_bracket` | The two MS1 retention times do not bracket the MS2. |
| `invalid_isolation_window` | The isolation window is invalid. |
| `no_peaks_in_window` | The window contains no positive peaks. |
| `no_target_isotope` | The numerator is 0. PIF is stored as 0 and the row cannot be selected. |

`03_pif_method.docx` is a method note and a count for this run. Quantification does not need it.

The mzML used for PIF must be the same file used in step 1. If scan numbers do not match, many rows become `ms2_not_found`.

### Step 4: highest PIF for each precursor

Only candidates with `PIF.Status = ok` are ranked. One row is kept for each `PrecursorRecordID`, in this order:

1. Higher `Precursor.PIF`.
2. Higher `PIF.MatchedIsotopeCount`.
3. Higher `PIF.TargetIsotopeIntensity`.
4. Smaller absolute `RT.Difference`.
5. Smaller `CandidateMS2.ScanID`.

Corrected intensities from step 2 are then joined on `PrecursorRecordID` + `CandidateMS2.ScanID`.

This table is **not** the final result. One MS2 may still be selected for more than one precursor. `SelectedMS2.IsShared` marks that case. A precursor with no `ok` candidate is absent.

`04_maxpif_master_QC.csv` reports the number of input precursors, how many had an `ok` candidate, how many had a tied maximum PIF, and how many selected rows lacked corrected intensities. The last count should be 0 when steps 1–3 come from the same candidate table.

### Step 5: one precursor per MS2

Step 5 uses `PrecursorRecordID` from the master table to decide which precursors enter assignment. Candidates and PIF values come from the full step-3 table. Corrected intensities come from step 2.

Each precursor proposes MS2 scans in its own priority order. Each MS2 keeps the higher-priority precursor. A displaced precursor proposes its next candidate. The priority order is:

1. Higher PIF.
2. Higher target-isotope intensity.
3. More matched isotopes.
4. Smaller absolute RT difference.
5. Smaller scan ID.

Items 2 and 3 are reversed relative to step 4. The final assignment places target-isotope intensity before the matched-isotope count. That is the order used by the locked method.

A precursor that exhausts its candidates without receiving an MS2 is omitted from the final table. `Assignment.UnassignedPrecursors` in the QC file is that count.

Quantification files in `05_final/`:

| File | Contents |
|---|---|
| `precursor_quantification_allSpecies_strictOneToOneReassigned_uncalibrated10ppm_officialIsotopeCorrected_18channel.csv` | Final table: 18 corrected channels, PIF, MS1 fields, and whether the assignment was changed by competition. Use this file for any other experimental design. |
| The matching `_12channel_minimal.csv` | Only the 12 channels used in this study. The channel list is fixed. It is not the 12 most intense channels. |
| `precursor_quantification_strictOneToOneReassigned_uncalibrated10ppm_officialIsotopeCorrected_QC.csv` | Whether assignment is one-to-one, whether every PIF is in [0, 1], and how many precursors were reassigned or left unassigned. |

In the final table each `PrecursorRecordID` and each `CandidateMS2.ScanID` occurs once. `StrictMS2.WasReassigned = True` means the precursor did not keep its first-choice scan and received a later one after competition. `StrictMS2.CompetitionCount` is how many precursors listed that MS2 as an `ok` candidate.

The 12 channels written to the minimal table are:

```text
127N  127C  128N  128C  129N  129C
130N  131N  132N  133N  134N  135N
```

The 18-channel table also contains `126`, `130C`, `131C`, `132C`, `133C`, and `134C`. For a different labeling scheme, use the 18-channel table and ignore the 12-channel file.

## 5. E. coli ratio figure used only by this study

After both quantification tables are on disk, step 5 also:

1. Keeps rows whose `Protein.Names` contain the literal string `ECOLI`.
2. Keeps rows in which all 12 experimental channels are finite and greater than 0.
3. Calculates ratios for this study's mixing design and writes a PDF.

The ratios are hard-coded in the helper scripts. They are not read from a sample table:

| Spacing | Condition | Ratio | Numerator / denominator | Theoretical ratio |
|---|---|---|---|---|
| 0.6 Da | Mix | 2:8 | 127N / 128N | 0.25 |
| 0.6 Da | Mix | 4:8 | 127C / 128N | 0.50 |
| 0.6 Da | Pure | 2:8 | 128C / 129C | 0.25 |
| 0.6 Da | Pure | 4:8 | 129N / 129C | 0.50 |
| 1 Da | Mix | 2:8 | 130N / 132N | 0.25 |
| 1 Da | Mix | 4:8 | 131N / 132N | 0.50 |
| 1 Da | Pure | 2:8 | 133N / 135N | 0.25 |
| 1 Da | Pure | 4:8 | 134N / 135N | 0.50 |

`05_final/` therefore also contains:

- `ecoli_reporter_ratios_..._long.csv`
- `ecoli_reporter_ratios_..._summary.csv`
- `ecoli_reporter_ratio_..._transparent80_smooth_8x5.pdf`

If no protein name contains `ECOLI`, or if no row has all 12 channels above 0, the quantification CSVs are already written, but the process then stops with `IndexError` on the empty ratio table. The QC file and the PDF are not written. The 18-channel CSV is still the valid quantification result. Do not interpret the ratio tables for a different channel design.

## 6. Optional PIF and intensity filter

Quantification ends at step 5. The MA plots in this study applied two further filters, and only to E. coli rows whose 12 experimental channels were all above 0:

1. Keep `Precursor.PIF >= 0.70`.
2. Remove the lowest 10% by the minimum corrected intensity across those 12 channels.

Steps 3 and 5 do not do this. The earlier PIF ≥ 0.75 rule, and taking the per-channel maximum across several MS2 scans of one precursor, are not this result.

To reproduce the plotted subset from this study, filter the 12-channel table yourself. Change the channel list for any other labeling scheme:

```python
import pandas as pd

channels = [
    "Corrected_TMT_127N", "Corrected_TMT_127C", "Corrected_TMT_128N",
    "Corrected_TMT_128C", "Corrected_TMT_129N", "Corrected_TMT_129C",
    "Corrected_TMT_130N", "Corrected_TMT_131N", "Corrected_TMT_132N",
    "Corrected_TMT_133N", "Corrected_TMT_134N", "Corrected_TMT_135N",
]
table = pd.read_csv(
    "05_final/precursor_quantification_allSpecies_"
    "strictOneToOneReassigned_uncalibrated10ppm_"
    "officialIsotopeCorrected_12channel_minimal.csv"
)
ecoli = table[table["Protein.Names"].fillna("").str.contains("ECOLI", regex=False)].copy()
values = ecoli[channels].apply(pd.to_numeric, errors="coerce")
ecoli = ecoli[values.notna().all(axis=1) & values.gt(0).all(axis=1)]
ecoli = ecoli[ecoli["Precursor.PIF"] >= 0.70]
minimum = ecoli[channels].min(axis=1)
ecoli = ecoli[minimum > minimum.quantile(0.10)]
```

## 7. Check the result

When step 5 finishes normally, the QC file should show:

- `Assignment.UniquePrecursors` equal to `Assignment.StrictOneToOneRows`
- `Assignment.UniqueMS2` equal to the number of rows
- `Assignment.MissingCorrectedReporterMatch` equal to 0
- `Validation.PIFMinimum` ≥ 0 and `Validation.PIFMaximum` ≤ 1
- `PIF.Status` equal to `ok` on every row

The final table can have fewer rows than step 4. The missing precursors lost the MS2 competition. The table was not truncated.

## 8. Common failures

| What you see | Cause |
|---|---|
| Step 1 reports no valid precursors for the run | `--run` does not match the `Run` column, or the parquet is not the DIA-NN main report. |
| Reporter-channel order error | `Channel` order or spelling in the CSV differs from the table in section 2. |
| Step 1 matches almost nothing | RT units differ, or the mzML does not contain isolation windows. |
| Many `ms2_not_found` rows | Step 3 used a different mzML from step 1, or spectrum ids have no stable scan number. |
| Missing corrected reporters, or a non-zero corrected-match count in QC | Steps 2, 3, and 4 were not all produced from the same step-1 table. |
| Step 5: `ModuleNotFoundError: TMTpro_MAplot_analysis` | The stub in section 1 is not in the working directory. |
| Step 5: `ModuleNotFoundError` for `make_strict_one_to_one_quant_and_plot` or `select_max_pif_unique_ms2_and_plot_ratios` | Those two helper scripts are not in the same directory, or an old copy in `mZML/Python/` two levels up is being imported instead. |
| Step 5: `IndexError` after the CSVs are written | No row has `ECOLI` in the protein name and all 12 channels above 0. The 18-channel quantification table is already valid. The ratio figure is specific to this study. |
| Out of memory | Steps 4 and 5 load the full candidate table. Check the row count from step 1, and do not merge several runs into one table. |
