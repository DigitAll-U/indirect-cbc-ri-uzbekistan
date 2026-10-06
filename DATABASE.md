# Analytic dataset

To run the analysis on your own data, provide one CSV file,
**`data/cbc_cohort.csv`**, with one row per blood count, and start at
step 2 (`validate_panel.py`).

## Columns

| Column | Content | Unit |
|---|---|---|
| `ID` | Any unique row number | |
| `SEX` | `Man` or `Woman` | |
| `AGE` | Age at sampling | years |
| `HGB` | Haemoglobin | g/L |
| `RBC` | Red blood cell count | 10¹²/L |
| `HCT` | Haematocrit | % |
| `MCV` | Mean cell volume | fL |
| `MCH` | Mean cell haemoglobin | pg |
| `MCHC` | Mean cell haemoglobin concentration | g/L |
| `RDW_CV`, `RDW_SD` | Red cell distribution width | %, fL |
| `PLT` | Platelet count | 10⁹/L |
| `MPV` | Mean platelet volume | fL |
| `PCT` | Plateletcrit | % |
| `PDW` | Platelet distribution width | as reported |
| `WBC` | White blood cell count | 10⁹/L |
| `NEU_ABS`, `LYMPH_ABS`, `MON_ABS`, `EOS_ABS`, `BAS_ABS` | Differential, absolute | 10⁹/L |
| `NEU_PCT`, `LYMPH_PCT`, `MON_PCT`, `EOS_PCT`, `BAS_PCT` | Differential, relative | % of WBC |
| `IMG_ABS`, `IMG_PCT`, `NRBC_ABS` | Immature granulocytes, nucleated red cells | 10⁹/L, % |

`SEX`, `AGE` and `HGB` are required. Any other parameter may be missing; a parameter that
is missing gets no reference interval, and the checks that depend on it are skipped.
Any adult whose `SEX` is not exactly `Man` is placed in the women's partition.

## Example

| ID | SEX | AGE | HGB | RBC | HCT | MCV | MCH | MCHC | PLT | MPV | PCT | WBC |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Woman | 34 | 128 | 4.41 | 38.6 | 87.5 | 29.0 | 332 | 268 | 9.8 | 0.263 | 6.42 |
| 2 | Man | 52 | 151 | 5.02 | 44.9 | 89.4 | 30.1 | 336 | 231 | 10.2 | 0.236 | 7.15 |
| 3 | Woman | 9 | 124 | 4.63 | 37.2 | 80.3 | 26.8 | 333 | | | | 8.03 |

Record 3 has no platelet results; those cells are left empty.

## Rules

- **One row per blood count.** If your laboratory system exports one row per test,
  pivot the data first.
- **Use the units above.** Step 2 recomputes MCV, MCH, MCHC and plateletcrit from the
  other values, so data in other units is rejected. For example, haemoglobin in
  g/dL × 10 = g/L, and haematocrit as a fraction × 100 = %.
- **Numbers only.** Leave a cell empty if the result is missing or failed; do not
  enter flags, text or `0` for a failed measurement.
- **Apply the step 1 exclusions.** Remove records aged under 1 or over 122 years and
  duplicate records. Empty any MCV outside 44.5–136.9 fL, platelet count outside
  1–1703 × 10⁹/L, or haemoglobin below 14 g/L. See `prepare_cohort.py`.
- **Keep abnormal results.** The method separates healthy from pathological results
  on its own; removing abnormal values biases the intervals.
- **CSV format:** UTF-8, comma-separated, with `.` as the decimal point.
