# Independent R reanalysis: anemia and stroke in women with benign uterine disease

A reproducible R project that (1) audits the data, (2) reproduces the original Python analysis, and (3) runs a revised, pre-specified primary analysis of **incident** stroke. The analysis plan (`SAP.md`) was committed before any model was fitted. Every number in the reports comes from `output/`.

## Requirements
- **R** 4.3.3 (tested). Packages: data.table, readxl, sandwich, lmtest, survival, splines, mice, ggplot2, jsonlite, openxlsx, brglm2, Hmisc, officer, flextable. Exact versions are in `output/logs/sessionInfo.txt`.
  - On Ubuntu: `apt-get install r-base-core r-cran-data.table r-cran-survival r-cran-rms r-cran-sandwich r-cran-lmtest r-cran-mice r-cran-ggplot2 r-cran-readxl r-cran-jsonlite r-cran-openxlsx r-cran-brglm2 r-cran-officer r-cran-flextable`
- **Python 3** with pandas: used only for `tools/export_original_results.py`, which reads the original pipeline's saved results so the comparison table never uses retyped numbers.
- **Source data:** in `../data/` (override with the environment variable `STROKE_DATA_DIR`):
  - `stroke_analysis_master.csv`
  - `stroke.xlsx`
  - `labs/`, `dx/` (including `_cache/` for Diagnosis_31), `enc/`, `meds/`

## Run
```bash
cd reanalysis
python3 tools/export_original_results.py ../outputs/_state.pkl output/original_reported.csv
Rscript run_all.R          # about 30-40 min; seed 20261022
```

## Structure

| Path | Content |
|---|---|
| `SAP.md` | Pre-specified analysis plan |
| `R/00_setup.R` | Paths, packages, logging, helpers (Rubin's rules, E-values, formatting) |
| `R/01_import.R` | Import, structural checks, variable dictionary |
| `R/02_labs.R` | Hb, MCV, MCH/MCHC, RDW, platelets re-derived from the raw extracts; validation against the Codebook Hb |
| `R/03_comorbid.R` | Dated comorbidity flags from the diagnosis extracts |
| `R/04_cohort.R` | Cohort flow, outcome audit, analysis datasets |
| `R/05_missing.R` | Missingness; women with vs without a baseline Hb |
| `R/06_reproduce.R` | Reproduction of the original models and comparison table |
| `R/model_fns.R` | Shared model functions (Poisson/Cox, MI, pooling) |
| `R/07_primary.R` | Revised primary analysis (incident stroke), splines, rates, E-values |
| `R/08_secondary.R` | Secondary incident outcomes (ischemic, stroke/TIA, TIA, subtypes) |
| `R/09_timeupdated.R` | Time-updated exposure, lags, testing frequency, landmark |
| `R/10_sensitivity.R` | Sensitivity analyses S2–S11 |
| `R/11_tables_figures.R` | Figures 1–3, S1; `reanalysis_tables.xlsx` / `.docx` |
| `output/tables/*.csv` | Every table, machine-generated |
| `output/figures/` | PNG (300 dpi) and PDF |
| `output/results.json` | Machine-readable key results |
| `output/logs/` | `run_log.md`, `deviations.md`, `model_warnings.md`, `sessionInfo.txt` |
| `output/dag.md` | Causal framework used to choose covariates |
| `derived/` | Patient-level intermediates (git-ignored; never shared) |

## Confidentiality
MRN and date of birth are dropped on import. Records are linked by clinic number in memory only, and outputs contain aggregate results only.
