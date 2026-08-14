# Stroke prevalence in non-cancerous uterine disease

R pipeline that reads the abstraction workbook, restricts to eligible patients,
estimates stroke prevalence, and runs subgroup / adjusted analyses.

## Run

1. Install packages (first run only):

```r
install.packages(c("readxl","janitor","lubridate","dplyr","tidyr","stringr",
                   "binom","broom","purrr","sandwich","lmtest","ggplot2",
                   "scales","forcats","openxlsx","readr","brglm2"))
```

2. Set `input_file` and `sheet_name` at the top of
   `stroke_prev_pipeline_no_mice.R`.
3. `Rscript stroke_prev_pipeline_no_mice.R`

Everything lands in `nomice/` — one Excel workbook with a sheet per result,
plus PNG forest/bar plots in `nomice/plots/`.

## Analysis population

`eligible == 1` **and** not active malignancy (`malignancy_active == 0` or
missing). Because active malignancy is an exclusion, `race` and
`malignancy_active` are dropped from every adjusted model (`adj_exclude`).

## Outcomes

- Primary: `stroke_any`
- Sensitivity: `stroke_any_incl_imaging` (counts imaging-only infarcts)

Overall prevalence is reported three ways so the missing-outcome assumption is
explicit: complete-case, including imaging, and treating missing outcome as
no stroke.

## What is produced

| Sheet | Contents |
|---|---|
| `flow_overview` | rows eligible, rows with observed outcome, events |
| `overall_prevalence_*` | Wilson 95% CI under each missingness assumption |
| `prevalence_by_<var>` | subgroup prevalence + Wilson CI for each risk factor |
| `summary_<var>`, `continuous_summary_all` | median / mean / SD / missingness |
| `adjusted_ORs_continuous_per_unit`, `..._per_SD` | one model per continuous exposure |
| `adjusted_ORs_categorical_models` | one model per categorical exposure |
| `adjusted_ORs_logistic_fullmodel_*` | all covariates at once |
| `adjusted_ORs_logistic_reduced_coremodel_*` | age, BMI, HTN, DM only |
| `model_complete_case_counts` | n and events actually used by each model |
| `audit_nonmissing_counts` | missingness audit per variable |
| `analysis_dataset_deidentified` | analysis frame with MRN and DOB dropped |

## Read `model_complete_case_counts` before trusting any adjusted OR

The per-exposure and full models adjust for ~35 covariates and use
`drop_na()`, so a patient missing a single lab (ferritin, D-dimer, HbA1c…) is
dropped from that model entirely. With labs typically missing in a large share
of a chart-abstracted cohort, the complete-case n can fall to a small fraction
of the eligible cohort, and stroke is a rare outcome — the full model can
easily end up with too few events per variable to be interpretable.

`model_complete_case_counts` and `full_model_complete_case_summary` exist to
make that visible. If the n or event count there is small relative to the
cohort, prefer the reduced core model (age, BMI, HTN, DM — few missing values)
and treat the full model as exploratory. The unadjusted prevalence tables are
unaffected by this, since they only require the outcome.

## No data in this repo

`.gitignore` blocks `*.xlsx`, `*.csv`, and `data/`. Keep PHI on the network
share; only code is versioned here.
