# Stroke prevalence in non-cancerous uterine pathology

Cross-sectional analysis of the `Data Entry` sheet: overall and subgroup stroke
prevalence with Wilson 95% CIs, plus unadjusted and adjusted odds ratios.

## Run

```r
install.packages(c("readxl","janitor","lubridate","dplyr","tidyr","stringr",
                   "binom","broom","purrr","ggplot2","scales","forcats",
                   "openxlsx","rlang","brglm2","car"))
```

Set `input_file` at the top of `stroke_prevalence.R`, then:

```
Rscript stroke_prevalence.R
```

Output: `results/stroke_prevalence_results.xlsx` (one sheet per table) and
`results/plots/`.

Two sheets are the write-up-ready tables; everything else is supporting detail:

- **`TABLE_1_characteristics`** — cohort characteristics overall and by stroke
  status. Continuous as median (IQR) with Wilcoxon p; categorical as n (%) with
  chi-square / Fisher p.
- **`TABLE_2_main_results`** — one row per exposure level: N, strokes,
  prevalence (95% CI), unadjusted OR, adjusted OR, both p-values, and the model
  denominator. Reference levels are printed as `1.00 (reference)`, and a `Note`
  column flags any estimate resting on fewer than 10 events per variable.

## Analysis population

`eligible == 1` and no active malignancy — **39,807 patients**. Active
malignancy is an exclusion criterion (`exclusion_reason = 1`), not a covariate,
so it never enters a model.

## Headline result

Stroke prevalence (`stroke_any`): **1,538 / 39,807 = 3.86% (95% CI 3.68–4.06)**.
Including imaging-only infarcts: 1,543 / 39,807 = 3.88% (3.69–4.07).

## How the data departs from the Codebook

Three fields do not match their codebook definition. Each is handled explicitly
and the raw values are dumped to an audit sheet so you can check the decision.

| Field | Codebook says | Data actually contains | Handling |
|---|---|---|---|
| `fibroid_count` | coded 0=None, 1=Single, 2=2–4, 3=≥5, 9=Unknown | raw counts 0–8 mixed with `"2+"` and `"Multiple"` | rebuilt as ordinal None / Single / 2-4 / ≥5 / Multiple (unspecified); see `audit_fibroid_count_raw` |
| `race` | coded 1–9 | free text, ~20 levels (`White`, `African American`, `Asian Indian`, `Choose Not to Disclose`, …) | collapsed to standard groups for reporting; see `audit_race_raw`. Not used in any model |
| `migraine`, `surgical_tx`, `antithrombotic`, `hormonal_type` | — | multi-level, **not** binary | kept as labelled factors, not forced to 0/1 |

`9` means "unknown" for coded fields only. It is **never** treated as missing on
continuous labs, where 99 is an ordinary HDL, LDL, platelet or triglyceride
value.

## Adjustment sets

Chosen from observed missingness in this dataset, not from a wish list.

| Set | Covariates | Cost |
|---|---|---|
| `core` | age, BMI | BMI 9% missing |
| **`clinical`** (primary) | age, BMI, HTN, DM, dyslipidaemia, CAD, AF, CHF, smoking, migraine, VTE, thrombophilia, antithrombotic, bleeding, hormonal Tx, surgical Tx, uterine dx group | all complete except BMI |
| `clinical_no_treatment` | clinical minus hormonal Tx, surgical Tx, antithrombotic | temporality sensitivity |
| `full_with_labs` | clinical + Hgb, MCV, platelets, LDL | exploratory |

### Collinearity rule (do not relax this)

`uterine_dx_group` is a deterministic recoding of `fibroids`, `adenomyosis` and
`endometriosis` — knowing the three flags fixes the group exactly. Putting the
group in the covariate set while one of the flags is the exposure enters the
same information twice, inflating the variance without changing confounding
control.

`build_covs()` enforces the split, and the two branches are mutually exclusive:

- exposure ∈ {`fibroids`, `adenomyosis`, `endometriosis`, `uterine_dx_group`,
  `fibroid_count_cat`, `fibroid_max_cm`} → drop `uterine_dx_group` **and** all
  three flags from the covariates;
- any other exposure → keep `uterine_dx_group`, never add the flags.

It `stop()`s if an exposure ever appears in its own covariate vector. Max VIF
and the design condition number are written to `09_model_diagnostics` for every
fitted model, and `15_before_after_collinearity_fix` refits the old
specification next to the new one so the change is visible, not asserted.

### Temporality

44.9% of strokes (691/1,538) are timed *before* index_date, so `hormonal_tx`,
`surgical_tx` and `antithrombotic` — all ascertained at or after index — are
post-outcome for those patients. Two guards:

- `clinical_no_treatment` drops all three covariates;
- an incident-stroke population (`13_sensitivity_population`) removes the 691
  patients whose stroke predates index, leaving 847 strokes in 39,116 patients.
  Strokes with blank (189) or unknown (13) timing are **retained** — they are
  not positively known to precede index.

`TABLE_2_main_results` carries a `Sensitivity OR` column combining both guards.

### Sparse levels

Factor levels with fewer than 10 patients or fewer than 5 events are set to `NA`
before fitting and logged to `14_dropped_sparse_levels`. `fibroid_count_cat`
level "Unknown" (n = 3) previously produced OR 0.00 with an upper CI of 2.5e11 —
complete separation, not a finding — and is now coded `NA`, since an explicitly
unknown count carries no information about fibroid burden.

Labs with heavy missingness are deliberately **not** in the primary model:
ferritin 58%, HbA1c 57%, INR 69%, D-dimer 85%, fibroid_count 57%. Requiring them
simultaneously would leave a small, selected subsample. They are still analysed
as exposures under `core` adjustment.

Smoking is 73% "Unknown" — that is retained as an explicit factor level rather
than converted to `NA`, which would otherwise drop three quarters of the cohort
from every adjusted model.

`09_model_diagnostics` reports n, events and events-per-variable for every model;
the `stable` column flags anything below 10 events per variable.

## Sanity checks

Computed independently from the workbook — the R run should reproduce these:

| Subgroup | n | Prevalence | Crude OR |
|---|---|---|---|
| Fibroids yes vs no | 24,201 / 15,606 | 4.39% vs 3.05% | 1.46 |
| Adenomyosis yes vs no | 9,815 / 29,992 | 4.30% vs 3.72% | 1.16 |
| Endometriosis yes vs no | 21,090 / 18,717 | 3.32% vs 4.47% | 0.73 |
| Hypertension yes vs no | 5,942 / 33,865 | 10.05% vs 2.78% | 3.91 |
| Diabetes yes vs no | 6,992 / 32,815 | 8.17% vs 2.95% | 2.93 |
| Atrial fibrillation yes vs no | 511 / 39,296 | 26.03% vs 3.58% | 9.49 |

## Two internal inconsistencies, resolved

Both turned out to be reporting artefacts rather than data errors:

- **`prev_age_group` showed a "60-69" band (n = 433).** Eligible ages run
  exactly 18–60, and all 433 patients are aged exactly 60. The old bins
  (`right = TRUE`, break at 59) put age 60 into a band labelled "60-69",
  implying ineligible patients had leaked through the filter. They had not —
  only the label was wrong. Bands now stop at the eligibility ceiling.
- **`audit_out_of_range` said 360 ages set to NA, `12_missingness` said 0.**
  Different denominators: the audit counted all 53,060 file rows, missingness
  counted the 39,807 eligible. All 360 out-of-range ages (345 under 18, 15 over
  100) belong to *ineligible* rows, so age is genuinely 0% missing in the
  analysis population. The audit now reports both denominators.

## Interpretation caveats worth stating in the paper

- `stroke_timing` includes events **before** index. Prevalence is directionless,
  but the timing breakdown is in `04_strokes_by_stroke_timing` — reviewers will
  ask, and the "before index" share bears on any causal reading.
- The codebook flags `antithrombotic` as containing no antiplatelet agents
  (categories 1 and 3 are empty), `hormonal_type` as undercapturing combined OC,
  and `stroke_etiology` as imaging-derived rather than a true TOAST assignment.
- `stroke_date` is the earliest coded diagnosis date, not a chart-verified
  event date.

## No data in this repo

`.gitignore` blocks `*.xlsx`/`*.csv`. Keep the workbook outside version control.
