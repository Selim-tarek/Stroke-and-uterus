# Stroke prevalence in non-cancerous uterine pathology

Cross-sectional analysis of the `Data Entry` sheet: overall and subgroup stroke
prevalence with Wilson 95% CIs, plus unadjusted and adjusted odds ratios.

## Run

```r
install.packages(c("readxl","janitor","lubridate","dplyr","tidyr","stringr",
                   "binom","broom","purrr","ggplot2","scales","forcats",
                   "openxlsx","rlang","brglm2","car"))
```

`janitor` and `binom` are optional — `compat_shims.R` supplies the two functions
used from them (`make_clean_names()`, Wilson `binom.confint()`) when either
package is unavailable.

Set `input_file` at the top of `stroke_prevalence.R`, or pass the path through
the environment, then:

```
STROKE_XLSX=/path/to/stroke.xlsx Rscript stroke_prevalence.R
```

Output: `results/stroke_prevalence_results.xlsx` (one sheet per table) and
`results/plots/`. The cleaned analysis frame is cached to
`results/prepared_data.rds`; `STROKE_PREP_ONLY=1` stops the run after that step,
which is useful when only the modelling code changed.

Four sheets are the write-up-ready tables; everything else is supporting detail:

- **`TABLE_1_characteristics`** — cohort characteristics overall and by stroke
  status. Continuous as median (IQR) with Wilcoxon p; categorical as n (%) with
  chi-square / Fisher p.
- **`TABLE_2_main_results`** — one row per exposure level: N, strokes,
  prevalence (95% CI), unadjusted OR, adjusted OR, both p-values, and the model
  denominator. Reference levels are printed as `1.00 (reference)`, and a `Note`
  column flags any estimate resting on fewer than 10 events per variable.
- **`TABLE_3_ischemic_multivariable`** — the dedicated ischaemic-stroke model
  (below): every covariate in one model, unadjusted and adjusted OR side by
  side.
- **`TABLE_3B_ischemic_with_hgb`** — the same model with haemoglobin added, on
  the patients who have one.

## Dedicated ischaemic-stroke model

`TABLE_2_main_results` fits one exposure at a time against `stroke_any`.
`TABLE_3_ischemic_multivariable` is a different object: a **single multivariable
logistic regression** in which every covariate is entered simultaneously and an
adjusted OR is reported for each, with **ischaemic stroke** as the outcome, in
women aged 18–60 with non-cancerous uterine pathology.

- **Outcome**: `stroke_any = 1` **and** `stroke_type = "Ischaemic stroke"`.
- **Competing outcomes**: patients whose only stroke was haemorrhagic, SAH, TIA
  or CVT are *not* controls — they are removed from the primary model. Strokes
  with unknown or blank type are removed for the same reason (they cannot be
  classified). Both counts are in `16_ischemic_flow`.
- **Population**: the analysis population with the 18–60 age band re-applied
  explicitly, so the model matches its stated definition regardless of the
  `eligible` flag.
- **Covariates**: age, BMI, the three uterine pathology flags, uterine bleeding,
  HTN, DM, dyslipidaemia, CAD, AF, CHF, smoking (Never/Ever/Unknown), migraine,
  VTE, thrombophilia, antithrombotic, hormonal Tx, surgical Tx.
  `uterine_dx_group` is excluded — it is a deterministic recoding of the three
  flags, and the collinearity rule below applies here too.
- **Sensitivity models** (`16_ischemic_all_models`): four-level smoking; Hgb
  added; the Hgb-complete population without Hgb; treatment covariates dropped;
  other stroke types kept as non-cases; incident strokes only.
- `16_ischemic_diagnostics` carries n, events, events-per-variable, C-statistic,
  AIC, max VIF and the design condition number for all seven models.

### Smoking collapsed to Never / Ever

`smoking_ever` merges Current and Former; the four-level variable is unchanged
everywhere else in the pipeline. Worth knowing before the collapsed estimate is
quoted: in the four-level fit the two merged levels point in **opposite**
directions — Current 0.84 (0.59–1.18), Former 1.39 (1.15–1.68), both vs Never —
so Ever 1.29 (1.07–1.55) averages a null and a raised estimate rather than
sharpening either. Discrimination is unchanged (C 0.799 → 0.798). "Unknown"
stays its own level, as it does throughout: it is 71% of the cohort.

### Haemoglobin

Hgb is **not** in the primary model, for a population reason rather than a
modelling one. It is 7.1% missing overall but only 0.6% missing among strokes:
requiring it drops ~1,660 stroke-free patients and 7 cases, so controls without
a recorded Hgb are largely those who never had bloods drawn. Adding it to the
primary model would change the control group and the covariate list at the same
time, leaving any shift uninterpretable.

`TABLE_3B_ischemic_with_hgb` therefore reports the Hgb model separately
(n = 34,214, 1,172 events, C = 0.795), and `16_ischemic_all_models` carries a
third fit — the same covariates on the same Hgb-complete patients but with Hgb
left out — so the effect of the variable can be separated from the effect of the
restriction. **Hgb: adjusted OR 0.94 per 1 g/dL (0.91–0.97), p < 0.001.** Every
other estimate is stable across the pair (fibroids 1.22 → 1.20; the population
restriction, not Hgb, accounts for most of that).

## Analysis population

`eligible == 1` and no active malignancy — **39,807 patients**. Active
malignancy is an exclusion criterion (`exclusion_reason = 1`), not a covariate,
so it never enters a model.

## Headline result

Stroke prevalence (`stroke_any`): **1,538 / 39,807 = 3.86% (95% CI 3.68–4.06)**.
Including imaging-only infarcts: 1,543 / 39,807 = 3.88% (3.69–4.07).

Ischaemic strokes: **1,189 of the 1,538 events**; 349 were another stroke type
and none had a missing type. The dedicated multivariable model runs on 35,878
patients with 1,179 ischaemic strokes (3.29%) after complete-case restriction
(BMI is the binding covariate), C-statistic 0.798, 45.3 events per variable, max
VIF 1.60.

Adjusted ORs from that model, largest first: CHF 2.26 (1.68–3.03), AF 2.23
(1.68–2.93), CAD 2.19 (1.68–2.84), thrombophilia 2.85 (2.20–3.65), migraine with
aura 3.17 (2.63–3.81), anticoagulant use 2.66 (2.28–3.09), VTE 2.07 (1.72–2.49),
hypertension 2.05 (1.78–2.36), dyslipidaemia 1.42 (1.24–1.64), **fibroids 1.22
(1.03–1.45, p = 0.020)**. Adenomyosis (0.97) and endometriosis (0.88) are null
after adjustment, and diabetes attenuates from a crude 3.02 to 1.13 (p = 0.107).

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
  Different denominators: the audit counted all 52,827 file rows, missingness
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
