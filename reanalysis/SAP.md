# Statistical analysis plan — independent R reanalysis

**Status:** pre-specified on 2026-10-22, before any model in this R project was fitted. It was committed to git before Phases 4–7 ran, so later changes show up in the history. Any deviation is logged in `output/logs/deviations.md`.

**Decision recorded 2026-10-22 (PI):** the primary estimand is **incident stroke after the index date** (recommendation accepted). The original cross-sectional analysis is reproduced (Phase 4) and reported as secondary.

## 1. Data sources
These files are used as delivered and are never modified:
- `data/stroke_analysis_master.csv` (one row per woman)
- `data/stroke.xlsx` (Codebook, `stroke_outcome`, Hb/ferritin cleaning sheets)
- raw lab extracts (Lab Tests 2, 3, 4, 7, 8)
- diagnosis extracts 22–31
- encounter extracts 6, 8, 13
- medications administered 9

Records are linked by clinic number in memory only. Outputs carry `study_id` only. No MRN or date of birth is written anywhere.

## 2. Exposure
- **Hb source (primary):** laboratory whole-blood CBC hemoglobin from Lab Tests 2. Rules:
  - Included descriptions are listed explicitly in `R/02_labs.R`: in-house and external CBC hemoglobin.
  - Excluded: urine, qualitative, HbA1c, Hb fractions, reticulocyte Hb, plasma Hb, oxy/met/carboxy-Hb, blood-gas and co-oximetry Hb, point-of-care (POC) Hb, and comment rows.
  - Unit must be g/dL (any spelling) or blank.
  - Value must be 3–20 g/dL. The sentinel 9999999 is dropped, as are dates before 1990.
  - Several results on one day are reduced to that day's median.
- POC Hb is excluded from the primary definition because it is drawn mainly in acute care and is less accurate. Including it is sensitivity analysis S4.
- **Baseline Hb (revised primary):** the most recent Hb in [index − 1095 d, index]. Values measured after the index date are never used to define baseline exposure.
- **Grades:** anemia threshold Hb < 12 g/dL (WHO, non-pregnant women). Severity bands as originally specified (CTCAE-style; not WHO bands): mild 10.0–11.9, moderate 8.0–9.9, severe < 8.0. Reference: Hb ≥ 12.
- **Continuous Hb:**
  - restricted cubic spline, 4 knots at the 5th/35th/65th/95th percentiles of the analysis sample, reference 13 g/dL
  - piecewise-linear terms below and above 13 g/dL
  - grade score 0–3 for trend
- **Codebook `hgb`** (closest to index within ±3 years in either direction, POC included) is used only to reproduce the original analysis (Phase 4).

## 3. Outcomes
The source is `stroke_any`, `stroke_type`, `stroke_date` and `stroke_timing` (Codebook definitions, unchanged), plus `stroke_outcome.source_kinds`.
- **Primary (incident):** first stroke dated after the index date (`stroke_timing` = 3) with `stroke_type` ≠ 4 (TIA). This covers ischemic stroke, ICH, SAH, CVT and stroke of unknown type.
  - `stroke_type` is the most specific type ever coded, so a woman coded TIA had no stroke code and is followed without an event.
- **Secondary:**
  - incident ischemic stroke
  - incident stroke or TIA (the Codebook `stroke_any` definition)
  - incident TIA
  - individual hemorrhagic subtypes and CVT, analyzed only where there are at least 10 incident events; otherwise counts only
- **Original outcome** (`stroke_any` at any date, cross-sectional): Phase 4 only.

## 4. Cohort for the revised primary analysis
- Start from eligible women (eligible = 1; eligibility re-derived from age and malignancy and checked against the flag).
- Exclude:
  - women with `stroke_any` = 1 whose event is before index, at the same encounter, or of unknown timing (prevalent or unclassifiable)
  - women with no baseline Hb
  - women with no follow-up time
- **Time origin:** index date.
- **End of follow-up:** the earliest of stroke date, death, last contact (Encounters 8 Create Date, capped at 2026-09-27), or administrative end (2026-09-25).

## 5. Covariates (DAG-based; rationale in `output/dag.md`)
- **Primary adjustment set:**
  - age (natural spline, 3 df), race (4 groups), BMI (continuous; multiple imputation), smoking (never / ever / unknown)
  - hypertension, diabetes, dyslipidemia, CAD (dated on or before index), migraine (any), hormonal therapy type (sparse levels merged into "Other/multiple" when < 5 events)
  - uterine diagnosis group
  - conditions dated on or before index from the diagnosis extracts: CKD, chronic liver disease, alcohol use disorder, IBD, HIV, hemoglobinopathy, malabsorption, GI bleeding in the year before index, pregnancy in the year before index
- **Not in the primary model:**
  - uterine bleeding: a determinant of the exposure, undated, with no clear separate path to stroke (S6)
  - thrombophilia, atrial fibrillation, heart failure, VTE: undated, and often recorded after a stroke work-up (S7)
  - antithrombotics: medication window extends after index, and they mark hospitalization
  - Hb testing frequency: a proxy for healthcare contact and partly a consequence of anemia (S8)

## 6. Models
- **Primary:** Poisson GLM with a log person-years offset and robust (sandwich) SE, giving adjusted rate ratios (RR).
  - Confirmed with a Cox model (hazard ratios, HR; proportional hazards checked by Schoenfeld residuals). RR and HR are never relabeled as each other.
  - Results are pooled over 20 imputations (Rubin's rules; seed 20261022). The imputation model includes the outcome indicator and the Nelson–Aalen cumulative hazard.
- **Events per variable:** if events per model degree of freedom fall below 10, a pre-declared reduced set is used instead: age, race, BMI, hypertension, diabetes, dyslipidemia, smoking, CKD. This is reported whenever it happens.
- **Sparse cells:** categories with fewer than 5 events are reported as counts only. Sparse cross-sectional subtype models use Firth (bias-reduced) logistic regression via `brglm2`, labeled as secondary.

## 7. Sensitivity analyses (each addresses one named bias)

| ID | Analysis | Bias addressed |
|---|---|---|
| S1 | Complete-case for BMI | missing-data method |
| S2 | Follow-up starts at day 31 (early events excluded) | reverse causation / detection around index |
| S3 | Baseline Hb within 365 d before index | stale exposure |
| S4 | POC Hb included | exposure definition |
| S5 | Events limited to women with ≥ 1 acute stroke code | outcome misdating (history/sequela codes) |
| S6 | + uterine bleeding | confounding vs exposure determinant |
| S7 | + thrombophilia, AF, HF, VTE (undated) | over-adjustment check |
| S8 | + number of Hb tests in the 2 y before index | healthcare contact / surveillance |
| S9 | Excluding women with any pre-index anemia-causing condition | residual confounding by illness |
| S10 | Fine–Gray with death as competing event (subdistribution HR) | competing risk |
| S11 | Inverse-probability weights for having a baseline Hb | selection into the analyzed cohort |
| S12 | E-values | unmeasured confounding |

## 8. Time-updated analysis (secondary)
- **Cohort:** women stroke-free at index.
- **Person-time:** starts at index and is split at every CBC Hb date.
- **Exposure in each interval:** grade of the most recent Hb measured at least L days before the interval. A value is never effective on its own measurement day. Carry-forward window W:

  | Model | Window W | Lag L |
  |---|---|---|
  | Primary | 365 days | 0 |
  | Sensitivity | 365 days | 30 |
  | Sensitivity | 365 days | 90 |
  | Original reproduction | 1095 days | 0 |

- Person-time with no Hb inside the window is a separate "no recent Hb" level. It is never dropped.
- **Model:** Poisson with a person-time offset and SE clustered by woman, same covariates as the primary model.
- **Additional models:** adjustment for the number of Hb tests in the preceding 365 d; reproduction of the 12-month landmark analysis.

## 9. Reproduction of the original analysis (Phase 4)
- The original Python models are re-fitted in R on the same Codebook variables: logistic regression with HC1 SE; time-to-event from the later of index and Hb date; time-updated; 30-day lag; testing frequency.
- Each reported estimate is shown with the R estimate, the difference, and an explanation.

## 10. Multiplicity and reporting
- One primary estimand: grade categories plus trend.
- Secondary outcomes and sensitivity analyses are labeled exploratory. Benjamini–Hochberg q-values are reported across the secondary-outcome table.
- Estimates with 95% CIs are reported. Null and unfavorable results are kept. No causal language is used.
