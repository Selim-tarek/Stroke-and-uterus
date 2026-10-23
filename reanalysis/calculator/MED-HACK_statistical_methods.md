# MED-HACK score: statistical basis of the risk calculator

**What it predicts:** the 5-year risk of a first stroke in women aged 18–60 with uterine fibroids, adenomyosis or endometriosis and no previous stroke. It uses information available at the gynecologic diagnosis.
**Status:** development and internal validation only (TRIPOD type 1b/2a). Research use only; not externally validated.
**Calculator:** https://selim-tarek.github.io/Med-Hack-score/
**Code:** `reanalysis/R/13_anchor_design.R`, `15_risk_model.R`, `16_calculator_model.R`. The analysis plan was pre-specified (`reanalysis/SAP.md`, Addendum 2) and deviations are logged (`reanalysis/output/logs/deviations.md`).

---

## 1. Data and cohort

**Source:** electronic health records from one US academic medical center. Uterine diagnoses are coded from 2018 onward (the IRB-approved window). Laboratory, diagnosis, procedure and encounter records come from the full record.

**Time zero (anchor):** the earliest evidence of uterine disease. This is the first coded diagnosis, or an earlier dated myomectomy, uterine artery embolization, hysterectomy or endometrial ablation.
- 382 of 39,807 eligible women had such an earlier procedure, moving the anchor a median 3.7 years earlier.
- Procedure dates implying an age under 18 were ignored.

| Step | n |
|---|---|
| Eligible women (18–60 y, no active malignancy) | 39,807 |
| Excluded: stroke/TIA on or before the anchor, or undated | 983 |
| Excluded: no follow-up after the anchor | 491 |
| Stroke-free at the anchor | 38,332 |
| With a laboratory CBC hemoglobin in the 3 years up to the anchor | 24,886 |
| **Development sample** (also complete BMI, age ≥18 at anchor) | **24,006** |

**Outcome:** first stroke after the anchor (ischemic, intracerebral or subarachnoid hemorrhage, cerebral venous thrombosis, or unspecified; **TIA excluded**).
- Identified from diagnosis/problem-list codes and imaging reports, dated by the earliest code. Not adjudicated by chart review.
- **345 strokes.** Median follow-up 2.8 years.

**Censoring:** death, last recorded contact, or 25 September 2026.

---

## 2. Predictors

**Candidate predictors** were fixed in advance and are all known at the diagnosis visit:
- age, BMI, hypertension, diabetes, dyslipidemia, smoking
- coronary artery disease, migraine, chronic kidney disease, hemoglobinopathy
- endometriosis, adenomyosis, estrogen-containing hormonal therapy (combined oral contraceptive or menopausal HT)
- hemoglobin

**How they were coded:**
- **Hemoglobin:** a laboratory CBC value, entered as a "deficit" = max(13 − Hb, 0) g/dL. Each 1 g/dL below 13 counts; values above 13 count as 0.
- **Race:** deliberately excluded on equity grounds. Adding it did not improve the model (likelihood-ratio P = 0.88).
- **No data-driven selection.** No stepwise or p-value-based selection was used.

**Dropped after the points were set.** Five predictors received **0 points**: BMI, dyslipidemia, smoking, endometriosis and adenomyosis. They were then removed from the calculator, and the exact-risk model was refitted with the nine MED-HACK items (§4). This was a PI decision; the full 14-predictor model is reported for comparison.

**Smoking:**
- Status was undocumented for 65% of women.
- Only 3.9% were recorded as current smokers.
- Some statuses came from notes written after the stroke.
- Under every coding tried (ever/never; never/current/former; unknown kept separate or set to 0), the association was null or implausible.

So smoking **could not be assessed**. That is not evidence of no effect, and risk in smokers is probably underestimated.

---

## 3. Model

**Cox proportional-hazards regression** for time to first stroke.
- **Proportional hazards:** global Schoenfeld test P = 0.95.
- **Age linearity:** non-linear term P = 0.13, so age is linear.

**5-year risk:**

> risk = 1 − S₀(5) ^ exp( s × Σ βᵢ (xᵢ − cᵢ) )

- **S₀(5):** baseline 5-year stroke-free survival at the centring values.
- **cᵢ:** centring values. These are the cohort means for age and the hemoglobin deficit, and 0 for every yes/no item (the convention the survival software uses).
- **s:** uniform shrinkage factor = the mean bootstrap calibration slope (200 resamples). It corrects for overfitting.

---

## 4. The calculator's exact-risk model (nine MED-HACK items)

S₀(5) = 0.99028; shrinkage s = 0.980; centring: age 41.113 y, Hb deficit 0.671 g/dL.

| Letter | Predictor | β | Hazard ratio (95% CI) |
|---|---|---|---|
| M | Migraine | 0.6759 | 1.97 (1.59–2.43) |
| E | Estrogen-containing therapy | 0.6926 | 2.00 (1.17–3.42) |
| D | Diabetes | 0.2019 | 1.22 (0.97–1.54) |
| H | Hypertension | 0.7628 | 2.14 (1.70–2.71) |
| H | Hemoglobinopathy (sickle cell disease/trait, thalassemia) | 1.0627 | 2.89 (1.68–4.98) |
| A | Anemia: per 1 g/dL hemoglobin below 13 | 0.1328 | 1.14 (1.06–1.23) |
| A | Age, per year | −0.0107 | 0.99 (0.98–1.00) |
| C | Coronary artery disease | 1.4302 | 4.18 (2.92–5.98) |
| K | Chronic kidney disease | 0.6803 | 1.97 (1.35–2.88) |

**About age:** in this cohort, stroke rates were highest in women under 30 (6.8 per 1,000 person-years) and flat from 30 to 60. This held whether age was taken at diagnosis or at the time of stroke (an attained-age time-scale analysis, which left every other hazard ratio unchanged). The age effect reflects this cohort, most likely testing and referral patterns, and should not be read as a protective effect of age.

---

## 5. The points score (CHA₂DS₂-VASc style)

**How points were assigned (Sullivan method).**
- Points come from the shrunken coefficients of the pre-specified 14-predictor model.
- The plan set 1 point = the effect of 5 years of age. Age had no independent effect, so instead **1 point = half the hypertension coefficient**, making hypertension 2 points as in CHA₂DS₂-VASc. This is a logged deviation.
- **Age and hemoglobin** used banded values: age at the band midpoint; Hb bands ≥12 / 10–11.9 / 8–9.9 / <8 g/dL.
- **Age display:** the age effect produced "age ≥45 = −1 point". It is presented as "age <45 = +1" instead. Every total shifts by one, and rankings and risks are unchanged.

| Item | Points |
|---|---|
| **M** Migraine | 2 |
| **E** Estrogen-containing hormonal therapy | 2 |
| **D** Diabetes | 1 |
| **H** Hypertension | 2 |
| **H** Hemoglobinopathy | 3 |
| **A** Anemia: Hb 8.0–11.9 g/dL / <8.0 g/dL | 1 / 2 |
| **A** Age under 45 | 1 |
| **C** Coronary artery disease | 4 |
| **K** Chronic kidney disease | 2 |

**Range** 0–19. Each point corresponds to a hazard ratio of 1.40 (about 40% higher stroke rate).

**Points to risk:** a Cox model with the total score as its only predictor gives

> risk = 1 − 0.98980 ^ exp(0.33504 × (points − 1))

| Points | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5-year risk, % | 0.7 | 1.0 | 1.4 | 2.0 | 2.8 | 3.8 | 5.3 | 7.4 | 10.1 | 13.9 | 18.9 | 25.3 | 33.5 | 43.5 |
| Women in cohort | 3,195 | 7,439 | 4,129 | 4,622 | 2,270 | 1,123 | 642 | 285 | 144 | 66 | 56 | 23 | 9 | 3 |
| Strokes observed | 18 | 44 | 35 | 76 | 56 | 42 | 18 | 18 | 15 | 11 | 5 | 4 | 2 | 1 |

Risks above about 9 points rest on very few women and should be read as imprecise.

---

## 6. Performance

### Discrimination (C-statistic; 0.5 = chance, 1.0 = perfect)

| Model | Apparent | Optimism-corrected (200 bootstraps) | Temporal validation (3 y) |
|---|---|---|---|
| **9-item MED-HACK exact model (calculator)** | 0.705 | **0.701** | **0.695** |
| Points score | 0.697 | — | — |
| Full 14-predictor model | 0.706 | 0.696 | 0.665 |

**Temporal validation:** the model was developed on women with an anchor before 2023 (11,396 women; 226 strokes) and tested on women with an anchor in 2023–2026 (12,610 women; 119 strokes), at a 3-year horizon.

### Calibration

**Development data, 5 years (9-item model), predicted vs observed risk by decile:**

| Decile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Predicted % | 0.80 | 0.89 | 0.95 | 1.03 | 1.14 | 1.48 | 1.79 | 2.03 | 2.52 | 6.02 |
| Observed % (Kaplan–Meier) | 0.93 | 0.77 | 0.48 | 0.88 | 1.19 | 1.45 | 1.82 | 1.50 | 3.91 | 6.03 |

Calibration slope 0.98.

**Points score, 5 years:** predicted vs observed was 0.93 vs 0.89%, 1.42 vs 1.20%, 1.98 vs 1.94%, 2.76 vs 3.31% and 6.07 vs 6.16% across score groups.

**Temporal validation:**
- Calibration slope 0.90.
- **Observed/expected 3-year risk = 1.56.** In women diagnosed from 2023 on, actual stroke risk was about 1.5 times the predicted risk.
- The calculator therefore probably **underestimates** absolute risk in newer patients and other settings. Recalibration (updating S₀) is recommended before use elsewhere.

### Decision-curve analysis (full model, 5 years)
Net benefit was above "treat all" and "treat none" at every risk threshold from 1% to 10%. The absolute net benefit is small, because 5-year stroke risk is low for most women (cohort average 2.0%; median predicted risk about 1.2%).

---

## 7. Limitations

1. **Single center; no external validation.** This is the most important next step.
2. **Outcome:** strokes are coded, not chart-validated, and dated by the earliest code.
3. **Follow-up:** median follow-up is 2.8 years, so the 5-year risks partly extrapolate.
4. **Selection:** only women with a laboratory hemoglobin (65% of stroke-free women) could be included. Tested women had more recorded comorbidity.
5. **Missing risk factors:** smoking could not be assessed. Oral contraceptives were under-recorded (administered-medication source).
6. **Age:** the age effect is cohort-specific (see §4).
7. **Hemoglobin is a marker:** in the wider analyses, its association with stroke behaved like a marker of general illness. The score does not imply that treating anemia lowers stroke risk.

---

## 8. Reproducibility

Everything is generated by code from the source data: `Rscript run_all.R` in `reanalysis/` (seed 20261022).

| File | Content |
|---|---|
| `output/calculator_model_spec.json` | Calculator exact model (S₀, shrinkage, β, centring) |
| `output/risk_model_spec.json` | Full model and points |
| `output/tables/R1–R9` | Coefficients, validation, calibration, decision curve, points, score-to-risk |

The calculator's arithmetic was checked against R for test profiles and matches exactly.
