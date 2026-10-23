# MED-HACK stroke risk model — statistical audit

**Date:** 2026-10-23
**Scope:** an audit of the completed MED-HACK development (`R/13`, `R/15`, `R/16`, calculator). It is not a new model; nothing under `output/`, `derived/` or `calculator/` was changed.
**Inputs:** the de-identified patient-level data were available (`../data`, `derived/*.rds`). Every number below comes from the audit scripts in `audit/R/` (run with `Rscript audit/run_audit.R`). Each script writes one Markdown file of tables to `audit/output/` (CSV copies are git-ignored).

**Status codes:**
- **VERIFIED:** supported directly by existing code, data or outputs.
- **REPRODUCED:** independently recalculated in this audit and matched.
- **NEW:** a new analysis run during this audit, on the existing data.
- **NEW ANALYSIS REQUIRED:** the data exist but the work is not done.
- **DATA REQUIRED:** the information is not in the project.
- **UNRESOLVED:** the evidence is insufficient or contradictory.

---

## Summary of what the audit found

**What holds up**
- The arithmetic reproduces exactly:
  - The coefficients, baseline survival and shrinkage reproduce to within 1e-10.
  - The calculator's JavaScript and the public web page reproduce R to display rounding.
  - The points and the score-to-risk table reproduce exactly.
- The cohort flow reconciles to the woman.
- In the development data the calculator is well calibrated: 5-year O/E 1.02 (0.91–1.14), calibration slope 1.03, mean absolute calibration error 0.16 percentage points.
- Discrimination is moderate and fairly stable: C about 0.70.
- The hazard ratios barely move when the outcome is restricted to acute-coded or ischemic events, when death is treated as a competing event, or when Hb and BMI are imputed.

**What could invalidate the calculator (Priority 1)**
1. **Several predictors are not measured at time zero.**
   - Estrogen exposure is defined from 365 days before the index date *onward*, with no end date.
   - Hemoglobinopathy is taken from codes of *any date*: 67 of 263 positive women were first coded after the anchor, including 7 of the 14 positive women with a stroke.
   - Smoking uses the earliest statement *after* the index date when none exists before it. For 43 of the 180 stroke cases with a documented status, the statement is dated on or after the stroke.
   - Migraine and diabetes have no dated source in the project. An indirect check suggests part of their recording happened after baseline.
   - A woman using the calculator at a clinic visit cannot supply future information, so some predictor effects may be overstated.
2. **The outcome is code-based and has not been validated.**
   - Ten of the 345 events (2.9%) are imaging-only, nine of them chronic or old infarcts dated by the scan.
   - 42 events (12.2%) fall within 30 days of the anchor.
   - For 33% of events, brain imaging was reviewed and showed no stroke.
   - No positive predictive value is known.
3. **Risk is underpredicted in recent patients, and not just through the baseline.**
   - In women anchored in 2023 or later, observed 3-year risk was 1.44 (1.19–1.73) times the expected risk.
   - The ratio was similar for each anchor year: 2023, 2024, and 2025–26.
   - It was concentrated in women with at least one risk factor (1.57, versus 1.09 with none).
   - A test of full coefficient revision was significant (P < 0.001). Updating the baseline alone would not fix it.
4. **The 5-year horizon has never been tested on later patients.** Maximum follow-up in the 2023+ set is 3.72 years, and only 28% of the development sample was observed for 5 years.
5. **Two predictors depend on a handful of events.** The estrogen association loses statistical significance after removing 2 influential strokes, hemoglobinopathy after removing 5.
6. **The calculator presents two different risks.** The headline points risk and the exact-model risk differ, for example 3.8% vs 5.3% for the worked example. Combinations never observed extrapolate to risks of up to 99.99%.

**Conclusion (Deliverable D):** the evidence supports **research use only**. It does not yet support prospective clinical evaluation, and not routine use.

---

## Part 0 — What was inspected

| Item | Files | Notes |
|---|---|---|
| SAP and addenda | `SAP.md` (main, Addendum 1, Addendum 2) | Addendum 2 (risk calculator) was committed at 03:52 on 2026-10-23, **6 minutes after** the E1 prediction results on the same cohort (`git log`: f298953 at 03:46). E1 already showed C 0.706–0.717 for a similar predictor set, so Addendum 2 is not fully blind. |
| Deviations log | `output/logs/deviations.md` | Two calculator deviations are logged: the points constant, and the removal of zero-point predictors (PI decision). The second cites temporal C values as justification, so the decision was made after seeing the temporal results. |
| Scripts | `R/00`, `R/01`, `R/03`, `R/04`, `R/13`, `R/15`, `R/16`, `calculator/bud_stroke_calculator.html`, `calculator/standalone/index.html` | Read in full. |
| Outputs | `output/tables/A1`, `A3`, `R1`–`R9`, `SM1`, `SM1a`; `output/calculator_model_spec.json`; `output/risk_model_spec.json`; `output/logs/model_warnings.md` | All traced to their generating lines. |
| Data | `stroke_analysis_master.csv`, `stroke.xlsx` (Codebook, stroke_outcome, smoking_source, hormonal_meds, migraine_terms, cad_terms, dyslipidemia_update), `dx/` extracts 22–31, `enc/` 6, 8 and 13, `meds/` file 9 | Analysis-ready derived files present. |

**Completed before this audit:** anchored cohort (`R/13`); 14-predictor model with 200-resample bootstrap, temporal validation, decision curve and points (`R/15`); 9-item model with bootstrap and temporal validation (`R/16`); calculator.

**Attempted and failed or abandoned:**
- The pre-specified points constant (5 years of age) was unusable because the age coefficient was negative. This is a logged deviation.
- Earlier E1 and E4 errors were fixed.
- The E1–E5 tables predate the anchor-age fix.

**Never run before this audit:**
- outcome validation
- confidence intervals for temporal performance
- flexible calibration and Brier score
- decision curve for the actual calculator and the points
- bootstrap of the full development process
- competing-risk comparison for the calculator
- subgroup performance
- multiple imputation for the calculator model
- conventional-risk-factor comparator
- external validation

**Files missing (DATA REQUIRED):**
- the raw stroke, hypertension, diabetes, dyslipidemia, CAD and migraine diagnosis extracts with dates (`dx/` holds only files 22–31: hemoglobinopathy, CKD, liver, alcohol, GI bleeding, HIV, menopause, anemia, IBD, pregnancy)
- the administered-medication extract behind `hormonal_type` (`meds/` holds file 9 only, which contains iron and spironolactone)
- BMI measurement dates
- clinical notes or charts for outcome adjudication

---

## Deliverable A — Answers by question

### Part I — Cohort and design

**Q3. Cohort reconstruction** — `audit/output/A01_cohort_design.md` (A01a)

- **Finding:** every count reconciles.

  | Step | n |
  |---|---|
  | Eligible women | 39,807 |
  | Excluded: stroke, TIA or undated event on or before the anchor | 983 (682 stroke, 99 TIA, 202 undated) |
  | Excluded: hormonal type unknown (not shown in A3) | 1 |
  | Excluded: no follow-up after the anchor | 491 (437 last contact on the anchor date, 53 before it, 1 death) |
  | Stroke-free at the anchor | 38,332 |
  | Excluded: no CBC Hb in the 3 years before the anchor | 13,446 |
  | With CBC Hb | 24,886 |
  | Excluded: BMI missing | 880 |
  | Excluded: age under 18 | 0 |
  | Development sample | **24,006**, with **345** strokes |

  An independent re-derivation of the anchor matched for all 38,332 women. `R/15` lacks the age ≥18 filter, but its sample is identical.
- **Status:** REPRODUCED. The one undisclosed exclusion is a single woman with hormonal type 9.
- **Implication:** none for validity; the flow table should list the hormonal-type exclusion.
- **Action:** add that line to the flow diagram.

**Q3.1–3.3. Target population** — A01b, A01c

- **Finding:** this is not a fibroid cohort.
  - The development sample has fibroids in 15,113 women (63%), endometriosis in 11,911 (50%) and adenomyosis in 6,558 (27%).
  - 9,019 have fibroids only, 6,041 endometriosis only, 1,166 adenomyosis only, and 7,780 more than one condition.
  - Stroke rates range from 3.57 to 5.23 per 1,000 person-years across diagnosis combinations.
  - The adenomyosis code exists only from October 2022 (Codebook), so earlier adenomyosis is recorded as fibroids or missed.
- **Status:** VERIFIED.
- **Implication:** the intended population ("benign uterine disease") matches the data. A fibroid-specific claim does not.
- **Action:** describe the population as women with fibroids, adenomyosis or endometriosis. Report performance by diagnosis (Q18).

**Q3.4–3.5. Predictors measured before time zero?** — A02a–A02d

- **Finding:** no, for several predictors (see the ascertainment-window table A02a).

  | Predictor | Measured relative to the anchor | Evidence |
  |---|---|---|
  | Hb | at or before | — |
  | Age | at | — |
  | CKD | at or before | — |
  | Hypertension, CAD, dyslipidemia | at or before the **index** date | For the 336 procedure-anchored women this is a median of about 3.7 years after the anchor (A1). |
  | Estrogen | from index − 365 days **onward** | Codebook (`hormonal_tx`). |
  | Hemoglobinopathy | **any date** | `R/13_anchor_design.R:65`. |
  | BMI | within ±12 months of the index date | Codebook. |
  | Smoking | latest on or before the index date, otherwise **earliest after** | Codebook. |
  | Migraine, diabetes | timing not documented | No dated source in the project. |

  Quantified:
  - **Hemoglobinopathy:** 67 of 263 first coded after the anchor; 7 of 14 positive strokes; 3 coded within 30 days before the stroke or later. Restricting to pre-anchor codes changes the HR from 2.89 (1.68–4.98) to 2.15 (1.01–4.56), based on 7 strokes.
  - **Smoking:** of 8,391 documented statuses, 3,268 come from statements after the anchor. For 43 of the 180 stroke cases with a documented status, the statement is on or after the stroke date.
  - **Indirect look-ahead check (A02d):** among women without a stroke, within the same anchor year, more observed follow-up goes with higher prevalence. The odds ratio per doubling of follow-up is 1.56 for migraine and 1.72 for diabetes. For predictors dated at baseline it is 1.31 (hypertension), 1.20 (CAD), 1.30 (dyslipidemia) and 1.00 (CKD). Diabetes prevalence rises from 15.9% to 27.3% between women followed under 1 year and those followed 3 years or more.
  - Women who stayed in care were much more likely to carry migraine or diabetes codes: OR for loss to follow-up 0.39 and 0.29 (A04f).
- **Status:** VERIFIED for the definitions; NEW for the quantification; **DATA REQUIRED** to date migraine, diabetes, hypertension, CAD and estrogen exactly.
- **Implication:** this is the most serious threat to validity.
  - Predictors recorded after baseline, and especially after a stroke work-up, can inflate associations.
  - Hemoglobinopathy and anemia work-ups are common after a stroke.
  - The look-ahead pattern also means women with longer follow-up had more opportunity to be classified as exposed.
  - None of this can occur when the calculator is used at a visit, so performance at the point of use may be lower than reported.
- **Action:** obtain dated extracts for every predictor. Re-derive each one strictly at or before the anchor. Refit and re-validate. This is Priority 1.

**Q3.6. Outcome-free eligibility and follow-up alignment** — A01a, A01f

- **Finding:**
  - Strokes and TIAs on or before the anchor, or undated, are excluded: 983 women.
  - No event falls after the end of follow-up. No stroke falls after the last encounter or after death. No death falls before the anchor.
  - A TIA ends stroke follow-up: 85 development women were censored at a TIA, and any later stroke in them is not counted.
  - There is one record per woman, so the "first stroke" is the earliest code of any kind.
- **Status:** REPRODUCED.
- **Implication:**
  - Censoring at TIA treats a high-risk group as uninformative and may lower observed risk slightly.
  - The earliest-code rule can date a remote stroke recorded first as a history code after the anchor (Q5).
- **Action:** in a revision, keep women under follow-up after a TIA (TIA as a non-terminal event). This needs event-level dates (DATA REQUIRED).

**Q3.7. Bias from procedure dates** — A01d

- **Finding:**
  - 336 development women (1.4%) are anchored at an earlier procedure, contributing 1,641 person-years before the first coded diagnosis and 5 strokes in that window.
  - They had to be alive, aged 18–60, malignancy-free and in the system at a 2018+ diagnosis to be included, so that window is partly immortal for death.
  - Their covariates are measured at the index date, after the anchor.
  - The extract window shows no pile-up of first diagnoses in early 2018: 614 in January–March versus a monthly median of 205.
- **Status:** NEW.
- **Implication:** small. It affects 1.4% of women and 5 events.
- **Action:** sensitivity analysis excluding procedure-anchored women, or anchoring them at the index date.

**Q3.8. Censoring rules** — A01e, A01f

- **Finding:** follow-up ends at the earliest of stroke, TIA, death, last encounter (`enc` file 8 create date, capped at 2026-09-27) or 2026-09-25.

  | End of follow-up | Women |
  |---|---|
  | Last encounter | 20,996 |
  | Administrative end | 2,419 |
  | Stroke | 345 |
  | Death | 161 |
  | TIA | 85 |

  2,028 women had no encounter in the last year, so they may have been lost to follow-up.
- **Status:** REPRODUCED.
- **Implication:** deaths outside the institution are probably missed, since only 1,471 deaths are recorded across the whole extract. Censoring depends on the predictors (Q6).
- **Action:** link to death records if possible (DATA REQUIRED).

### Q4. Selection and missing data — A02e–A02j, A11

- **Exclusion reasons:** see Q3. The main exclusion is no CBC Hb in the 3 years before the anchor: 13,446 women (35.1% of the anchored cohort).
- **Missingness (anchored cohort, A02g):**
  - Hb 35.1%, BMI 8.4%, smoking unknown 73.5%; 65.0% in the development sample.
  - Migraine code 9 (aura not stated, counted as migraine): 9.9%.
  - Hypertension, diabetes, dyslipidemia, CAD, CKD, hemoglobinopathy, endometriosis and adenomyosis are coded 0 when absent, so true missingness is undetectable.
- **Differences between included and excluded women (A02h, A02i):**
  - Large differences: standardized difference 0.59 for dyslipidemia, 0.54 for smoking documented, 0.45 for diabetes, 0.43 for migraine, 0.31 for hypertension, 0.27 for follow-up, and −0.23 for endometriosis.
  - Having a CBC was strongly predicted by recorded comorbidity: OR 3.98 for CKD, 2.50 for dyslipidemia, 2.31 for migraine.
  - Crude stroke rate was 4.18 (3.75–4.64) per 1,000 person-years in the included women versus 2.93 (2.40–3.53) in the excluded.
  - **After adjustment, the stroke hazard did not differ: HR 1.00 (0.80–1.26).**
  - Exploratory check: applying the model to excluded women with an assumed Hb deficit of 0 gives O/E 0.94 at 5 years and 1.14 at 3 years.
- **Multiple imputation (A11a):** imputing Hb and BMI in all 38,332 women (m = 10, predictive mean matching) gave similar HRs. The largest change was migraine, from 1.97 to 2.18; hb_def went from 1.14 to 1.12, with a fraction of missing information of 0.30. Two caveats:
  - Missing-at-random is doubtful: not being tested is informative, and predictive mean matching can only borrow values from tested women.
  - mice dropped the Nelson–Aalen term as collinear.

  The design for imputation inside resampling is written out in A11, but it was **not executed**.
- **Smoking (A02e, A02f):** the 65.0% unknown and 3.9% current rates are **REPRODUCED**. Adding smoking to the 9-item model:

  | Coding | HR (95% CI) |
  |---|---|
  | Ever vs never or unknown | 1.15 (0.91–1.45) |
  | Current vs never | 0.59 (0.33–1.06) |
  | Former vs never | 1.03 (0.76–1.40) |
  | Unknown vs never | 0.75 (0.57–0.98) |
  | Known status only, ever vs never | 0.94 (0.70–1.27) |

  The "unknown" HR below 1 shows that documentation itself tracks contact with care.
- **Status:** NEW for the analyses; VERIFIED for smoking.
- **Implication:**
  - Complete-case selection makes the model a model **for women who had a CBC**. Coefficients look robust to that selection.
  - The smoking effect **cannot be estimated** from these data. That is not evidence of no effect, and risk in smokers is probably underestimated.
  - Effects of current smoking, oral contraceptives (under-captured) and any variable coded 0 when absent are not reliably estimable.
- **Action:**
  - State the CBC requirement in the intended-use statement.
  - Run imputation inside resampling only if validation in all women is needed.
  - Do not impute smoking: the missingness is clearly not at random.

### Part II — Outcome and follow-up

**Q5. Outcome validity** — `audit/output/A03_outcome.md`

- **Algorithm (A03a):**
  - Problem list and diagnosis codes from any date, plus positive brain-imaging reports.
  - Imaging-only cases are included, even chronic, incidental or hedged ones.
  - The event type is the most specific code; the date is the earliest code. TIA is excluded from the outcome but censors follow-up.
  - There is one record per woman, so recurrent strokes cannot be counted and duplicates are collapsed by design.
  - The code lists are in the workbook's stroke_outcome sheet, but the raw stroke extracts are not in the project.
- **Event provenance, 345 strokes (A03b):**

  | Category | Events |
  |---|---|
  | Ischemic | 301 |
  | Intracerebral hemorrhage | 24 |
  | Subarachnoid hemorrhage | 11 |
  | Unspecified | 7 |
  | Cerebral venous thrombosis | 2 |
  | Has an acute-event code | 328 (95.1%) |
  | History or sequela codes only | 7 |
  | Imaging-only (9 chronic or remote) | 10 |
  | Imaging reviewed and positive | 147 (42.6%) |
  | Imaging reviewed and negative | 115 (33.3%) |
  | No imaging | 83 |
  | Within 30 days of the anchor | 42 (12.2%) |
  | 31–90 days after the anchor | 26 |
  | More than one event type coded | 23 |

  - Imaging only supports a code-based diagnosis, except for the 10 imaging-only cases.
  - Events within 30 days of the anchor may be strokes during which the uterine diagnosis was coded, or misdated earlier strokes.
- **Robustness (A03d):** the HRs barely change when the outcome requires an acute code (328 events), an acute code more than 30 days after the anchor (287), or ischemic stroke only (255). C ranges from 0.701 to 0.712.
- **Subtypes (A03c):** only ischemic stroke supports a separate model: 301 events, 33 per parameter. ICH has 24, SAH 11, CVT 2. Subtype-specific models are **not feasible** and were not fitted.
- **Chart validation:**
  - Not possible, because there are no charts or notes in the project.
  - A sampling frame has been prepared in `derived/audit/chart_review_sample_events.csv` (git-ignored). It holds all 58 non-acute, imaging-only or within-30-days events plus 42 randomly sampled events, weighted for estimating PPV. The precision of a PPV near 0.85 with 100 reviews is about ±0.07.
  - Sensitivity would need thousands of non-event reviews.
- **Status:** VERIFIED for the algorithm; NEW for provenance and robustness; **DATA REQUIRED** for PPV.
- **Implication:**
  - What the code-based outcome supports: a model for *coded* stroke that is insensitive to the main misclassification scenarios that could be tested.
  - What it does not support: an absolute risk of clinically confirmed stroke. The share of false positives is unknown, and chronic imaging-only infarcts are not incident clinical events.
- **Action:** run the 100-event blinded chart review (Priority 1). Exclude imaging-only chronic infarcts from the primary outcome.

**Q6. Follow-up and competing risks** — `audit/output/A04_followup_competing.md`

- **Follow-up (A04a):**
  - Median observed follow-up is **2.95 years** (IQR 1.25–5.33); reverse-KM potential follow-up 3.01 years; maximum 25.41 years (procedure anchors); 82,607 person-years.
  - **Correction:** the methods summary's "2.8 years" was the earlier index-date cohort, not the development sample.
- **By horizon (A04b):**

  | Horizon | Still observed | Strokes | Deaths | Censored without event | KM risk |
  |---|---|---|---|---|---|
  | 1 year | 19,129 (79.7%) | 140 | 58 | 4,679 | 0.63% |
  | 3 years | 11,866 (49.4%) | 247 | 100 | 11,793 | 1.34% |
  | 5 years | 6,753 (28.1%) | 306 | 130 | 16,817 | 1.99% (1.75–2.23) |

  The 5-year risk rests on the 28% observed that long. It is a Kaplan–Meier extrapolation that assumes censoring is non-informative.
- **Competing death (A04c–A04e):**
  - 161 deaths without stroke.
  - Aalen–Johansen 5-year stroke incidence is 1.984%, versus 1.992% by Kaplan–Meier.
  - Fine–Gray subdistribution HRs equal the cause-specific HRs to two decimals.
  - Individual 5-year risks differ by a median ratio of 1.002; the 99th percentile of the difference is 0.23 percentage points.
- **Censoring (A04f):** loss to follow-up depends strongly on the predictors (OR 0.39 for migraine, 0.29 for diabetes, 1.14 per g/dL of Hb deficit). This is compatible with non-informative censoring *given* the predictors, but cannot test it.
- **Status:** NEW.
- **Implication:**
  - The published model is a cause-specific Cox model whose risk ignores competing death. With so few deaths this makes no material difference here.
  - The real limitation is horizon extrapolation.
- **Action:**
  - Report 3-year risk alongside 5-year risk, or as the primary output, until longer follow-up exists.
  - Keep the cause-specific model, but state the estimand.

### Part III — Model specification

**Q7. Reproduction** — `audit/output/A05_reproduce_calculator.md`

- **Finding:**
  - Refitting with the original seed reproduces the stored values to within 4e-11: S0(5) 0.99027729, shrinkage 0.97961740, centring values age 41.112745 and Hb deficit 0.67095101, and all nine β.
  - Apparent C is 0.7054, optimism-corrected 0.7008.
  - Fourteen test profiles were compared across sources:

    | Comparison | Maximum discrepancy |
    |---|---|
    | Independent R vs JSON spec | 9e-9 percentage points |
    | Spec vs calculator JavaScript | 5e-6 percentage points |
    | JavaScript vs public page | 0 at display precision |
    | Points, R vs JavaScript vs local page vs public page | 0 |

  - The public page was byte-identical to `calculator/standalone/index.html` on 2026-10-23. It was fetched with curl and driven in headless Chromium.
  - Unit conversion (g/L, mmol/L) is correct.
- **Coefficient convention:**
  - The spec, R8 and the methods summary report **unshrunk** coefficients and HRs. The calculator applies shrinkage to the whole linear predictor.
  - The baseline survival is taken at coxph's centring values, which are 0 for binary predictors, matching the calculator.
- **Input handling:**
  - Every item is required, so unknown values are not silently treated as absent. But there is no "unknown" option either, so a user must answer "No" when a factor is unknown.
  - Age 18–60 and Hb 3–20 g/dL are enforced.
  - Risk is not capped: an all-positive profile gives 99.99%, and the page only warns when risk exceeds the 99th percentile (9.78%).
- **Two risks shown:** the headline risk comes from the points table (a Cox model on the 14-predictor points total); the exact risk comes from the nine-item model. For the worked example they are 3.84% and 5.30%.
- **Status:** REPRODUCED.
- **Action:**
  - Show one risk.
  - Grey out estimates beyond the observed range: published points above 13, or exact risk above roughly 20%.
  - Add "unknown" handling, or a statement that unknown means not documented.

**Q8. Assumptions** — `audit/output/A06_diagnostics_age.md`

- **Proportional hazards (A06a):** no single-predictor test is below P 0.09 (age is lowest at 0.092); global P 0.85; plots are in `figures/A06a_schoenfeld.png`.
- **Functional form (A06b, A06c):**
  - Age: spline versus linear P 0.087. The spline shows a raised hazard below about 30 (HR 1.82 at age 20 versus 41) and a flat curve from 30 to 60.
  - Hb: the hinge at 13 had a lower AIC than both a 4-knot spline and linear Hb. The spline gives HR 1.17 at 12 g/dL and 2.14 at 7 g/dL.
  - BMI in the 14-predictor model: linear is adequate (P 0.62).
- **Influence (A06d, A06k):**
  - Refitting without the 20 most influential women moves estrogen from 2.00 to 0.15, and hemoglobinopathy from 2.89 to 1.71.
  - Fragility: the 95% CI includes 1 after removing **2** estrogen-positive strokes, 5 hemoglobinopathy-positive strokes, 8 CKD-positive strokes or 22 CAD-positive strokes.
- **Collinearity:** all VIFs are 1.0–1.1.
- **Sparse predictors:**
  - Estrogen: 338 women, 14 strokes.
  - Hemoglobinopathy: 263 women, 14 strokes.
  - CAD: 378 women, 37 strokes.
  - Rare combinations hold 0–10 strokes (A06f).
- **Interaction:** migraine × estrogen, the clinically motivated check, gives P 0.041 and HR 3.44 (0.94–12.57), from only 11 strokes among 130 women with both. It was not pre-specified in the SAP and should not be added.
- **Status:** NEW.
- **Implication:**
  - The estrogen and hemoglobinopathy points rest on a handful of events whose timing is uncertain (Q3.4).
  - Non-significant global tests do not establish the assumptions. The per-term tests and plots show no gross violation.
- **Action:**
  - Report the fragility.
  - Consider pooling or penalizing the rare binary predictors, for example ridge or lasso shrinkage with predictor-specific penalties.
  - Re-estimate on correctly dated predictors.

**Q9. The age association** — A06g–A06j

- **Crude rates by age at the anchor, per 1,000 person-years:**

  | Age | Rate (95% CI) |
  |---|---|
  | 18–29 | 6.35 (4.85–8.17) |
  | 30–34 | 2.73 |
  | 35–59 (bands) | 3.30–5.11 |

  By attained age the rate is 6.83 below 30, then 2.9–4.7 (A06h).
- **The 60 strokes in women under 30:**
  - 52 ischemic, 3 ICH, 3 SAH, 2 unspecified
  - 10 within 30 days of the anchor
  - 4 without an acute code
  - 19 with migraine
  - 2 with a pregnancy code in the year before
- **HR per 10 years:**

  | Model | HR per 10 years (95% CI) |
  |---|---|
  | Crude | 0.94 (0.84–1.05) |
  | Nine-item | 0.90 (0.81–1.00) |
  | Plus a utilization proxy | 0.90 (0.80–1.00) |
  | Women without endometriosis | 0.83 (0.72–0.96) |
  | Women without fibroids | 0.95 (0.77–1.17) |
  | Acute-coded events only | 0.91 (0.81–1.01) |
  | Ischemic only | 0.91 (0.81–1.02) |
  | All anchored women, no Hb requirement | **0.97 (0.88–1.07)** |

  On the attained-age time scale the other HRs are unchanged (A06j).
- **Status:** NEW. The earlier attained-age summary is REPRODUCED.
- **Implication:**
  - The inverse age association is weak and depends on the model and the sample: it disappears without the Hb requirement.
  - It is driven by a small excess below age 30.
  - A linear term spreads that excess into a monotone decline that the data do not show, and the calculator turns it into "+1 point under 45". The spline gives HR about 1.0 at 35–44, the same as 45–54.
  - This is not a biological protective effect. Plausible explanations, none proven, include different reasons for presenting young, testing patterns, and young-onset stroke mechanisms.
- **Action:**
  - Either remove age from the points, or replace it with a pre-specified flexible or under-30 term and re-validate.
  - Do not describe age as protective.

### Part IV — Points and simplification

**Q10. Points audit** — `audit/output/A07_points.md`

- **Derivation (A07a):**
  - Points are the shrunken 14-predictor coefficients divided by B = β_HTN/2 = 0.3730. All raw values and rounded points reproduce exactly.
  - Hb bands use midpoints 13.2, 11, 9 and 7 against a reference deficit of 0 at Hb 13. Their raw values are 0.70, 1.40 and 2.10, which round to 1, 1 and 2.
  - Age bands use midpoints 30, 40, 50 and 57.5, with raw values 0, −0.28, −0.56 and −0.77, rounding to 0, 0, −1 and −1.
  - Diabetes (raw 0.533) and age 45–54 (−0.562) both lie within 0.1 of a rounding boundary, so slightly different data would change their points.
- **Shift:** development score = published score − 1 for every woman. The ranking (C 0.6967) and the score-to-risk mapping risk = 1 − 0.98980^exp(0.33504 × (points − 1)) are unchanged. All 14 rows of R6 reproduce, and the calculator's constants match.
- **Calibration of the points (A07b):** observed 5-year risk (Kaplan–Meier) tracks the table up to about 9 points. At 9 points or more, only 157 women and 23 strokes inform the table; for example, at 10 points the table says 18.9% and the observed risk is 6.4% (0.0–13.3).
- **All 1,024 combinations of items (A07d):**
  - Points rank them almost identically to the exact model: Spearman 0.992, 0.8% of pairs disagree.
  - Within one points total the exact risk varies widely; at 7 points it ranges from 5.8% to 11.1%.
  - The table extends to 19 points (98.6%), but the highest score observed in the cohort is 13.
  - A 30-year-old with no risk factors gets 1 point and a 50-year-old 0, while the exact model gives 1.00% and 0.81%.
- **Meaning of a point (A07c):** one point is a hazard ratio of 1.40 (1.34–1.46) only within the score model. Its implied HRs differ from the exact model's, for example diabetes 1.40 versus 1.22, and CAD 3.82 versus 4.06. In absolute terms a point adds more risk the higher the starting risk: going from 0 to 1 point adds 0.3 percentage points, from 9 to 10 adds 5.
- **Cost of using points:** apparent C falls by 0.0089 (95% CI 0.0014 to 0.0164). The two risks differ by more than 1 percentage point for 652 women (2.7%). Brier-based accuracy (IPA) is 1.82% for points versus 2.02% for the exact model.
- **Status:** REPRODUCED and NEW.
- **Implication:** the arithmetic is internally consistent, and the loss of information is small on average. The problems are the duplicated risk display and the extrapolated top of the table.
- **Action:** pick one output. If points are kept, report risk for points 0–9 and "≥10" as one group with observed risk.

**Q11. Removing five predictors** — A08a–A08c, A09e–A09f

- **Was it pre-specified?** No. It was a PI decision logged as a deviation, and it cites temporal C values, so it was made after the temporal results were known.
- **Evidence:**
  - Likelihood-ratio test of the five dropped predictors together: χ² 2.1 on 5 df, P 0.84.
  - Development difference in apparent C (nine minus fourteen): −0.0012 (−0.0068 to 0.0040). Optimism-corrected C: 0.701 versus 0.696.
  - IPA at 5 years: 2.02% versus 2.04%. Decision curves are nearly identical (A08g).
  - Temporal comparison, with models refitted on pre-2023 data:

    | Measure | Nine-item | Fourteen-item |
    |---|---|---|
    | C | 0.695 | 0.665 |
    | Calibration slope | 0.90 | 0.82 |
    | O/E | 1.44 | 1.42 |

    The difference in C is 0.030 (0.013 to 0.048).
  - **But** when the whole process, including predictor choice, is re-run on pre-2023 data only, it selects ten predictors (BMI enters) and gets temporal C 0.676 (0.623–0.728) (A09b). The nine-item advantage partly reflects having chosen the nine with all the data.
- **Zero points do not mean no prognostic value.** Smoking had a raw value of 0.37 points, about HR 1.15, and is badly measured. Dropping it makes risk in smokers depend on correlated items.
- **Status:** NEW.
- **Action:**
  - Do not claim the simplified model is superior.
  - Pre-register the predictor set before any external validation, and report the full-process estimates.

### Part V — Validation

**Q12. Internal validation** — A08a, A08b

- **Finding:**
  - Resampling was at the patient level (`R/16:13`), with 200 resamples.
  - The bootstrap refit only the **final nine-item model**. The 14-predictor fit, the points, the dropping of predictors and the shrinkage were **not** repeated, so the optimism of predictor selection was not captured.
  - Repeating the whole process gives optimism 0.0085 versus 0.0046, an optimism-corrected **C of 0.697** versus 0.701, and a calibration slope of 0.968 versus 0.980.
  - Only **13.5%** of resamples selected exactly the published nine predictors. Selection frequencies: diabetes 59%, age 73%, estrogen 91.5%, smoking 37.5%, BMI 28%.
  - Baseline survival and calibration-in-the-large were not bootstrapped; only the calibration slope was.
  - Monte-Carlo standard errors with 200 resamples: 0.0011 for optimism and 0.0035 for shrinkage, adequate for reporting C to 2–3 decimals.
  - No confidence intervals were reported originally. Apparent C is 0.705 (0.674–0.737).
- **Status:** NEW. (The process bootstrap uses heuristic shrinkage inside each resample, an approximation.)
- **Implication:** a bootstrap-corrected C is not evidence of transportability, and the published 0.701 slightly overstates it.
- **Action:** report 0.697 (full process), the selection-stability table and confidence intervals.

**Q13. Temporal validation** — `audit/output/A09_temporal_validation.md`

- **Code and split:** `R/16:22–28`. Development is anchors before 2023 (11,396 women, 226 strokes); validation is 2023 or later (12,610 women, 119 strokes). The model is the nine-item model refitted on pre-2023 data, **unshrunk**, with a 3-year horizon.
- **Reproduced with confidence intervals:**

  | Measure | Value (95% CI) |
  |---|---|
  | C at 3 years | 0.695 (0.644–0.746) |
  | Calibration slope | 0.90 (0.72–1.09) |
  | O/E, KM method | 1.56 (1.23–1.89) |
  | O/E, expected counts (handles censoring better) | 1.44 (1.19–1.73) |
  | Mean predicted 3-year risk | 1.06% |
  | Observed 3-year risk | 1.66% (1.30–2.01) |

- **Independence (A09b):** the validation data were **not** independent of development decisions. The predictor set, the points and the shrinkage were all derived using the full sample, and the temporal check validates the refitting procedure, not the published equation. The full-process temporal C is 0.676.
- **The 5-year horizon:** never temporally validated. Maximum follow-up in the validation set is 3.72 years.
- **Is the O/E meaningful?** Yes, statistically and clinically: the CI excludes 1, and absolute risk is about 45% higher than predicted.
  - It is not an artefact of short follow-up. The excess is in years 1–3 after the anchor: 4.92 versus 2.93 strokes per 1,000 person-years. The rate in the first 30 days is similar in both periods: 22.5 versus 20.4 (A09c).
  - It holds in every anchor year: 2023 1.42, 2024 1.33, 2025–26 1.57 (A09d).
  - Development-period apparent O/E is 1.00 (0.83–1.19).
- **Updating (A09h, analysis only):**
  - Recalibrating the slope adds nothing (P 0.30), but revising the coefficients does (χ² 28.8 on 9 df, P < 0.001).
  - In the 2023+ data the HRs are larger than in development for migraine (ratio 1.80, 1.25–2.59), smaller for diabetes (0.55) and higher for age (1.02 per year).
  - O/E is 1.09 (0.70–1.63) in women with no binary risk factor and 1.57 (1.27–1.92) in women with at least one.
- **Status:** REPRODUCED and NEW.
- **Implication:** updating the baseline alone would not fix calibration, because the predictor effects themselves changed over calendar time. That pattern fits changes in coding or ascertainment, including the look-ahead problem in Q3.4, at least as well as real change.
- **Action:**
  - Do not recalibrate the live calculator in place.
  - First fix predictor timing.
  - Then pre-specify an updating plan: develop on 2018–2022, update on 2023–2024, evaluate on 2025 and later or in an external cohort.

### Part VI — Calibration, utility, comparators

**Q14. Calibration** — A08d, A08e, `figures/A08_calibration_dev.png`, A09g, `figures/A09_calibration_temporal.png`

- **Development (apparent):**
  - O/E 1.02 (0.91–1.14) at 5 years by expected counts, 1.07 by KM over mean predicted.
  - Slope on cloglog(predicted risk) 1.03 (0.91–1.15).
  - Mean absolute error (ICI) 0.157 percentage points; E50 0.150; E90 0.271.
  - The decile table was coarse: deciles 1–4 span predicted risks of only 0.80–1.03%, and decile 9 was an outlier (2.52 predicted, 3.91 observed).
  - By the calculator's own bands, predicted versus observed 5-year risk:

    | Band | Women | Predicted | Observed (95% CI) |
    |---|---|---|---|
    | Under 1% | 7,580 | 0.89% | 0.72% |
    | 1–2.5% | 13,032 | 1.60% | 1.66% |
    | 2.5–5% | 2,472 | 3.54% | 4.17% |
    | 5–10% | 699 | 6.74% | 6.95% |
    | 10% or more | 223 (33 strokes) | 17.0% | 16.2% (10.3–21.8) |

- **Temporal:**
  - Underprediction in every band from 0.5% to 5%. For example, the 1–2% band was predicted at 1.33% and observed at 2.31% (1.52–3.11).
  - ICI 0.556 percentage points.
- **Status:** NEW.
- **Implication:**
  - A slope of 0.98 summarized the development pattern correctly, but it said nothing about calibration-in-the-large in new patients, which is where the model fails.
  - KM-based observed risks treat deaths as censored. Given the Q6 results, this makes no material difference.
- **Action:** report calibration-in-the-large with confidence intervals, flexible curves and band-level tables in every validation.

**Q15. Decision curves** — A08g, A09i

- **Original analysis (`R/15:224`):**
  - It handled censoring correctly, using KM risk within the high-risk group (Vickers et al. 2008).
  - The horizon was 5 years.
  - It covered the 14-predictor model only.
  - It had no confidence intervals and defined no clinical decision.
- **Repeated for the actual calculator** (5 years, development, apparent):
  - Net benefit beats both treat-all and treat-none from 1.5% to 10%. The difference over the better default is 0.0056 (0.0041–0.0068) at 2% and 0.0017 (0.0009–0.0027) at 5%.
  - At 1% the CI includes 0.
  - The points perform slightly worse and the 14-predictor model about the same.
  - Temporally, at 3 years, net benefit beats treat-all from 1% upward.
- **Absolute size:** 0.0056 means about 6 net true positives per 1,000 women at a 2% threshold.
- **Status:** NEW.
- **Implication:** no evidence-based action is linked to any threshold. Possible uses would be blood-pressure, migraine or estrogen counseling, or avoiding estrogen-containing therapy. None is tested, so **clinical utility is unestablished**.
- **Action:** define a use case and its threshold before any further decision-curve analysis.

**Q16. Comparators** — A08c, A09e, A09f

- **Why CHA₂DS₂-VASc was not used:** it is designed for atrial fibrillation, uses vascular history and age 65 or older, and is **not an appropriate comparator** for first stroke in women aged 18–60.
- **Why Framingham Stroke Risk or the Pooled Cohort Equations were not used:** they need systolic blood pressure, treatment status and cholesterol. Blood pressure is not in the data, and the LDL field is flagged as implausible. The Pooled Cohort Equations also start at age 40. They cannot be implemented validly here.
- **The comparator used instead:** a model of conventional risk factors fitted in the same data (age, hypertension, diabetes, dyslipidemia, smoking, BMI, CAD, CKD).

  | Measure | Conventional | Nine-item |
  |---|---|---|
  | Optimism-corrected C (development) | 0.640 | 0.701 |
  | Temporal C | 0.580 | 0.695 |
  | Temporal IPA | −0.09% | 1.04% |

  - Adding migraine, estrogen, hemoglobinopathy and Hb deficit to the conventional model: likelihood-ratio χ² 65.9 on 4 df, P < 0.001; ΔC 0.058 (0.033–0.081) in development and 0.091 (0.051–0.136) temporally.
- **Status:** NEW.
- **Implication:** the added predictors improve ranking substantially compared with conventional risk factors. But overall accuracy is low: IPA is 2% at 5 years. Migraine is the predictor most exposed to the timing problem.
- **Action:** repeat after the predictors are re-dated. This comparison is the most important test of whether the MED-HACK items add anything.

### Part VII — Interpretation

**Q17. What the model establishes**

| Claim | Verdict | Reason |
|---|---|---|
| 1. Predicts first stroke in women with uterine disease | **Partially supported** | Moderate discrimination internally and temporally, in tested women (with a CBC), for coded stroke. Underpredicts in recent years. Not externally validated. Some predictors are not baseline-only. |
| 2. Predicts stroke specifically in women with fibroids | **Not established** | 37% of the cohort has no fibroids. Performance among women with fibroids is similar (dev C 0.71, temporal C 0.70, O/E 1.34) but is a subgroup result. |
| 3. Fibroids independently increase stroke risk | **Not established** | There is no comparison group without uterine disease, and the design is prediction. Within the cohort, rates by diagnosis combination are 3.6–5.2 per 1,000 person-years, with overlapping CIs. |
| 4. Anemia increases stroke risk | **Not established (association only)** | Hb deficit is associated with stroke: HR 1.14 per g/dL, robust to outcome definition and imputation. Earlier analyses suggest it behaves like a marker of illness. Association is not causation. |
| 5. Treating anemia reduces stroke risk | **Not supported** | No intervention data. A positive coefficient is not a treatment effect. |
| 6. Estrogen-containing therapy increases risk in this population | **Not established** | 14 strokes; fragile (2 events); the exposure window extends after baseline; misclassified items (F-18 fluoroestradiol, IV conjugated estrogens); oral contraceptives under-captured. |
| 7. The score can guide treatment or screening | **Not supported** | No defined action or threshold; net benefit small; no external or prospective evaluation. |

The mixed endpoint (87% ischemic, plus ICH, SAH, CVT and chronic imaging-only infarcts) means each claim is about "any coded stroke". The predictors probably relate differently to ischemic and hemorrhagic stroke. Without external validation, every claim is limited to this institution and period.

**Q18. Subgroups and fairness** — `audit/output/A10_subgroups.md`

- **Development (apparent, 5 years):** C is 0.61–0.83 and O/E 0.82–1.16 across age, diagnosis, Hb, hypertension, diabetes, estrogen, migraine and race groups. All O/E confidence intervals include 1.
- **Temporal (3 years):** underprediction is concentrated in:

  | Subgroup | O/E (95% CI) |
  |---|---|
  | Migraine | 1.94 (1.50–2.47) |
  | Age 45–60 | 1.87 (1.40–2.45) |
  | Endometriosis only | 1.85 (1.28–2.59) |
  | No diabetes | 1.68 (1.34–2.08) |
  | Ages 18–34 | 1.00 (0.65–1.47) |
  | Diabetes | 1.05 (0.72–1.48) |

- **Race:** White 1.46 (1.18–1.78); Black 1.25 (0.60–2.30, 10 strokes); Asian 1.17 (5 strokes); other 1.79 (9 strokes). The non-White groups are too small for reliable estimates.
- **Measurement by race:** CBC availability (63.8–66.1%), smoking documentation, BMI recording and acute-code share were similar across groups. Imaging-positive events were less frequent in Black women: 30% versus 43% in White women, based on 33 strokes.
- **Hemoglobinopathy** is coded in 0.8% of Black and 0.8% of White women. That is far below the general prevalence of sickle cell trait in Black Americans (an external figure), so the flag reflects documented disease, not carrier status.
- **Race LR test:** P 0.877 shows only that race adds no discrimination. It does not establish fairness.
- **Status:** NEW.
- **Action:** report subgroup calibration with confidence intervals in any external validation; pre-specify the groups. Do not add race to improve C.

---

## Deliverable B — Audit outputs

All audit outputs are in `reanalysis/audit/` (scripts) and `reanalysis/audit/output/` (results):

| Script | Output | Content |
|---|---|---|
| `R/00_audit_setup.R`, `R/val_fns.R` | — | Redirects outputs; development-sample definition; validation helpers (O/E by expected counts, IPCW Brier, flexible calibration, net benefit, full-process selection) |
| `R/A01_cohort_design.R` | `A01_cohort_design.md` | Cohort reconciliation, diagnoses, time zero, censoring checks |
| `R/A02_predictors_missingness.R` | `A02_predictors_missingness.md` | Ascertainment windows, hemoglobinopathy timing, look-ahead check, smoking, missingness, included vs excluded, selection model, model in excluded women |
| `R/A03_outcome.R` | `A03_outcome.md` | Outcome algorithm, provenance, subtypes, strict-outcome refits, chart-review frame |
| `R/A04_followup_competing.R` | `A04_followup_competing.md` | Follow-up, horizons, Aalen–Johansen, Fine–Gray, censoring |
| `R/A05_reproduce_calculator.R`, `js/calc_test.js` | `A05_reproduce_calculator.md` | Refit and bootstrap reproduction; 14 profiles across R, spec, JavaScript, local page and public page |
| `R/A06_diagnostics_age.R` | `A06_diagnostics_age.md`, `figures/A06a_schoenfeld.png`, `figures/A06b_age_hb_curves.png` | PH, functional form, influence, fragility, collinearity, interaction, age |
| `R/A07_points.R` | `A07_points.md` | Points derivation, shift, score-to-risk, all combinations, information loss |
| `R/A08_internal_validation.R` | `A08_internal_validation.md`, `figures/A08_calibration_dev.png` | Full-process bootstrap, selection stability, model comparison, calibration, Brier, decision curves |
| `R/A09_temporal_validation.R` | `A09_temporal_validation.md`, `figures/A09_calibration_temporal.png` | Temporal validation with CIs, independence, rates by period, model comparison, updating test, decision curve |
| `R/A10_subgroups.R` | `A10_subgroups.md` | Subgroup C and O/E, differential measurement |
| `R/A11_multiple_imputation.R` | `A11_multiple_imputation.md` | MI vs complete case; design for imputation inside resampling (not executed) |
| `run_audit.R` | `audit_run_log.md`, `sessionInfo.txt` | Runs everything in order; log of tables per script |

**Not executed:**
- imputation inside resampling (designed in A11)
- chart review (sampling frame only)
- external validation
- re-derivation with dated predictors (data not in the project)

---

## Deliverable C — Prioritized action plan

| # | Task | Needs |
|---|---|---|
| **Priority 1 — potentially invalidating** | | |
| 1.1 | Re-derive every predictor strictly at or before the anchor (estrogen, hemoglobinopathy, smoking, migraine, diabetes, BMI; hypertension and CAD for procedure-anchored women); refit and re-validate | Dated raw extracts (DATA REQUIRED) plus patient-level data |
| 1.2 | Blinded chart review of 100 coded events (frame ready) to estimate PPV; exclude chronic imaging-only infarcts; decide whether events within 30 days are incident | Chart access (chart review) |
| 1.3 | External validation in an independent cohort (other institution or system): discrimination, calibration-in-the-large, flexible calibration, subgroups | Independent external cohort |
| 1.4 | Replace the 5-year output with a horizon the follow-up supports (3 years), or label 5-year risk as extrapolated | Existing files |
| 1.5 | Calculator: show one risk; mark estimates beyond the observed range; add unknown handling; remove or re-specify the age point | Existing files (no model change needed for the display) |
| **Priority 2 — before clinical interpretation** | | |
| 2.1 | Report the full-process optimism (C 0.697), selection stability, and confidence intervals for every performance measure | Done in this audit (A08, A09); needs incorporating into the report |
| 2.2 | Pre-specified updating plan for calendar drift (develop, update, evaluate), only after 1.1 | Patient-level data |
| 2.3 | Penalized estimation or pooling for the rare predictors (estrogen, hemoglobinopathy) and report their fragility | Patient-level data |
| 2.4 | Imputation inside resampling, if the intended use includes untested women | Patient-level data (several hours) |
| 2.5 | Define the clinical use and threshold; then decision-curve analysis on external data | Clinical input plus external cohort |
| 2.6 | TIA as a non-terminal event (keep women under follow-up after TIA) | Event-level dates (DATA REQUIRED) |
| **Priority 3 — strengthening** | | |
| 3.1 | Ischemic-only model as a pre-specified secondary analysis (301 events) | Patient-level data |
| 3.2 | Sensitivity analyses: exclude procedure-anchored women; exclude events within 30 days | Patient-level data |
| 3.3 | Subgroup calibration with pre-specified groups in the external cohort | External cohort |
| 3.4 | Correct the methods summary (below); TRIPOD checklist; labeling | Existing files |

---

## Deliverable D — Statistical-review conclusion

**Strongest evidence:**
- A large single-center cohort: 24,006 women and 345 coded strokes over 82,607 person-years.
- Moderate discrimination, about 0.70, that holds up on later patients (temporal C 0.695, 0.676 after full-process re-selection).
- Good apparent calibration.
- Hazard ratios that are stable under stricter outcome definitions, competing-risk modeling, imputation and the attained-age time scale.
- A clear gain over conventional risk factors, ΔC about 0.06–0.09 (the conventional model itself is weak).
- An exact, reproducible implementation of the calculator.

**Most serious limitations:**
1. Several predictors, including three of the most influential (migraine, estrogen, hemoglobinopathy), were partly recorded after time zero, in some cases after the stroke.
2. The outcome is unvalidated and code-based, and it includes chronic imaging-only infarcts.
3. Absolute risk is underpredicted by about 45% in recent patients, with evidence that the predictor effects themselves are unstable over time.
4. The 5-year horizon is beyond the observed follow-up of most women and was never temporally validated.
5. The validation sample was not independent of predictor selection.
6. The estrogen and hemoglobinopathy associations rest on a few events.
7. There is no defined clinical action, so utility is unestablished.

**Verdict:** **research use only.** The current evidence does not support prospective clinical evaluation or routine clinical use.

**Minimum changes before publication:**
- Re-derive the predictors with correct timing, or else present the model explicitly as a model using information available up to the end of follow-up. That framing would not be valid for prediction.
- Report the full-process optimism and temporal performance with confidence intervals, including O/E 1.44 (1.19–1.73).
- Present 3-year risk, or label 5-year risk as extrapolated.
- Report the outcome's code provenance.
- Use one risk display.
- Describe the population as fibroids, adenomyosis or endometriosis, in women with a CBC.
- Describe the age term honestly.

**Minimum evidence before clinical deployment:**
- external validation in an independent population with chart-validated outcomes and correctly dated predictors, showing adequate calibration-in-the-large and in subgroups
- a defined clinical decision with evidence that acting on the risk changes outcomes, or at least a prospective impact study

---

## Corrections needed in `calculator/MED-HACK_statistical_methods.md` (not changed by this audit)

1. "Candidate predictors … are all known at the diagnosis visit": **incorrect** for estrogen, hemoglobinopathy, smoking, BMI, and possibly migraine and diabetes (Q3.4).
2. "Median follow-up 2.8 years": the development sample's median is **2.95 years**. The 2.8 came from the earlier index-date cohort.
3. "TRIPOD type 1b/2a": temporal splitting is **type 2b** (non-random split), so the correct label is 1b and 2b.
4. The optimism-corrected C of 0.701 does not include predictor-selection optimism; the full-process value is **0.697**.
5. Temporal validation should be described as validating the refitting procedure, with the predictor set chosen on all data. O/E should be given as 1.44 (1.19–1.73) by expected counts, or 1.56 (1.23–1.89) by the KM method.
6. The outcome description should state that imaging-only chronic infarcts are included and that 33% of events had negative brain imaging.
7. The hazard ratios shown are unshrunk; the calculator applies shrinkage of 0.980.

*Data note:* the audit tables contain only aggregate counts, but some cells are small (for example 2–5 strokes). The repository is public. Consider making it private, or suppressing cells under 11, before sharing further.
