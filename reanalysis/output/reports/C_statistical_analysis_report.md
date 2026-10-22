# Deliverable C — Statistical analysis report

All estimates are from the clean end-to-end run of `reanalysis/run_all.R` (R 4.3.3, seed 20261022; 35 min). Table and figure names refer to `output/tables/` and `output/figures/`. The analysis plan (`SAP.md`) was committed before any model was fitted. Deviations are listed in `output/logs/deviations.md`.

The report has three separate parts:
- **Part 1:** reproduction of the original analysis
- **Part 2:** the revised primary analysis
- **Part 3:** secondary, longitudinal and sensitivity analyses

---

## Part 1 — Reproduction of the original analysis [S2a_reproduction_comparison]

**Method.**
- The original Python models were re-fitted in R using the original Codebook variables and rules:
  - anemia grade from Codebook `hgb` (closest to index within ±3 y, either direction)
  - logistic regression with HC1 SEs
  - the same covariates and the hormonal-type sparse-merge rule
  - the original "2a" pre-Hb conditions, original time-to-event start, and original time-updated algorithm
- The originally reported values were read by script from the original results file, not retyped.

**Result: 73 of 80 estimates that have an original value were reproduced within 2%; 4 within 10%; 3 differed by more than 10%.**

| Original analysis | Reported | R | Status |
|---|---|---|---|
| Cross-sectional OR, any stroke: mild | 1.20 (1.03–1.38) | 1.20 (1.03–1.38) | Reproduced (N 31,472; 1,460 events) |
| Cross-sectional OR: moderate | 1.77 (1.41–2.21) | 1.77 (1.41–2.21) | Reproduced |
| Cross-sectional OR: severe | 2.01 (1.39–2.91) | 2.01 (1.38–2.91) | Reproduced |
| Ischemic / incident-vs-none ORs, 2a, Hb ≥30 d before, hemoglobinopathy exclusion, + bleeding/diagnosis group | — | — | All reproduced |
| Spline OR vs 13 g/dL, at 7–16 g/dL | — | — | All reproduced |
| Per 1 g/dL below 13 (2a) | 1.11 (1.06–1.16) | 1.11 (1.06–1.16) | Reproduced |
| Original time-to-event (31,400 women; 100,309 PY; 458 events): moderate, 2a | 1.57 (1.09–2.28) | 1.57 (1.09–2.28) | Reproduced |
| Time-updated (31,924 women; 465 events): moderate, 2a | 1.86 (1.29–2.69) | 1.80 (1.24–2.62) | Close |
| Time-updated: severe, 2a | 2.25 (0.99–5.13) | 2.63 (1.23–5.64) | **Not reproduced** (6 events) |
| Time-updated + testing: moderate | 1.54 (1.05–2.25) | 1.49 (1.01–2.19) | Close |
| Time-updated, 30-day lag: moderate | 1.30 (0.85–1.98) | 1.29 (0.85–1.98) | Reproduced |

**Explanation of the time-updated differences.**
- The original code kept every same-day Hb row, so the value in force was the last row of the day. R uses the daily median.
- On 11,326 woman-days, same-day results fall in different grades, and 5,946 of these involve Hb <8 g/dL.
- The "severe" time-updated estimate rests on 6 events, so it is unstable to this arbitrary choice. It should not be reported as a finding.

**Conclusion on reproducibility.** The original numbers are computationally correct. The problems are in **design and interpretation**, not arithmetic (Part 2).

**Issues found in the original specification:**
1. The primary outcome `stroke_any` includes **786 strokes before or at the index encounter (51%)** and 202 undated strokes. The ORs therefore mainly describe prevalent stroke.
2. The exposure Hb could be measured **after** the index date (13,661 women) and **after the stroke** (699 of the stroke cases with known order).
3. The original "pre-Hb" covariates had several definition problems:
   - GI bleeding and pregnancy were counted up to 30 days **after** the Hb.
   - HIV counted at any date.
   - CKD omitted women coded only with ESRD.
4. The original time-updated model:
   - dropped person-time with no Hb in the window, together with the 16 strokes that fell in it
   - let a Hb drawn in the days before a coded stroke define the exposure at that stroke
   - used non-clustered SEs. Clustering by woman changed the CIs negligibly (moderate 1.80, 1.23–2.63).
5. The previous briefing stated that adjusting for testing frequency "did not change the associations". The reproduced estimates (1.80 → 1.49) show it attenuated them.

---

## Part 2 — Revised primary analysis (PI decision 2026-10-22)

### Estimand
Among women with benign uterine disease who were free of recorded stroke/TIA at their first uterine diagnosis (index date):
- **Exposure:** anemia on the most recent laboratory CBC Hb in the 3 years up to and including the index date.
- **Outcome:** the rate of a first **stroke (excluding TIA)** coded after the index date.

### Cohort [Fig1_flow_counts; Table1_baseline_by_grade; Figure 1]
- **24,871 women** in the models, 82,184 person-years, median follow-up 2.8 years (IQR 1.2–5.2).
- **342 incident strokes:** 298 ischemic, 24 ICH, 11 SAH, 2 CVT, 7 unknown type.
- 85 TIAs (secondary outcome) and 162 deaths without stroke.
- **Baseline anemia: 5,471 (22.0%)** — mild 4,008; moderate 1,062; severe 401.
- **Distinctive features of anemic women** (Table 1):
  - more fibroids (56% of severe vs 36% of non-anemic) and heavy bleeding (85% vs 57%)
  - more hemoglobinopathy (4.0% vs 0.6%) and CKD (6.5–7.7% vs 1.9%)
  - more recorded heart failure, VTE and AF (undated)
  - less migraine and less dyslipidemia

### Model
- Poisson regression for the stroke rate (offset log person-years, robust SE), pooled over 20 imputations of BMI.
- Confirmed with a Cox model.
- **Adjusted for:** age (spline), race, BMI, smoking, hypertension, diabetes, dyslipidemia, CAD, migraine, hormonal therapy type, uterine diagnosis group, and conditions dated before index (hemoglobinopathy, CKD, liver disease, alcohol use disorder, IBD/HIV, malabsorption, GI bleeding and pregnancy in the prior year).
- **Model checks:**
  - 32 degrees of freedom; 10.7 events per parameter; all fits converged
  - proportional hazards: global Schoenfeld P = 0.40 [S2b]
- **Why these covariates** (`output/dag.md`):
  - Uterine bleeding is a cause of anemia, not a confounder, so it is not in the model; adding it changes nothing (S6).
  - Thrombophilia, AF, heart failure and VTE are undated and often recorded after a stroke, so they are excluded (S7).

### Results [Table2_primary_models; Table2a_rates_by_grade; Figure 3A; Figure S1]

| Baseline Hb | Women | Strokes | Rate /1,000 PY (95% CI) | 5-y cumulative incidence | Adjusted RR (95% CI) | Adjusted HR (95% CI) |
|---|---|---|---|---|---|---|
| ≥12 g/dL | 19,400 | 236 | 3.67 (3.22–4.17) | 1.8% (1.5–2.0) | 1.00 | 1.00 |
| 10.0–11.9 | 4,008 | 76 | 5.79 (4.56–7.25) | 2.6% (2.0–3.3) | 1.49 (1.14–1.96) | 1.48 (1.14–1.94) |
| 8.0–9.9 | 1,062 | 22 | 6.27 (3.93–9.49) | 2.9% (1.9–4.6) | 1.49 (0.96–2.31) | 1.48 (0.94–2.31) |
| <8.0 | 401 | 8 | 6.58 (2.84–12.96) | 3.1% (1.5–6.4) | 1.64 (0.79–3.40) | 1.60 (0.78–3.28) |
| Any anemia (<12) | 5,471 | 106 | — | — | **1.50 (1.18–1.91)** | — |
| Per grade (trend) | | | | | 1.26 (1.09–1.45) | 1.25 (1.08–1.45) |

- **Crude to adjusted:** the any-anemia RR was 1.62 (1.29–2.04) crude, 1.63 after age and race, and 1.50 fully adjusted. The confounders measured here explain little of the association.
- **E-value** for any anemia: 2.37 for the point estimate and 1.64 for the CI limit. An unmeasured confounder associated with both anemia and stroke by a rate ratio of about 1.6 could move the CI to include 1.
- **Continuous Hb** [Table2b; Figure 2]:
  - Below 13 g/dL, each 1 g/dL lower Hb: RR **1.14 (1.05–1.23)**. Above 13 g/dL: 1.02 (0.85–1.23).
  - Spline vs 13 g/dL: 12 g/dL 1.23 (1.03–1.47); 10 g/dL 1.52 (1.18–1.94); 8 g/dL 1.72 (1.16–2.56).
  - Overall P = 0.010; non-linearity P = 0.21.
  - The categorical and continuous analyses agree that the rate is higher below about 12–13 g/dL. Neither shows a well-estimated gradient beyond that.

**Interpretation.**
- Anemia measured before the uterine diagnosis is associated with about a 50% higher subsequent rate of recorded stroke.
- The categories do not show a clear severity gradient: point estimates are 1.49, 1.49 and 1.64. The moderate and severe estimates rest on 22 and 8 events and are compatible with no association.
- The per-grade "trend" (1.26) is a summary of this pattern, not evidence of a biological dose–response.
- **In absolute terms:** 5-year stroke incidence was 1.8% without anemia and 2.6–3.1% with anemia.

---

## Part 3 — Secondary, longitudinal and sensitivity analyses

### 3.1 Secondary incident outcomes [Table3; Table3a] (exploratory; BH q-values across trend tests)

| Outcome | Events | Any anemia RR | Per grade RR | q |
|---|---|---|---|---|
| Ischemic stroke¹ | 298 | 1.52 (1.18–1.95) | 1.26 (1.09–1.46) | 0.011 |
| Stroke or TIA (Codebook definition) | 427 | 1.44 (1.16–1.80) | 1.20 (1.05–1.37) | 0.019 |
| TIA² | 85 | 1.10 (0.66–1.81) | 0.84 (0.63–1.14) | 0.34 |
| ICH²,³ | 24 | 1.46 (0.61–3.49) | 1.59 (0.93–2.71) | 0.15 |
| SAH²,³ | 11 | 2.13 (0.61–7.40) | 1.22 (0.67–2.24) | 0.51 |
| CVT, unknown type | 2, 7 | counts only | — | — |

¹ Reduced pre-declared covariate set (events per variable would be <10 with the full set).
² Age + hypertension only; no TIA occurred in women with moderate or severe anemia.
³ Events per variable <10: unstable, exploratory.

- Ischemic stroke drives the primary result.
- **TIA shows no association.** Combining TIA with stroke (as the Codebook outcome does) dilutes the moderate-grade estimate to 1.21 (0.79–1.87). This supports analyzing them separately.
- Hemorrhagic subtypes are too sparse to interpret.

### 3.2 Sensitivity analyses [Table5_sensitivity; each row names the bias it addresses]

| | Analysis | Bias addressed | Any anemia | Per grade |
|---|---|---|---|---|
| P | Primary (MI) | — | 1.50 (1.18–1.91) | 1.26 (1.09–1.45) |
| S1 | Complete case (BMI) | Missing-data method | 1.49 (1.17–1.90) | 1.25 (1.08–1.45) |
| S2 | Follow-up from day 31 | Reverse causation / detection at index | 1.38 (1.06–1.80) | 1.20 (1.02–1.40) |
| S3 | Hb within 1 y before index | Stale exposure | 1.39 (1.07–1.80) | 1.21 (1.04–1.42) |
| S4 | Including point-of-care Hb | Exposure definition | 1.49 (1.16–1.90) | 1.25 (1.08–1.45) |
| S5 | Events with an acute stroke code only | Outcome misdating | 1.57 (1.23–2.01) | 1.28 (1.11–1.47) |
| S6 | + uterine bleeding | Exposure determinant | 1.52 (1.19–1.94) | 1.27 (1.10–1.47) |
| S7 | + thrombophilia/AF/HF/VTE | Over-adjustment | 1.46 (1.14–1.87) | 1.23 (1.06–1.42) |
| S8 | + Hb tests in 2 y before index | Healthcare contact | 1.45 (1.14–1.85) | 1.23 (1.07–1.43) |
| S9 | Excluding anemia-causing conditions | Confounding by illness | 1.46 (1.08–1.97) | 1.24 (1.03–1.48) |
| S10 | Fine–Gray (subdistribution HR) | Competing death | — | 1.24 (1.07–1.43) |
| S11 | IPW for having a baseline Hb | Selection | — | 1.24 (1.07–1.43) |

- The baseline association is stable: any-anemia RR 1.38–1.57 and per-grade 1.20–1.28 across every pre-specified sensitivity analysis.
- The moderate-anemia category alone ranges from 1.08 (S9) to 1.54 (S5), and its CI includes 1 in every analysis.

### 3.3 Time-updated anemia [Table4; Table4a; Figure 3C] (secondary)
**Design:**
- 38,326 women stroke-free at index, 119,305 PY, 451 strokes.
- The most recent CBC Hb carries forward for up to 1 year and never takes effect on its own measurement day.
- Time without a recent Hb is a separate category.
- 20 imputations; SEs clustered by woman.

| Model (stroke, excl. TIA) | Mild | Moderate | Severe | Per grade | No recent Hb |
|---|---|---|---|---|---|
| No lag | 1.15 (0.86–1.55) | **2.40 (1.63–3.52)** | 2.49 (1.01–6.13) | 1.40 (1.19–1.66) | **0.55 (0.44–0.69)** |
| + Hb tests in prior year | 1.03 (0.76–1.39) | 1.96 (1.32–2.92) | 2.09 (0.85–5.13) | 1.28 (1.08–1.52) | 0.60 (0.46–0.80) |
| 30-day lag | 1.13 (0.83–1.54) | 1.56 (0.98–2.48) | 1.60 (0.51–5.05) | 1.20 (1.00–1.44) | 0.82 (0.67–1.01) |
| 90-day lag | 0.97 (0.69–1.36) | 1.57 (0.96–2.57) | 2.28 (0.85–6.15) | 1.19 (0.97–1.45) | 1.03 (0.84–1.27) |
| 30-day lag + Hb tests | 0.96 (0.70–1.32) | **1.20 (0.75–1.92)** | 1.26 (0.40–3.96) | **1.05 (0.87–1.27)** | 1.00 (0.75–1.33) |

**Interpretation.** This is the most informative and least favorable set of results.
- **Without a lag, the no-recent-Hb category has about half the stroke rate.** Person-time without blood tests also has fewer recorded strokes, so testing and stroke coding happen at the same encounters.
- **The strong no-lag moderate estimate (2.40) halves on the log scale with a 30-day lag** (1.56). With both a lag and adjustment for testing frequency, it disappears (1.20, 0.75–1.92; per grade 1.05). Much of the time-updated association comes from Hb drawn shortly before a recorded stroke, and from how often women are tested.
- Testing frequency may partly lie downstream of anemia, so the fully adjusted model may over-adjust. **The data cannot separate surveillance from a true short-term association.**

### 3.4 Landmark analysis [Table4b]
Women anemic at baseline whose anemia had **resolved** by 12 months had a stroke rate similar to those whose anemia **persisted**:
- persistent vs resolved RR **0.93 (0.53–1.62)** at 12 months
- 1.07 (0.58–1.98) at 24 months

Both groups had higher point estimates than non-anemic women (12 months: resolved 1.47, 0.91–2.35; persistent 1.36, 0.89–2.08), but with wide CIs. If anemia raised stroke risk through low oxygen-carrying capacity, its resolution would be expected to lower the rate. The data do not show that, although they are too imprecise to exclude it.

### 3.5 Selection and missing data [S1a–S1d]
- **Selection into testing.** Women without a baseline Hb (35% of stroke-free women) had less recorded comorbidity, fewer recorded deaths and a lower crude stroke rate (2.94 vs 4.16 per 1,000 PY). After adjustment the rate ratio was 0.85 (0.69–1.07).
  - Weighting the tested cohort to resemble all stroke-free women on measured characteristics (S11) did not change the estimate.
  - Selection on unmeasured health-care intensity cannot be corrected, and multiple imputation does not address it.
- **Missing values.**
  - BMI (3.5%): MI and complete case agree.
  - Smoking is "unknown" for 66% and is modeled as a category. Its effect on confounding control is unknown and probably not ignorable, because missing smoking status signals less documentation.
- **Ferritin** is not used: it was measured in a selected subset and cannot be traced to source. The earlier "protective" ferritin finding should be removed or reported only as a selection artifact.

### 3.6 Multiplicity
- One primary estimand.
- Secondary-outcome trend tests: BH-adjusted q-values reported.
- Everything in Part 3 is exploratory.
- The original paper reported several dozen subgroup, morphology, pattern and subtype analyses without correction. These should be reduced, labeled exploratory, or moved to the supplement.

---

## Answers to the reviewer questions

| Question | Answer |
|---|---|
| Is the primary association reproducible? | **Yes, computationally** (73/80 within 2%). But it describes mostly prevalent stroke and Hb often measured after the stroke. |
| Does exposure precede outcome? | Original cross-sectional: **often not**. Revised primary: **yes** (Hb on or before index; events after index). Original time-to-event: yes, but its start was the later of index and Hb. Time-updated: yes, but values close to the event drive the association. |
| Is the primary model appropriately specified? | The original adjusted for an undated variable likely recorded after stroke (thrombophilia) and had mistimed "pre-Hb" flags. The revised model is DAG-based; results are robust to adding or removing the debatable covariates. |
| Does missingness or selective testing matter? | Selective testing is real and substantial. Measured selection does not change the estimate (S11); unmeasured surveillance remains the main threat (3.3). |
| Robust to reasonable alternatives? | The baseline any-anemia association: yes (1.38–1.57). The severity gradient: no. The time-updated association: no (it depends on lag and testing). |
| Are subtype analyses adequately powered? | Ischemic only. TIA null. Hemorrhagic subtypes, CVT and unknown type: no. |
| Do continuous and categorical analyses agree? | Yes: both show a higher rate below about 12–13 g/dL and no well-estimated gradient beyond that. |
| Clinically meaningful and precise? | Absolute 5-year difference of about 1 percentage point (1.8% vs 2.6–3.1%). The any-anemia RR is reasonably precise (1.18–1.91); grade-specific estimates are not. |
| Marker or cause? | The data favor anemia as a **marker** of underlying illness, health-care contact and vascular risk. This is supported by the time-updated attenuation, the resolved-vs-persistent landmark result, the absent gradient and the null TIA result. The data cannot fully separate these explanations, and they do not test whether correcting anemia changes risk. |
