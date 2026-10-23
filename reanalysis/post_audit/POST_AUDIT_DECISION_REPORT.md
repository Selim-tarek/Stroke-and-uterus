# MED-HACK: post-audit decision report

**Date:** 2026-10-23
**Scripts:** `post_audit/R/P01`–`P08` (runner: `post_audit/run_post_audit.R`). Outputs: `post_audit/output/`. Starting evidence: `audit/AUDIT_REPORT.md`.
**No original file was changed:** not the scripts, specifications, calculator, outputs or derived data.

## Direct answers

**1. Is the current model useless, potentially useful, or ready for clinical use?**
- **As specified, it cannot be defended, and it is not ready for clinical use.**
- Its distinctive advantage, about +0.05 in C over a model of dated baseline predictors, comes entirely from predictors that cannot be shown to be known at the prediction time: migraine, estrogen, diabetes, and hemoglobinopathy coded at any date.
- That advantage may be real signal, or it may be information leakage. These data cannot tell the two apart.
- **A modest predictive signal is recoverable.** A five-predictor model using only dated baseline information (age, Hb deficit, hypertension, CAD, CKD) is **potentially useful for research**, with optimism-corrected C 0.661 and temporal C 0.663. Its absolute risks are not yet trustworthy.

**2. What supports this?**
- **The datable predictors carry modest signal.** The baseline-only models beat a conventional comparator in development (ΔC +0.022 to +0.027, CIs excluding 0) and temporally (C4 +0.045, 0.004–0.093) (P05d, P06c).
- **Removing the undatable predictors costs about 0.05–0.07 in C** (P05d, P06c).
- **Every model underpredicts by about 45% in 2023+ patients,** and a baseline recalibration estimated on 2023–24 did not correct 2025+ (P06b, P06d).
- **Overall accuracy is low:** IPA about 1% at 3 years (P05c).
- **C4 discriminates poorly in women without fibroids:** temporal C 0.59 (P08a).

**3. What can be repaired with the current data?**

| Fix | Where |
|---|---|
| Restrict to dated predictors; re-date hemoglobinopathy | P02 |
| Drop the needless BMI requirement (+880 women) | P02a |
| Use an acute-coded outcome | P03 |
| Use a 3-year horizon | P04 |
| Pre-specify predictors, so no selection optimism | P05 |
| Remove the points system and the misleading age point | P07, CALCULATOR_RECOMMENDATIONS |
| Correct the population statement | P08 |
| Report uncertainty and temporal performance honestly | VALIDATION_REPORT |

**4. What needs missing data?** (See `MISSING_DATA_REQUEST.md`.)

| Problem | Needs |
|---|---|
| Whether migraine, diabetes and estrogen add genuine baseline information | Dated diagnosis and medication extracts |
| Outcome PPV | Chart access and a raw stroke extract |
| A complete conventional comparator | Blood pressure; dated diabetes |
| The cause of calibration drift | The above, plus continued follow-up |
| Transportability | An external cohort |

**5. Single most important next analysis:**
- Obtain **dated diagnosis and medication records** for migraine, diabetes and estrogen-containing therapy (request items 1–2).
- Re-derive them at or before the anchor.
- Re-run P05 and P06 with them added to C4.
- **This decides whether MED-HACK's distinctive predictors are real.**

**6. Should the calculator remain research-only?**
- **Yes, at most.** Add the warning in `CALCULATOR_RECOMMENDATIONS.md` §1 now, or take the page offline; that is the PI's decision.
- Its current numbers rest on undatable predictors, a points table and an age point the data contradict, and a 5-year horizon that is not validated.

**7. What would make me abandon this approach?** Any one of these:
- **(a) Dated predictors add nothing:** once dated, migraine, diabetes and estrogen add < 0.02 in temporal C over C4, *and* C4 itself stays below about 0.65 against a complete conventional comparator.
- **(b) The outcome is unreliable:** blinded chart review gives a PPV < 0.75 for acute-coded events.
- **(c) External failure:** in an external cohort, C < 0.60, or calibration-in-the-large outside 0.7–1.4 after one pre-specified recalibration on separate data.
- **(d) No use case:** no clinical decision can be defined for which a 3-year coded-stroke risk would change management (see below).

## Verdict

**Redesign, do not retain as is, do not yet abandon.**

The data support a small, genuine improvement over conventional risk factors using dated predictors. They do not yet support the MED-HACK-specific predictors, accurate absolute risk, or any clinical use.

## Intended-use statement (Part 9)

| Option | Supported? | Consequence |
|---|---|---|
| A. Stroke risk in women with fibroids | **No.** 37% of the cohort has no fibroids. Fibroid-only performance is a subgroup result: C4 temporal C 0.70 (0.63–0.77), O/E 1.30. | Would require a fibroid-only cohort and model |
| B. Stroke risk in women with fibroids, endometriosis or adenomyosis | **Partly.** It matches the cohort, but omits the CBC requirement, the coded outcome and the single center. C4 performs worse without fibroids (0.59). | Acceptable only with the qualifiers in C |
| **C. Risk of *coded* stroke in a selected single-center cohort of women with benign uterine disease who had a CBC** | **Yes** | Use it in the title, abstract, methods and calculator. Results must report the drop without fibroids. The discussion must not claim fibroids or anemia cause stroke. |

## Clinical utility and readiness (Part 8)

| Possible use | Decision | Horizon and threshold | Intervention | Evidence of benefit in this population | Can a DCA evaluate it now? |
|---|---|---|---|---|---|
| Choosing contraception or menopausal hormone therapy | Avoid estrogen-containing products in higher-risk women | 3 years; no justified threshold | Use a non-estrogen alternative | Guideline contraindications already rest on individual factors (e.g. migraine with aura, hypertension), not a risk score; no trial of score-guided choice | No. The score's estrogen input is unreliable, and no threshold exists. |
| Intensifying vascular risk-factor care | Blood pressure or lipid management, smoking cessation | 3 years; no justified threshold | Guideline management, which is indicated anyway by the individual factors | Benefits come from treating the factors themselves; no evidence that a score adds to this | No |
| Anemia work-up or treatment | Iron studies or treatment | — | Iron, or treating the bleeding | Indicated by the anemia itself; **no evidence that treating it lowers stroke risk** | No. Prediction is not causation. |
| Antiplatelet or anticoagulant therapy, or imaging | — | — | — | **None.** Not appropriate, and must not be suggested. | No |

**Conclusion: clinical utility is unestablished.** No defensible decision exists for which the score changes management.

**Readiness levels:**

| Level | Requirements | Status |
|---|---|---|
| 1. Retrospective research use | — | **Current status** (C4 or a successor) |
| 2. Prospective observational evaluation, not directing care | Dated predictors; PPV ≥ 0.85; a fixed model frozen before the prospective period; a pre-registered protocol | Not reached |
| 3. Clinical decision support | External validation with acceptable calibration; a defined decision and threshold; an impact study showing benefit | Not reached |

## Completed, failed, infeasible, unresolved

**Completed (run, outputs inspected):**

| Script | Content |
|---|---|
| P01 | Inventory, integrity check, git timeline |
| P02 | Provenance and baseline reconstruction; utilization test (inconclusive) |
| P03 | Outcome definitions |
| P04 | Horizon and censoring |
| P05 | Seven candidate models with full-process bootstrap |
| P06 | Temporal validation and drift |
| P07 | Calculator support |
| P08 | Fibroid subgroups |

**Failed or inconclusive:**
- The utilization-bias test (P02d) cannot separate look-ahead from sicker women being tested more, because the baseline-dated controls show the same pattern.
- The hypothesis that drift came from events in the first 30 days was rejected (P06h).

**Infeasible:**
- Dating migraine, diabetes and estrogen.
- PPV.
- A three-way temporal design at 3 years.
- 5-year temporal validation.
- External validation.
- A complete conventional comparator.

**Unresolved:**
- Whether the undatable predictors are real signal or leakage.
- The cause of calibration drift.
- A clinical decision.

**Statistical judgments (not empirical findings):**
- The abandonment criteria.
- Choosing C4 over C3. They perform equally; C4 is simpler and contains no noise predictors.
- 3 years as the primary horizon.

## Minimum next steps, in order

1. Add the research-prototype warning to the calculator, or take it offline (PI decision). Correct the methods summary (`MANUSCRIPT_CORRECTIONS.md`).
2. Request items 1–4 of `MISSING_DATA_REQUEST.md`: dated diagnoses and medications, the raw stroke extract, chart access.
3. Run the blinded 100-event chart review (protocol in `ISSUE_RESOLUTION_MATRIX.md`).
4. Re-derive the dated predictors. Pre-register C4 and C4 plus dated migraine, diabetes and estrogen. Re-run P05 and P06 unchanged.
5. If abandonment criteria (a) and (b) are not met, freeze one model and seek an external cohort.
6. Only then consider a new research calculator built to the `CALCULATOR_RECOMMENDATIONS.md` specification.

## Major remaining uncertainty

**Leakage versus signal in migraine, diabetes and estrogen.**
- If dated versions keep their association, the model regains about +0.05–0.07 in C and becomes a plausible research tool.
- If they lose it, what remains is essentially a conventional-risk-factor model plus hemoglobin, with modest discrimination and unstable calibration.

Second, the calibration drift. Its cause is unknown; if it reflects changing outcome coding, absolute risk may never be reportable from this source.
