# Calculator recommendations

**Status: proposal only.** The production calculator (`calculator/`, https://selim-tarek.github.io/Med-Hack-score/) was **not changed** in this task.

Evidence: `post_audit/output/P07_calculator_support.md`, `audit/output/A05_reproduce_calculator.md`, `audit/output/A07_points.md`, and P05/P06.

## 1. Should the current calculator stay online?

Only as **research-only, with an added warning**. Its exact-model risk relies on predictors that cannot be dated to the prediction time: migraine, estrogen, diabetes, and hemoglobinopathy from any date. In women anchored from 2023 on, it underpredicts by about 45–50% (P06b). Its 5-year risk has never been tested on later patients. The minimum immediate change is a text change:

> "Research prototype. Not for clinical decisions. Predicts *coded* stroke in one center's records; underestimates risk in recent patients; several inputs could not be confirmed as known at the time of prediction."

The alternative is to take the page offline until a re-specified model exists. That is the PI's decision.

## 2. Points system: do not retain

- **What it adds:** nothing beyond presentation. It loses discrimination: ΔC 0.0089 (0.0014–0.0164) (A07e). Its IPA is lower: 1.82% versus 2.02% (A08f). Its net benefit is equal or lower at every threshold (A08g).
- **It shows a second, conflicting risk:** 3.84% versus 5.30% for the worked example.
- **Its scale extends past the data:** the table runs to 19 points, but the highest observed score is 13, and only 91 women (12 strokes) scored 10 or more (P07a).
- **The age point is misleading:** "+1 under 45" goes against the observed data. Three-year risk is 1.05–1.09% at ages 30–44 (+1 point) and 1.59% at 45–54 (0 points) (P07b). It also suggests younger age is a risk factor, which the data do not support beyond a small excess below age 30.
- **Several items sit on rounding boundaries:** diabetes 0.53, age −0.56 (A07a).

**Recommendation:** remove the points and the points-to-risk table from any future version. If a paper-based tool is wanted later, derive it from the re-specified model and validate it separately.

## 3. Proposed output specification for a future research calculator

The future calculator should be built only after a re-specified model passes the steps in `POST_AUDIT_DECISION_REPORT.md`.

| Element | Specification |
|---|---|
| Model | The re-specified baseline-only model (predictors dated at or before the anchor; see `PREDICTOR_PROVENANCE.md`) |
| Horizon | **3-year** risk only. No 5-year output until follow-up supports it; 57% of 2018–2022 anchors had reached 5 years (P04a). |
| Number of risks shown | **One** risk, from the exact model |
| Format | Risk as a percentage with one decimal place, plus a **95% uncertainty interval** from the bootstrap distribution of the prediction for that profile; and a comparison with the average 3-year risk in the development cohort |
| Supported range | Show a numeric risk only when **both** hold: the predicted risk is at or below the development-sample 99th percentile; and the yes/no combination occurred in at least 20 development women (33 of 128 combinations for the current items, covering 99.2% of women, P07a). Otherwise show "Outside the range of the data used to build the model: estimate not shown." |
| Unknown status | Every yes/no item has three choices: **Yes / No / Not known**. "Not known" must never be treated as "No". If any item is "Not known", show the **range** of risks across the possible answers (minimum to maximum), with a note. The model was fitted on coded data where absence of a code means "no", so the help text must say "No = no documented diagnosis". |
| Hemoglobin | Required: a laboratory CBC value at or before the visit, within 3 years. Without it, show no estimate, because the model applies only to tested women. |
| Calibration notice | Until recalibrated and validated on separate data: "In patients diagnosed in 2023 or later, observed risk was about 1.5 times the predicted risk." |
| Intended-use statement | "Estimates the 3-year risk of a first *coded* stroke in women aged 18–60 with fibroids, adenomyosis or endometriosis who had a CBC, based on one US center's records. Research use only; not externally validated; not for treatment decisions." |
| Language | No causal wording. No "protective" factors. Avoid "risk factor" for items that are predictors only. |
| Prohibited outputs | No risk category tied to an action (no "treat", "screen" or "stop hormones") until a clinical decision and threshold are justified (`ISSUE_RESOLUTION_MATRIX.md` row 12) |
| Versioning | Model version, development dates, and a link to the validation report on the page; every change logged |
| Validation before release | Verify the JavaScript against R for pre-specified profiles (the method in `audit/R/A05`). The update data must not be the same data used to evaluate it. |

## 4. Safe communication of uncertainty

- Report the interval for the individual prediction, not just the point estimate.
- State the cohort average next to the individual risk. Do not use "×-fold" language for small absolute risks.
- Suppress estimates outside the supported range, rather than showing extrapolated values: the current calculator shows up to 99.99%.
