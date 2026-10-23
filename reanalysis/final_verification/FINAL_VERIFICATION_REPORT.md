# MED-HACK: final evidence verification

**Date:** 2026-10-23. **Starting commit:** `21e94f0`.

**Evidence base:**
- verification scripts `final_verification/R/V00`–`V06`, with outputs in `final_verification/output/`
- the claim ledger `CLAIM_LEDGER.csv` (38 claims: 27 verified, 10 corrected, 1 unresolved)
- the audit (`audit/`) and post-audit (`post_audit/`) outputs, traced line by line

**Preservation:** no original file changed. Checksums of all 243 protected files are identical before and after (REPRODUCTION_LOG).

**Labels used below:**
- **[V]** verified empirical result
- **[I]** interpretation
- **[J]** statistical judgment
- **[U]** unresolved
- **[B]** blocked by unavailable data

**Small cells:** counts of 1–10 are shown as "<11".

## 1. Executive verdict

- **[V]** The arithmetic of every reported figure reproduces:
  - The P05 bootstrap re-run from scratch is identical to the cache (max |Δ| = 0).
  - Every temporal C, O/E and IPA was recomputed independently, and matches.
- **[V]** Ten claims needed correction, mostly in interpretation and uncertainty rather than in the numbers themselves. The four that change conclusions:
  1. **The baseline model's gain over conventional risk factors is smaller and less certain than reported.** The conventional comparator (C1) was weakened by two noise predictors. Against age + hypertension + CAD + CKD, adding Hb deficit (C4) gains ΔC **+0.022 (−0.010 to 0.056)** temporally.
  2. **The development "CIs" for ΔC in P05d are not confidence intervals.** They show model-fitting variability only.
  3. **The recent-patient underprediction is a diagnosis-cohort effect, not a calendar-period effect.** Women diagnosed before 2023 and followed in 2023+ have O/E **0.59 (0.36–0.90)**; women diagnosed in 2023+ have **1.35 (1.11–1.62)**. A system-wide coding change is therefore unlikely to be the main cause.
  4. **The "C4 performs worse without fibroids" claim is not statistically established:** Δ 0.111 (−0.003 to 0.225), interaction P 0.40.
- **[I]** Hypertension, CAD and diabetes (master-file flags) decline 6–8% per diagnosis year after adjustment, while CKD (dated from raw extracts) does not. That pattern fits post-anchor code accumulation, which would make hypertension and CAD less "baseline-valid" than the post-audit assumed. It is not proven.
- **[J] Decision: pause model development; freeze C4 only as a provisional research candidate; do not publish or promote the nine-predictor model.** Restart when dated predictor records and chart validation are available.

## 2. Numerical reconciliation

All figures are at 3 years on the baseline cohort with outcome O2 (acute-coded), unless stated otherwise.

| Claim | Source A | Source B | Explanation | Verified value to report | Remaining uncertainty |
|---|---|---|---|---|---|
| Baseline-only temporal C 0.644 vs C4 0.663 | P06b, C3 row | P06b, C4 row | **Different models**, same patients, outcome, horizon and estimator. C3 = C4 + `hbp_pre` + dyslipidemia + smoking. | C3 0.644 (0.585–0.703); C4 0.663 (0.608–0.718). Paired C4 − C3 = 0.018 (−0.007 to 0.040) (V01c). | No evidence that either is better |
| Nine-predictor 0.711 vs "≈0.70" | P06b (all 2023+) | P08a / V04a (subgroups: 0.70 without fibroids, 0.72 with); P05c (development, corrected 0.704) | Different quantities: subgroup or optimism-corrected, not a different estimate of the same thing | Temporal 0.711 (0.662–0.760); development corrected 0.704 | — |
| C4 corrected 0.661 | P05c | V06a re-run | Identical | 0.661 (apparent 0.667, 0.628–0.706). No valid CI exists for an optimism-corrected C here. | — |
| Ridge corrected 0.652, ΔC +0.022 | P05c, P05d | V01c | Point values verified; the P05d interval [0.014–0.027] is a fitting-variability interval, not a CI | 0.652; ΔC vs C1 +0.022, no valid CI computed | Proper CI likely about ±0.025 (by analogy with C3: 0.028, 0.003–0.054) |
| C4 vs conventional temporal +0.045 (0.004–0.093) | P06c (500 resamples) | V01c (1,000 resamples) | Same estimand; the CI moves slightly with resample count | **+0.045 (0.003–0.089)**. But vs C1min (no noise predictors): **+0.022 (−0.010 to 0.056)**. | The comparator lacks diabetes and blood pressure |
| Nine minus baseline-only, "+0.05 to +0.07" | P05d, P06c (vs C3) | V01c (vs C3 and C4) | The range was computed against C3 | vs C4: dev 0.046 (0.020–0.078); temporal **0.048 (0.007–0.093)**. vs C3: 0.049 dev; 0.067 temporal. | Signal vs leakage [B] |
| IPA "≈1%" | P05c (development corrected) | V01b (temporal) | Development vs temporal | Development 0.99–1.00% (C3, C4); **temporal 0.40–0.59%** for baseline-valid models; C2 1.11% | — |
| O/E 1.44–1.56; "underpredicts by ≈45%" | A09a (audit: development sample, O1) | P06b / V01b (post-audit: baseline cohort, O2) | Two estimators: expected counts 1.44–1.49, KM over mean 1.56–1.65 | **Expected counts 1.45 (1.20–1.74) for C4.** Wording: observed ≈1.45 × predicted, i.e. predictions ≈31% too low, not "45% underprediction". | — |
| 345 vs 322 events | A01, A03 (development sample, O1) | P02a, P03a (baseline cohort, O2) | Different population **and** different outcome: 345 → minus procedure-anchored (<11) → plus BMI-missing index-anchored (<11) → 339 (O1) → minus imaging-only (<11) and history/sequela-only (<11) → 322 (V01a) | Both correct for their definitions | — |
| Median follow-up 2.95 | A04a (development sample) | V01d (baseline cohort 2.82) | Different populations | 2.95 (development sample); 2.82 (baseline cohort) | — |
| 28% to 5 years | A04b | P04a | Different populations | 28.1% / 26.6% | — |
| "4 of 1,400 fits" | P05 note | V06b | Verified, but from **2** resamples (primary outcome: 1; ischemic: 2) | 4 fits; quasi-separation; optimism change ≤0.0001 | — |

## 3. Verified model comparisons

**Population:** baseline cohort (24,550 women, no missing predictor values; verified by `complete.cases`). **Outcome:** O2. **Horizon:** 3 years.

**Methods:**
- Censoring at death, TIA, last encounter, and the administrative end of follow-up.
- Ties: Efron in `coxph`; Breslow baseline in P05.
- C is Harrell's, with follow-up truncated at 3 years and ties in the predictor counted as half.
- CIs are paired percentile bootstraps (1,000 resamples) unless marked.

| Comparison | Development, apparent paired CI | Temporal paired CI (2023+) | Defensibility |
|---|---|---|---|
| C2 nine − C4 | 0.046 (0.020–0.078) | 0.048 (0.007–0.093) | **[V]** gap real; whether signal or leakage is **[B]** |
| C4 − C1 conventional | 0.031 (0.005–0.057) | 0.045 (0.003–0.089) | Comparator handicapped by noise predictors |
| **C4 − C1min (age, hypertension, CAD, CKD) = increment from Hb deficit** | **0.026 (0.003–0.049)** | **0.022 (−0.010 to 0.056)** | **Most defensible comparison**: a nested, fair one. Exploratory: C1min was defined in this verification. |
| C3 − C1 | 0.028 (0.003–0.054) | 0.027 (−0.005 to 0.066) | — |
| C4 − C3 | 0.003 (−0.004 to 0.010) | 0.018 (−0.007 to 0.040) | No difference; preferring C4 is parsimony **[J]**, not performance |

**Answers:**
- **Strongest comparison:** C4 versus C1min, nested and temporal. The development gain is real but small. The temporal gain is not statistically established.
- **Exploratory:** C1min; everything involving C4 plus migraine; the subgroup contrasts; the fragility indices.
- **Potentially optimistic:**
  - All "development" figures, because the post-audit predictor sets were chosen after the audit had seen these data (P01e timeline). No data-driven selection was done inside the post-audit.
  - The nine-predictor model in particular, whose predictor choice used the full cohort.
- **Incremental value over a *complete* conventional model: not demonstrable.** Diabetes and blood pressure are absent or undatable **[B]**. Over the incomplete comparator the evidence is weak (above).
- **[V]** Overall accuracy is low: temporal IPA 0.4–0.6% for the baseline-valid models.

## 4. Sparse predictors

| Predictor | Women / events | Evidence | Pre-specified? | [J] Recommended handling |
|---|---|---|---|---|
| Estrogen (as coded) | 338 / 14 (development sample, O1) | Loses significance after 2 influential events (A06k); window extends after the anchor; some products misclassified | Influence analysis exploratory (added in the audit) | **Remove** until dated prescriptions exist. Instability does not show the absence of an association. |
| Hemoglobinopathy (any date) | 263 / 14 | Fragile at 5 events; 7 of 14 strokes first coded after the anchor | Exploratory | **Replace** by the pre-anchor definition |
| Hemoglobinopathy before the anchor | 197 / <11 (baseline cohort, O2) | HR 2.06 (0.92–4.64). 4 of 1,400 bootstrap fits diverged (quasi-separation in 2 resamples, V06b); excluding them changes optimism by ≤0.0001 (V06c). | Bootstrap pre-specified; diagnosis new | **Defer**: too few events to estimate. Could enter later with penalization if events accrue. |
| Smoking at baseline | 3,046 positive; 79% unknown | Same sign in 63% of resamples; HR 0.96 (0.69–1.33) | — | **Remove**: not measurable at baseline |
| Dyslipidemia | — | Same sign in 75% of resamples; HR 1.08 (0.84–1.37) | — | Remove (adds noise; C3 ≈ C4) |
| Age point (calculator) | — | Unadjusted 3-year KM: 1.05–1.09% at 30–44 (+1 point) vs 1.59% at 45–54 (0 points) (P07b). Adjusted spline HR ≈1.0 across 35–55 (A06c). | — | Contradiction **[V]** by both unadjusted and adjusted evidence. Remove the point. |

## 5. Temporal calibration drift

| Claim | Status | Evidence |
|---|---|---|
| All models underpredict in 2023+ | **[V]** | O/E 1.45–1.49 by expected counts, for C1–C4 (V01b) |
| O/E ≈1.44–1.56 | **[V]**, with the estimator and population made explicit | §2 |
| Excluding events within 30 days does not remove it | **[V]** | O/E 1.59–1.63 (P06h); 30-day rate stable (P06i) |
| A 2023–24 update does not fix 2025+ | **[V]** at 1 year | 41 events, wide CIs (P06d) |
| Hypertension and diabetes coding falls by diagnosis year | **[V]**, plus CAD | OR per year 0.937, 0.928, 0.917; CKD 1.006 (V03e) |
| Cause undeterminable | **Partly corrected** | See below |

**Calendar period versus diagnosis cohort** (V03a, V03b). The model was developed with follow-up censored at 2022-12-31, so it never sees calendar-2023+ outcomes. Expected counts were left-truncated, and both cohorts were compared over the same 0–3 years since the anchor.
- Women diagnosed before 2023, followed during calendar 2023+: O/E **0.59 (0.36–0.90)**.
- Women diagnosed in 2023+: O/E **1.35 (1.11–1.62)**.
- Poisson model: period RR 0.87 (0.62–1.23); diagnosis-cohort RR **1.58 (1.09–2.29)**.

**[V]** The excess is attached to recently diagnosed women, not to calendar time. Outcome composition by calendar period does not suggest degraded coding: corroborated events are 32% in 2023+ versus 23% before (V03d).

**[I]** The leading candidate mechanisms are:
- **(a) Fewer accumulated comorbidity codes in recent cohorts**, so their predictors understate risk. The decline in master-file flags alongside a flat, verified-dated CKD flag fits this.
- **(b) A real case-mix change not captured by the predictors:** recent cohorts are younger, with fewer fibroids.

**[U]** Neither is established. **[B]** Telling them apart needs dated diagnosis records.

**Why the 2023–24 recalibration failed for 2025+ [I]:** a cohort-specific excess is not transferable between cohorts. In addition, the 1-year 2023–24 O/E (≈1.0) and the 2025+ O/E (≈1.5) rest on 32 and 41 events. The result establishes neither a time-varying baseline hazard nor a change in predictor effects. The revision tests (P06g; P 0.046 for C3) are weak and rest on 117 events.

## 6. Outcome and horizon

- **[V] Outcome definitions by report:**
  - Audit and the published model: O1, all coded strokes excluding TIA (345 events in the development sample).
  - Post-audit P05, P06 and P08: O2, acute-coded (322 events in the baseline cohort).
  - Each report labels its definition.
- **[V] The "89 confirmed" events:** acute-coded events with a positive-imaging earliest date or a first NIHSS flowsheet record within ±30 days of the coded date.
  - This is **record-based corroboration, not adjudication.** "Corrected" in the ledger means the wording "confirmed" was wrong.
  - Positive imaging includes chronic findings, NIHSS is also recorded for mimics, and missing corroboration does not exclude a stroke.
- **[V]** 115 of 345 events had brain imaging reviewed and negative (262 had imaging reviewed).
- **[V] 3-year risk:**
  - The estimand is 1 − Kaplan–Meier at 3 years: cause-specific, with death and TIA censored.
  - KM risk is 1.27% (1.10–1.44); Aalen–Johansen with death as a competing event gives 1.266%.
  - Pattern-mixture analysis, assuming a 1.5–3× hazard after loss to follow-up, gives up to 1.36%.
  - The pattern-mixture assumption is untestable.
- **[V] 5-year risk is not validated.** No 2023+ woman reached 5 years, and only 57% of 2018–2022 anchors did. 5-year output should not be presented as validated, although the live calculator still shows it.

## 7. Fibroid subgroup (V04)

- **Definitions:** fibroids = Codebook flag `fibroids == 1` (any coded leiomyoma in the 2018+ extract); no fibroids = endometriosis and/or adenomyosis only.
- **Data:** temporal (2023+) set, 3 years, O2. 7,897 women with 68 events in the fibroid group; 5,258 women with 49 events in the no-fibroid group.

| Result | Value |
|---|---|
| C4 temporal C | 0.70 (0.63–0.77) with fibroids vs 0.59 (0.50–0.68) without |
| Difference | 0.111 (−0.003 to 0.225) |
| Linear-predictor × fibroid interaction | ratio 1.26 (0.74–2.15), P 0.40 |
| C2 temporal C | 0.72 vs 0.70; difference 0.013 (−0.087 to 0.111) |

- **[V]** The point estimates match the reports. **The difference is not statistically established;** no formal test had been run before this verification.
- **[I]** Adding as-coded migraine to C4 raises C in women without fibroids by 0.128 (0.054–0.199). This fits the migraine interpretation, but the variable is the undatable one. **[U]** whether that is signal or leakage.
- **[J]** These data support only a model for the **mixed, CBC-tested benign uterine disease cohort** (statement C), not a fibroid-specific model.

## 8. Calculator audit (V05b, A05, P07)

**Reproducibility [V]:**
- Every displayed risk is reproducible from the underlying model: R, the JSON spec, the JavaScript and the public page agree to display rounding (A05).
- Points and coefficients are internally consistent.
- The **material discrepancy is between the two displayed risks**: 3.84% from the points table versus 5.30% from the exact model for the worked example. They come from different models.

**Inputs and range [V]:**
- All 9 items are required. There is no "unknown" option, so a user whose information is missing must answer "No".
- Age is limited to 18–60 and Hb to 3–20 g/dL.
- There is **no cap or suppression**: the all-positive profile gives 99.99%. There is a warning above the 99th percentile (9.78%).
- 95 of 128 yes/no combinations have fewer than 20 development women or were never observed (P07a).
- The points table rows stop at 13, but the total can reach 19.

**Warnings [V]:**
- "Research use only · not externally validated" is visible as a tag (line 96) and in copied summaries (line 230).
- The Limitations tab states that risk is "likely underestimated in recent patients" (line 256). **The post-audit implication that no research-only warning existed is corrected.**
- **Absent:** that several predictors may postdate the prediction time; that 5-year risk is not validated; why two different risks are shown.

**Recommendation [J]:** a warning alone is not adequate. The headline number is an unvalidated 5-year risk from undatable predictors, presented without uncertainty. A tag does not prevent a reader from treating that number as a clinical estimate. I recommend **taking the page offline** until a re-specified model exists. If the PI keeps it online, the banner above the inputs needs the three missing caveats. **This decision rests with the PI.** No intervention (antiplatelet, anticoagulation, imaging, hormone changes) should be tied to the score.

## 9. Missing-data decision matrix

| Question | Existing evidence | Resolvable now? | Exact analysis or data required | Decision impact |
|---|---|---|---|---|
| Dated migraine | none (Codebook: any date) | **Requires external records** | Dated problem-list and diagnosis extract (request item 1) | Decides whether ≈0.05 C is real |
| Dated diabetes | none; prevalence falls with anchor year | Requires external records | Item 1, or dated HbA1c (item 7) | Comparator completeness; leakage |
| Dated estrogen exposure | window with no end; misclassified products | Requires external records | Prescriptions with dates (item 2) | Whether estrogen can enter |
| Hemoglobinopathy dating | dated extracts present | **Resolved** | — | `hbp_pre` built; sparse |
| Smoking dating | one chosen statement per woman | **Partially resolvable**; built, 79% unknown | All dated statements (item 6) | Smoking stays out |
| Hypertension and CAD dating | Codebook claim only; trend against the dated CKD control | **Requires external records** (new concern) | Item 1 | Whether C4 itself is baseline-valid |
| Blood pressure / complete comparator | absent | Requires external records | Vitals (item 5) | Real incremental value |
| Raw stroke extract; chart review | sample frame ready; no charts | Requires chart access | Items 3–4; protocol ready | Outcome validity (abandonment criterion b) |
| Cause of drift | cohort effect shown; mechanism unknown | **Partially resolvable** (done) | Items 1, 3, plus follow-up | Whether absolute risk is reportable |
| 5-year validation | maximum 3.7 years in 2023+ | **Not answerable now**; needs time (≈2028 for 2023 anchors) | Continued follow-up | Horizon |
| External validation | none | Requires an independent cohort | Item 10 | Any use beyond research |
| Decision threshold | none | **Not answerable from this dataset** (clinical) | Clinical use-case definition | Clinical utility |

## 10. Final recommendation and unresolved questions

**Evidence-based [V]:**
- The nine-predictor model's extra discrimination over C4 (≈0.05) depends entirely on undatable predictors.
- C4's temporal gain over a fair conventional comparator is not established (+0.022, CI includes 0).
- Calibration fails in recent diagnosis cohorts.
- 5-year output is unvalidated.
- No outcome PPV exists.

**Recommendations [J]:**
1. **Model to freeze:** C4 (age, Hb deficit, hypertension, CAD, CKD; O2; 3 years), as a *provisional research candidate*. It is preferred over C3 for parsimony and the absence of noise predictors, not for better performance; their difference is 0.003 in development and 0.018 temporally, both with CIs including 0.
2. **Reintroducing migraine, diabetes or estrogen** requires *all* of:
   - (a) dated records showing the predictor is recorded at or before the anchor
   - (b) a pre-registered refit
   - (c) temporal ΔC over C4 with a paired CI excluding 0, and stable HRs
   - (d) no worsening of calibration
3. **Abandonment criteria (post-audit):**
   - (a) dated predictors add <0.02 and C4 <0.65 against a complete comparator: **not evaluable [B]**
   - (b) PPV <0.75: **not evaluable [B]**
   - (c) external C <0.60 or O/E outside 0.7–1.4: **not evaluable [B]**
   - (d) no definable use case: **currently unmet but undecided**

   **None has been met.** Note that C4's temporal C (0.663) sits close to the 0.65 bound, and its gain over a fair comparator is uncertain.
4. **Before an external validation study:**
   - dated predictors re-derived
   - chart-review PPV ≥0.85
   - one frozen, pre-registered model with a 3-year horizon
   - a data-use agreement for an external cohort with the same variable definitions
5. **Before any prospective clinical evaluation:** all of item 4, plus successful external validation (calibration-in-the-large within about 0.8–1.25 after at most one pre-specified recalibration), plus a defined clinical decision and threshold.
6. **Manuscript:**
   - A nine-predictor development and validation paper is **not defensible now**, even with corrections.
   - Pausing model development pending outcome validation and dated extracts is warranted.
   - A corrected methods-and-limitations report would be defensible: Hb deficit as an associated predictor in a CBC-tested cohort, the leakage findings, and temporal drift. That is the PI's choice.

**Unresolved [U]:**
- leakage versus signal in migraine, diabetes and estrogen
- whether the hypertension and CAD flags are truly baseline
- the mechanism of the cohort-level underprediction
- outcome PPV
- any clinical decision
