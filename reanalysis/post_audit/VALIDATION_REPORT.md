# Validation report (post-audit)

**Validation levels (labeled throughout):**
- **Development:** apparent, in-sample.
- **Internal:** full-process bootstrap.
- **Temporal:** chronological split at the same center.
- **External:** none. No external validation exists.

## 1. Internal validation: full process

**Script:** `post_audit/R/P05_candidate_models.R` (P05c, P05d, P05e).

**Procedure:**
- 200 bootstrap resamples at the patient level (seed 20261022).
- For each pre-specified candidate, everything data-dependent is repeated in each resample: Cox fitting, spline knot placement, ridge λ by 5-fold cross-validation, and the Breslow baseline.
- No predictor selection exists in the post-audit candidates; predictor eligibility was fixed by date provenance, not by data. Missing-data handling is complete case on Hb; the BMI requirement was removed.

**Results:**

| Measure | Value |
|---|---|
| C optimism | 0.006–0.013 |
| Monte-Carlo SE of optimism | 0.0012–0.0014 |
| Corrected C: C4 | 0.661 |
| Corrected C: C3 | 0.652 |
| Corrected C: C2 (as coded) | 0.704 |
| Corrected C: C1 | 0.625 |
| Calibration slope (C3/C4) | 0.998 / 1.017 |
| Calibration-in-the-large optimism | negligible (corrected O/E 0.997–0.999) |

**Trade-off:** 200 resamples × 7 models took about 25 minutes, and the Monte-Carlo SE is small relative to the differences between models. More resamples would not change any conclusion.

**For the original nine-item model:** with predictor selection repeated, corrected C is 0.697 (A08a). Only 13.5% of resamples reproduced the same nine predictors.

## 2. Temporal validation

**Script:** `post_audit/R/P06_temporal_drift.R`.

**Sizes (verified, P06a):**

| Anchor period | Women | Events | Events within 3 years | Maximum follow-up |
|---|---|---|---|---|
| 2018–2022 | 11,395 | 205 | 118 | 8.7 years |
| 2023–2024 | 6,354 | 75 | 75 | 3.7 years |
| 2025–2026 | 6,801 | 42 | 42 | 1.7 years |

**Primary (pre-specified):** develop on 2018–2022, evaluate on 2023+ at 3 years.
- Predictor sets were fixed by provenance before this split was examined. The split itself was used in the audit, so it is not new to the analyst. That is acknowledged; no cut-off was tuned.
- Results are in `CANDIDATE_MODEL_COMPARISON.md`.
- C4: C 0.663 (0.608–0.718); slope 0.88 (0.66–1.09); O/E 1.45 (1.20–1.74); ICI 0.61 percentage points.

**Secondary (pre-specified three-way design):** develop 2018–2022, recalibrate the baseline only on 2023–2024, evaluate on 2025+.
- This can only be evaluated at **1 year**, because 2025+ anchors have at most 1.7 years of follow-up. The proposed three-way design at 3 years is **not feasible** with current data.
- Results (P06d):

  | Measure | C3 | C2 |
  |---|---|---|
  | O/E in 2023–24 | 0.98 | 1.00 |
  | O/E in 2025+ before the update | 1.51 (1.08–2.04) | 1.57 |
  | O/E in 2025+ after the update | 1.54 | 1.57 |
  | 1-year C in 2025+ (41 events) | 0.57 | 0.69 |

- **A baseline update did not help:** the drift is not stable from one year to the next.

## 3. Calibration drift: decomposition

| Candidate source | Evidence | Conclusion |
|---|---|---|
| Baseline-risk drift | O/E about 1.45 for **all** models, including conventional (P06b). First-year rate by anchor year: 3.2–10.0 per 1,000 person-years, non-monotone (P06e). | Present, and not monotone |
| Early events / time zero | O/E rises to 1.59–1.63 when events within 30 days are excluded; the 30-day rate is stable (18–23 per 1,000 person-years) (P06h, P06i) | Not the explanation |
| Predictor-effect change | Revision test P 0.046 (C3) and 0.004 (C2); HR ratios: CAD 0.49–0.52, hypertension 1.47, hemoglobinopathy 1.85 (P06g) | Some evidence. Effects are estimated from about 117 validation events, so imprecise. |
| Predictor measurement | Hypertension prevalence 22.3% → 11.7% and diabetes 28.0% → 16.1% by anchor year; estrogen 2.2% → 0.8–1.0% (P06e) | Fits ascertainment time accumulating after the anchor, or case-mix change; not distinguishable |
| Outcome coding | Acute-code share 94.0% vs 96.7%; ischemic 85.3% vs 90.1%; corroborated 30.3% vs 27.3% (P06f) | Small shifts |
| Case mix and utilization | Mean age 42.4 → 40.0; fibroids 69.4% → 52.9% (2026); CBC testing rate 62–70%; median CBC days before anchor 1 → 2 (P06e) | Moderate shift towards younger women and endometriosis |

**Conclusion:** the drift is real and affects every model. Its cause cannot be identified from these data. The calculator must **not** be recalibrated with the data used to evaluate that recalibration (none was).

## 4. Prediction horizon

From P04a (baseline cohort):

| Horizon | At risk | Events | KM risk | Pattern-mixture range (δ = 1–3) |
|---|---|---|---|---|
| 1 year | 19,184 (78.1%) | 134 | 0.59% | 0.59–0.61% |
| 3 years | 11,724 (47.8%); 90.2% of 2018–2022 anchors | 235 | 1.27% (1.10–1.44) | 1.27–1.36% |
| 5 years | 6,538 (26.6%); 57.4% of 2018–2022 and 0% of 2023+ anchors | 287 | 1.86% | 1.86–2.03% |

- Competing death changes these by less than 0.01 percentage points (P04c).
- **Recommended primary horizon: 3 years.** It is directly observed in 90% of the development period, and it can be temporally evaluated, though only 11% of 2023+ women had reached 3 years.
- **5 years is not validated and should not be presented.**
- **1 year** is the only horizon testable in 2025+, but it is dominated by early events.

## 5. What has not been done

- **External validation:** none (MISSING_DATA_REQUEST item 10).
- **Validation with dated migraine, diabetes and estrogen:** not possible (items 1, 2).
- **MI-in-bootstrap:** designed but not run (A11).
