# Candidate-model comparison

**Design (identical for every model unless labeled otherwise):**
- **Cohort:** the baseline cohort, 24,550 CBC-tested women aged 18 or over who were anchored at their first coded diagnosis.
- **Outcome:** acute-coded stroke (O2), 322 events, 235 within 3 years.
- **Estimand:** the 3-year risk of a first acute-coded stroke after the gynecologic diagnosis, treating death and TIA as censoring (cause-specific). Competing death changes 3-year risk by less than 0.01 percentage points (P04c).
- **Model:** Cox; baseline hazard by Breslow.
- **Pre-specification:** all predictor sets were fixed before any post-audit model was fitted (`post_audit/R/P00_setup.R`, `CAND`), and no data-driven selection was used. The same cohort, outcome and horizon apply throughout. **C6** uses a different outcome (ischemic) and is not directly comparable.

**Sources:** `post_audit/output/P05_candidate_models.md` (CSV copies in `output/tables/`, git-ignored) and `P06_temporal_drift.md`.

## Candidates

| Model | Predictors | Parameters | Events per parameter | Rationale |
|---|---|---|---|---|
| C1 conventional (baseline-valid) | age, hypertension, CAD, CKD, dyslipidemia, smoking at baseline | 6 | 53.7 | Comparator. Lacks diabetes (undatable) and blood pressure (absent). |
| C2 current nine (as coded) | migraine, estrogen, diabetes, hypertension, hemoglobinopathy (any date), Hb deficit, age, CAD, CKD | 9 | 35.8 | Benchmark. **Not defensible**, because 4 predictors cannot be dated. |
| C3 baseline-only | age, Hb deficit, hypertension, CAD, CKD, hemoglobinopathy before the anchor, dyslipidemia, smoking at baseline | 8 | 40.2 | Every predictor that is dated at or before the anchor |
| C3s | C3 with an age spline (3 knots, Harrell default) | 9 | 35.8 | Pre-specified functional-form check |
| C4 reduced baseline | age, Hb deficit, hypertension, CAD, CKD | 5 | 64.4 | Clinically justified: established vascular factors plus Hb, all dated |
| C5 ridge | C3, ridge-penalized (λ by cross-validation) | 8 | — | Addresses sparse predictors |
| C6 ischemic (secondary) | C3 | 8 | 35.6 (285 events) | Ischemic outcome only |

## Results (3 years)

| Model | Apparent C (95% CI) | **Optimism-corrected C** | Temporal C, 2023+ (95% CI) | Temporal O/E (95% CI) | Temporal slope | Corrected IPA | ΔC vs C1, development | ΔC vs C1, temporal |
|---|---|---|---|---|---|---|---|---|
| C1 | 0.637 (0.595–0.678) | 0.625 | 0.617 (0.554–0.680) | 1.45 (1.20–1.74) | 0.86 | 0.92% | — | — |
| C2 | 0.713 (0.678–0.747) | **0.704** | 0.711 (0.662–0.760) | 1.49 (1.23–1.79) | 0.95 | 1.14% | +0.072 (0.062–0.084) | +0.094 (0.042–0.145) |
| C3 | 0.664 (0.624–0.703) | 0.652 | 0.644 (0.585–0.703) | 1.45 (1.20–1.74) | 0.85 | 0.99% | +0.022 (0.012–0.028) | +0.028 (−0.008 to 0.066) |
| C3s | 0.664 | 0.651 | — | — | — | 0.97% | ΔC vs C3: −0.001 | — |
| **C4** | 0.667 (0.628–0.706) | **0.661** | **0.663 (0.608–0.718)** | 1.45 (1.20–1.74) | 0.88 | 1.00% | +0.027 (0.017–0.037) | +0.045 (0.004–0.093) |
| C5 | 0.663 | 0.652 | — | — | — | 1.01% | +0.022 (0.014–0.027) | — |
| C6 (ischemic) | 0.678 | 0.665 | — | — | — | 0.92% | not comparable | — |

**How the columns were computed:**
- Optimism-corrected values come from a full-process bootstrap (200 resamples). Every resample repeats the fitting, the spline knots and the ridge cross-validation. Monte-Carlo SE of the optimism is 0.0012–0.0014.
- Calibration slopes in development: 0.99–1.06 (bootstrap intervals about 0.89–1.20), so no shrinkage is needed (P05c).
- Temporal results: models developed on 2018–2022 anchors and evaluated on 2023+ (P06b, P06c).
- Development ΔC uses bootstrap-fitted models evaluated on the original data (P05d).

## Coefficients (P05a; unshrunk HR, 95% CI)

| Predictor | C3 | C4 |
|---|---|---|
| CAD | 4.35 (2.98–6.36) | 4.39 (3.02–6.38) |
| Hypertension | 2.17 (1.70–2.78) | 2.21 (1.73–2.81) |
| CKD | 2.02 (1.37–2.99) | 2.04 (1.38–3.02) |
| Hb deficit, per g/dL below 13 | 1.13 (1.05–1.22) | 1.14 (1.05–1.22) |
| Hemoglobinopathy before the anchor | 2.06 (0.92–4.64) | — |
| Age, per year | 0.99 (0.98–1.00) | 0.99 (0.98–1.00) |
| Dyslipidemia | 1.08 (0.84–1.37); same sign in 75% of resamples | — |
| Smoking at baseline | 0.96 (0.69–1.33); same sign in 63% of resamples | — |

## Which measures matter most for the decision

1. **Temporal calibration-in-the-large.** Every model, including the conventional one, underpredicts by about 45% in 2023+. A baseline-only recalibration estimated on 2023–24 did not fix 2025+ (P06d). Absolute risks from **any** candidate are therefore not yet reportable as accurate.
2. **Temporal ΔC over the conventional comparator.**
   - The only baseline-valid gain with a CI excluding 0 is C4's +0.045 (0.004–0.093). (This audit did not decompose which predictor drives it.)
   - The comparator is incomplete: it lacks diabetes and blood pressure.
3. **Overall accuracy is low for every model:** IPA about 1% at 3 years. Because the outcome is rare, these models rank women moderately well but explain little of the individual uncertainty.

## Interpretation

- **The current nine-item model's advantage comes from predictors that cannot be dated.** Its ΔC over the baseline-only model is 0.050 (0.041–0.063) in development and 0.066 (0.025–0.116) temporally. That advantage is either real prognostic information (migraine, estrogen, diabetes) or leakage; with the present data the two cannot be separated. *Interpretation, not a computed decomposition:* migraine is probably the main contributor. It is the most prevalent undatable item with a large HR, and the gap between C2 and C4 is widest in women without fibroids (P08a).
- **Among defensible models,** C4 (5 parameters) performs as well as or better than C3 (8 parameters). Smoking and dyslipidemia add only noise, and hemoglobinopathy before the anchor is too sparse to estimate reliably.
- **Ridge penalization changes nothing material:** λ is small, and 40–64 events per parameter is adequate.
- **The spline for age adds nothing.**
- **Recommended next model:** C4 (age, Hb deficit, hypertension, CAD, CKD), plus migraine, diabetes and estrogen **if and only if** they can be dated at or before the anchor.
- **Subgroups:** C4 discriminates poorly in women without fibroids (temporal C 0.59, 0.50–0.68) (P08a).
