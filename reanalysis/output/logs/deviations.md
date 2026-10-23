# Deviations from the original analysis and from the SAP

| Item | Reason |
|---|---|
| Time-updated reproduction: same-day multiple Hb results | R uses the daily median (Codebook rule); the original kept every row, so the value in force differs on 11,326 woman-days whose same-day results span grades (5,946 involve Hb <8). This explains the severe-grade discrepancy (6 events). |
| Primary cohort | 1 woman with hormonal_type = 9 (missing) excluded from models (complete case for this variable). |
| Sparse binary covariates | HIV before index had 0 incident strokes and IBD 4; they are pooled into one 'IBD or HIV' indicator (analogous to the pre-specified hormonal-type rule) instead of being entered separately. |
| Imputation model | mice removed the Nelson-Aalen cumulative hazard as collinear (logged event); BMI imputation used the event indicator and all covariates. BMI missing 3.5%; complete-case results (S1) are reported alongside. |
| Uterine diagnosis group in the primary model | Included as a confounder per the DAG (output/dag.md): it affects anemia through fibroid bleeding, and endometriosis may itself be associated with stroke. The original revision dropped it at the co-author's request; adding it changed the original estimates by <2% (reproduction table). Heavy bleeding stays out of the model (exposure determinant, S6). |
| Primary estimand | Changed from prevalent-or-incident stroke/TIA (cross-sectional, Hb +-3 y in either direction) to incident stroke excluding TIA after index, with Hb on/before index (PI decision 2026-10-22). |
| Hb source | Primary exposure uses laboratory CBC Hb only; the Codebook hgb also included point-of-care Hb (S4). |
| Sensitivity analyses S2-S5, S10, S11 | Complete case for BMI (3.5% missing) instead of MI, to limit computation; S1 shows MI and complete-case estimates agree. S6-S9 reuse the primary imputations. |
| Anchored design (PI request 2026-10-23) | Anchor = earliest evidence of uterine pathology: first coded diagnosis (2018 onward, IRB window) or an earlier dated myomectomy, uterine artery embolisation, hysterectomy or endometrial ablation. 'Other gynaecological surgery' (tubal, ovarian, D&C, IUD removal) not used as evidence. TIA severe-grade cells in the nested case-control analysis have 0 cases (not estimable). |
| Points-score constant | Pre-specified constant (5 years of age) unusable: adjusted age coefficient -0.0104 per year (HR 0.99). Constant set to half the hypertension coefficient (hypertension = 2 points). |
| Points-score constant | Pre-specified constant (5 years of age) unusable: adjusted age coefficient -0.0104 per year (HR 0.99). Constant set to half the hypertension coefficient (hypertension = 2 points). |
| Anchor plausibility (2026-10-23) | Procedure dates implying age <18 at the anchor (2 women, no strokes) were ignored; the anchor falls back to the first coded diagnosis. |
| Points-score constant | Pre-specified constant (5 years of age) unusable: adjusted age coefficient -0.0105 per year (HR 0.99). Constant set to half the hypertension coefficient (hypertension = 2 points). |
| Calculator model (PI decision 2026-10-23) | Predictors that carry no points (BMI, dyslipidemia, smoking, endometriosis, adenomyosis) removed from the calculator; exact-risk model refitted with the nine MED-HACK items (R/16_calculator_model.R). Optimism-corrected C 0.701 vs 0.696 for the 14-predictor model; temporal C 0.695 vs 0.665. |
