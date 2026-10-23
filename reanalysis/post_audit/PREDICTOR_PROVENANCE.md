# Predictor provenance

**Dates:**
- **Anchor:** the prediction time, the earliest evidence of uterine disease (`R/13_anchor_design.R:26`).
- **Index date:** the first coded uterine diagnosis from 2018 onward. It equals the anchor for 24,550 of the 24,886 CBC-tested women; the other 336 have an earlier procedure anchor.
- **Stroke date:** the earliest code of any kind (Codebook `stroke_date`).
- **Ascertainment date:** the date of the record that defines the predictor.

**Reconstructability classes:**

| Class | Meaning |
|---|---|
| 1 | Fully reconstructable from dated data in the project |
| 2 | Partially reconstructable, with a stated limitation |
| 3 | Not reconstructable without data that are not in the project |

Counts are from `post_audit/output/P02_baseline_reconstruction.md`, `P06_temporal_drift.md` (P06e) and `audit/output/A02_predictors_missingness.md`.

## Model predictors

### Age
- **Coding:** years at the anchor, linear (`R/13:29`).
- **Source:** `age_index` (completed years at index_date, master file), shifted by the anchor offset.
- **Date and timing:** the anchor itself, so it **cannot follow the stroke**.
- **Missing:** 0.
- **Type:** demographic.
- **Class:** 1. Precision is ±1 year, because the master file gives completed years only; date of birth is deliberately not used.
- **Handling:** keep. Its functional form is unresolved (P05: linear versus spline).

### Hb deficit
- **Coding:** max(13 − Hb, 0) g/dL.
- **Source:** Lab Tests (2), CBC hemoglobin; the latest value in the 1,095 days up to and including the anchor (`R/13:70`).
- **Date and timing:** at or before the anchor, so it cannot follow the stroke.
- **Missing:** 35.1% of the anchored cohort have no CBC (selection).
- **Type:** laboratory value. Testing reflects clinical indication.
- **Class:** 1.
- **Handling:** keep. State that the model applies to women with a CBC.

### Hypertension
- **Coding:** yes/no; master `htn`.
- **Source:** the Codebook says "documented on or before index_date". The source extract is **not in the project**.
- **Date and timing:** at or before the index date. That equals the anchor for index-anchored women, but is after the anchor for the 336 procedure-anchored women.
- **Can it follow the stroke?** Only in procedure-anchored women.
- **Missing:** none recorded; absence of a code is coded "no".
- **Type:** pre-existing condition.
- **Class:** 2. The date rule is asserted by the Codebook but can't be verified. Prevalence falls from 22.3% (2018 anchors) to 11.7% (2026 anchors), which fits either case-mix change or ascertainment over time.
- **Handling:** keep in index-anchored women; verify with a dated extract.

### CAD
- **Coding:** yes/no; master `cad`.
- **Source:** Codebook: dated CAD diagnosis and problem-list extracts, on or before index_date. The extract is not in the project.
- **Date and timing:** same as hypertension.
- **Missing:** as for hypertension.
- **Type:** pre-existing condition.
- **Class:** 2.
- **Handling:** keep in index-anchored women; verify.

### CKD
- **Coding:** yes/no; `ckd_pre`.
- **Source:** `dx/` extracts 22–31 (CKD and ESRD code rules, `R/03`), first code on or before the anchor (`R/13:66`).
- **Date and timing:** at or before the anchor, so it cannot follow the stroke.
- **Missing:** as for hypertension.
- **Type:** pre-existing condition.
- **Class:** 1.
- **Handling:** keep.

### Diabetes
- **Coding:** yes/no; master `dm`.
- **Source:** the Codebook gives **no timing**. The source extract is not in the project.
- **Date and timing:** unknown, so it **may follow the anchor and the stroke**.
- **Missing:** as for hypertension.
- **Type:** condition; possibly partly ascertained later.
- **Class:** 3.
  - The look-ahead check gives OR 1.72 per doubling of follow-up, versus 1.20–1.31 for the date-verified negative controls.
  - Prevalence falls from 28.0% to 16.1% by anchor year.
  - Prevalence is 23% in a cohort with a mean age of 41.
- **Handling:** **remove** until a dated extract is available.

### Migraine
- **Coding:** yes/no; master `migraine` (codes 1, 2 and 9 count as yes).
- **Source:** Codebook: problem list plus migraine diagnosis extracts, with **no date restriction**. Not in the project.
- **Date and timing:** unknown, so it **may follow the anchor and the stroke** (post-stroke headache codes, migrainous infarction G43.6).
- **Missing:** code 9 (aura not stated) in 9.9% of women, counted as yes.
- **Type:** condition; possibly partly ascertained later.
- **Class:** 3. Look-ahead OR 1.56.
- **Handling:** **remove** until a dated extract is available. This is the single most influential undatable predictor (P05, P06).

### Estrogen-containing therapy
- **Coding:** yes/no; `hormonal_type` 1 (combined OC) or 5 (menopausal HT).
- **Source:** Codebook: Administered Medications, window **index_date − 365 days onward, with no end date**. Not in the project; `meds/` holds file 9 only (iron and spironolactone).
- **Date and timing:** the window extends after the anchor by definition, so it **can follow the stroke**.
- **Missing:** under-captured, because only administered medications are recorded.
- **Type:** treatment exposure.
- **Class:** 3.
  - The medication list includes F-18 fluoroestradiol (a PET tracer) and IV conjugated estrogens (acute bleeding treatment).
  - 338 exposed women, 14 strokes; the association is fragile (2 events).
- **Handling:** **remove** until dated prescriptions are available.

### Hemoglobinopathy
- **Coding:** yes/no; any sickle cell, thalassemia or other hemoglobinopathy code.
- **Source:** `dx/` extracts 22–31.
- **Date and timing:** the published model uses codes of **any date** (`R/13:65`), so it **can follow the stroke**: 7 of 14 strokes in the development sample were first coded after the anchor. The reconstructed `hbp_pre` uses only codes on or before the anchor.
- **Missing:** as for hypertension.
- **Type:** congenital condition, but the *diagnosis* is time-dependent.
- **Class:** **1**. In the baseline cohort, 197 women are positive versus 262 as coded; 65 are reclassified, 6 of them with strokes.
- **Handling:** replace with `hbp_pre`. Prevalence is 0.8% in both Black and White women, so this flags diagnosed disease, not trait.

## Other candidates (14-predictor model)

### Smoking
- **Source:** `smoking_source` sheet: a single chosen statement, the latest on or before the index date, otherwise the **earliest after** it.
- **Date and timing:** the statement date; 3,268 of 8,391 documented statuses were stated after the anchor.
- **Can it follow the stroke?** Yes: 43 of 180 cases with a documented status.
- **Missing:** unknown for 65.0% as coded, and 79.2% at baseline.
- **Type:** behavior; documentation tracks contact with care.
- **Class:** 2. A baseline status can be built (`smoke_b`), but an earlier statement may exist that the sheet doesn't contain.
- **Handling:** not estimable. Leave it out of the model and report it as a known missing risk factor.

### BMI
- **Source:** master `bmi`, the closest value within ±12 months of the index date. No date in the project.
- **Can it follow the stroke?** Possibly, if the stroke happened within 12 months of the anchor.
- **Missing:** 8.4%.
- **Class:** 3.
- **Handling:** not used. Dropping the BMI requirement removes 880 needless exclusions.

### Dyslipidemia
- **Source:** `dyslipidemia_update` sheet: a code or an abnormal lipid value on or before the index date.
- **Date and timing:** at or before the index date.
- **Missing:** as for hypertension.
- **Class:** 2.
- **Handling:** used only in C1 and C3. HR about 1.08.

### Endometriosis and adenomyosis
- **Source:** Codebook flags from diagnosis codes, 2018–2026 extract.
- **Date and timing:** any date in the extract, so they can follow the anchor.
- **Missing:** none.
- **Class:** 3 for timing.
- **Handling:** describe the population; not used as predictors.

## Recommended predictor set for the next model

**Use:** age, Hb deficit, hypertension, CAD, CKD, and hemoglobinopathy coded on or before the anchor (`hbp_pre`). This is C4 plus `hbp_pre`, or the pre-specified C3, depending on the P05 and P06 results. See `CANDIDATE_MODEL_COMPARISON.md`.

**Restrict to:** women anchored at the first coded diagnosis.

**Remove until dated records exist:** migraine, diabetes, estrogen-containing therapy, and hemoglobinopathy coded on any date.

**Do not use:** smoking and BMI. Smoking cannot be estimated from these data; BMI has no date.
