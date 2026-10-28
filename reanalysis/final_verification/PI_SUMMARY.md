# MED-HACK stroke risk score: summary for the PI

**Bottom line:** the MED-HACK score is **not ready for publication or clinical use**. I recommend **pausing** model development until two kinds of data arrive: dated patient records and chart review. The project is not dead: a modest predictive signal remains. But the published version cannot be defended as it stands.

## What we checked

Three rounds of review of the stroke risk model built from 24,006 women with fibroids, adenomyosis or endometriosis (345 coded strokes):
1. a full statistical audit
2. a repair attempt using only information known at diagnosis
3. an independent re-check of every key number

All the calculations reproduce exactly. The problems are in the data and the design, not the arithmetic.

## The main problems

1. **Some predictors may come from after the diagnosis, in some cases after the stroke.**
   - Migraine, diabetes and estrogen use have no dates in our data.
   - Hemoglobinopathy and smoking were sometimes recorded after the stroke.
   - A model that "knows the future" looks better than it really is.
2. **The extra accuracy of MED-HACK depends on those undated predictors.** With only predictors known at diagnosis (age, hemoglobin, hypertension, CAD, CKD), discrimination drops from about 0.71 to about 0.66.
   - The remaining gain over standard risk factors is small, about +0.02, and not statistically certain in later patients.
3. **It underestimates risk in recently diagnosed women.** In women diagnosed from 2023 on, observed strokes were about 1.45 times the predicted number. The cause is unknown; it is not simply a coding change across the hospital system.
4. **The stroke outcome is not chart-verified.** Strokes were identified from codes. We don't know what share are true strokes.
5. **The 5-year risk is not validated.** Follow-up is too short; 3 years is the supported horizon.
6. **The online calculator has issues.** It shows two different risk numbers, can give near-100% risks for unusual combinations, and its age point ("+1 under 45") contradicts the data. It does say "Research use only", but its other caveats sit on a back tab.

## What this means

- Fibroids, anemia or estrogen should **not** be described as causing stroke. This is a prediction study, not a causal one.
- The model applies only to coded strokes, in one center, in women who had a blood count.
- No clinical decision (treatment, imaging, stopping hormones) should be based on the score.

## Decisions needed from you

1. **Calculator:** I recommend taking it offline. If it stays, add a prominent warning at the top.
2. **Manuscript:** pause submission of the current nine-predictor paper.
3. **Data request:** ask the data team for **dated** diagnosis and prescription records (migraine, diabetes, hypertension, CAD, estrogen). This is the single most important step: it will show whether MED-HACK's distinctive predictors are real.
4. **Chart review:** approve chart access for a blinded review of 100 stroke cases (a sample is ready) to confirm the outcome.
5. **Repository:** make the GitHub repository private before sharing further. It is public and contains small counts.

## What happens next

Once the dated records and chart review are available, I will re-run the pre-planned analysis. Then we can decide whether to:
- continue to external validation
- publish a smaller, honest model
- stop

**Detailed reports:**
- `reanalysis/final_verification/FINAL_VERIFICATION_REPORT.md`
- `reanalysis/final_verification/PI_ACTION_LIST.md`
- `reanalysis/post_audit/MISSING_DATA_REQUEST.md`
