# Deliverable A — Executive audit report

**Study:** Anemia and stroke in women with benign uterine disease (fibroids, adenomyosis, endometriosis). Single US academic center, EHR, index years 2018–2026.

**Declared limitation of this audit.**
- The original Python analysis, manuscript and briefing were produced earlier in this same project by the same assistant that performed this audit. The audit is therefore **not fully independent**.
- Mitigation:
  - Every variable was re-derived in R from the raw source files.
  - The plan was pre-registered in git before any model was fitted.
  - Original results were compared by script, not by memory.
- A human statistician should still review the work before submission.

## 1. What was received and analyzed
- **Patient-level data:**
  - analysis master file (52,827 women)
  - Codebook workbook
  - raw dated laboratory extracts (Hb, MCV, MCH/MCHC, RDW, platelets; one duplicate file not used)
  - dated diagnosis extracts
  - one-row-per-woman encounter files (last contact, death)
  - administered medications
- **All of it was used**, except: ferritin (no raw extract), medications (not needed for the revised estimand) and the duplicate lab file.
- **Not available:**
  - stroke-code extraction queries
  - chart validation of stroke codes
  - admission histories
  - transfusion data
  - outpatient prescriptions

## 2. Could the original analysis be reproduced?
**Yes.**
- 73 of 80 original estimates were reproduced within 2%, with identical sample sizes and event counts:
  - cross-sectional ORs (N 31,472; 1,460 strokes)
  - original time-to-event cohort (31,400 / 458)
  - original time-updated cohort (31,924 / 465)
- The few differences are confined to time-updated estimates resting on 6 events. They are explained by how same-day repeat Hb results are combined.
- **The arithmetic is right. The design is the problem.**

## 3. Most consequential data-quality problems
1. **Half the "primary" strokes happened before the uterine diagnosis.** Of 1,538 strokes/TIAs, 786 (51%) occurred before or at the index encounter and 202 (13%) are undated.
2. **The exposure often followed the outcome.** The Codebook Hb was the value closest to the index date in either direction:
   - measured after index in 13,661 women
   - after the stroke in 699 cases
3. **Stroke dates are earliest code dates, not verified event dates.** 10 incident events rest on sequela codes only and 13 on imaging only.
4. **Testing is far from random.** 35% of stroke-free women had no baseline Hb, and they had much less recorded health care. Without blood tests, fewer strokes are also recorded.
5. **Several covariates are undated or recorded after index:**
   - thrombophilia, AF, heart failure, VTE, migraine (any date)
   - hormonal therapy and antithrombotics (window extends after index)
   - smoking unknown in 66–73%
6. **Ferritin cannot be traced to source.**

## 4. Most consequential statistical problems in the original
1. The headline ORs (moderate 1.77, severe 2.01) come from a **cross-sectional comparison dominated by prevalent stroke**, with Hb frequently measured after the event.
2. **The "severity gradient" does not survive a design with correct time ordering.**
3. **The time-updated estimates depend on Hb measured just before a recorded stroke, and on testing frequency:**
   - no lag: moderate 2.40
   - 30-day lag: 1.56 (0.98–2.48)
   - lag + testing frequency: 1.20 (0.75–1.92)
4. Mistimed "pre-Hb" covariates (GI bleeding and pregnancy counted up to 30 days after Hb; HIV at any date). Adjustment for thrombophilia, which is likely recorded after stroke.
5. Many exploratory analyses without multiplicity control. A ferritin "inverse association" open to misreading.

## 5. Original findings: confirmed, changed, or unverifiable

| Original finding | Status after reanalysis |
|---|---|
| Anemia is associated with stroke | **Confirmed, in weaker form**: baseline anemia → incident stroke RR 1.50 (1.18–1.91); robust (1.38–1.57) across 11 sensitivity analyses |
| Moderate anemia OR 1.77; severe OR 2.01 | **Changed**: in the incident design, moderate 1.49 (0.96–2.31) and severe 1.64 (0.79–3.40) are imprecise, with no clear gradient |
| Per-g/dL association below 13 g/dL | **Confirmed**: RR 1.14 (1.05–1.23) per 1 g/dL lower |
| Association mainly with ischemic stroke | **Confirmed**: ischemic RR 1.52 (1.18–1.95). TIA null (no TIA in moderate or severe anemia) |
| Time-updated moderate RR 1.86, persisting after testing adjustment | **Changed**: attenuated to the null with a 30-day lag plus testing adjustment |
| Persistent anemia carries more risk than resolved | **Not supported**: persistent vs resolved 0.93 (0.53–1.62) |
| Morphology / RDW patterns, ferritin, treatment, mediation | **Not re-analyzed** (exploratory; ferritin untraceable). They should not appear in the main paper. |

## 6. Is the manuscript's primary conclusion supported?
**Partly.**
- **Supported:** "Women with benign uterine disease who are anemic at diagnosis have a modestly higher subsequent rate of recorded (mainly ischemic) stroke." In absolute terms, 5-year incidence was 1.8% without anemia vs 2.6–3.1% with anemia.
- **Not supported:**
  - a severity gradient
  - "moderate and severe anemia were associated with higher risk" as a headline
  - any implication that anemia itself raises stroke risk
  - suggestions that Hb at diagnosis should trigger stroke-risk assessment, which rests on a marker relationship of uncertain specificity
- **Most consistent interpretation:** the data favor anemia as a **marker of underlying illness, health-care contact and vascular risk**. They cannot rule out a modest direct association.

---

# Deliverable G — Statistical reviewer decision

## Recommendation: **Proceed after substantial revision**

The study is potentially publishable. The data are large, the code is reproducible, and the baseline anemia → incident stroke association is consistent across all pre-specified sensitivity analyses. However, the submitted analysis answers the wrong question (prevalent stroke, exposure after outcome). Its headline numbers and severity claims do not hold in a correctly time-ordered design.

## The three most important actions before journal submission
1. **Replace the primary analysis with the incident design and rewrite the abstract and results around it.**
   - New primary: any anemia RR 1.50 (1.18–1.91); per grade 1.26; continuous 1.14 per g/dL below 13.
   - Absolute 5-year incidence 1.8% vs 2.6–3.1%.
   - The cross-sectional ORs move to the supplement, labeled "prevalent or incident stroke".
   - The AAN abstract numbers (1.77 / 2.01 / 1.57 / 1.86 / 1.54) should not be reused in a journal submission.
2. **Report the time-updated lag and testing-frequency results and the landmark result in the main text, and frame anemia as a marker.**
   - Remove the severity-gradient language, the "time-updated association persisted" statement and the clinical-screening implication.
   - State plainly that the data cannot show whether treating anemia changes stroke risk.
   - Keep TIA separate and report its null result.
3. **Strengthen outcome and exposure documentation before submission:**
   - Describe exactly how `stroke_any` and its dates were derived (code lists, problem list, imaging). If at all possible, report the positive predictive value of the stroke codes from a small chart-review sample.
   - State the CBC-only Hb definition.
   - Resolve the co-author disagreement on adjusting for uterine diagnosis group. Results are unchanged either way (<2%), but the paper must state one justified choice.

## Lesser items
- Remove or confine to the supplement: morphology/RDW patterns, ferritin, treatment, mediation and subgroup analyses. Label them exploratory.
- Report smoking "unknown" (66%) and the testing-selection analysis (S11) as limitations.
- Single center; left-truncated index (diagnosis extract begins in 2018).
- Use STROBE. Report rate ratios, hazard ratios and odds ratios with their correct names; never interchange them.
