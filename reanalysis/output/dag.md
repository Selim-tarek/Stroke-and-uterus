# Causal framework (written DAG) for the revised primary model

**Estimand:** the association between baseline anemia (most recent laboratory Hb in the 3 years up to the index date) and the rate of a first stroke after the index date, among women with benign uterine disease who were free of recorded stroke at index.

The design is observational, so this framework says which variables to adjust for. It does not make the estimate causal.

## Diagram

```
                 age, race, BMI, smoking, hypertension, diabetes, dyslipidemia, CAD
                 (vascular risk factors)  ───────────────────────────────┐
                        │ (CKD, diabetes, inflammation → anemia)          │
                        ▼                                                 ▼
 uterine dx group ─► heavy bleeding ─► ANEMIA (baseline Hb) ─ ─ ─ ? ─ ─ ─► STROKE (after index)
        │                 ▲                 ▲    │                          ▲
        │                 │                 │    └─► more Hb testing ◄───── │ ─ healthcare contact
        │       hormonal therapy ──────────(treats bleeding)                │   (finds both anemia
        │                 └──────────────────────────── estrogen ──────────►│    and stroke codes)
        └── endometriosis / inflammation ─────────────────────────────────►│
 chronic illness: CKD, liver disease, alcohol, IBD/HIV, hemoglobinopathy (sickle cell), malabsorption,
 GI bleeding, recent pregnancy ─────► ANEMIA   and ─────► STROKE
 recorded only after stroke, undated: thrombophilia (work-up), AF, heart failure, VTE, anticoagulants
```

## Classification

| Variable | Role | Decision |
|---|---|---|
| Age, race, BMI, smoking, hypertension, diabetes, dyslipidemia, CAD before index | Confounders: affect stroke; linked to anemia via CKD, diabetes, inflammation and nutrition | **Adjust** |
| Migraine | Possible confounder: iron deficiency is linked to migraine; migraine with aura is linked to stroke | **Adjust** (Codebook timing: any date) |
| Hormonal therapy type | Confounder: treats heavy bleeding (reduces anemia); estrogen raises stroke risk | **Adjust** (medication window starts 1 year before index; partly after index) |
| Uterine diagnosis group | Confounder: affects anemia via fibroid bleeding; endometriosis has been associated with stroke [ref needed] | **Adjust**. This differs from the co-author request to drop it; see the deviation log. |
| Heavy/abnormal uterine bleeding | Cause of the exposure. Its only clear route to stroke runs through anemia, so it acts like an instrument; adjusting adds variance and can amplify bias. It is undated. | **Not adjusted** (S6) |
| CKD, liver disease, alcohol use disorder, IBD/HIV, hemoglobinopathy, malabsorption, GI bleeding and pregnancy in the year before index | Confounders: cause anemia and are associated with stroke | **Adjust** (dated on or before index) |
| Thrombophilia, AF, heart failure, VTE | Undated. Often diagnosed during a stroke work-up, so they may be consequences of the outcome. Thrombophilia is not a plausible cause of anemia. | **Not adjusted** (S7, over-adjustment check) |
| Antithrombotics | Medication window extends after index; marks hospitalization | Not adjusted |
| Number of Hb tests before index | Proxy for healthcare contact and illness, and partly a consequence of known anemia (descendant of the exposure) | **Not adjusted** (S8): adjusting may remove part of a real association or reduce detection bias. The data cannot tell these apart. |
| Being tested at all (having a baseline Hb) | Selection variable | Weighted in S11 |

## What this framework cannot fix
- Unmeasured general ill health and inflammation, which drive both anemia and stroke. The E-values describe how strong such a confounder would need to be.
- Socioeconomic status: not available.
- Outcome detection: women with more contact have more chances to be coded with stroke, and also more Hb tests.
