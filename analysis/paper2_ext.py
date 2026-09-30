"""Paper 2 extension: temporality, additional confounders from the supplementary
diagnosis extracts (Diagnosis_22-26), anaemia morphology, stroke subtypes,
subgroups and E-values.

Pre-specified (fixed before fitting)
------------------------------------
Exposure: anaemia grade (none / mild / moderate / severe), plus per-grade trend.
Base covariates: as Paper 2 (age, race, BMI, HTN, DM, dyslipidaemia, smoking,
migraine, thrombophilia, hormonal type, uterine bleeding, dx group).

Comorbidity timing (from the supplementary extracts, linked in memory by MRN):
- hereditary (sickle disease, sickle trait, thalassaemia, other
  haemoglobinopathy): any record, any date;
- chronic acquired (CKD, ESRD/dialysis/transplant, cirrhosis/portal
  hypertension, chronic liver disease, alcohol use disorder, malabsorption/
  coeliac/bariatric, menopause): first record on or before the Hb date (index
  date if no Hb) - i.e. documented before the Hb was measured;
- recent GI bleeding: any record from 365 days before to 30 days after the Hb
  date;
- pregnancy/delivery: any record from 365 days before to 30 days after the Hb
  date (ICD-9 only in the extract - coverage after 2015 is minimal; reported).

Analyses
1. Temporality: cases restricted to strokes whose Hb was measured >= 30 days
   BEFORE stroke_date (controls: all women without stroke). Repeated for
   ischaemic and incident strokes. Also: excluding Hb within 30 days of stroke.
2. Extended adjustment.
   2a (primary): conditions documented BEFORE the Hb measurement - CAD
   (baseline per Codebook), haemoglobinopathy, CKD, chronic liver disease,
   alcohol use disorder, GI bleeding in the year before Hb, malabsorption, HIV,
   pregnancy in the year before Hb.
   2b (over-adjustment check): 2a + undated conditions and treatments that may
   follow the stroke (heart failure, AF, VTE, malignancy ever, anticoagulant
   use - anticoagulants are largely started after stroke).
3. Restricted cohort: excluding (a) any haemoglobinopathy; (b) all anaemia-
   causing conditions listed in 2 (haemoglobinopathy, CKD, ESRD, cirrhosis,
   chronic liver disease, alcohol disorder, recent GI bleed, malabsorption).
4. Anaemia morphology: no anaemia / microcytic (MCV <80) / normocytic /
   macrocytic (>100) anaemia.
5. Stroke subtypes vs no stroke: ischaemic, TIA, ICH, SAH, CVT, unspecified
   (per-grade trend; categorical ORs only where every level has >= 5 events).
6. Subgroups with interaction tests: age (<45 / >=45), race (Black / not),
   uterine group (any fibroids / no fibroids).
7. E-values (VanderWeele & Ding) for the main ORs (outcome prevalence <15%, so
   OR approximates RR).
"""
import numpy as np
import pandas as pd

from .comorbid_codes import classify_rows, load_extracts
from .paper2 import COV, EXPO, cov_for
from .utils import (RESULTS, Term, add_table, fit_logit, fmt_or, fmt_p, get_or, log, register_fit, wald)

GRADES = ["Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"]
# conditions documented before the Hb (adjustment 2a); shared with the longitudinal module
PRE_VARS = ["cad", "haemoglobinopathy_any", "ckd", "liver_any", "alcohol_disorder", "gi_bleed_recent",
                "malabsorption", "ibd", "hiv", "pregnancy_recent"]
HEREDITARY = ["sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy"]
CHRONIC = ["ckd", "ckd_esrd", "cirrhosis", "chronic_liver", "alcohol_disorder", "malabsorption", "ibd", "menopause"]
ANAEMIA_CODED = ["anaemia_iron", "anaemia_b12_folate", "anaemia_nutritional", "anaemia_blood_loss_acute",
                 "anaemia_neoplastic_chemo", "anaemia_chronic_disease", "anaemia_haemolytic_aplastic",
                 "anaemia_pregnancy", "anaemia_unspecified"]
LABELS = {"sickle_disease": "Sickle cell disease", "sickle_trait": "Sickle cell trait",
          "thal_minor": "Thalassaemia minor/trait (incl. HPFH)", "thal_other": "Thalassaemia (other/unspecified)",
          "other_haemoglobinopathy": "Other haemoglobinopathy / hereditary haemolytic",
          "haemoglobinopathy_any": "Any haemoglobinopathy", "ckd": "Chronic kidney disease (any)",
          "ckd_esrd": "ESRD / dialysis / kidney transplant", "cirrhosis": "Cirrhosis / portal hypertension",
          "chronic_liver": "Chronic liver disease (incl. NAFLD/NASH)", "liver_any": "Any chronic liver disease",
          "alcohol_disorder": "Alcohol use disorder", "gi_bleed_recent": "GI bleeding within 1 y before Hb",
          "gi_bleed_ever": "GI bleeding (ever)", "malabsorption": "Malabsorption / coeliac / bariatric",
          "pregnancy_recent": "Pregnancy/delivery/postpartum within 1 y before Hb", "hiv": "HIV infection",
          "ibd": "Inflammatory bowel disease",
          "postmenopausal": "Postmenopausal at Hb/index (coded)",
          "anaemia_iron_nearhb": "Coded iron-deficiency anaemia (±1 y of Hb)",
          "anaemia_b12_folate_nearhb": "Coded B12/folate anaemia (±1 y of Hb)",
          "anaemia_blood_loss_acute_nearhb": "Coded acute blood-loss anaemia (±1 y of Hb)",
          "anaemia_chronic_disease_nearhb": "Coded anaemia of chronic disease/CKD (±1 y of Hb)",
          "anaemia_neoplastic_chemo_nearhb": "Coded neoplastic/chemotherapy anaemia (±1 y of Hb)",
          "anaemia_pregnancy_nearhb": "Coded anaemia of pregnancy/puerperium (±1 y of Hb)",
          "any_anaemia_cause": "Any anaemia-causing condition (exclusion set)"}


def build_flags(elig):
    x, files = load_extracts()
    c = classify_rows(x)
    RESULTS["p2x_files"] = [f.rsplit("/", 1)[-1] for f in files]
    unm = c[c.flag.isin(["UNMATCHED", "EXCLUDED"])].groupby(["flag", "code", "desc"]).mrn.nunique() \
        .reset_index(name="patients").sort_values("patients", ascending=False)
    add_table("P2X_codes_not_counted", unm, "Supplementary-extract terms that were NOT counted as a condition "
                                            "(non-qualifying by design, or unmatched). Patients = all records.")
    kept = c[~c.flag.isin(["UNMATCHED", "EXCLUDED"])]
    summ = kept.groupby("flag").agg(patients=("mrn", "nunique"), rows=("mrn", "size"),
                                    first=("date", "min"), last=("date", "max")).reset_index()
    add_table("P2X_codes_counted", summ, "Supplementary-extract conditions: patients and date range (all records).")
    e = elig[["mrn", "index_date", "hgb_date"]].copy()
    e["idx"] = pd.to_datetime(e.index_date)
    e["hbd"] = pd.to_datetime(e.hgb_date, errors="coerce")
    e["ref"] = e["hbd"].fillna(e["idx"])
    k = kept.merge(e, on="mrn", how="inner")
    out = pd.DataFrame(index=elig.index)
    mrn_index = pd.Series(elig.index, index=elig.mrn)
    for fl in HEREDITARY:
        ids = k.loc[k.flag == fl, "mrn"].unique()
        out[fl] = elig.mrn.isin(ids).astype(float)
    for fl in CHRONIC:
        s = k[(k.flag == fl) & (k.date <= k.ref)]
        out[fl] = elig.mrn.isin(s.mrn.unique()).astype(float)
    out["hiv"] = elig.mrn.isin(k.loc[k.flag == "hiv", "mrn"].unique()).astype(float)
    out["postmenopausal"] = out.pop("menopause")
    for a in ANAEMIA_CODED:
        s_ = k[(k.flag == a) & k.hbd.notna() & ((k.date - k.hbd).dt.days.abs() <= 365)]
        out[f"{a}_nearhb"] = elig.mrn.isin(s_.mrn.unique()).astype(float)
    gb = k[k.flag == "gi_bleed"]
    out["gi_bleed_ever"] = elig.mrn.isin(gb.mrn.unique()).astype(float)
    rec = gb[gb.hbd.notna() & (gb.date >= gb.hbd - pd.Timedelta(days=365)) & (gb.date <= gb.hbd + pd.Timedelta(days=30))]
    out["gi_bleed_recent"] = elig.mrn.isin(rec.mrn.unique()).astype(float)
    pg = k[k.flag == "pregnancy"]
    prec = pg[pg.hbd.notna() & (pg.date >= pg.hbd - pd.Timedelta(days=365)) & (pg.date <= pg.hbd + pd.Timedelta(days=30))]
    out["pregnancy_recent"] = ((elig.mrn.isin(prec.mrn.unique())) | (out["anaemia_pregnancy_nearhb"] == 1)).astype(float)
    out["haemoglobinopathy_any"] = out[HEREDITARY].max(axis=1)
    out["liver_any"] = out[["cirrhosis", "chronic_liver"]].max(axis=1)
    out["any_anaemia_cause"] = out[["haemoglobinopathy_any", "ckd", "ckd_esrd", "liver_any", "alcohol_disorder",
                                    "gi_bleed_recent", "malabsorption", "ibd", "hiv", "pregnancy_recent",
                                    "anaemia_neoplastic_chemo_nearhb"]].max(axis=1)
    RESULTS["p2x_preg_recent"] = int(out.pregnancy_recent.sum())
    RESULTS["p2x_preg_icd10_patients"] = int(kept[(kept.flag == "pregnancy") &
                                                  kept.code.str.upper().str.startswith("O", na=False)].mrn.nunique())
    log("Paper 2 ext", f"Supplementary extracts {RESULTS['p2x_files']}; eligible-cohort flags: " +
        ", ".join(f"{f} {int(out[f].sum())}" for f in out.columns))
    return out


def evalue(o, lo, hi):
    def ev(r):
        r = 1 / r if r < 1 else r
        return r + np.sqrt(r * (r - 1))
    e_pt = ev(o)
    if lo <= 1 <= hi:
        e_ci = 1.0
    else:
        e_ci = ev(lo if o > 1 else hi)
    return e_pt, e_ci


def grade_row(f, label, extra=None):
    row = {"Analysis": label, "N": f.n, "Events": f.events, "EPV": round(f.epv, 1), "Flags": "; ".join(f.flags)}
    for g in GRADES:
        c = f"anemia_cat={g}"
        if c in f.params.index:
            r = get_or(f, c)
            sparse = any(s.startswith(c + " ") for s in f.sparse)
            row[g.split(" ")[0]] = "not estimated (sparse)" if sparse else r["txt"]
            row[g.split(" ")[0] + " OR"] = np.nan if sparse else r["OR"]
            row[g.split(" ")[0] + " lo"] = np.nan if sparse else r["lo"]
            row[g.split(" ")[0] + " hi"] = np.nan if sparse else r["hi"]
    if extra:
        row.update(extra)
    return row


def run(elig):
    d = elig.copy()
    flags = build_flags(d)
    for cname in flags.columns:
        d[cname] = flags[cname]
    d["anticoag"] = (d.antithrombotic == 2).astype(float)
    d["y_isch"] = np.where(d.stroke_any == 0, 0.0, np.where(d.stroke_type == 1, 1.0, np.nan))
    d["y_incident"] = np.where(d.stroke_any == 0, 0.0, np.where(d.stroke_timing == 3, 1.0, np.nan))
    hbd = pd.to_datetime(d.hgb_date, errors="coerce")
    sd = pd.to_datetime(d.stroke_date, errors="coerce")
    gap = (sd - hbd).dt.days  # >0: Hb before stroke
    d["hb_timing"] = np.select([d.stroke_any == 0, sd.isna() | hbd.isna(), gap >= 30, gap.abs() < 30, gap <= -30],
                               ["no stroke", "unknown", "Hb >=30 d before stroke", "Hb within 30 d of stroke",
                                "Hb >=30 d after stroke"], "unknown")
    RESULTS["p2x_hb_timing"] = d.loc[(d.stroke_any == 1) & d.hgb.notna(), "hb_timing"].value_counts().to_dict()

    # ---------------- comorbidity prevalence table
    prev = []
    for fl in list(LABELS):
        if fl not in d:
            continue
        row = {"Condition": LABELS[fl], "All eligible n": int(d[fl].sum()),
               "All %": round(100 * d[fl].mean(), 2)}
        for g in ["None (Hb ≥12)"] + GRADES:
            s = d[d.anemia_cat == g]
            row[f"{g.split(' ')[0]} anaemia %"] = round(100 * s[fl].mean(), 2)
        for v, lab in [(1, "Stroke %"), (0, "No stroke %")]:
            row[lab] = round(100 * d.loc[d.stroke_any == v, fl].mean(), 2)
        prev.append(row)
    add_table("P2X_comorbidity_prev", pd.DataFrame(prev), "Prevalence of anaemia-related conditions by anaemia grade "
                                                          "and stroke status (eligible cohort).")
    RESULTS["p2x_prev"] = pd.DataFrame(prev)

    rows = []
    base_terms = lambda o, dd=None, tag="": [EXPO["anemia_cat"]] + cov_for(d if dd is None else dd, o, ["anemia_cat"], tag)  # noqa: E731
    trend = lambda f_terms: [Term("anemia", "cont") if t.var == "anemia_cat" else t for t in f_terms]  # noqa: E731

    def fit_pair(data, o, terms, label, section):
        f = register_fit(fit_logit(data, o, terms, name=label), "P2X", label)
        ft = fit_logit(data, o, trend(terms), compute_vif=False)
        g = get_or(ft, "anemia")
        e_pt, e_ci = evalue(g["OR"], g["lo"], g["hi"])
        rows.append(grade_row(f, label, {"Section": section, "Per grade (trend)": g["txt"], "p-trend": fmt_p(g["p"]),
                                         "trend OR": g["OR"], "trend lo": g["lo"], "trend hi": g["hi"],
                                         "trend p": g["p"]}))
        return f, g

    # 0. reference: primary models
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke"), ("y_incident", "Incident stroke")]:
        fit_pair(d, o, base_terms(o), f"Primary: {lab}", "0. Primary (for comparison)")

    # 1. temporality
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke"), ("y_incident", "Incident stroke")]:
        dd = d.copy()
        keep = (dd[o] == 0) | (dd.hb_timing == "Hb >=30 d before stroke")
        dd.loc[~keep, o] = np.nan
        fit_pair(dd, o, base_terms(o, dd, "HbBefore"), f"Hb ≥30 d before stroke: {lab}", "1. Temporality")
    dd = d.copy()
    dd.loc[(dd.stroke_any == 1) & (dd.hb_timing == "Hb within 30 d of stroke"), "stroke_any"] = np.nan
    fit_pair(dd, "stroke_any", base_terms("stroke_any", dd, "noAcute"), "Excluding Hb within 30 d of stroke: Any stroke",
             "1. Temporality")

    # 2. extended adjustment
    pre_vars = PRE_VARS
    post_vars = ["chf", "afib", "vte_history", "malignancy_ever", "anticoag"]
    ext_a = [Term(v, "bin", label=LABELS.get(v, v)) for v in pre_vars]
    ext_b = ext_a + [Term(v, "bin", label=LABELS.get(v, v)) for v in post_vars]
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke"), ("y_incident", "Incident stroke")]:
        f, _ = fit_pair(d, o, base_terms(o) + ext_a, f"2a Pre-Hb conditions: {lab}", "2. Extended adjustment")
        if o == "stroke_any":
            from .utils import or_table
            add_table("P2X_extended_full_model", or_table(f, model_label="Any stroke, adjustment 2a (pre-Hb conditions)"))
            RESULTS["p2x_ext_fit"] = f
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke"), ("y_incident", "Incident stroke")]:
        fit_pair(d, o, base_terms(o) + ext_b, f"2b Over-adjustment check (+ undated/post-stroke): {lab}",
                 "2. Extended adjustment")
    # original primary model (with uterine bleeding and diagnosis group) as a sensitivity analysis
    from .paper2 import UTERINE_TERMS
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        fit_pair(d, o, base_terms(o) + UTERINE_TERMS, f"+ uterine bleeding and diagnosis group: {lab}",
                 "2. Extended adjustment")
    # temporality + extended together (strictest defensible)
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        dd = d.copy()
        dd.loc[~((dd[o] == 0) | (dd.hb_timing == "Hb >=30 d before stroke")), o] = np.nan
        fit_pair(dd, o, base_terms(o, dd, "HbBeforeExt") + ext_a,
                 f"Hb ≥30 d before stroke + 2a: {lab}", "2. Extended adjustment")
        fit_pair(dd, o, base_terms(o, dd, "HbBeforeExtB") + ext_b,
                 f"Hb ≥30 d before stroke + 2b: {lab}", "2. Extended adjustment")

    # 3. restricted cohorts
    for excl, lab in [("haemoglobinopathy_any", "Excluding haemoglobinopathies"),
                      ("any_anaemia_cause", "Excluding all anaemia-causing conditions")]:
        dd = d[d[excl] == 0].copy()
        RESULTS[f"p2x_n_{excl}"] = int(d[excl].sum())
        for o, olab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
            fit_pair(dd, o, base_terms(o, dd, excl), f"{lab}: {olab}", "3. Restricted cohort")

    # 3b. markers of general illness: early death and serious chronic illness (adjustment 2a)
    from .followup import encounters as _enc
    _lc, _dth = _enc()
    dd_dt = pd.to_datetime(d.mrn.map(_dth))
    d["death_1y_hb"] = ((dd_dt.notna()) & ((dd_dt - hbd).dt.days <= 365)).astype(float)
    d["serious_illness"] = d[["malignancy_ever", "chf", "ckd", "liver_any", "hiv"]].max(axis=1)
    RESULTS["p2x_n_death1y"] = int(d.loc[d.hgb.notna(), "death_1y_hb"].sum())
    RESULTS["p2x_n_serious"] = int(d.loc[d.hgb.notna(), "serious_illness"].sum())
    for excl, lab in [("death_1y_hb", "Excluding deaths within 1 y of Hb (+2a)"),
                      ("serious_illness", "Excluding cancer, heart failure, CKD, liver disease, HIV (+2a)")]:
        dd = d[d[excl] == 0].copy()
        ea = [x for x in ext_a if dd[x.var].nunique() > 1]
        for o, olab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
            fit_pair(dd, o, base_terms(o, dd, excl) + ea, f"{lab}: {olab}", "3. Restricted cohort")

    # 4. anaemia morphology
    d["anaemia_type"] = np.select(
        [d.anemia_any == 0, (d.anemia_any == 1) & (d.mcv < 80), (d.anemia_any == 1) & (d.mcv <= 100),
         (d.anemia_any == 1) & (d.mcv > 100)],
        ["No anaemia", "Microcytic anaemia", "Normocytic anaemia", "Macrocytic anaemia"], default="")
    d["anaemia_type"] = d["anaemia_type"].replace("", np.nan)
    morph = []
    at = Term("anaemia_type", "cat", ref="No anaemia",
              levels=["No anaemia", "Microcytic anaemia", "Normocytic anaemia", "Macrocytic anaemia"], label="Anaemia type")
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        f = register_fit(fit_logit(d, o, [at] + cov_for(d, o, ["anaemia_type"], "morph")), "P2X", f"Morphology: {lab}")
        for lv in at.levels:
            s = d[d.anaemia_type == lv]
            r = {"Outcome": lab, "Anaemia type": lv, "n": len(s), "Events": int(s[o].sum())}
            if lv == at.ref:
                r["Adjusted OR (95% CI)"] = "1.00 (reference)"
            else:
                c = f"anaemia_type={lv}"
                sp = any(x.startswith(c + " ") for x in f.sparse)
                g = get_or(f, c)
                r.update({"Adjusted OR (95% CI)": "not estimated (sparse)" if sp else g["txt"],
                          "OR": np.nan if sp else g["OR"], "CI low": np.nan if sp else g["lo"],
                          "CI high": np.nan if sp else g["hi"], "p": np.nan if sp else g["p"]})
            morph.append(r)
        chi2, df_, p = wald(f, f.colmap["anaemia_type"][1:]) if len(f.colmap["anaemia_type"]) > 1 else (np.nan, 0, np.nan)
        # heterogeneity between micro/normo/macro: test equality of the anaemia-type coefficients
        cols = f.colmap["anaemia_type"]
        if len(cols) >= 2:
            L = np.zeros((len(cols) - 1, len(f.params)))
            idx = list(f.params.index)
            for i in range(len(cols) - 1):
                L[i, idx.index(cols[0])] = 1
                L[i, idx.index(cols[i + 1])] = -1
            b = L @ f.params.values
            V = L @ f.cov.values @ L.T
            from scipy import stats as _st
            chi = float(b @ np.linalg.pinv(V) @ b)
            morph.append({"Outcome": lab, "Anaemia type": f"Heterogeneity across types χ²({len(cols) - 1})",
                          "Adjusted OR (95% CI)": f"{chi:.2f}", "p": float(_st.chi2.sf(chi, len(cols) - 1))})
    # 4b. laboratory anaemia pattern (MCV + RDW-CV from the dated lab extracts, within ±30 d of Hb;
    # Bessman classification; RDW-CV > 14.5% = high)
    PAT = ["No anaemia", "Microcytic, high RDW (iron-deficiency pattern)",
           "Microcytic, normal RDW (thalassaemia-trait pattern)", "Normocytic, normal RDW",
           "Normocytic, high RDW (mixed / early iron deficiency)", "Macrocytic"]
    hi_rdw = d.rdw > 14.5
    have = d.mcv.notna() & d.rdw.notna()
    an = d.anemia_any == 1
    d["anaemia_pattern"] = np.select(
        [d.anemia_any == 0, an & have & (d.mcv < 80) & hi_rdw, an & have & (d.mcv < 80) & ~hi_rdw,
         an & have & d.mcv.between(80, 100) & ~hi_rdw, an & have & d.mcv.between(80, 100) & hi_rdw,
         an & have & (d.mcv > 100)], PAT, default="")
    d["anaemia_pattern"] = d["anaemia_pattern"].replace("", np.nan)
    RESULTS["p2x_pattern_missing"] = int((an & ~have).sum())
    pt = Term("anaemia_pattern", "cat", ref="No anaemia", levels=PAT, label="Laboratory anaemia pattern")
    pat = []
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        for adj, extra in [("Paper 2 covariates", []), ("+ pre-Hb conditions (2a)", ext_a)]:
            f = register_fit(fit_logit(d, o, [pt] + cov_for(d, o, ["anaemia_pattern"], f"pat{o}") + extra), "P2X",
                             f"Anaemia pattern: {lab} ({adj})")
            for lv in PAT:
                s_ = d[d.anaemia_pattern == lv]
                r = {"Outcome": lab, "Adjustment": adj, "Pattern": lv, "Women": len(s_), "Strokes": int(s_[o].sum()),
                     "Median MCV": round(float(s_.mcv.median()), 1) if s_.mcv.notna().any() else np.nan,
                     "Median RDW-CV": round(float(s_.rdw.median()), 1) if s_.rdw.notna().any() else np.nan}
                if lv == "No anaemia":
                    r["Adjusted OR (95% CI)"] = "1.00 (reference)"
                else:
                    c = f"anaemia_pattern={lv}"
                    sp = any(x.startswith(c + " ") for x in f.sparse)
                    g = get_or(f, c)
                    r.update({"Adjusted OR (95% CI)": "not estimated (sparse)" if sp else g["txt"],
                              "OR": np.nan if sp else g["OR"], "CI low": np.nan if sp else g["lo"],
                              "CI high": np.nan if sp else g["hi"], "p": np.nan if sp else g["p"]})
                pat.append(r)
            cols = f.colmap["anaemia_pattern"]
            L = np.zeros((len(cols) - 1, len(f.params)))
            idx = list(f.params.index)
            for i in range(len(cols) - 1):
                L[i, idx.index(cols[0])] = 1
                L[i, idx.index(cols[i + 1])] = -1
            b = L @ f.params.values
            V = L @ f.cov.values @ L.T
            from scipy import stats as _stp
            chi = float(b @ np.linalg.pinv(V) @ b)
            pat.append({"Outcome": lab, "Adjustment": adj, "Pattern": f"Heterogeneity across patterns χ²({len(cols) - 1})",
                        "Adjusted OR (95% CI)": f"{chi:.2f}", "p": float(_stp.chi2.sf(chi, len(cols) - 1))})
    pat = pd.DataFrame(pat)
    pat["p (text)"] = pat["p"].map(lambda v: fmt_p(v) if pd.notna(v) else "")
    add_table("P2X_anaemia_pattern", pat, "Laboratory anaemia pattern (MCV and RDW-CV from the dated laboratory "
              "record nearest the Hb, within 30 days; RDW-CV >14.5% = high) vs no anaemia. Anaemic women without both "
              "indices are excluded.")
    RESULTS["p2x_pattern"] = pat

    # coded aetiology among anaemic women: iron-deficiency coded vs anaemia without an iron-deficiency code
    d["anaemia_coded"] = np.select([d.anemia_any == 0, (d.anemia_any == 1) & (d.anaemia_iron_nearhb == 1),
                                    d.anemia_any == 1], ["No anaemia", "Anaemia, iron deficiency coded",
                                                         "Anaemia, no iron-deficiency code"], default="")
    d["anaemia_coded"] = d["anaemia_coded"].replace("", np.nan)
    ac = Term("anaemia_coded", "cat", ref="No anaemia",
              levels=["No anaemia", "Anaemia, iron deficiency coded", "Anaemia, no iron-deficiency code"])
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        f = register_fit(fit_logit(d, o, [ac] + cov_for(d, o, ["anaemia_coded"], "coded")), "P2X", f"Coded type: {lab}")
        for lv in ac.levels[1:]:
            s = d[d.anaemia_coded == lv]
            g = get_or(f, f"anaemia_coded={lv}")
            morph.append({"Outcome": lab, "Anaemia type": lv, "n": len(s), "Events": int(s[o].sum()),
                          "Adjusted OR (95% CI)": g["txt"], "OR": g["OR"], "CI low": g["lo"], "CI high": g["hi"],
                          "p": g["p"]})
    morph = pd.DataFrame(morph)
    morph["p (text)"] = morph["p"].map(fmt_p)
    add_table("P2X_anaemia_morphology", morph, "Anaemia type by MCV (vs no anaemia), adjusted as Paper 2.")
    RESULTS["p2x_morph"] = morph

    # 5. stroke subtypes
    sub = []
    types = {1: "Ischaemic stroke", 4: "TIA", 2: "Intracerebral haemorrhage", 3: "Subarachnoid haemorrhage",
             5: "Cerebral venous thrombosis", 9: "Unspecified"}
    for code, lab in types.items():
        o = f"y_type{code}"
        d[o] = np.where(d.stroke_any == 0, 0.0, np.where(d.stroke_type == code, 1.0, np.nan))
        ev = int((d[o] == 1).sum())
        ev_hb = int(((d[o] == 1) & d.anemia_cat.notna()).sum())
        r = {"Subtype": lab, "Events (all)": ev, "Events with Hb": ev_hb}
        if ev_hb >= 20:
            terms = [Term("anemia", "cont")] + (cov_for(d, o, ["anemia"], f"sub{code}") if ev_hb >= 100 else
                                                [t for t in COV if t.var in ("age_index", "htn")])
            f = register_fit(fit_logit(d, o, terms, name=f"subtype {lab}"), "P2X", f"Subtype: {lab} (trend)")
            g = get_or(f, "anemia")
            f2 = register_fit(fit_logit(d, o, [Term("anemia_any", "bin")] + terms[1:], name=f"subtype {lab} any"),
                              "P2X", f"Subtype: {lab} (any anaemia)")
            g2 = get_or(f2, "anemia_any")
            r.update({"Adjustment": "full Paper 2 set" if ev_hb >= 100 else "age + hypertension (few events)",
                      "Any anaemia OR (95% CI)": g2["txt"], "Per grade OR (95% CI)": g["txt"],
                      "OR": g["OR"], "CI low": g["lo"], "CI high": g["hi"], "p-trend": fmt_p(g["p"]),
                      "Flags": "; ".join(f.flags)})
        else:
            r.update({"Adjustment": "not estimated (<20 events with Hb)"})
        sub.append(r)
    sub = pd.DataFrame(sub)
    add_table("P2X_stroke_subtypes", sub, "Anaemia and stroke subtypes (each subtype vs no stroke).")
    RESULTS["p2x_sub"] = sub

    # 6. subgroups
    sg = []
    d["age45"] = (d.age_index >= 45).astype(float)
    d["black"] = (d.race4 == "Black").astype(float)
    d["fib_any"] = d.fibroids.astype(float)
    for var, lab, levels in [("age45", "Age", {0: "<45", 1: "≥45"}), ("black", "Race", {0: "Not Black", 1: "Black"}),
                             ("fib_any", "Fibroids", {0: "No fibroids (adenomyosis/endometriosis only)",
                                                      1: "Any fibroids"}),
                             ("uterine_bleeding", "Heavy/abnormal bleeding", {0: "No", 1: "Yes"}),
                             ("postmenopausal", "Menopausal status (coded)", {0: "Not coded postmenopausal",
                                                                             1: "Postmenopausal"})]:
        for v, vl in levels.items():
            dd = d[d[var] == v].copy()
            cv = [t for t in cov_for(dd, "stroke_any", ["anemia"], f"sg{var}{v}")
                  if not (var == "black" and t.var == "race4") and not (var == "fib_any" and t.var == "dxgrp")
                  and not (var == "uterine_bleeding" and t.var == "uterine_bleeding")]
            f = register_fit(fit_logit(dd, "stroke_any", [Term("anemia", "cont")] + cv, name=f"sg {var}={v}"), "P2X",
                             f"Subgroup {lab} {vl}")
            g = get_or(f, "anemia")
            sg.append({"Subgroup": lab, "Level": vl, "N": f.n, "Events": f.events, "Per grade OR (95% CI)": g["txt"],
                       "OR": g["OR"], "CI low": g["lo"], "CI high": g["hi"], "Flags": "; ".join(f.flags)})
        dd = d.copy()
        dd["ix"] = dd["anemia"] * dd[var]
        cv = [t for t in cov_for(dd, "stroke_any", ["anemia"], f"sgint{var}")
              if not (var == "black" and t.var == "race4") and not (var == "fib_any" and t.var == "dxgrp")
              and not (var == "uterine_bleeding" and t.var == "uterine_bleeding")]
        fi = fit_logit(dd, "stroke_any", [Term("anemia", "cont"), Term(var, "bin"), Term("ix", "bin")] + cv,
                       compute_vif=False)
        gi = get_or(fi, "ix")
        sg.append({"Subgroup": lab, "Level": "Interaction (ratio of per-grade ORs)", "Per grade OR (95% CI)": gi["txt"],
                   "p interaction": gi["p"], "p (text)": fmt_p(gi["p"])})
    sg = pd.DataFrame(sg)
    add_table("P2X_subgroups", sg, "Per-grade anaemia OR for any stroke within subgroups; interaction = product term.")
    RESULTS["p2x_sg"] = sg

    # 8. Time-to-event from the Hb measurement (prospective ordering guaranteed)
    import statsmodels.api as sm
    from .followup import encounters, treatments
    from .utils import design
    lc_, dth_ = encounters()
    t = d[d.anemia_cat.notna()].copy()
    t = t[~t.stroke_timing.isin([1, 2])]
    t = t[~((t.stroke_any == 1) & ~t.stroke_timing.isin([3]))]
    t["start"] = pd.concat([pd.to_datetime(t.index_date), pd.to_datetime(t.hgb_date, errors="coerce")], axis=1).max(axis=1)
    sdt = pd.to_datetime(t.stroke_date, errors="coerce")
    n_before_start = int(((t.stroke_any == 1) & (sdt <= t.start)).sum())
    t = t[~((t.stroke_any == 1) & (sdt <= t.start))]
    sdt = pd.to_datetime(t.stroke_date, errors="coerce")
    lastc = t.mrn.map(lc_)
    dth = t.mrn.map(dth_)
    cens = pd.concat([lastc, dth], axis=1).min(axis=1)
    t["ev_any"] = (t.stroke_any == 1).astype(int)
    t["ev_isch"] = ((t.stroke_any == 1) & (t.stroke_type == 1)).astype(int)
    t["end"] = np.where(t.ev_any == 1, sdt, cens)
    t["end"] = pd.to_datetime(t["end"])
    t["died"] = ((t.ev_any == 0) & dth.notna() & (pd.to_datetime(dth) <= t["end"])).astype(int)
    t = t[t.end.notna() & (t.end > t.start)]
    t["py"] = (t.end - t.start).dt.days / 365.25
    RESULTS["p2x_tte_cohort"] = dict(n=len(t), py=float(t.py.sum()), ev=int(t.ev_any.sum()), ev_isch=int(t.ev_isch.sum()),
                                     excl_before_start=n_before_start, fu_median=float(t.py.median()))

    def pois(data, oc, terms):
        data = data.dropna(subset=[x.var for x in terms])
        X, cm = design(data, terms)
        X = X.loc[:, (X != 0).any(axis=0)]
        r = sm.GLM(data[oc].values, X, family=sm.families.Poisson(), offset=np.log(data.py.values)).fit(cov_type="HC1")
        return r, len(data), int(data[oc].sum()), data[oc].sum() / (X.shape[1] - 1)

    def rr(r, c):
        b, se = r.params[c], r.bse[c]
        from scipy import stats as _st
        return dict(txt=fmt_or(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)), RR=np.exp(b),
                    lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=2 * _st.norm.sf(abs(b / se)))
    tte = []
    for oc, lab in [("ev_any", "Any stroke"), ("ev_isch", "Ischaemic stroke")]:
        g = t.groupby("anemia_cat").agg(ev=(oc, "sum"), py=("py", "sum"), n=("mrn", "size"))
        cvb = cov_for(t, oc, ["anemia_cat"], f"tte{oc}")
        ext_a2 = [Term(v, "bin") for v in pre_vars]
        fits_ = {}
        for nm, terms in [("Crude", [EXPO["anemia_cat"]]), ("Adjusted (Paper 2 covariates)", [EXPO["anemia_cat"]] + cvb),
                          ("Adjusted + pre-Hb conditions (2a)", [EXPO["anemia_cat"]] + cvb + ext_a2)]:
            fits_[nm] = pois(t, oc, terms)
            tr_terms = [Term("anemia", "cont") if x.var == "anemia_cat" else x for x in terms]
            rt, _, _, _ = pois(t, oc, tr_terms)
            fits_[nm + " trend"] = rr(rt, "anemia")
        for lv in ["None (Hb ≥12)"] + GRADES:
            row = {"Outcome": lab, "Anaemia grade": lv, "Women": int(g.n.get(lv, 0)), "Events": int(g.ev.get(lv, 0)),
                   "Person-years": round(float(g.py.get(lv, 0)), 1),
                   "Rate /1,000 PY": round(1000 * g.ev.get(lv, 0) / g.py.get(lv, 1), 2)}
            for nm in ["Crude", "Adjusted (Paper 2 covariates)", "Adjusted + pre-Hb conditions (2a)"]:
                r_, n_, e_, epv_ = fits_[nm]
                if lv == "None (Hb ≥12)":
                    row[nm + " RR"] = "1.00 (reference)"
                elif row["Events"] < 5:
                    row[nm + " RR"] = "not estimated (<5 events)"
                else:
                    q = rr(r_, f"anemia_cat={lv}")
                    row[nm + " RR"] = q["txt"]
                    if nm.startswith("Adjusted + pre"):
                        row.update({"RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "EPV": round(epv_, 1)})
            tte.append(row)
        tte.append({"Outcome": lab, "Anaemia grade": "Per grade (trend)",
                    **{nm + " RR": fits_[nm + " trend"]["txt"] + f"; p={fmt_p(fits_[nm + ' trend']['p'])}"
                       for nm in ["Crude", "Adjusted (Paper 2 covariates)", "Adjusted + pre-Hb conditions (2a)"]},
                    "RR": fits_["Adjusted + pre-Hb conditions (2a) trend"]["RR"],
                    "CI low": fits_["Adjusted + pre-Hb conditions (2a) trend"]["lo"],
                    "CI high": fits_["Adjusted + pre-Hb conditions (2a) trend"]["hi"],
                    "p": fits_["Adjusted + pre-Hb conditions (2a) trend"]["p"]})
    tte = pd.DataFrame(tte)
    add_table("P2X_time_to_event", tte, "Stroke rates after the Hb measurement by anaemia grade. Follow-up from the "
                                        "later of index and Hb date to stroke, death or last encounter; women with "
                                        "a stroke before or at that start excluded. Poisson rate ratios, HC1 SEs.")
    RESULTS["p2x_tte"] = tte
    log("Paper 2 ext", f"Time-to-event cohort: {RESULTS['p2x_tte_cohort']}")

    # 9. Anaemia treatment (exploratory landmark analysis)
    tx, audit = treatments()
    add_table("P2X_medication_classes", audit, "Medications Administered extract: classification of products "
                                              "(excluded = multivitamins/prenatal/OC iron placebo/spironolactone).")
    hbd_t = pd.to_datetime(t.hgb_date, errors="coerce")
    txm = tx.merge(pd.DataFrame({"mrn": t.mrn, "hbd": hbd_t}), on="mrn")
    win = txm[(txm.date >= txm.hbd - pd.Timedelta(days=30)) & (txm.date <= txm.hbd + pd.Timedelta(days=90))]
    for c_ in ["iv_iron", "oral_iron", "esa"]:
        t[f"tx_{c_}"] = t.mrn.isin(win.loc[win.cls == c_, "mrn"].unique()).astype(int)
    t["tx_any"] = t[["tx_iv_iron", "tx_oral_iron", "tx_esa"]].max(axis=1)
    lm = t.copy()
    lm["lstart"] = lm.start + pd.Timedelta(days=90)
    lm = lm[lm.end > lm.lstart]
    lm["py"] = (lm.end - lm.lstart).dt.days / 365.25
    lm["tx_group"] = np.select([lm.anemia_any == 0, lm.tx_any == 1], ["No anaemia", "Anaemia, treated (iron/ESA)"],
                               "Anaemia, no administered iron/ESA")
    tg = Term("tx_group", "cat", ref="No anaemia",
              levels=["No anaemia", "Anaemia, no administered iron/ESA", "Anaemia, treated (iron/ESA)"])
    txr = []
    for oc, lab in [("ev_any", "Any stroke"), ("ev_isch", "Ischaemic stroke")]:
        cvb = cov_for(lm, oc, ["tx_group"], f"lm{oc}")
        r_, n_, e_, epv_ = pois(lm, oc, [tg, Term("anemia", "cont")] + cvb)
        r2, _, _, _ = pois(lm, oc, [tg] + cvb)
        g = lm.groupby("tx_group").agg(ev=(oc, "sum"), py=("py", "sum"), n=("mrn", "size"))
        for lv in tg.levels:
            row = {"Outcome": lab, "Group": lv, "Women": int(g.n.get(lv, 0)), "Events": int(g.ev.get(lv, 0)),
                   "Rate /1,000 PY": round(1000 * g.ev.get(lv, 0) / g.py.get(lv, 1), 2)}
            if lv == tg.ref:
                row["Adjusted RR"] = "1.00 (reference)"
            elif row["Events"] < 5:
                row["Adjusted RR"] = "not estimated (<5 events)"
            else:
                q = rr(r2, f"tx_group={lv}")
                row.update({"Adjusted RR": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"]})
            txr.append(row)
        if all(f"tx_group={lv}" in r_.params.index for lv in tg.levels[1:]) and \
                int(g.ev.get("Anaemia, treated (iron/ESA)", 0)) >= 5:
            q = rr(r_, "tx_group=Anaemia, treated (iron/ESA)")
            q0 = rr(r_, "tx_group=Anaemia, no administered iron/ESA")
            L = np.zeros(len(r_.params))
            idx = list(r_.params.index)
            L[idx.index("tx_group=Anaemia, treated (iron/ESA)")] = 1
            L[idx.index("tx_group=Anaemia, no administered iron/ESA")] = -1
            b = float(L @ r_.params.values)
            se = float(np.sqrt(L @ r_.cov_params().values @ L))
            from scipy import stats as _st
            txr.append({"Outcome": lab, "Group": "Treated vs untreated anaemia, additionally adjusted for grade",
                        "Adjusted RR": fmt_or(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)),
                        "RR": np.exp(b), "CI low": np.exp(b - 1.96 * se), "CI high": np.exp(b + 1.96 * se),
                        "p": 2 * _st.norm.sf(abs(b / se))})
    txr = pd.DataFrame(txr)
    RESULTS["p2x_tx"] = txr
    RESULTS["p2x_tx_counts"] = {c_: int(lm.loc[lm.anemia_any == 1, f"tx_{c_}"].sum()) for c_ in ["iv_iron", "oral_iron", "esa"]}
    RESULTS["p2x_tx_n_anaemic"] = int((lm.anemia_any == 1).sum())
    add_table("P2X_treatment_landmark", txr, "EXPLORATORY. Iron/ESA administered from 30 d before to 90 d after the Hb "
                                             "(facility-administered only; outpatient oral iron not captured). "
                                             "Landmark at 90 d after follow-up start; strokes after the landmark. "
                                             "Poisson, HC1 SEs, Paper 2 covariates.")

    # 10. Continuous Hb: risk per 1 g/dL lower Hb (whole range; piecewise below/above 13 g/dL)
    from scipy import stats as _st2
    for dd_ in (d, t):
        dd_["hb_drop"] = -dd_["hgb"]
        dd_["hb_below13"] = np.where(dd_.hgb.notna(), np.maximum(13 - dd_.hgb, 0), np.nan)
        dd_["hb_above13"] = np.where(dd_.hgb.notna(), np.maximum(dd_.hgb - 13, 0), np.nan)
    ext_a3 = [Term(v, "bin") for v in pre_vars]
    per = []
    for o, lab in [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke")]:
        cv = cov_for(d, o, ["hgb"], "perhb")
        for adj, extra in [("Paper 2 covariates", []), ("+ pre-Hb conditions (2a)", ext_a3)]:
            f = register_fit(fit_logit(d, o, [Term("hb_drop")] + cv + extra, name=f"per g {o}"), "P2X",
                             f"Per 1 g/dL lower Hb: {lab} ({adj})")
            g = get_or(f, "hb_drop")
            per.append({"Design": "Cross-sectional (odds ratio)", "Outcome": lab, "Hb range": "Whole range",
                        "Adjustment": adj, "Estimate (95% CI)": g["txt"], "Est": g["OR"], "lo": g["lo"], "hi": g["hi"],
                        "p": g["p"], "N": f.n, "Events": f.events})
            f = register_fit(fit_logit(d, o, [Term("hb_below13"), Term("hb_above13")] + cv + extra, name="pw"), "P2X",
                             f"Per 1 g/dL, piecewise: {lab} ({adj})")
            for col, rng in [("hb_below13", "Below 13 g/dL (per 1 g/dL lower)"),
                             ("hb_above13", "Above 13 g/dL (per 1 g/dL higher)")]:
                g = get_or(f, col)
                per.append({"Design": "Cross-sectional (odds ratio)", "Outcome": lab, "Hb range": rng,
                            "Adjustment": adj, "Estimate (95% CI)": g["txt"], "Est": g["OR"], "lo": g["lo"],
                            "hi": g["hi"], "p": g["p"], "N": f.n, "Events": f.events})
    for oc, lab in [("ev_any", "Any stroke"), ("ev_isch", "Ischaemic stroke")]:
        cv = cov_for(t, oc, ["hgb"], f"perhbtte{oc}")
        for adj, extra in [("Paper 2 covariates", []), ("+ pre-Hb conditions (2a)", ext_a3)]:
            r_, n_, e_, _ = pois(t, oc, [Term("hb_drop")] + cv + extra)
            q = rr(r_, "hb_drop")
            per.append({"Design": "After the Hb measurement (rate ratio)", "Outcome": lab, "Hb range": "Whole range",
                        "Adjustment": adj, "Estimate (95% CI)": q["txt"], "Est": q["RR"], "lo": q["lo"],
                        "hi": q["hi"], "p": q["p"], "N": n_, "Events": e_})
            r_, n_, e_, _ = pois(t, oc, [Term("hb_below13"), Term("hb_above13")] + cv + extra)
            for col, rng in [("hb_below13", "Below 13 g/dL (per 1 g/dL lower)"),
                             ("hb_above13", "Above 13 g/dL (per 1 g/dL higher)")]:
                q = rr(r_, col)
                per.append({"Design": "After the Hb measurement (rate ratio)", "Outcome": lab, "Hb range": rng,
                            "Adjustment": adj, "Estimate (95% CI)": q["txt"], "Est": q["RR"], "lo": q["lo"],
                            "hi": q["hi"], "p": q["p"], "N": n_, "Events": e_})
    per = pd.DataFrame(per)
    per["p (text)"] = per["p"].map(fmt_p)
    add_table("P2X_per_g_dL", per, "Stroke risk per 1 g/dL lower haemoglobin. Piecewise-linear models allow separate "
                                   "slopes below and above 13 g/dL (the spline showed no association above 13 g/dL).")
    RESULTS["p2x_per"] = per
    # absolute rates by grade with exact Poisson CIs and crude 5-year risk (constant-rate approximation)
    ab = []
    for oc, lab in [("ev_any", "Any stroke"), ("ev_isch", "Ischaemic stroke")]:
        for lv in ["None (Hb ≥12)"] + GRADES:
            s_ = t[t.anemia_cat == lv]
            k, py = int(s_[oc].sum()), float(s_.py.sum())
            lo_ = _st2.chi2.ppf(0.025, 2 * k) / 2 / py if k > 0 else 0.0
            hi_ = _st2.chi2.ppf(0.975, 2 * k + 2) / 2 / py
            rt = k / py
            dk = int(s_.died.sum())
            ab.append({"Outcome": lab, "Anaemia grade": lv, "Women": len(s_), "Events": k, "Person-years": round(py, 1),
                       "Deaths without stroke": dk, "Death rate /1,000 PY": round(1000 * dk / py, 2),
                       "Median follow-up (y)": round(float(s_.py.median()), 2),
                       "Rate /1,000 PY": 1000 * rt, "Rate low": 1000 * lo_, "Rate high": 1000 * hi_,
                       "Rate (95% CI)": f"{1000 * rt:.2f} ({1000 * lo_:.2f}–{1000 * hi_:.2f})",
                       "Crude 5-y risk %": 100 * (1 - np.exp(-5 * rt)),
                       "Crude 5-y risk (95% CI)": f"{100 * (1 - np.exp(-5 * rt)):.1f}% "
                                                  f"({100 * (1 - np.exp(-5 * lo_)):.1f}–{100 * (1 - np.exp(-5 * hi_)):.1f})"})
    ab = pd.DataFrame(ab)
    add_table("P2X_absolute_rates", ab, "Crude stroke rates after the Hb measurement by anaemia grade (exact Poisson "
                                        "95% CI). 5-year risk = 1 − exp(−5 × rate), assuming a constant rate; unadjusted. Deaths "
                                        "without stroke end follow-up (competing risk).")
    RESULTS["p2x_abs"] = ab

    # 10b. Hb curves (RCS, same knots as the primary spline), adjusted for covariates + pre-Hb conditions (2a):
    # cross-sectional OR and rate ratio after the Hb measurement, both vs Hb 13 g/dL
    from .utils import rcs_basis
    kn = RESULTS["p2_spline_knots"]
    ref13 = rcs_basis(np.array([13.0]), kn)
    curves, cpts, cinfo = [], [], {}
    for design_, data_ in [("Cross-sectional (odds ratio)", d), ("After the Hb measurement (rate ratio)", t)]:
        oc = "stroke_any" if design_.startswith("Cross") else "ev_any"
        cv = cov_for(data_, oc, ["hgb"], f"curve{oc}")
        terms = [Term("hgb", "rcs", knots=kn)] + cv + ext_a3
        if design_.startswith("Cross"):
            f = register_fit(fit_logit(d, oc, terms, name="Hb curve 2a"), "P2X", "Hb spline, adjustment 2a")
            cols = [f"hgb__rcs{i}" for i in range(len(kn) - 1)]
            beta, V, n_, e_ = f.params[cols].values, f.cov.loc[cols, cols].values, f.n, f.events
        else:
            r_, n_, e_, _ = pois(t, oc, terms)
            cols = [f"hgb__rcs{i}" for i in range(len(kn) - 1)]
            beta, V = r_.params[cols].values, r_.cov_params().loc[cols, cols].values
        wov = float(beta @ np.linalg.solve(V, beta))
        bn, Vn = beta[1:], V[1:, 1:]
        wnl = float(bn @ np.linalg.solve(Vn, bn))
        p_ov, p_nl = _st2.chi2.sf(wov, len(beta)), _st2.chi2.sf(wnl, len(bn))
        src = d if design_.startswith("Cross") else t
        lo_, hi_ = np.nanpercentile(src.hgb, [1, 99])
        grid = np.linspace(lo_, hi_, 200)
        B = rcs_basis(grid, kn) - ref13
        lp = B @ beta
        se = np.sqrt(np.einsum("ij,jk,ik->i", B, V, B))
        curves.append(pd.DataFrame({"Design": design_, "Hb (g/dL)": grid, "Est": np.exp(lp),
                                    "lo": np.exp(lp - 1.96 * se), "hi": np.exp(lp + 1.96 * se)}))
        cinfo[design_] = dict(n=n_, events=e_, p_overall=p_ov, p_nonlin=p_nl)
        for h in [7, 8, 9, 10, 11, 12, 14, 15, 16]:
            b = (rcs_basis(np.array([float(h)]), kn) - ref13)[0]
            l, s_ = float(b @ beta), float(np.sqrt(b @ V @ b))
            cpts.append({"Design": design_, "Hb (g/dL)": h,
                         "Estimate vs 13 g/dL (95% CI)": fmt_or(np.exp(l), np.exp(l - 1.96 * s_), np.exp(l + 1.96 * s_)),
                         "Est": np.exp(l), "lo": np.exp(l - 1.96 * s_), "hi": np.exp(l + 1.96 * s_),
                         "N": n_, "Events": e_, "p overall": fmt_p(p_ov), "p non-linearity": fmt_p(p_nl)})
    cpts = pd.DataFrame(cpts)
    add_table("P2X_Hb_curves", cpts, "Any stroke by Hb (restricted cubic spline, knots "
              + ", ".join(f"{k:.1f}" for k in kn) + "), vs 13 g/dL, adjusted for Paper 2 covariates and pre-Hb "
              "conditions (2a). Cross-sectional: logistic OR; after the Hb measurement: Poisson rate ratio.")
    RESULTS["p2x_curves"] = pd.concat(curves, ignore_index=True)
    RESULTS["p2x_curve_pts"] = cpts
    RESULTS["p2x_curve_info"] = cinfo

    # 10c. time-to-event with the illness-marker exclusions (adjustment 2a)
    t["death_1y"] = ((pd.to_datetime(t.mrn.map(_dth)).notna()) &
                     ((pd.to_datetime(t.mrn.map(_dth)) - t.start).dt.days <= 365)).astype(float)
    t["serious_illness"] = t[["malignancy_ever", "chf", "ckd", "liver_any", "hiv"]].max(axis=1)
    ill = []
    for excl, lab in [(None, "All women (as Table 2)"), ("death_1y", "Excluding deaths within 1 y of Hb"),
                      ("serious_illness", "Excluding cancer, heart failure, CKD, liver disease, HIV")]:
        tt_ = t if excl is None else t[t[excl] == 0]
        ea = [x for x in ext_a3 if tt_[x.var].nunique() > 1]
        cvb = cov_for(tt_, "ev_any", ["anemia_cat"], f"ill{excl}")
        r_, n_, e_, epv_ = pois(tt_, "ev_any", [EXPO["anemia_cat"]] + cvb + ea)
        rt_, _, _, _ = pois(tt_, "ev_any", [Term("anemia", "cont")] + cvb + ea)
        m_, pg_ = rr(r_, "anemia_cat=Moderate (8–9.9)"), rr(rt_, "anemia")
        sv_ = rr(r_, "anemia_cat=Severe (<8)")
        ill.append({"Analysis": lab, "Women": n_, "Strokes": e_, "EPV": round(float(epv_), 1),
                    "Moderate RR (95% CI)": m_["txt"], "Severe RR (95% CI)": sv_["txt"], "Per grade RR (95% CI)": pg_["txt"],
                    "Mod RR": m_["RR"], "Mod lo": m_["lo"], "Mod hi": m_["hi"],
                    "PG RR": pg_["RR"], "PG lo": pg_["lo"], "PG hi": pg_["hi"], "PG p": pg_["p"]})
    ill = pd.DataFrame(ill)
    add_table("P2X_illness_marker_tte", ill, "Rate ratios for any stroke after the Hb measurement (Poisson, adjusted "
              "for Paper 2 covariates and pre-Hb conditions 2a) after excluding women who died within 1 year of the Hb "
              "or who had serious chronic illness recorded at any time. Healthcare-use counts are not available "
              "(one encounter row per patient).")
    RESULTS["p2x_ill"] = ill

    # 10d. time-to-event by laboratory anaemia pattern (adjustment 2a)
    ptt = []
    for oc, lab in [("ev_any", "Any stroke"), ("ev_isch", "Ischaemic stroke")]:
        cvb = cov_for(t, oc, ["anaemia_pattern"], f"pattte{oc}")
        tt_ = t.dropna(subset=["anaemia_pattern"])
        r_, n_, e_, epv_ = pois(tt_, oc, [pt] + cvb + ext_a3)
        for lv in PAT:
            s_ = tt_[tt_.anaemia_pattern == lv]
            k, py = int(s_[oc].sum()), float(s_.py.sum())
            row = {"Outcome": lab, "Pattern": lv, "Women": len(s_), "Strokes": k, "Person-years": round(py, 1),
                   "Rate /1,000 PY": round(1000 * k / py, 2) if py > 0 else np.nan}
            if lv == "No anaemia":
                row["Adjusted RR (95% CI)"] = "1.00 (reference)"
            elif k < 5:
                row["Adjusted RR (95% CI)"] = f"not estimated (<5 strokes)"
            else:
                q = rr(r_, f"anaemia_pattern={lv}")
                row.update({"Adjusted RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"],
                            "p": q["p"]})
            ptt.append(row)
    ptt = pd.DataFrame(ptt)
    ptt["p (text)"] = ptt["p"].map(lambda v: fmt_p(v) if pd.notna(v) else "") if "p" in ptt else ""
    add_table("P2X_pattern_tte", ptt, "Stroke rate after the Hb measurement by laboratory anaemia pattern (Poisson, "
              "adjusted for Paper 2 covariates and pre-Hb conditions 2a). Patterns with <5 strokes not estimated.")
    RESULTS["p2x_pattern_tte"] = ptt

    # 7. E-values
    tab = pd.DataFrame(rows)
    evr = []
    for _, r in tab[tab.Section.str.startswith(("0.", "1.", "2."))].iterrows():
        for key, o, lo, hi in [("Moderate", r.get("Moderate OR"), r.get("Moderate lo"), r.get("Moderate hi")),
                               ("Severe", r.get("Severe OR"), r.get("Severe lo"), r.get("Severe hi")),
                               ("Per grade", r["trend OR"], r["trend lo"], r["trend hi"])]:
            if pd.isna(o):
                continue
            e_pt, e_ci = evalue(o, lo, hi)
            evr.append({"Analysis": r["Analysis"], "Contrast": key, "OR (95% CI)": fmt_or(o, lo, hi),
                        "E-value (point)": round(e_pt, 2), "E-value (CI limit)": round(e_ci, 2)})
    evr = pd.DataFrame(evr)
    add_table("P2X_evalues", evr, "E-values: minimum strength of association (risk-ratio scale) an unmeasured "
                                  "confounder would need with both anaemia and stroke to explain away the estimate.")
    RESULTS["p2x_ev"] = evr
    cols = ["Section", "Analysis", "N", "Events", "Mild", "Moderate", "Severe", "Per grade (trend)", "p-trend", "EPV",
            "Flags", "Mild OR", "Mild lo", "Mild hi", "Moderate OR", "Moderate lo", "Moderate hi", "Severe OR",
            "Severe lo", "Severe hi", "trend OR", "trend lo", "trend hi", "trend p"]
    tab = tab[[c for c in cols if c in tab.columns]]
    add_table("P2X_sensitivity_summary", tab, "Anaemia grade ORs (vs Hb ≥12) across temporality, extended-adjustment "
                                              "and restricted-cohort analyses. Logistic regression, HC1 SEs.")
    RESULTS["p2x_tab"] = tab
    RESULTS["_p2x_t"] = t
    return d
