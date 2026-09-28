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
HEREDITARY = ["sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy"]
CHRONIC = ["ckd", "ckd_esrd", "cirrhosis", "chronic_liver", "alcohol_disorder", "malabsorption", "menopause"]
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
          "pregnancy_recent": "Pregnancy/delivery within 1 y before Hb", "hiv": "HIV infection",
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
                                    "gi_bleed_recent", "malabsorption", "hiv", "pregnancy_recent",
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
    pre_vars = ["cad", "haemoglobinopathy_any", "ckd", "liver_any", "alcohol_disorder", "gi_bleed_recent",
                "malabsorption", "hiv", "pregnancy_recent"]
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
    return d
