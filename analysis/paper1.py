"""Paper 1 - Covert (silent) brain infarcts among women with benign uterine
pathology who had brain imaging reviewed.
"""
import re

import numpy as np
import pandas as pd

from .utils import (LINEARITY, RESULTS, WORKBOOK, Term, add_table, fit_logit, fmt_or, fmt_p, fmt_prev,
                    get_or, linearity_rows, log, mice, or_table, pool_fits, register_fit, table1, wilson)

# ---------------------------------------------------------------------------
# Imaging-only flags built from the notes column
# ---------------------------------------------------------------------------
TAGS = {
    "tag_acute_features": r"no qualifying diagnosis code \(acute features\)",
    "tag_chronic_only": r"no qualifying diagnosis code \(chronic/incidental only\)",
    "tag_found_on_imaging": r"CASE FOUND ON IMAGING ONLY",
    "tag_revoked": r"REVOKED imaging-only",
    "tag_excluded": r"imaging-only case EXCLUDED",
    "tag_undecided": r"imaging-only case UNDECIDED",
    "tag_incl_D_acute": r"INCLUDED \(D\. Acute",
    "tag_incl_E_hedged": r"INCLUDED \(E\. Chronic infarct, hedged",
    "tag_incl_C_haem": r"INCLUDED \(C - spontaneous ICH/SAH\)",
    "cc1_acute_new": r"imaging-only case: acute or new infarct documented",
    "cc3_symptoms": r"imaging-only case: infarct on imaging with stroke symptoms documented",
    "cc4_multi_no_sx": r"imaging-only case: multiple or territorial infarcts on imaging without documented symptoms",
    "cc5_incidental": r"imaging-only case: chronic/incidental brain infarct found on imaging ordered for another indication",
}
IO_ANY = "|".join(TAGS.values()) + r"|radiological ischaemia"

SUBCATS = [  # (name, covert?)
    ("Imaging-only ICH/SAH", False),
    ("Imaging-only acute/new infarct", False),
    ("Imaging-only infarct with stroke symptoms", False),
    ("Covert: incidental, scan for another indication", True),
    ("Covert: multiple/territorial, no symptoms", True),
    ("Covert: chronic infarct, hedged wording", True),
]


# PI decision 2026-09-28: imaging-only cases that were revoked, excluded or left
# undecided at adjudication are NOT revoked for Paper 1 when the imaging text
# describes a brain lesion - they count as covert infarcts. The remaining
# revoked/excluded/undecided patients (non-brain organ, truncated/negated
# template, or organ not stated) count as imaged without infarct. All of them
# had imaging reviewed, so all enter the Paper 1 denominator. stroke_any is
# not changed (Papers 2-3 unaffected).
RESTORED = "Covert: restored after revocation/exclusion (brain wording)"
BRAIN_RE = (r"cerebr|brain|lacun|subcortical|cortical infarct|cerebell|white matter|encephalomalacia|temporal|"
            r"frontal lobe|demyelinat")
REVIEWED_NEG = ["revoked", "excluded", "undecided (stroke_any_incl_imaging only)"]


def classify(df):
    d = df.copy()
    notes = d["notes"].fillna("")
    for k, p in TAGS.items():
        d[k] = notes.str.contains(p, regex=True).astype(int)
    d["io_any"] = notes.str.contains(IO_ANY, regex=True).astype(int)
    d["io_included"] = ((d["io_any"] == 1) & (d["stroke_any"] == 1)).astype(int)
    sub = np.select(
        [d.tag_incl_C_haem == 1,
         (d.tag_incl_D_acute == 1) | (d.cc1_acute_new == 1) | (d.tag_acute_features == 1),
         d.cc3_symptoms == 1,
         d.cc5_incidental == 1,
         d.cc4_multi_no_sx == 1,
         d.tag_incl_E_hedged == 1],
        [s for s, _ in SUBCATS], default="UNCLASSIFIED")
    d["io_subcat"] = np.where(d["io_included"] == 1, sub, "")
    covert_names = [s for s, c in SUBCATS if c]
    d["covert"] = d["io_subcat"].isin(covert_names).astype(int)
    d["covert_strict"] = (d["io_subcat"] == "Covert: incidental, scan for another indication").astype(int)
    d["io_status"] = np.select(
        [d.io_included == 1, d.tag_undecided == 1, d.tag_revoked == 1, d.tag_excluded == 1, d.io_any == 1],
        ["included in stroke_any", "undecided (stroke_any_incl_imaging only)", "revoked", "excluded",
         "flagged, other"], default="")
    ev = notes.str.extract(r"imaging: (.*?)(?: \| |$)")[0].fillna("").str.lower()
    d["io_reviewed_neg"] = d["io_status"].isin(REVIEWED_NEG).astype(int)
    d["io_restored"] = ((d["io_reviewed_neg"] == 1) & ev.str.contains(BRAIN_RE, regex=True)).astype(int)
    d.loc[d.io_restored == 1, "io_subcat"] = RESTORED
    d["covert"] = ((d["covert"] == 1) | (d["io_restored"] == 1)).astype(int)
    d["imaged"] = (d["stroke_confirmed_imaging"].isin([0, 1]) | (d["io_reviewed_neg"] == 1)).astype(int)
    d["p1_group"] = np.select(
        [d.imaged == 0,
         d.covert == 1,
         (d.stroke_any == 0),
         (d.stroke_any == 1) & (d.io_included == 0),
         d.io_included == 1],
        ["not imaged", "Covert infarct", "Imaged, no infarct", "Clinical stroke",
         "Imaging-only acute/symptomatic/haemorrhagic"], default="?")
    return d


def read_sheet(name):
    x = pd.read_excel(WORKBOOK, sheet_name=name)
    x.columns = [str(c).lower() for c in x.columns]
    return x


# ---------------------------------------------------------------------------
# Pre-specified model terms (fixed BEFORE fitting; see brief, Paper 1)
# ---------------------------------------------------------------------------
# Full pre-specified covariate set for covert infarct vs imaged-negative:
#   age (linear, per year), race (White ref / Black / Asian / Other-unknown),
#   BMI (linear, per kg/m2), hypertension, diabetes, dyslipidaemia,
#   smoking (Never ref / Ever / Unknown - 3 level, not imputed),
#   migraine (any diagnosis, codes 1/2/9 = yes; PI decision - type/aura not used),
#   thrombophilia, atrial fibrillation,
#   uterine_dx_group (Fibroids only ref / Adenomyosis only / Endometriosis only / >1),
#   anaemia grade (None ref / mild / moderate / severe).
FULL_TERMS = [
    Term("age_index", "cont", label="Age (per year)"),
    Term("race4", "cat", ref="White", levels=["White", "Black", "Asian", "Other/unknown"], label="Race"),
    Term("bmi", "cont", label="BMI (per kg/m²)"),
    Term("htn", "bin", label="Hypertension"),
    Term("dm", "bin", label="Diabetes"),
    Term("dyslipidemia", "bin", label="Dyslipidaemia"),
    Term("smoking3", "cat", ref="Never", levels=["Never", "Ever", "Unknown"], label="Smoking"),
    Term("migraine_any", "bin", label="Migraine (any type)"),
    Term("thrombophilia", "bin", label="Thrombophilia"),
    Term("afib", "bin", label="Atrial fibrillation"),
    Term("dxgrp", "cat", ref="Fibroids only",
         levels=["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"],
         label="Uterine diagnosis group"),
    Term("anemia_cat", "cat", ref="None (Hb ≥12)",
         levels=["None (Hb ≥12)", "Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"], label="Anaemia"),
]
# Pre-specified simplification if the full model has EPV < 10:
#   CORE = age, hypertension, diabetes, dyslipidaemia, smoking (3-level), migraine (3-level)
#   (8 df), then each remaining pre-specified factor added ONE AT A TIME to CORE.
CORE_VARS = ["age_index", "htn", "dm", "dyslipidemia", "smoking3", "migraine_any"]
ADDON_VARS = ["race4", "bmi", "thrombophilia", "afib", "dxgrp", "anemia_cat"]
# Pre-specified second-level simplification: if a CORE+factor model has EPV<10
# or a level with <5 events, a categorical factor is refitted in collapsed form
# (race: White/Black/Asian-Other-unknown; dx group: fibroids only /
# endometriosis only / adenomyosis only or >1; anaemia: none/mild/moderate-
# severe). A binary factor with <5 events in either cell is not estimated.
COLLAPSE = {
    "race4": Term("race3", "cat", ref="White", levels=["White", "Black", "Asian/Other/unknown"], label="Race (collapsed)"),
    "dxgrp": Term("dxgrp3", "cat", ref="Fibroids only",
                  levels=["Fibroids only", "Endometriosis only", "Adenomyosis only or >1"],
                  label="Uterine diagnosis group (collapsed)"),
    "anemia_cat": Term("anemia3", "cat", ref="None (Hb ≥12)",
                       levels=["None (Hb ≥12)", "Mild (10–11.9)", "Moderate/severe (<10)"], label="Anaemia (collapsed)"),
}


def terms_for(vars_):
    return [t for t in FULL_TERMS if t.var in vars_]


def core_plus_one(df, outcome, label, paper="P1", force=()):
    """CORE model, then CORE + each add-on factor. Returns tidy OR table
    (CORE estimates from the CORE model; each add-on estimate from its own
    CORE+factor model) and the list of fits."""
    rows = []
    fits = {}
    f = register_fit(fit_logit(df, outcome, terms_for(CORE_VARS), name=f"{label}: CORE"), paper, f"{label}: CORE")
    fits["CORE"] = f
    t = or_table(f, model_label=f"{label}: CORE")
    t["Estimate from"] = "CORE model"
    rows.append(t)
    for v in ADDON_VARS:
        f2 = register_fit(fit_logit(df, outcome, terms_for(CORE_VARS + [v]), name=f"{label}: CORE + {v}"),
                          paper, f"{label}: CORE + {v}")
        fits[v] = f2
        t2 = or_table(f2, only_vars=[v], model_label=f"{label}: CORE + {v}")
        t2["Estimate from"] = f"CORE + {v}"
        unstable = f2.epv < 10 or bool(f2.sparse) or not f2.converged or v in force
        if unstable and v in COLLAPSE:
            ct = COLLAPSE[v]
            f3 = register_fit(fit_logit(df, outcome, terms_for(CORE_VARS) + [ct],
                                        name=f"{label}: CORE + {ct.var}"), paper, f"{label}: CORE + {ct.var} (collapsed)")
            fits[v + "_collapsed"] = f3
            log("Paper 1", f"{label}: CORE + {v} unstable ({'; '.join(f2.flags)}) -> collapsed to {ct.var}; "
                           f"new flags: {f3.flags or 'none'}")
            t2 = or_table(f3, only_vars=[ct.var], model_label=f"{label}: CORE + {ct.var}")
            t2["Estimate from"] = f"CORE + {ct.var} (collapsed because CORE + {v}: {'; '.join(f2.flags)})"
        elif unstable and v not in COLLAPSE:
            log("Paper 1", f"{label}: CORE + {v} unstable ({'; '.join(f2.flags)}) -> estimate not reported")
            t2 = pd.DataFrame([{"Model": f"{label}: CORE + {v}", "Term": t2["Term"].iloc[0],
                                "Level": "yes vs no", "OR (95% CI)": "not estimated (sparse)",
                                "n in level": t2["n in level"].iloc[0], "events in level": t2["events in level"].iloc[0],
                                "N": f2.n, "Events": f2.events, "EPV": round(f2.epv, 1), "Flags": "; ".join(f2.flags)}])
            t2["Estimate from"] = "not reported"
        rows.append(t2)
    return pd.concat(rows, ignore_index=True), fits


def run(df_all, elig):
    d = classify(elig)
    dall = classify(df_all)

    # ------------------------------------------------------------------ flags
    RESULTS["p1_io_status"] = d.loc[d.io_any == 1, "io_status"].value_counts().to_dict()
    RESULTS["p1_io_subcat"] = d.loc[d.io_included == 1, "io_subcat"].value_counts().to_dict()
    n_uncl = int((d.io_subcat == "UNCLASSIFIED").sum())
    RESULTS["p1_unclassified"] = n_uncl
    log("Paper 1", f"Imaging-only notes flags (eligible): {RESULTS['p1_io_status']}")
    log("Paper 1", f"Included imaging-only subcategories (eligible): {RESULTS['p1_io_subcat']}; unclassified {n_uncl}")
    assert (d.loc[d.io_included == 1, "stroke_confirmed_imaging"] == 1).all()

    # ------------------------------------------------------ reconciliation
    rec = []
    sheets = {"imaging_only_candidates": read_sheet("imaging_only_candidates"),
              "imaging_only_review": read_sheet("imaging_only_review"),
              "deferred_for_review": read_sheet("deferred_for_review")}
    key = dall.set_index("mrn")  # internal linkage only; never written out
    for nm, sh in sheets.items():
        m = key.loc[key.index.intersection(sh["mrn"])]
        rec.append({"Source": nm, "Scope": "all records", "Rows in sheet": len(sh),
                    "Linked to analysis_master": len(m),
                    "…with any imaging-only note": int(m.io_any.sum()),
                    "…included in stroke_any": int(m.io_included.sum()),
                    "…eligible": int((m.eligible == 1).sum()),
                    "…eligible & included": int(((m.eligible == 1) & (m.io_included == 1)).sum()),
                    "…eligible & covert (primary)": int(((m.eligible == 1) & (m.covert == 1)).sum())})
    in_any_sheet = dall["mrn"].isin(pd.concat([s["mrn"] for s in sheets.values()]))
    for scope, x, sel in [("all records", dall, in_any_sheet), ("eligible", d, in_any_sheet[d.index])]:
        rec.append({"Source": "notes column (imaging-only flag)", "Scope": scope,
                    "Rows in sheet": int(x.io_any.sum()), "Linked to analysis_master": int(x.io_any.sum()),
                    "…with any imaging-only note": int(x.io_any.sum()),
                    "…included in stroke_any": int(x.io_included.sum()),
                    "…eligible": int(((x.io_any == 1) & (x.eligible == 1)).sum()),
                    "…eligible & included": int(((x.eligible == 1) & (x.io_included == 1)).sum()),
                    "…eligible & covert (primary)": int(((x.eligible == 1) & (x.covert == 1)).sum()),
                    "Flagged in notes but in none of the 3 sheets": int((x.io_any.astype(bool) & ~sel).sum())})
    rec = pd.DataFrame(rec)
    add_table("P1_reconciliation", rec, "Imaging-only flag built from notes vs the three review sheets.")
    RESULTS["p1_rec"] = rec

    # cross-tab: review category x notes classification (eligible)
    rv = sheets["imaging_only_review"].merge(dall[["mrn", "eligible", "io_status", "io_subcat"]], on="mrn", how="left")
    rv = rv[rv.eligible == 1]
    ct1 = pd.crosstab(rv["category"], rv["io_status"] + " | " + rv["io_subcat"].replace("", "-")).reset_index()
    add_table("P1_recon_review_x_notes", ct1, "imaging_only_review category vs notes-based status (eligible).")
    df_ = sheets["deferred_for_review"].merge(dall[["mrn", "eligible", "io_subcat", "io_status"]], on="mrn", how="left")
    df_ = df_[df_.eligible == 1]
    ct2 = pd.crosstab(df_["clinical_category"], df_["io_status"] + " | " + df_["io_subcat"].replace("", "-")).reset_index()
    add_table("P1_recon_deferred_x_notes", ct2, "deferred_for_review clinical category vs notes-based status (eligible).")
    ct3 = pd.crosstab(df_["scan_indication_class"], df_["io_subcat"].replace("", "(not included)")).reset_index()
    add_table("P1_recon_indication", ct3, "Scan indication class (deferred_for_review) by covert subcategory (eligible).")

    # revoked/excluded/undecided imaging-only: blank stroke_confirmed_imaging although imaged
    n_rev_blank = int(((d.io_reviewed_neg == 1) & d.stroke_confirmed_imaging.isna()).sum())
    RESULTS["p1_rev_excl_not_in_denominator"] = n_rev_blank
    RESULTS["p1_n_reviewed_neg"] = int(d.io_reviewed_neg.sum())
    RESULTS["p1_n_restored"] = int(d.io_restored.sum())
    RESULTS["p1_restored_by_status"] = d.loc[d.io_restored == 1, "io_status"].value_counts().to_dict()
    log("Data problem", f"{n_rev_blank} eligible patients whose imaging-only infarct was revoked/excluded/undecided "
                        f"have blank stroke_confirmed_imaging although their imaging was reviewed.")
    log("PI decision", f"All {RESULTS['p1_n_reviewed_neg']} revoked/excluded/undecided patients added to the Paper 1 "
                       f"imaged denominator; {RESULTS['p1_n_restored']} whose imaging text describes a brain lesion "
                       f"({RESULTS['p1_restored_by_status']}) counted as covert infarcts, the rest as imaged without "
                       f"infarct. stroke_any unchanged.")
    rv_ = d[d.io_reviewed_neg == 1]
    rest_tab = pd.crosstab(rv_["io_status"], np.where(rv_["io_restored"] == 1, "counted as covert infarct (brain wording)",
                                                      "counted as imaged, no infarct")).reset_index()
    add_table("P1_restored_cases", rest_tab, "Revoked/excluded/undecided imaging-only patients (eligible): PI decision - "
                                             "brain-lesion wording = covert infarct; all enter the imaged denominator.")
    n_desc_blank = int(((d.io_included == 1) & d.location.isna()).sum())
    RESULTS["p1_io_desc_blank"] = n_desc_blank
    RESULTS["p1_io_date_blank"] = int(((d.io_included == 1) & d.stroke_date.isna()).sum())
    log("Data problem", f"{n_desc_blank} of {int(d.io_included.sum())} eligible imaging-only cases have blank location/"
                        f"vascular_territory/num_infarct in analysis_master (and {RESULTS['p1_io_date_blank']} have blank "
                        f"stroke_date, hence blank stroke_timing). Descriptors are taken from imaging_only_candidates "
                        f"for the covert description; master fields reported alongside.")

    # ------------------------------------------------------------ population
    im = d[d.imaged == 1].copy()
    RESULTS["p1_n_imaged"] = len(im)
    grp = im["p1_group"].value_counts().to_dict()
    RESULTS["p1_groups"] = grp
    log("Paper 1", f"Imaged denominator {len(im):,}; groups {grp}")
    n_cov = grp.get("Covert infarct", 0)
    RESULTS["p1_prev_imaged"] = wilson(n_cov, len(im))
    RESULTS["p1_prev_imaged_txt"] = fmt_prev(n_cov, len(im))
    n_nonclin = n_cov + grp.get("Imaged, no infarct", 0)
    RESULTS["p1_prev_nonclin_txt"] = fmt_prev(n_cov, n_nonclin)
    n_strict = int(im.covert_strict.sum())
    RESULTS["p1_n_strict"] = n_strict
    RESULTS["p1_prev_strict_txt"] = fmt_prev(n_strict, len(im))
    RESULTS["p1_prev_cohort_txt"] = fmt_prev(n_cov, len(d))
    prev_rows = [
        {"Definition": "Covert infarct (primary)", "Denominator": "all imaged", "Text": RESULTS["p1_prev_imaged_txt"]},
        {"Definition": "Covert infarct (primary)", "Denominator": "imaged, excluding clinical and other imaging-only strokes",
         "Text": RESULTS["p1_prev_nonclin_txt"]},
        {"Definition": "Covert infarct (strict: incidental on scan for another indication)", "Denominator": "all imaged",
         "Text": RESULTS["p1_prev_strict_txt"]},
        {"Definition": "Covert infarct (primary)", "Denominator": "whole eligible cohort (lower bound; most never imaged)",
         "Text": RESULTS["p1_prev_cohort_txt"]},
    ]
    for g in ["18–39", "40–60"]:
        s = im[im.age2 == g]
        prev_rows.append({"Definition": f"Covert infarct (primary), age {g}", "Denominator": "all imaged",
                          "Text": fmt_prev(int(s.covert.sum()), len(s))})
    for g in [0.0, 1.0]:
        s = im[im.migraine_any == g]
        prev_rows.append({"Definition": f"Covert infarct (primary), migraine {'yes' if g else 'no'}", "Denominator": "all imaged",
                          "Text": fmt_prev(int(s.covert.sum()), len(s))})
    add_table("P1_prevalence", pd.DataFrame(prev_rows), "Covert infarct prevalence with Wilson 95% CI.")
    sub_tab = im[im.io_included == 1]["io_subcat"].value_counts().rename_axis("Imaging-only subcategory").reset_index(name="n")
    sub_tab["Covert (primary definition)"] = sub_tab["Imaging-only subcategory"].str.startswith("Covert")
    add_table("P1_imaging_only_subcats", sub_tab, "Imaging-only included cases among imaged eligible patients.")

    # ------------------------------------------------------------ Table 1
    groups = ["Imaged, no infarct", "Covert infarct", "Clinical stroke"]
    spec = [("age_index", "Age at index, years", "cont"), ("race4", "Race", "cat"), ("bmi", "BMI, kg/m²", "cont"),
            ("obesity", "Obesity (BMI ≥30)", "bin"), ("htn", "Hypertension", "bin"), ("dm", "Diabetes", "bin"),
            ("dyslipidemia", "Dyslipidaemia", "bin"), ("cad", "Coronary artery disease", "bin"),
            ("afib", "Atrial fibrillation", "bin"), ("smoking3", "Smoking", "cat"),
            ("migraine_any", "Migraine (any type)", "bin"), ("thrombophilia", "Thrombophilia", "bin"),
            ("vte_history", "Prior VTE", "bin"), ("dxgrp", "Uterine diagnosis group", "cat"),
            ("anemia_cat", "Anaemia grade", "cat"), ("uterine_bleeding", "Heavy/abnormal uterine bleeding", "bin"),
            ("hormonal_tx", "Hormonal therapy", "bin")]
    t1 = table1(im, "p1_group", groups, spec)
    add_table("P1_Table1", t1, "Characteristics of imaged patients by group. Imaging-only acute/symptomatic/"
                               "haemorrhagic cases are excluded from the three columns (listed in P1_imaging_only_subcats).")

    # ------------------------------------------------------------ Imaging selection
    sel_rows = []
    for var, lab in [("migraine_any", "Migraine"), ("age2", "Age"), ("htn", "Hypertension"), ("race4", "Race"),
                     ("dxgrp", "Uterine dx group")]:
        for lv, s in d.groupby(var):
            sel_rows.append({"Factor": lab, "Level": lv, "n": len(s), "Imaged n": int(s.imaged.sum()),
                             "Imaged %": round(100 * s.imaged.mean(), 1)})
    add_table("P1_imaging_selection", pd.DataFrame(sel_rows), "Proportion of the eligible cohort with brain imaging reviewed.")
    fsel = register_fit(fit_logit(d, "imaged", [Term("age_index"), FULL_TERMS[1], FULL_TERMS[7]],
                                  name="P(imaged) ~ age + race + migraine"), "P1", "Selection: P(imaged)")
    sel_or = or_table(fsel, model_label="Selection model: odds of being imaged")
    add_table("P1_selection_model", sel_or)
    RESULTS["p1_sel_mig"] = get_or(fsel, "migraine_any")
    RESULTS["p1_imaged_pct_by_mig"] = {r["Level"]: r["Imaged %"] for r in sel_rows if r["Factor"] == "Migraine"}

    # ------------------------------------------------------------ Models
    ana = im[im.p1_group.isin(["Imaged, no infarct", "Covert infarct"])].copy()
    ana["y_covert"] = ana["covert"].astype(float)
    ffull = register_fit(fit_logit(ana, "y_covert", FULL_TERMS, name="Full pre-specified"), "P1",
                         "Main: full pre-specified (covert vs imaged-negative)")
    RESULTS["p1_full_epv"] = ffull.epv
    RESULTS["p1_full_n"], RESULTS["p1_full_events"], RESULTS["p1_full_df"] = ffull.n, ffull.events, ffull.df_model
    RESULTS["p1_full_flags"] = ffull.flags
    full_tab = or_table(ffull, model_label="Full pre-specified model (UNSTABLE if EPV<10 - not for reporting)")
    log("Paper 1", f"Full pre-specified model: n={ffull.n}, events={ffull.events}, df={ffull.df_model}, "
                   f"EPV={ffull.epv:.1f}; flags: {ffull.flags}")
    use_core = ffull.epv < 10 or not ffull.converged or bool(ffull.sparse)
    RESULTS["p1_use_core"] = use_core
    main_tab, main_fits = core_plus_one(ana, "y_covert", "Main (covert vs imaged-negative)")
    add_table("P1_main_model", main_tab,
              "Primary reported model. CORE = age, HTN, DM, dyslipidaemia, smoking(3), migraine(3); other "
              "pre-specified factors added one at a time. Logistic regression, HC1 robust SEs.")
    add_table("P1_full_model_unstable", full_tab,
              "Full pre-specified model shown for transparency only; see Flags/EPV.")
    RESULTS["p1_core"] = main_fits["CORE"]
    RESULTS["p1_main_fits"] = main_fits

    # linearity for continuous terms
    LINEARITY.extend(linearity_rows(ana, "y_covert", terms_for(CORE_VARS + ["bmi"]), ["age_index", "bmi"],
                                    "P1", "CORE + BMI"))

    # ---- Sensitivity A: covert vs whole eligible cohort without any stroke
    whole = d[(d.covert == 1) | (d.stroke_any == 0)].copy()
    whole["y_covert"] = whole["covert"].astype(float)
    # same simplifications as the main model so the two comparators are like-for-like
    forced = [v for v in ADDON_VARS if v + "_collapsed" in main_fits or
              (v not in COLLAPSE and (main_fits[v].sparse or main_fits[v].epv < 10))]
    whole_tab, whole_fits = core_plus_one(whole, "y_covert", "Sens A (covert vs all eligible non-stroke)", force=forced)
    RESULTS["p1_whole_fits"] = whole_fits
    RESULTS["p1_whole_n"] = len(whole)
    # comparison table (same terms) main vs whole
    comp = []
    for key_ in ["CORE"] + ADDON_VARS:
        kk = key_ + "_collapsed" if key_ + "_collapsed" in main_fits else key_
        fm, fw = main_fits[kk], whole_fits[kk]
        cols = [c for c in fm.params.index if c != "const"]
        if key_ != "CORE":
            cols = fm.colmap[COLLAPSE[key_].var if kk.endswith("_collapsed") else key_]
        for c in cols:
            a, b = get_or(fm, c), get_or(fw, c) if c in fw.params.index else None
            if key_ != "CORE" and key_ not in COLLAPSE and (fm.sparse or fw.sparse):
                comp.append({"Term": c, "Model": f"CORE + {kk}", "Note": "; ".join(fm.flags),
                             "OR vs imaged-negative": "not estimated (sparse)",
                             "OR vs whole eligible cohort": "not estimated (sparse)"})
                continue
            comp.append({"Term": c, "Model": "CORE" if key_ == "CORE" else f"CORE + {kk}",
                         "Note": "; ".join(fm.flags),
                         "OR vs imaged-negative": a["txt"], "OR vs whole eligible cohort": b["txt"] if b else "",
                         "OR_imaged": a["OR"], "OR_whole": b["OR"] if b else np.nan,
                         "lo_imaged": a["lo"], "hi_imaged": a["hi"],
                         "lo_whole": b["lo"] if b else np.nan, "hi_whole": b["hi"] if b else np.nan,
                         "Ratio of ORs (whole/imaged)": (b["OR"] / a["OR"]) if b else np.nan})
    comp = pd.DataFrame(comp)
    add_table("P1_sens_selection_bias", comp,
              "How estimates change when the comparator is the whole eligible cohort (never-imaged included).")
    add_table("P1_sensA_model", whole_tab)
    RESULTS["p1_comp"] = comp

    # ---- Sensitivity B: strict covert definition (incidental on scan for another indication only)
    strict = im[(im.p1_group == "Imaged, no infarct") | (im.covert_strict == 1)].copy()
    strict["y"] = strict["covert_strict"].astype(float)
    fs = register_fit(fit_logit(strict, "y", terms_for(["age_index", "htn", "migraine_any"]),
                                name="Sens B strict covert: age+HTN+migraine"), "P1", "Sens B strict covert")
    add_table("P1_sensB_strict", or_table(fs, model_label=f"Strict covert (n events={fs.events}); reduced to age, "
                                                             f"HTN, migraine because of EPV"))
    RESULTS["p1_strict_fit"] = fs

    # ---- Sensitivity C: incident-only (stroke_timing == 3), excluding timing 1/2
    inc = im[(im.p1_group == "Imaged, no infarct") | ((im.covert == 1) & (im.stroke_timing == 3))].copy()
    inc = inc[~inc.stroke_timing.isin([1, 2])]
    inc["y"] = inc["covert"].astype(float)
    n_inc = int(inc.y.sum())
    RESULTS["p1_incident_events"] = n_inc
    RESULTS["p1_covert_timing"] = im.loc[im.covert == 1, "stroke_timing"].fillna(-1).map(
        {-1: "blank", 1: "before", 2: "same", 3: "after", 9: "unknown"}).value_counts().to_dict()
    inc_rows = []
    if n_inc >= 10:
        fi = register_fit(fit_logit(inc, "y", terms_for(["age_index"]), name="Sens C incident: age"), "P1",
                          "Sens C incident-only (age-adjusted)")
        inc_rows.append(or_table(fi, model_label="Incident-only covert vs imaged-negative, age-adjusted"))
    inc_note = (f"Incident-only covert infarcts (stroke_timing == 3): {n_inc}. Covert-case timing distribution: "
                f"{RESULTS['p1_covert_timing']}. Timing is blank for most imaging-only cases because stroke_date "
                f"is blank; a multivariable incident-only model is not estimable")
    if n_inc < 10:
        inc_note += " (<10 events): only counts reported."
    log("Paper 1", inc_note)
    add_table("P1_sensC_incident", pd.concat(inc_rows) if inc_rows else pd.DataFrame([{"Note": inc_note}]), inc_note)

    # ---- Sensitivity D: MICE. Brief: run MI where the exposure or a key covariate
    # has >10% missing. With migraine code 9 kept as a level, nothing in the
    # imaged analysis sample reaches 10%, so MI is not run (logged).
    core_terms = terms_for(CORE_VARS)
    RESULTS["p1_anemia_missing_imaged_pct"] = round(100 * ana.anemia_cat.isna().mean(), 1)
    RESULTS["p1_bmi_missing_imaged_pct"] = round(100 * ana.bmi.isna().mean(), 1)
    RESULTS["p1_max_missing_pct"] = max(round(100 * ana[t.var].isna().mean(), 1) for t in FULL_TERMS)
    log("Paper 1", f"Max missingness of any pre-specified term in the imaged analysis sample: "
                   f"{RESULTS['p1_max_missing_pct']}% (anaemia {RESULTS['p1_anemia_missing_imaged_pct']}%, BMI "
                   f"{RESULTS['p1_bmi_missing_imaged_pct']}%) -> MICE not required.")

    # ------------------------------------------------------------ Covert descriptors
    cand = sheets["imaging_only_candidates"]
    cv = im[im.covert == 1][["mrn", "location", "vascular_territory", "num_infarct", "stroke_etiology", "io_subcat"]]
    cvm = cv.merge(cand.rename(columns={"location": "loc_sheet", "vascular_territory": "vt_sheet",
                                        "num_infarct": "ni_sheet", "stroke_etiology": "et_sheet"}),
                   on="mrn", how="left")
    RESULTS["p1_covert_in_candidates"] = int(cvm.evidence.notna().sum())
    RESULTS["p1_n_covert"] = len(cvm)
    vt = {1: "Anterior", 2: "Posterior", 3: "Both", 9: "Unknown"}
    ni = {1: "One", 2: "Two", 3: "Three or more", 9: "Unknown"}
    et = {1: "Cardioembolic pattern", 2: "Large-artery", 3: "Small-vessel (lacunar)", 4: "Intracranial athero",
          5: "Cryptogenic", 8: "Other determined", 9: "Unknown"}
    desc = []
    for lab, col_m, col_s, mp in [("Vascular territory", "vascular_territory", "vt_sheet", vt),
                                  ("Number of infarcts", "num_infarct", "ni_sheet", ni),
                                  ("Imaging-inferred mechanism (stroke_etiology; NOT TOAST)", "stroke_etiology",
                                   "et_sheet", et)]:
        s = cvm[col_s].map(mp).fillna("not in sheet / blank")
        sm_ = cvm[col_m].map(mp).fillna("blank")
        for lv in sorted(set(s) | set(sm_)):
            desc.append({"Descriptor": lab, "Level": lv,
                         "n (imaging_only_candidates)": int((s == lv).sum()),
                         "% (candidates)": round(100 * (s == lv).mean(), 1),
                         "n (analysis_master field)": int((sm_ == lv).sum())})
    loc = cvm["loc_sheet"].fillna("").str.lower()
    regions = {"Basal ganglia / internal capsule": r"basal ganglia|putamen|caudate|lentiform|pallid|internal capsule",
               "Thalamus": r"thalam", "Cerebellum": r"cerebell|pica|aica|sca\b",
               "Brainstem": r"pons|pontine|midbrain|medulla|brainstem|brain stem",
               "Deep white matter (corona radiata / centrum semiovale / periventricular)":
                   r"corona radiata|centrum|periventricular|white matter",
               "Cortical / lobar": r"frontal|parietal|temporal|occipital|insula|cortical|cortex",
               "Location text available": r".+"}
    for lab, p in regions.items():
        desc.append({"Descriptor": "Location (any mention; categories overlap)", "Level": lab,
                     "n (imaging_only_candidates)": int(loc.str.contains(p, regex=True).sum()),
                     "% (candidates)": round(100 * loc.str.contains(p, regex=True).mean(), 1),
                     "n (analysis_master field)": int(cvm["location"].fillna("").str.lower()
                                                      .str.contains(p, regex=True).sum())})
    for lab, p in {"Left": r"\bleft\b", "Right": r"\bright\b", "Bilateral": r"bilateral|left/right|right/left"}.items():
        desc.append({"Descriptor": "Side (any mention)", "Level": lab,
                     "n (imaging_only_candidates)": int(loc.str.contains(p, regex=True).sum()),
                     "% (candidates)": round(100 * loc.str.contains(p, regex=True).mean(), 1),
                     "n (analysis_master field)": ""})
    desc = pd.DataFrame(desc)
    add_table("P1_covert_descriptors", desc,
              "Covert infarct characteristics. Denominator = covert cases; % of covert cases.")
    RESULTS["p1_desc"] = desc
    return d
