"""Paper 2 - Iron-deficiency red-cell indices, reactive thrombocytosis and stroke."""
import numpy as np
import pandas as pd
from scipy import stats

from .utils import (LINEARITY, RESULTS, Term, add_table, design, fit_logit, fmt_or, fmt_p, get_or,
                    linearity_rows, log, mice, or_table, pool_fits, rcs_basis, register_fit, wald)

# ---------------------------------------------------------------------------
# Pre-specified terms (fixed before fitting)
# ---------------------------------------------------------------------------
# Covariates (brief, Paper 2): age (linear), race (4-level, White ref), BMI (linear),
# hypertension, diabetes, dyslipidaemia, smoking (Never/Ever/Unknown),
# migraine (any diagnosis, yes/no; PI decision), thrombophilia,
# hormonal_type (None ref; Codebook levels 1-6; code 9 -> missing),
# uterine_bleeding, uterine_dx_group (Fibroids only ref).
# Pre-specified sparse rule: a hormonal_type level with <5 outcome events is
# merged into "Other/multiple" (applied per outcome, logged).
COV = [
    Term("age_index", "cont", label="Age (per year)"),
    Term("race4", "cat", ref="White", levels=["White", "Black", "Asian", "Other/unknown"], label="Race"),
    Term("bmi", "cont", label="BMI (per kg/m²)"),
    Term("htn", "bin", label="Hypertension"),
    Term("dm", "bin", label="Diabetes"),
    Term("dyslipidemia", "bin", label="Dyslipidaemia"),
    Term("smoking3", "cat", ref="Never", levels=["Never", "Ever", "Unknown"], label="Smoking"),
    Term("migraine_any", "bin", label="Migraine (any type)"),
    Term("thrombophilia", "bin", label="Thrombophilia"),
    Term("hormonal_type_cat", "cat", ref="None", label="Hormonal therapy type"),
]
# Uterine bleeding and diagnosis group lie on the pathway from the uterine condition to anaemia, so adjusting for
# them can remove part of the effect of interest (reviewer comment, 2026-09-30). They are no longer in the primary
# model; the original model that included them is kept as a sensitivity analysis (UTERINE_TERMS).
UTERINE_TERMS = [
    Term("uterine_bleeding", "bin", label="Heavy/abnormal uterine bleeding"),
    Term("dxgrp", "cat", ref="Fibroids only",
         levels=["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"],
         label="Uterine diagnosis group"),
]
EXPO = {
    "anemia_cat": Term("anemia_cat", "cat", ref="None (Hb ≥12)",
                       levels=["None (Hb ≥12)", "Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"], label="Anaemia grade"),
    "mcv_cat": Term("mcv_cat", "cat", ref="Normal (80–100)",
                    levels=["Normal (80–100)", "Microcytic (<80)", "Macrocytic (>100)"], label="MCV"),
    "plt_cat": Term("plt_cat", "cat", ref="Normal (150–400)",
                    levels=["Normal (150–400)", "Low (<150)", "Thrombocytosis (>400)"], label="Platelets"),
    "iron_deficiency": Term("iron_deficiency", "bin", label="Iron deficiency (ferritin <30)"),
}
HT_LEVELS = ["None", "Combined OC", "Progestin-only", "LNG-IUD", "GnRH agonist/antagonist", "Menopausal HT",
             "Other/multiple"]


def cov_for(df, outcome, extra=(), tag=""):
    """Covariate list with the pre-specified sparse rule for hormonal_type,
    evaluated in the complete-case sample of the model it will be used in
    (outcome + all covariates + `extra` exposure columns)."""
    need = [outcome] + [t.var for t in COV] + list(extra)
    d = df.dropna(subset=need)
    ev = d.groupby("hormonal_type_cat")[outcome].sum()
    merge = [lv for lv in HT_LEVELS[1:-1] if ev.get(lv, 0) < 5]
    col = "hormonal_type_cat"
    if merge:
        col = f"ht_{outcome}_{abs(hash((tuple(extra), tag, len(d)))) % 10**6}"
        df[col] = df["hormonal_type_cat"].replace({m: "Other/multiple" for m in merge})
        log("Paper 2", f"[{outcome}{' ' + tag if tag else ''}{' + ' + ','.join(extra) if extra else ''}] hormonal_type "
                       f"levels with <5 events merged into Other/multiple: {merge}")
    levels = [lv for lv in HT_LEVELS if lv not in merge]
    out = []
    for t in COV:
        if t.var == "hormonal_type_cat":
            out.append(Term(col, "cat", ref="None", levels=levels, label="Hormonal therapy type"))
        else:
            out.append(t)
    return out


def _irls(X, y, iters=25):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = X @ b
        p = 1 / (1 + np.exp(-eta))
        w = p * (1 - p)
        z = eta + (y - p) / np.maximum(w, 1e-10)
        XtW = X.T * w
        b_new = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-8:
            b = b_new
            break
        b = b_new
    return b


def run(elig):
    d = elig.copy()
    RESULTS["p2_n"] = len(d)
    RESULTS["p2_events"] = int(d.stroke_any.sum())
    RESULTS["p2_isch_events"] = int(d.ischaemic.sum())
    # ischaemic-only analysis: ischaemic cases vs no stroke (other stroke types excluded)
    d["y_isch"] = np.where(d.stroke_any == 0, 0.0, np.where(d.stroke_type == 1, 1.0, np.nan))
    # incident-only: timing==3 cases vs no stroke; timing 1/2 excluded; timing 9/blank neither case nor control
    d["y_incident"] = np.where(d.stroke_any == 0, 0.0, np.where(d.stroke_timing == 3, 1.0, np.nan))
    RESULTS["p2_incident_events"] = int((d.y_incident == 1).sum())
    RESULTS["p2_incident_excluded_prior"] = int(d.stroke_timing.isin([1, 2]).sum())
    RESULTS["p2_incident_excluded_unknown"] = int(((d.stroke_any == 1) & ~d.stroke_timing.isin([1, 2, 3])).sum())

    outcomes = {"stroke_any": "Any stroke (primary)", "y_isch": "Ischaemic stroke vs no stroke",
                "y_incident": "Incident stroke (after index) vs no stroke"}
    covs = {o: cov_for(d, o) for o in outcomes}

    def C(o, extra=(), data=None, tag=""):
        return cov_for(d if data is None else data, o, extra, tag)

    # Pre-specified fallback when an adjusted model has EPV < 10: reduced
    # covariate set (age, race, BMI, HTN, DM, dyslipidaemia, smoking).
    REDUCED = ["age_index", "race4", "bmi", "htn", "dm", "dyslipidemia", "smoking3"]

    # -------------------------------------------------- descriptive (crude prevalence by exposure)
    desc = []
    for ev, t in EXPO.items():
        for lv in (t.levels or [0.0, 1.0]):
            s = d[d[ev] == lv]
            desc.append({"Exposure": t.label, "Level": str(lv if t.kind == "cat" else ("yes" if lv == 1 else "no")),
                         "n": len(s), "Any stroke n": int(s.stroke_any.sum()),
                         "Any stroke %": round(100 * s.stroke_any.mean(), 2) if len(s) else np.nan,
                         "Ischaemic n": int(s.ischaemic.sum())})
        desc.append({"Exposure": t.label, "Level": "missing", "n": int(d[ev].isna().sum()),
                     "Any stroke n": int(d.loc[d[ev].isna(), "stroke_any"].sum()),
                     "Any stroke %": round(100 * d.loc[d[ev].isna(), "stroke_any"].mean(), 2)})
    add_table("P2_Table1_exposures", pd.DataFrame(desc), "Stroke prevalence by exposure level (crude).")

    # Table 1 by anaemia grade
    from .utils import table1
    spec = [("age_index", "Age at index, years", "cont"), ("race4", "Race", "cat"), ("bmi", "BMI, kg/m²", "cont"),
            ("htn", "Hypertension", "bin"), ("dm", "Diabetes", "bin"), ("dyslipidemia", "Dyslipidaemia", "bin"),
            ("smoking3", "Smoking", "cat"), ("migraine_any", "Migraine (any type)", "bin"),
            ("thrombophilia", "Thrombophilia", "bin"), ("hormonal_type_cat", "Hormonal therapy type", "cat"),
            ("uterine_bleeding", "Heavy/abnormal uterine bleeding", "bin"), ("dxgrp", "Uterine diagnosis group", "cat"),
            ("hgb", "Haemoglobin, g/dL", "cont"), ("mcv", "MCV, fL", "cont"), ("platelets", "Platelets, ×10³/µL", "cont"),
            ("ferritin", "Ferritin, ng/mL", "cont"), ("iron_deficiency", "Iron deficiency (ferritin <30)", "bin"),
            ("stroke_any", "Any stroke", "bin"), ("ischaemic", "Ischaemic stroke", "bin")]
    t1 = table1(d, "anemia_cat", EXPO["anemia_cat"].levels, spec)
    add_table("P2_Table1", t1, "Characteristics by anaemia grade (patients with Hb within ±3 years of index). "
                               f"{int(d.anemia_cat.isna().sum()):,} patients without Hb are not in the grade columns.")

    # -------------------------------------------------- main models
    main_rows, forest = [], []
    fits = {}
    for o, olab in outcomes.items():
        for ev in ["anemia_cat", "mcv_cat", "plt_cat"]:
            et = EXPO[ev]
            fc = register_fit(fit_logit(d, o, [et], name=f"{olab}: {ev} crude"), "P2", f"{olab}: {ev} crude")
            fa = register_fit(fit_logit(d, o, [et] + C(o, [ev]), name=f"{olab}: {ev} adjusted"), "P2",
                              f"{olab}: {ev} adjusted")
            fits[(o, ev)] = fa
            tc = or_table(fc, only_vars=[ev])
            ta = or_table(fa, only_vars=[ev])
            m = ta.rename(columns={"OR (95% CI)": "Adjusted OR (95% CI)", "p (text)": "Adjusted p"})
            m["Crude OR (95% CI)"] = tc["OR (95% CI)"].values
            m["Crude p"] = tc["p (text)"].values
            m["Outcome"] = olab
            m["Model"] = "separate"
            main_rows.append(m)
            for _, r in ta.dropna(subset=["OR"]).iterrows():
                forest.append({"Outcome": olab, "Exposure": et.label, "Level": r["Level"], "OR": r["OR"],
                               "lo": r["CI low"], "hi": r["CI high"], "model": "separate"})
        # joint model
        fj = register_fit(fit_logit(d, o, [EXPO["anemia_cat"], EXPO["mcv_cat"], EXPO["plt_cat"]] +
                                    C(o, ["anemia_cat", "mcv_cat", "plt_cat"]),
                                    name=f"{olab}: joint"), "P2", f"{olab}: joint (anaemia+MCV+platelets)")
        fits[(o, "joint")] = fj
        tj = or_table(fj, only_vars=["anemia_cat", "mcv_cat", "plt_cat"])
        tj = tj.rename(columns={"OR (95% CI)": "Adjusted OR (95% CI)", "p (text)": "Adjusted p"})
        tj["Outcome"] = olab
        tj["Model"] = "joint"
        main_rows.append(tj)
        # p-trend across anaemia grades (0-3 as a score), adjusted
        ft = fit_logit(d, o, [Term("anemia", "cont", label="Anaemia grade (score 0-3)")] + C(o, ["anemia"]),
                       name=f"{olab}: trend")
        g = get_or(ft, "anemia")
        RESULTS[f"p2_trend_{o}"] = g
        main_rows.append(pd.DataFrame([{"Outcome": olab, "Model": "trend", "Term": "Anaemia grade (score 0–3)",
                                        "Level": "per grade (p-trend)", "Adjusted OR (95% CI)": g["txt"],
                                        "OR": g["OR"], "CI low": g["lo"], "CI high": g["hi"], "p": g["p"],
                                        "Adjusted p": fmt_p(g["p"]), "N": ft.n, "Events": ft.events}]))
    main = pd.concat(main_rows, ignore_index=True)
    first = ["Outcome", "Model", "Term", "Level", "Crude OR (95% CI)", "Crude p", "Adjusted OR (95% CI)",
             "Adjusted p", "OR", "CI low", "CI high", "p", "n in level", "events in level", "N", "Events", "EPV", "Flags"]
    main = main[[c for c in first if c in main.columns]]
    add_table("P2_main_models", main, "Logistic regression, HC1 robust SEs. Separate = one exposure + covariates; "
                                      "joint = anaemia + MCV + platelets + covariates.")
    RESULTS["p2_fits"] = fits
    RESULTS["p2_forest"] = pd.DataFrame(forest)

    # full covariate table for the primary anaemia model (supplement)
    add_table("P2_anemia_full_model", or_table(fits[("stroke_any", "anemia_cat")],
                                               model_label="Any stroke ~ anaemia + covariates"))
    LINEARITY.extend(linearity_rows(d, "stroke_any", [EXPO["anemia_cat"]] + C("stroke_any", ["anemia_cat"]),
                                    ["age_index", "bmi"], "P2", "Any stroke ~ anaemia + covariates"))

    # -------------------------------------------------- iron deficiency (secondary, complete-case subset)
    idr = []
    for o, olab in outcomes.items():
        fc = register_fit(fit_logit(d, o, [EXPO["iron_deficiency"]], name="ID crude"), "P2", f"{olab}: iron def crude")
        fa = register_fit(fit_logit(d, o, [EXPO["iron_deficiency"]] + C(o, ["iron_deficiency"]), name="ID adj"), "P2",
                          f"{olab}: iron def adjusted")
        fa2 = register_fit(fit_logit(d, o, [EXPO["iron_deficiency"], EXPO["anemia_cat"]] +
                                     C(o, ["iron_deficiency", "anemia_cat"]), name="ID+anaemia"),
                           "P2", f"{olab}: iron def + anaemia adjusted")
        mods = [("crude", fc), ("adjusted", fa), ("adjusted + anaemia grade", fa2)]
        for lab, f in list(mods[1:]):
            if f.epv < 10:
                red = [t for t in COV if t.var in REDUCED]
                ex = [EXPO["iron_deficiency"]] + ([EXPO["anemia_cat"]] if "anaemia" in lab else [])
                fr = register_fit(fit_logit(d, o, ex + red, name="ID reduced"), "P2", f"{olab}: iron def {lab} REDUCED")
                log("Paper 2", f"{olab}, iron deficiency {lab}: EPV {f.epv:.1f} < 10 -> reduced covariate set; "
                               f"EPV now {fr.epv:.1f}; flags {fr.flags or 'none'}")
                mods[[m[0] for m in mods].index(lab)] = (lab + " [UNSTABLE, EPV<10 - not for reporting]", f)
                mods.append((lab + " (reduced covariate set, EPV fallback)", fr))
        for lab, f in mods:
            g = get_or(f, "iron_deficiency")
            idr.append({"Outcome": olab, "Model": lab, "OR (95% CI)": g["txt"], "OR": g["OR"], "CI low": g["lo"],
                        "CI high": g["hi"], "p": g["p"], "p (text)": fmt_p(g["p"]), "N": f.n, "Events": f.events,
                        "EPV": round(f.epv, 1), "Flags": "; ".join(f.flags)})
    idr = pd.DataFrame(idr)
    add_table("P2_iron_deficiency", idr, "Secondary: ferritin available for "
              f"{int(d.iron_deficiency.notna().sum()):,}/{len(d):,} "
              f"({100 * d.iron_deficiency.notna().mean():.1f}%); complete-case, selected subset.")
    RESULTS["p2_id"] = idr
    RESULTS["p2_ferritin_n"] = int(d.iron_deficiency.notna().sum())
    RESULTS["p2_ferritin_pct"] = 100 * d.iron_deficiency.notna().mean()
    # who has ferritin measured? (selection)
    sel = d.groupby(d.ferritin.notna())[["stroke_any", "uterine_bleeding", "anemia_any"]].mean().T
    sel.columns = ["ferritin not measured", "ferritin measured"]
    sel = (100 * sel).round(1).reset_index().rename(columns={"index": "% with"})
    add_table("P2_ferritin_selection", sel, "Characteristics of patients with vs without a ferritin result.")
    RESULTS["p2_ferr_sel"] = sel

    # -------------------------------------------------- Hb spline
    hs = d.dropna(subset=["hgb"])
    ht = Term("hgb", "rcs", label="Haemoglobin (RCS, 4 knots)")
    fsp = register_fit(fit_logit(hs, "stroke_any", [ht] + C("stroke_any", ["hgb"]), name="Hb spline"), "P2",
                       "Any stroke ~ RCS(Hb) + covariates")
    cols = fsp.colmap["hgb"]
    chi2, dfn, pnl = wald(fsp, cols[1:])
    chi2o, dfo, pov = wald(fsp, cols)
    RESULTS["p2_spline_p_nonlin"] = pnl
    RESULTS["p2_spline_p_overall"] = pov
    RESULTS["p2_spline_knots"] = list(ht.knots)
    lo_, hi_ = np.nanpercentile(hs.hgb, [1, 99])
    grid = np.linspace(lo_, hi_, 200)
    B = rcs_basis(grid, ht.knots) - rcs_basis(np.array([13.0]), ht.knots)
    beta = fsp.params[cols].values
    V = fsp.cov.loc[cols, cols].values
    lp = B @ beta
    se = np.sqrt(np.einsum("ij,jk,ik->i", B, V, B))
    curve = pd.DataFrame({"Hb (g/dL)": grid, "OR": np.exp(lp), "CI low": np.exp(lp - 1.96 * se),
                          "CI high": np.exp(lp + 1.96 * se)})
    RESULTS["p2_curve"] = curve
    RESULTS["p2_spline_n"], RESULTS["p2_spline_events"] = fsp.n, fsp.events
    pts = []
    for h in [7, 8, 9, 10, 11, 12, 14, 15, 16]:
        b = (rcs_basis(np.array([float(h)]), ht.knots) - rcs_basis(np.array([13.0]), ht.knots))[0]
        l, s_ = float(b @ beta), float(np.sqrt(b @ V @ b))
        pts.append({"Hb (g/dL)": h, "OR vs 13 g/dL (95% CI)": fmt_or(np.exp(l), np.exp(l - 1.96 * s_), np.exp(l + 1.96 * s_)),
                    "OR": np.exp(l), "CI low": np.exp(l - 1.96 * s_), "CI high": np.exp(l + 1.96 * s_)})
    pts = pd.DataFrame(pts)
    RESULTS["p2_spline_pts"] = pts
    pts.insert(0, "Model", f"RCS Hb, knots {', '.join(f'{k:.1f}' for k in ht.knots)}; n={fsp.n}, events={fsp.events}; "
                           f"overall p={fmt_p(pov)}, non-linearity p={fmt_p(pnl)}")
    add_table("P2_Hb_spline", pts, "Adjusted OR for any stroke at selected Hb values vs 13 g/dL.")

    # -------------------------------------------------- interactions
    ints = []
    # microcytosis x thrombocytosis
    for o, olab in [("stroke_any", outcomes["stroke_any"]), ("y_isch", outcomes["y_isch"])]:
        x = d.dropna(subset=["microcytosis", "thrombocytosis"]).copy()
        x["micro_x_thromb"] = x["microcytosis"] * x["thrombocytosis"]
        terms = [Term("microcytosis", "bin", label="Microcytosis (MCV<80)"),
                 Term("thrombocytosis", "bin", label="Thrombocytosis (>400)"),
                 Term("micro_x_thromb", "bin", label="Microcytosis × thrombocytosis")] + \
            C(o, ["microcytosis", "thrombocytosis"])
        fi = register_fit(fit_logit(x, o, terms, name="micro x thromb"), "P2", f"{olab}: microcytosis × thrombocytosis")
        g = get_or(fi, "micro_x_thromb")
        cells = []
        for mi in [0, 1]:
            for th in [0, 1]:
                s = x[(x.microcytosis == mi) & (x.thrombocytosis == th)].dropna(subset=[o])
                cells.append(f"micro={mi}/thromb={th}: {int(s[o].sum())}/{len(s)}")
        # joint-category ORs vs neither (from the interaction model)
        comb = {}
        for lab_, v in [("microcytosis only", ["microcytosis"]), ("thrombocytosis only", ["thrombocytosis"]),
                        ("both", ["microcytosis", "thrombocytosis", "micro_x_thromb"])]:
            L = np.zeros(len(fi.params))
            for c in v:
                L[list(fi.params.index).index(c)] = 1
            est = float(L @ fi.params.values)
            se_ = float(np.sqrt(L @ fi.cov.values @ L))
            comb[lab_] = fmt_or(np.exp(est), np.exp(est - 1.96 * se_), np.exp(est + 1.96 * se_))
        ints.append({"Outcome": olab, "Interaction": "microcytosis × thrombocytosis",
                     "Interaction OR (ratio of ORs)": g["txt"], "p interaction": g["p"], "p (text)": fmt_p(g["p"]),
                     "Cells (events/n)": "; ".join(cells),
                     "OR microcytosis only vs neither": comb["microcytosis only"],
                     "OR thrombocytosis only vs neither": comb["thrombocytosis only"],
                     "OR both vs neither": comb["both"], "N": fi.n, "Events": fi.events, "Flags": "; ".join(fi.flags)})
    # anaemia x uterine bleeding
    for o, olab in [("stroke_any", outcomes["stroke_any"]), ("y_isch", outcomes["y_isch"])]:
        x = d.dropna(subset=["anemia_cat", "uterine_bleeding"]).copy()
        ic = []
        for lv in EXPO["anemia_cat"].levels[1:]:
            c = f"ax_{lv}"
            x[c] = (x.anemia_cat == lv).astype(float) * x.uterine_bleeding
            ic.append(c)
        terms = [EXPO["anemia_cat"]] + C(o, ["anemia_cat"]) + [Term(c, "bin", label=c) for c in ic]
        fi = register_fit(fit_logit(x, o, terms, name="anaemia x bleeding"), "P2", f"{olab}: anaemia × bleeding")
        chi2, dfi, pint = wald(fi, [c for c in ic if c in fi.params.index])
        strata = []
        for bl in [0, 1]:
            xs = x[x.uterine_bleeding == bl].copy()
            cv = [t for t in C(o, ["anemia_cat"], data=xs, tag=f"bleeding={bl}") if t.var != "uterine_bleeding"]
            fs = register_fit(fit_logit(xs, o, [EXPO["anemia_cat"]] + cv, name=f"anaemia | bleeding={bl}"), "P2",
                              f"{olab}: anaemia | bleeding={bl}")
            parts = [f"{c.split('=')[1]}: {get_or(fs, c)['txt']}" for c in fs.colmap["anemia_cat"]]
            strata.append(f"bleeding={'yes' if bl else 'no'} (n={fs.n}, events={fs.events}"
                          f"{'; ' + '; '.join(fs.flags) if fs.flags else ''}): " + "; ".join(parts))
        ints.append({"Outcome": olab, "Interaction": "anaemia grade × uterine bleeding",
                     "Interaction OR (ratio of ORs)": f"Wald χ²({dfi})={chi2:.2f}", "p interaction": pint,
                     "p (text)": fmt_p(pint), "Cells (events/n)": "",
                     "OR microcytosis only vs neither": "", "OR thrombocytosis only vs neither": "",
                     "OR both vs neither": "", "Stratum-specific anaemia ORs": " || ".join(strata),
                     "N": fi.n, "Events": fi.events, "Flags": "; ".join(fi.flags)})
    ints = pd.DataFrame(ints)
    add_table("P2_interactions", ints, "Product terms added to the adjusted model; robust Wald tests.")
    RESULTS["p2_ints"] = ints

    # -------------------------------------------------- lab timing: Hb within ±1 year
    lt = []
    near = d[d.hgb_days.abs() <= 365].copy()
    RESULTS["p2_hb1y_n"] = len(near)
    for o, olab in outcomes.items():
        f1 = register_fit(fit_logit(near, o, [EXPO["anemia_cat"]] + C(o, ["anemia_cat"], data=near, tag="Hb±1y"),
                                    name="Hb ±1y"), "P2",
                          f"{olab}: anaemia, Hb within ±1 y")
        f3 = fits[(o, "anemia_cat")]
        for c in f1.colmap["anemia_cat"]:
            a, b = get_or(f3, c), get_or(f1, c)
            lt.append({"Outcome": olab, "Level": c.split("=")[1], "Hb within ±3 y (primary)": a["txt"],
                       "Hb within ±1 y": b["txt"], "OR ±1y": b["OR"], "CI low ±1y": b["lo"], "CI high ±1y": b["hi"],
                       "p ±1y": b["p"], "N ±1y": f1.n, "Events ±1y": f1.events, "Flags ±1y": "; ".join(f1.flags)})
        ft1 = fit_logit(near, o, [Term("anemia", "cont")] + C(o, ["anemia"], data=near, tag="Hb±1y"))
        g = get_or(ft1, "anemia")
        lt.append({"Outcome": olab, "Level": "p-trend (per grade)", "Hb within ±3 y (primary)": RESULTS[f"p2_trend_{o}"]["txt"],
                   "Hb within ±1 y": g["txt"], "OR ±1y": g["OR"], "p ±1y": g["p"], "N ±1y": ft1.n, "Events ±1y": ft1.events})
    # iron deficiency within ±1 y
    nf = d[d.ferritin_days.abs() <= 365].copy()
    fid = register_fit(fit_logit(nf, "stroke_any", [EXPO["iron_deficiency"]] +
                                 C("stroke_any", ["iron_deficiency"], data=nf, tag="ferritin±1y"), name="ID ±1y"),
                       "P2", "Any stroke: iron def, ferritin within ±1 y")
    g = get_or(fid, "iron_deficiency")
    lt.append({"Outcome": outcomes["stroke_any"], "Level": "Iron deficiency (ferritin within ±1 y)",
               "Hb within ±3 y (primary)": RESULTS["p2_id"].query("Model == 'adjusted'").iloc[0]["OR (95% CI)"] + " (ferritin ±3 y)",
               "Hb within ±1 y": g["txt"], "OR ±1y": g["OR"], "p ±1y": g["p"], "N ±1y": fid.n, "Events ±1y": fid.events})
    lt = pd.DataFrame(lt)
    add_table("P2_sens_lab_timing", lt, "Restricting to labs within ±1 year of index (hgb_days / ferritin_days). "
                                        "MCV/platelet dates are not available (lab_date unreliable; see data_flags).")
    RESULTS["p2_lab_timing"] = lt

    # -------------------------------------------------- MICE sensitivity (BMI, migraine) for anaemia & joint models
    base = d.dropna(subset=["anemia_cat"]).copy()
    for c in ["mcv_cat", "plt_cat"]:
        base[c + "_aux"] = base[c].fillna("not measured")
    preds = ["stroke_any", "stroke_type_aux", "age_index", "race4", "htn", "dm", "dyslipidemia", "smoking3",
             "thrombophilia", "hormonal_type_aux", "uterine_bleeding", "dxgrp", "anemia_cat", "mcv_cat_aux",
             "plt_cat_aux", "afib", "cad", "vte_history"]
    base["stroke_type_aux"] = base["stroke_type"].fillna(0).astype(int).astype(str)
    base["hormonal_type_aux"] = base["hormonal_type_cat"].fillna("None")
    imps = mice(base, {"bmi": "pmm"}, preds + ["migraine_any"], m=20, iters=10, seed=202)
    RESULTS["p2_mi_n"] = len(base)
    mi_rows = []
    for o in ["stroke_any", "y_isch"]:
        cv = C(o, ["anemia_cat", "mcv_cat", "plt_cat"])
        for x in imps:
            for c in cv:
                if c.var not in x:
                    x[c.var] = d.loc[x.index, c.var]
        fl = [fit_logit(x, o, [EXPO["anemia_cat"]] + cv, compute_vif=False) for x in imps]
        cols = fl[0].colmap["anemia_cat"]
        pr = pool_fits(fl, cols)
        cc = fits[(o, "anemia_cat")]
        pr.insert(0, "Outcome", outcomes[o])
        pr.insert(2, "Complete-case OR (95% CI)", [get_or(cc, c)["txt"] for c in cols])
        pr["n (MI)"] = int(imps[0][o].notna().sum())
        pr["n (complete case)"] = cc.n
        mi_rows.append(pr)
        fj = [fit_logit(x, o, [EXPO["anemia_cat"], EXPO["mcv_cat"], EXPO["plt_cat"]] + cv, compute_vif=False)
              for x in imps]
        cols = [c for c in fj[0].params.index if c.startswith(("anemia_cat", "mcv_cat", "plt_cat"))]
        prj = pool_fits(fj, cols)
        prj.insert(0, "Outcome", outcomes[o] + " (joint model)")
        cj = fits[(o, "joint")]
        prj.insert(2, "Complete-case OR (95% CI)", [get_or(cj, c)["txt"] for c in cols])
        prj["n (MI)"] = int(fj[0].n)
        prj["n (complete case)"] = cj.n
        mi_rows.append(prj)
    mi = pd.concat(mi_rows, ignore_index=True)
    add_table("P2_sens_MICE", mi, "MICE m=20: BMI (PMM) imputed among patients with "
                                  "Hb; outcome and stroke type in the imputation model; exposures not imputed.")
    RESULTS["p2_mi"] = mi

    # -------------------------------------------------- Mediation (exploratory, difference method)
    # Exposure: fibroids (1 vs 0). Mediator: anaemia grade. Covariates: as above
    # EXCEPT uterine_bleeding (a potential intermediate between fibroids and
    # anaemia) and uterine_dx_group (derived from the fibroid flag); adenomyosis
    # and endometriosis flags are adjusted for instead.
    med_cov = [t for t in C("stroke_any", ["anemia_cat", "fibroids"]) if t.var not in ("uterine_bleeding", "dxgrp")] + \
              [Term("adenomyosis", "bin", label="Adenomyosis"), Term("endometriosis", "bin", label="Endometriosis")]
    fib = Term("fibroids", "bin", label="Fibroids")
    md = d.dropna(subset=["stroke_any", "fibroids", "anemia_cat"] + [t.var for t in med_cov]).copy()
    ftot = register_fit(fit_logit(md, "stroke_any", [fib] + med_cov, name="med total"), "P2", "Mediation: total")
    fdir = register_fit(fit_logit(md, "stroke_any", [fib, EXPO["anemia_cat"]] + med_cov, name="med direct"), "P2",
                        "Mediation: direct (+anaemia)")
    bt, bd = ftot.params["fibroids"], fdir.params["fibroids"]
    pm = (bt - bd) / bt
    Xt, _ = design(md, [fib] + med_cov)
    Xd, _ = design(md, [fib, EXPO["anemia_cat"]] + med_cov)
    Xt, Xd = Xt[[c for c in ftot.params.index]].values, Xd[[c for c in fdir.params.index]].values
    y = md["stroke_any"].values.astype(float)
    it = list(ftot.params.index).index("fibroids")
    idd = list(fdir.params.index).index("fibroids")
    rng = np.random.default_rng(303)
    boots = []
    for _ in range(1000):
        ix = rng.integers(0, len(y), len(y))
        b1 = _irls(Xt[ix], y[ix])[it]
        b2 = _irls(Xd[ix], y[ix])[idd]
        boots.append(((b1 - b2) / b1, b1, b2))
    boots = np.array(boots)
    ci = np.nanpercentile(boots[:, 0], [2.5, 97.5])
    ci_diff = np.nanpercentile(boots[:, 1] - boots[:, 2], [2.5, 97.5])
    gt, gd = get_or(ftot, "fibroids"), get_or(fdir, "fibroids")
    med = pd.DataFrame([{
        "Analysis": "EXPLORATORY, cross-sectional difference method (log-odds scale)",
        "N": ftot.n, "Events": ftot.events,
        "Fibroids OR, without anaemia (total)": gt["txt"], "p total": fmt_p(gt["p"]),
        "Fibroids OR, + anaemia grade (direct)": gd["txt"], "p direct": fmt_p(gd["p"]),
        "Difference in log-OR (95% bootstrap CI)": f"{bt - bd:.4f} ({ci_diff[0]:.4f} to {ci_diff[1]:.4f})",
        "Proportion mediated": round(pm, 4),
        "Proportion mediated 95% bootstrap CI (1,000 reps, percentile)": f"{ci[0]:.3f} to {ci[1]:.3f}",
        "Bootstrap reps with total log-OR ≤ 0": int((boots[:, 1] <= 0).sum()),
        "Note": "Proportion mediated is unstable when the total association is near null; the difference method "
                "on the logit scale is affected by non-collapsibility. No temporal ordering between Hb, fibroids "
                "and stroke can be established in these data."}])
    add_table("P2_mediation", med)
    RESULTS["p2_med"] = dict(pm=pm, ci=ci, gt=gt, gd=gd, n=ftot.n, ev=ftot.events, diff=bt - bd, ci_diff=ci_diff,
                             nonpos=int((boots[:, 1] <= 0).sum()))
    log("Paper 2", f"Mediation: total {gt['txt']}, direct {gd['txt']}, PM {pm:.3f} (95% CI {ci[0]:.3f} to {ci[1]:.3f})")
    return d
