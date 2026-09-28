"""Paper 3 - Vascular risk-factor burden across fibroids, adenomyosis and endometriosis."""
import numpy as np
import pandas as pd

from .utils import (LINEARITY, RESULTS, Term, add_table, bh, fit_logit, fmt_p, get_or, linearity_rows, log,
                    or_table, register_fit, table1, wald)

GROUPS = ["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"]
OUTCOMES = {
    "htn": "Hypertension", "dm": "Diabetes", "dyslipidemia": "Dyslipidaemia", "obesity": "Obesity (BMI ≥30)",
    "smoking_ever": "Ever smoking", "migraine_any": "Migraine (any)", "migraine_aura": "Migraine with aura",
    "afib": "Atrial fibrillation", "cad": "Coronary artery disease", "vte_history": "Prior VTE",
    "thrombophilia": "Thrombophilia",
}

# ---------------------------------------------------------------------------
# Pre-specified models (fixed before fitting)
# ---------------------------------------------------------------------------
# Exposure: uterine_dx_group, reference = fibroids only.
# M1 (minimal): exposure + age (linear) + race (White/Black/Asian/Other-unknown).
# M2 (primary): M1 + BMI (linear), except when obesity is the outcome (M2 = M1).
# Outcome definitions: smoking ever = current/former vs never (code 9 excluded);
# migraine any = with or without aura vs none (code 9 excluded); migraine with
# aura = with aura vs (none or without aura) (code 9 excluded); obesity = BMI >= 30.
# Pre-specified simplification when a model is unstable:
#   - sparse covariate level (<5 events)  -> race collapsed to 3 levels, then dropped;
#   - EPV < 10                           -> exposure + age only;
#   - sparse exposure level (<5 events)  -> that contrast reported as "not estimated".
EXPOSURE = Term("dxgrp", "cat", ref="Fibroids only", levels=GROUPS, label="Uterine diagnosis group")
AGE = Term("age_index", "cont", label="Age (per year)")
RACE = Term("race4", "cat", ref="White", levels=["White", "Black", "Asian", "Other/unknown"], label="Race")
RACE3 = Term("race3", "cat", ref="White", levels=["White", "Black", "Asian/Other/unknown"], label="Race (collapsed)")
BMI = Term("bmi", "cont", label="BMI (per kg/m²)")


def model_terms(outcome, full=True):
    t = [EXPOSURE, AGE, RACE]
    if full and outcome != "obesity":
        t.append(BMI)
    return t


def stable_fit(d, outcome, terms, label):
    f = fit_logit(d, outcome, terms, name=label)
    steps = []
    exp_cols = set(f.colmap.get("dxgrp", [])) if f.res is not None else set()

    def cov_sparse(ff):
        return [s for s in ff.sparse if not s.startswith("dxgrp=")]

    if f.res is not None and cov_sparse(f) and any(t.var == "race4" for t in terms):
        terms = [RACE3 if t.var == "race4" else t for t in terms]
        f = fit_logit(d, outcome, terms, name=label)
        steps.append("race collapsed (sparse)")
        if cov_sparse(f):
            terms = [t for t in terms if t.var != "race3"]
            f = fit_logit(d, outcome, terms, name=label)
            steps.append("race dropped (sparse)")
    if f.res is None or f.epv < 10 or not f.converged:
        terms = [EXPOSURE, AGE]
        f = fit_logit(d, outcome, terms, name=label)
        steps.append("EPV<10/non-convergence: exposure + age only")
    f.simplification = "; ".join(steps)
    return f


def contrasts(f):
    out = {}
    for g in GROUPS[1:]:
        c = f"dxgrp={g}"
        if f.res is None or c not in f.params.index:
            out[g] = None
            continue
        sparse = any(s.startswith(c + " ") for s in f.sparse)
        r = get_or(f, c)
        r["sparse"] = sparse
        out[g] = r
    return out


def std_prev(d, outcome, std_weights):
    """Direct age standardisation (5-year bands) to the whole eligible cohort."""
    rows = {}
    x = d.dropna(subset=[outcome])
    for g in GROUPS:
        s = x[x.dxgrp == g]
        est, var = 0.0, 0.0
        empty = []
        for band, w in std_weights.items():
            sb = s[s.age_band == band]
            if len(sb) == 0:
                empty.append(band)
                continue
            p = sb[outcome].mean()
            est += w * p
            var += w ** 2 * p * (1 - p) / len(sb)
        rows[g] = (est, est - 1.96 * np.sqrt(var), est + 1.96 * np.sqrt(var), empty,
                   int(s[outcome].sum()), len(s))
    return rows


def run(elig):
    d = elig.copy()
    RESULTS["p3_n"] = len(d)
    RESULTS["p3_group_n"] = d.dxgrp.value_counts().to_dict()
    RESULTS["p3_age_by_group"] = d.groupby("dxgrp").age_index.median().to_dict()

    # Table 1 by condition
    spec = [("age_index", "Age at index, years", "cont"), ("age2", "Age group", "cat"), ("race4", "Race", "cat"),
            ("bmi", "BMI, kg/m²", "cont"), ("uterine_bleeding", "Heavy/abnormal uterine bleeding", "bin"),
            ("hormonal_tx", "Hormonal therapy", "bin")] + \
           [(o, lab + (" (code 9 excluded)" if o in ("smoking_ever", "migraine_any", "migraine_aura") else ""), "bin")
            for o, lab in OUTCOMES.items()] + [("stroke_any", "Any stroke", "bin")]
    t1 = table1(d, "dxgrp", GROUPS, spec)
    add_table("P3_Table1", t1, "Characteristics by uterine diagnosis group (crude).")

    # standard population = whole eligible cohort, 5-year bands
    w = d.age_band.value_counts(normalize=True).sort_index()
    std_w = w.to_dict()
    RESULTS["p3_std_bands"] = {str(k): round(v, 4) for k, v in std_w.items()}

    prev_rows, main_rows, forest, strata_rows, inc_rows = [], [], [], [], []
    pvals = []
    for o, olab in OUTCOMES.items():
        # prevalences
        sp = std_prev(d, o, std_w)
        for g in GROUPS:
            est, lo, hi, empty, k, n = sp[g]
            prev_rows.append({"Outcome": olab, "Group": g, "n with outcome": k, "N (outcome known)": n,
                              "Crude %": round(100 * k / n, 2) if n else np.nan,
                              "Age-standardised % (95% CI)": f"{100 * est:.2f} ({100 * max(lo, 0):.2f}–{100 * hi:.2f})",
                              "Age-std %": 100 * est, "Age-std low": 100 * max(lo, 0), "Age-std high": 100 * hi,
                              "Empty age bands": ", ".join(map(str, empty))})
        # models
        f1 = register_fit(stable_fit(d, o, model_terms(o, full=False), f"{olab} M1"), "P3", f"{olab}: M1 age+race")
        f2 = register_fit(stable_fit(d, o, model_terms(o, full=True), f"{olab} M2"), "P3", f"{olab}: M2 age+race+BMI")
        c1, c2 = contrasts(f1), contrasts(f2)
        _, _, gp1 = wald(f1, f1.colmap["dxgrp"])
        _, _, gp2 = wald(f2, f2.colmap["dxgrp"])
        for g in GROUPS[1:]:
            a, b = c1[g], c2[g]
            row = {"Outcome": olab, "Group vs fibroids only": g,
                   "M1 OR (95% CI) age+race": ("not estimated (sparse)" if a and a["sparse"] else (a["txt"] if a else "")),
                   "M1 p": a["p"] if a and not a["sparse"] else np.nan,
                   "M2 OR (95% CI) age+race+BMI": ("not estimated (sparse)" if b and b["sparse"] else (b["txt"] if b else "")),
                   "OR": b["OR"] if b and not b["sparse"] else np.nan,
                   "CI low": b["lo"] if b and not b["sparse"] else np.nan,
                   "CI high": b["hi"] if b and not b["sparse"] else np.nan,
                   "M2 p (raw)": b["p"] if b and not b["sparse"] else np.nan,
                   "Global p M2 (3 df)": gp2, "N M2": f2.n, "Events M2": f2.events, "EPV M2": round(f2.epv, 1),
                   "M2 simplification": f2.simplification, "M2 flags": "; ".join(f2.flags),
                   "M2 note": "BMI not adjusted (outcome)" if o == "obesity" else ""}
            main_rows.append(row)
            if b and not b["sparse"]:
                forest.append({"Outcome": olab, "Group": g, "OR": b["OR"], "lo": b["lo"], "hi": b["hi"]})
        pvals.append((o, gp2))
        # age strata (M2)
        for ag in ["18–39", "40–60"]:
            ds = d[d.age2 == ag]
            fs = register_fit(stable_fit(ds, o, model_terms(o, True), f"{olab} M2 {ag}"), "P3", f"{olab}: M2 age {ag}")
            cs = contrasts(fs)
            for g in GROUPS[1:]:
                b = cs[g]
                n_g = int((ds.dxgrp == g).sum())
                ev_g = int(ds.loc[ds.dxgrp == g, o].sum())
                strata_rows.append({"Outcome": olab, "Age": ag, "Group vs fibroids only": g,
                                    "n in group": n_g, "events in group": ev_g,
                                    "OR (95% CI)": "not estimated (sparse)" if b and b["sparse"] else (b["txt"] if b else ""),
                                    "OR": b["OR"] if b and not b["sparse"] else np.nan,
                                    "p": b["p"] if b and not b["sparse"] else np.nan,
                                    "N": fs.n, "Events": fs.events, "Simplification": fs.simplification,
                                    "Flags": "; ".join(fs.flags)})
        # incident-only sensitivity: exclude patients whose stroke preceded / coincided with index
        di = d[~d.stroke_timing.isin([1, 2])]
        fi = register_fit(stable_fit(di, o, model_terms(o, True), f"{olab} M2 excl prior stroke"), "P3",
                          f"{olab}: M2 excluding prior stroke")
        ci_ = contrasts(fi)
        for g in GROUPS[1:]:
            b = ci_[g]
            inc_rows.append({"Outcome": olab, "Group vs fibroids only": g,
                             "Primary M2 OR": c2[g]["txt"] if c2[g] and not c2[g]["sparse"] else "",
                             "Excluding stroke before/at index": ("not estimated (sparse)" if b and b["sparse"]
                                                                  else (b["txt"] if b else "")),
                             "N": fi.n, "Events": fi.events, "Flags": "; ".join(fi.flags)})
        # linearity
        LINEARITY.extend(linearity_rows(d, o, model_terms(o, True), ["age_index"] + (["bmi"] if o != "obesity" else []),
                                        "P3", f"{olab} M2"))

    prev = pd.DataFrame(prev_rows)
    add_table("P3_prevalence", prev, "Crude and age-standardised prevalence (direct standardisation, 5-year bands, "
                                     "standard = whole eligible cohort).")
    main = pd.DataFrame(main_rows)
    main["M2 p (BH-adjusted, 33 contrasts)"] = bh(main["M2 p (raw)"].values)
    main["M1 p (BH-adjusted, 33 contrasts)"] = bh(main["M1 p"].values)
    gl = pd.DataFrame(pvals, columns=["o", "p"])
    gl["p_bh"] = bh(gl["p"].values)
    main["Global p M2 (BH-adjusted, 11 outcomes)"] = main["Outcome"].map(
        dict(zip([OUTCOMES[o] for o in gl.o], gl.p_bh)))
    for c in ["M1 p", "M2 p (raw)", "M2 p (BH-adjusted, 33 contrasts)", "M1 p (BH-adjusted, 33 contrasts)",
              "Global p M2 (3 df)", "Global p M2 (BH-adjusted, 11 outcomes)"]:
        main[c + " (text)"] = main[c].map(fmt_p)
    add_table("P3_main_models", main, "Logistic regression, HC1 robust SEs; reference = fibroids only. "
                                      "Benjamini-Hochberg correction across 33 contrasts (11 outcomes × 3 groups) "
                                      "and across the 11 global tests.")
    RESULTS["p3_main"] = main
    RESULTS["p3_prev"] = prev
    RESULTS["p3_forest"] = pd.DataFrame(forest)
    st = pd.DataFrame(strata_rows)
    add_table("P3_age_strata", st, "M2 stratified by age at index.")
    RESULTS["p3_strata"] = st
    # effect-modification by age: interaction test (M2 + dxgrp x age2)
    em = []
    for o, olab in OUTCOMES.items():
        x = d.copy()
        x["old"] = (x.age2 == "40–60").astype(float)
        ic = []
        for g in GROUPS[1:]:
            c = f"old_x_{g}"
            x[c] = x["old"] * (x.dxgrp == g)
            ic.append(c)
        terms = model_terms(o, True) + [Term("old", "bin")] + [Term(c, "bin") for c in ic]
        f = fit_logit(x, o, terms, compute_vif=False)
        chi2, dfi, p = wald(f, [c for c in ic if c in f.params.index])
        em.append({"Outcome": olab, "Wald χ² (condition × age group)": round(chi2, 2), "df": dfi, "p": p,
                   "p (text)": fmt_p(p), "Flags": "; ".join(f.flags)})
    em = pd.DataFrame(em)
    em["p BH (11 outcomes)"] = bh(em["p"].values)
    em["p BH (text)"] = em["p BH (11 outcomes)"].map(fmt_p)
    add_table("P3_age_interaction", em, "Test of condition × age-group (18–39 vs 40–60) interaction, M2.")
    RESULTS["p3_em"] = em
    add_table("P3_sens_excl_prior_stroke", pd.DataFrame(inc_rows),
              "Incident-only adaptation for Paper 3: patients with stroke before/at index (stroke_timing 1 or 2) "
              "excluded, since a prior stroke can prompt risk-factor ascertainment.")
    RESULTS["p3_inc"] = pd.DataFrame(inc_rows)
    RESULTS["p3_excl_prior_n"] = int(d.stroke_timing.isin([1, 2]).sum())

    # sensitivity: migraine code 9 counted as migraine (aura unspecified) for migraine_any
    d["migraine_any_incl9"] = d["migraine"].map({0: 0.0, 1: 1.0, 2: 1.0, 9: 1.0})
    f9 = register_fit(stable_fit(d, "migraine_any_incl9", model_terms("migraine_any", True), "mig incl 9"), "P3",
                      "Migraine any incl code 9 (sensitivity)")
    c9 = contrasts(f9)
    rows9 = [{"Group vs fibroids only": g, "Primary (code 9 excluded)":
              main.query("Outcome == 'Migraine (any)' and `Group vs fibroids only` == @g")["M2 OR (95% CI) age+race+BMI"].iloc[0],
              "Code 9 counted as migraine": c9[g]["txt"], "N": f9.n, "Events": f9.events} for g in GROUPS[1:]]
    add_table("P3_sens_migraine9", pd.DataFrame(rows9), "Sensitivity only - not a redefinition (see analysis_log).")
    RESULTS["p3_mig9"] = pd.DataFrame(rows9)
    RESULTS["p3_smoking_known"] = int(d.smoking_ever.notna().sum())
    return d
