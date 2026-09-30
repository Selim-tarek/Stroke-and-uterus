"""Paper 4 - Incident stroke by uterine diagnosis group (fibroids only, adenomyosis only, endometriosis only,
>1 condition), and the role of migraine and cardiometabolic risk factors.

Cohort: eligible women without a stroke before or at the index date and without an undatable stroke; follow-up
from the index date to the first of stroke (stroke_timing = 3, dated by stroke_date), death or last encounter
(Encounters extracts). Poisson regression with log person-time offset and HC1 robust SEs:
  M1: group + age + race
  M2: M1 + BMI, hypertension, diabetes, dyslipidaemia, smoking, atrial fibrillation, coronary artery disease
  M3: M2 + migraine (migraine as a possible pathway)
  M4: M3 + hormonal therapy + heavy/abnormal bleeding
Migraine contribution: % change in log RR from M2 to M3 with a percentile bootstrap CI (difference method;
exploratory, cross-sectional ascertainment of migraine). Effect modification by age group (18-39 / 40-60).
Stroke-type descriptives (Paper 4b) and completeness of stroke detail fields are also tabulated.
Sparse rule: exposure levels with < 5 events are reported as not estimated; EPV < 10 -> M1 covariates only.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from .surgery import DATE_COLS
from .utils import RESULTS, Term, add_table, design, fmt_or, fmt_p, log

GROUPS = ["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"]
EXPO = Term("dxgrp", "cat", ref="Fibroids only", levels=GROUPS)
M1 = [Term("age_index", "cont"), Term("race4", "cat", ref="White", levels=["White", "Black", "Asian", "Other/unknown"])]
M2 = M1 + [Term("bmi", "cont"), Term("htn", "bin"), Term("dm", "bin"), Term("dyslipidemia", "bin"),
           Term("smoking3", "cat", ref="Never", levels=["Never", "Ever", "Unknown"]), Term("afib", "bin"),
           Term("cad", "bin")]
M3 = M2 + [Term("migraine_any", "bin")]
M4 = M3 + [Term("hormonal_tx", "bin"), Term("uterine_bleeding", "bin")]
MODELS = [("M1 age + race", M1), ("M2 + cardiometabolic", M2), ("M3 + migraine", M3),
          ("M4 + hormonal therapy, bleeding", M4)]
STYPE = {1: "Ischemic stroke", 2: "Intracerebral hemorrhage", 3: "Subarachnoid hemorrhage", 4: "TIA",
         5: "Cerebral venous thrombosis", 9: "Unknown type"}


def cohort(elig):
    from .followup import encounters
    end = max(pd.to_datetime(elig[c], errors="coerce").max() for c in DATE_COLS)
    f = elig.copy()
    n0 = len(f)
    n_prior = int(f.stroke_timing.isin([1, 2]).sum())
    n_unk = int(((f.stroke_any == 1) & ~f.stroke_timing.isin([1, 2, 3])).sum())
    f = f[~f.stroke_timing.isin([1, 2])]
    f = f[~((f.stroke_any == 1) & ~f.stroke_timing.isin([3]))]
    f["t0"] = pd.to_datetime(f.index_date)
    sd = pd.to_datetime(f.stroke_date, errors="coerce")
    f["ev_any"] = ((f.stroke_any == 1) & (f.stroke_timing == 3)).astype(int)
    f["ev_isch"] = (f.ev_any.astype(bool) & (f.stroke_type == 1)).astype(int)
    f["ev_stroke_no_tia"] = (f.ev_any.astype(bool) & (f.stroke_type != 4)).astype(int)
    f["t1"] = pd.to_datetime(np.where(f.ev_any == 1, sd, end))
    lc, dth = encounters()
    lastc, dd = f.mrn.map(lc), f.mrn.map(dth)
    m = (f.ev_any == 0) & lastc.notna()
    f.loc[m, "t1"] = np.minimum(f.loc[m, "t1"], lastc[m])
    m = (f.ev_any == 0) & dd.notna() & (dd < f.t1)
    f.loc[m, "t1"] = dd[m]
    n_zero = int((f.t1 <= f.t0).sum())
    f = f[f.t1 > f.t0].copy()
    f["py"] = (f.t1 - f.t0).dt.days / 365.25
    info = dict(n_eligible=n0, n_prior=n_prior, n_unknown=n_unk, n_zero=n_zero, n=len(f), py=float(f.py.sum()),
                fu_median=float(f.py.median()), ev=int(f.ev_any.sum()), ev_isch=int(f.ev_isch.sum()),
                ev_no_tia=int(f.ev_stroke_no_tia.sum()), end=end.date().isoformat())
    return f, info


def fit(x, oc, terms):
    x = x.dropna(subset=[EXPO.var] + [t.var for t in terms])
    X, cmap = design(x, [EXPO] + terms)
    X = X.loc[:, (X != 0).any(axis=0)]
    r = sm.GLM(x[oc].values, X, family=sm.families.Poisson(), offset=np.log(x.py.values)).fit(cov_type="HC1")
    return r, len(x), int(x[oc].sum()), x[oc].sum() / (X.shape[1] - 1)


def rr(r, c):
    b, se = r.params[c], r.bse[c]
    return dict(RR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=2 * stats.norm.sf(abs(b / se)),
                txt=fmt_or(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)), b=b)


def run(elig):
    f, info = cohort(elig)
    RESULTS["p4_cohort"] = info
    log("Paper 4", f"Incidence cohort {info['n']:,} women, {info['py']:,.0f} person-years, {info['ev']} strokes/TIA "
                   f"({info['ev_isch']} ischemic); excluded {info['n_prior']} with stroke before/at index, "
                   f"{info['n_unknown']} undatable, {info['n_zero']} with no follow-up time.")
    # rates by group
    rates = []
    for oc, lab in [("ev_any", "Stroke or TIA"), ("ev_isch", "Ischemic stroke"), ("ev_stroke_no_tia", "Stroke excluding TIA")]:
        for g in GROUPS:
            s = f[f.dxgrp == g]
            k, py = int(s[oc].sum()), float(s.py.sum())
            lo = stats.chi2.ppf(0.025, 2 * k) / 2 / py if k else 0.0
            hi = stats.chi2.ppf(0.975, 2 * k + 2) / 2 / py
            rates.append({"Outcome": lab, "Group": g, "Women": len(s), "Events": k, "Person-years": round(py, 1),
                          "Median age": float(s.age_index.median()),
                          "Rate /1,000 PY (95% CI)": f"{1000 * k / py:.2f} ({1000 * lo:.2f}–{1000 * hi:.2f})",
                          "Rate": 1000 * k / py, "Rate low": 1000 * lo, "Rate high": 1000 * hi})
    rates = pd.DataFrame(rates)
    add_table("P4_rates", rates, "Crude incident stroke rates by uterine diagnosis group (exact Poisson 95% CI).")
    RESULTS["p4_rates"] = rates
    # models
    rows = []
    for oc, lab in [("ev_any", "Stroke or TIA"), ("ev_isch", "Ischemic stroke"), ("ev_stroke_no_tia", "Stroke excluding TIA")]:
        for mname, terms in MODELS:
            r, n, e, epv = fit(f, oc, terms)
            note = ""
            if epv < 10:
                r, n, e, epv = fit(f, oc, M1)
                note = "EPV<10: age and race only"
            for g in GROUPS[1:]:
                ev_g = int(f.loc[(f.dxgrp == g), oc].sum())
                row = {"Outcome": lab, "Model": mname, "Group vs fibroids only": g, "Events in group": ev_g,
                       "N": n, "Events": e, "EPV": round(float(epv), 1), "Note": note}
                if ev_g < 5:
                    row["RR (95% CI)"] = "not estimated (<5 events)"
                else:
                    q = rr(r, f"dxgrp={g}")
                    row.update({"RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"],
                                "p": q["p"], "p (text)": fmt_p(q["p"])})
                rows.append(row)
    mod = pd.DataFrame(rows)
    add_table("P4_models", mod, "Incident stroke rate ratios vs fibroids only (Poisson, log person-time offset, HC1 SEs).")
    RESULTS["p4_models"] = mod
    # migraine contribution (M2 -> M3), bootstrap
    rng = np.random.default_rng(20260930)
    x = f.dropna(subset=[EXPO.var] + [t.var for t in M3]).reset_index(drop=True)
    r2, _, _, _ = fit(x, "ev_any", M2)
    r3, _, _, _ = fit(x, "ev_any", M3)
    med = []
    boots = {g: [] for g in GROUPS[1:]}
    for _ in range(300):
        b = x.iloc[rng.integers(0, len(x), len(x))]
        try:
            a2, _, _, _ = fit(b, "ev_any", M2)
            a3, _, _, _ = fit(b, "ev_any", M3)
        except Exception:
            continue
        for g in GROUPS[1:]:
            c = f"dxgrp={g}"
            if c in a2.params and c in a3.params and abs(a2.params[c]) > 1e-9:
                boots[g].append(100 * (a2.params[c] - a3.params[c]) / a2.params[c])
    for g in GROUPS[1:]:
        c = f"dxgrp={g}"
        pc = 100 * (r2.params[c] - r3.params[c]) / r2.params[c]
        bl = np.array(boots[g])
        med.append({"Group vs fibroids only": g, "RR M2": rr(r2, c)["txt"], "RR M3 (+ migraine)": rr(r3, c)["txt"],
                    "% of log RR explained by migraine": round(pc, 1),
                    "Bootstrap 95% CI": f"{np.percentile(bl, 2.5):.1f} to {np.percentile(bl, 97.5):.1f}" if len(bl) else "",
                    "Bootstrap reps": len(bl)})
    med = pd.DataFrame(med)
    add_table("P4_migraine_contribution", med, "Change in the log rate ratio (stroke or TIA) when migraine is added to "
              "M2 (difference method; percentile bootstrap, 300 resamples). Exploratory: migraine is recorded at any "
              "time, so its ordering relative to stroke is not established; unstable when the M2 RR is near 1.")
    RESULTS["p4_med"] = med
    # age strata (M2) + interaction
    st = []
    for ag in ["18–39", "40–60"]:
        s = f[f.age2 == ag]
        r, n, e, epv = fit(s, "ev_any", M2)
        for g in GROUPS[1:]:
            ev_g = int(s.loc[s.dxgrp == g, "ev_any"].sum())
            row = {"Age": ag, "Group vs fibroids only": g, "Events in group": ev_g, "N": n, "Events": e,
                   "EPV": round(float(epv), 1)}
            if ev_g < 5:
                row["RR (95% CI)"] = "not estimated (<5 events)"
            else:
                q = rr(r, f"dxgrp={g}")
                row.update({"RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"]})
            st.append(row)
    x = f.copy()
    x["old"] = (x.age2 == "40–60").astype(float)
    ic = []
    for g in GROUPS[1:]:
        x[f"ox_{g}"] = x.old * (x.dxgrp == g)
        ic.append(f"ox_{g}")
    r, _, _, _ = fit(x, "ev_any", M2 + [Term("old", "bin")] + [Term(c, "bin") for c in ic])
    bb = r.params[ic].values
    V = r.cov_params().loc[ic, ic].values
    chi = float(bb @ np.linalg.pinv(V) @ bb)
    p_int = float(stats.chi2.sf(chi, len(ic)))
    st = pd.DataFrame(st)
    add_table("P4_age_strata", st, f"M2 by age group; group × age interaction χ²({len(ic)}) = {chi:.2f}, "
                                   f"P = {fmt_p(p_int)}.")
    RESULTS["p4_strata"] = st
    RESULTS["p4_age_int_p"] = p_int
    # 4b: stroke type and detail-field completeness among all strokes in eligible women
    s = elig[elig.stroke_any == 1]
    typ = []
    for g in GROUPS + ["All"]:
        x = s if g == "All" else s[s.dxgrp == g]
        r_ = {"Group": g, "Strokes/TIA": len(x)}
        for k, v in STYPE.items():
            r_[v] = f"{int((x.stroke_type == k).sum())} ({100 * (x.stroke_type == k).mean():.1f}%)" if len(x) else ""
        typ.append(r_)
    add_table("P4_stroke_type", pd.DataFrame(typ), "Stroke type by uterine diagnosis group (all strokes/TIA in "
                                                   "eligible women, any timing).")
    RESULTS["p4_type"] = pd.DataFrame(typ)
    comp = []
    for c, lab, unk in [("stroke_type", "Stroke type", 9), ("stroke_timing", "Timing vs index", 9),
                        ("stroke_confirmed_imaging", "Imaging confirmation", None),
                        ("stroke_etiology", "Mechanism (imaging-derived, not TOAST)", 9),
                        ("vascular_territory", "Vascular territory", 9), ("num_infarct", "Number of lesions", 9),
                        ("pfo", "PFO", 9), ("cardiac_thrombus", "Cardiac thrombus", 9), ("nihss", "NIHSS", None)]:
        v = s[c]
        known = v.notna() & (v != unk) if unk is not None else v.notna()
        comp.append({"Field": lab, "Known": int(known.sum()), "Known %": round(100 * known.mean(), 1),
                     "Strokes/TIA": len(s)})
    add_table("P4_field_completeness", pd.DataFrame(comp), "Completeness of stroke detail fields among strokes/TIA in "
                                                           "eligible women (unknown codes counted as missing).")
    RESULTS["p4_comp"] = pd.DataFrame(comp)
    # migraine among TIA vs ischemic stroke (possible migraine mimics)
    mt = []
    for k, lab in [(1, "Ischemic stroke"), (4, "TIA"), (None, "Other/unknown type")]:
        x = s[s.stroke_type == k] if k is not None else s[~s.stroke_type.isin([1, 4])]
        mt.append({"Event type": lab, "Events": len(x), "With migraine": int(x.migraine_any.sum()),
                   "Migraine %": round(100 * x.migraine_any.mean(), 1)})
    a_, b_ = s[s.stroke_type == 4].migraine_any, s[s.stroke_type == 1].migraine_any
    tab = [[int(a_.sum()), int(len(a_) - a_.sum())], [int(b_.sum()), int(len(b_) - b_.sum())]]
    p_mt = float(stats.chi2_contingency(tab)[1])
    mt = pd.DataFrame(mt)
    add_table("P4_migraine_by_event", mt, f"Migraine among women with TIA vs ischemic stroke (all strokes/TIA in "
                                          f"eligible women); chi-square P = {fmt_p(p_mt)} (TIA vs ischemic).")
    RESULTS["p4_mig_event"] = mt
    RESULTS["p4_mig_event_p"] = p_mt
    return f
