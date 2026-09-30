"""Paper 2, longitudinal extension: repeated CBC hemoglobin measurements (Lab_Tests_2, cleaned in labs.py).

Uses every dated CBC hemoglobin for each eligible woman (median 7 per woman) instead of the single Codebook value:

1. Measurement pattern (per woman, per year; before/after index).
2. Time-updated anemia exposure: follow-up from the index date (or the first Hb if none within the previous
   3 years) to stroke/TIA, death or last encounter, split at every Hb measurement; each interval carries the most
   recent Hb (last observation carried forward, valid for 3 years). Poisson regression with log person-time offset,
   HC1 SEs, Paper 2 covariates and pre-Hb conditions (2a). Sensitivity: 30-day lag (a Hb becomes effective 30 days
   after measurement, so values drawn during a stroke admission cannot define exposure), 1-year carry-forward cap,
   no cap, and adjustment for the number of Hb tests in the previous year (healthcare-contact proxy).
3. Cumulative anemia burden: years with Hb <12 accrued so far (time-varying), alone and with the current grade.
4. Persistence: among women anemic at baseline (Codebook Hb), the last Hb within 6, 12 and 24 months; landmark
   analysis at 12 months (primary) and 24 months comparing resolved vs persistent anemia with women not anemic at
   baseline, followed forward from the landmark.
5. Healthcare-contact proxy for the cross-sectional models: number of CBC Hb tests in the 2 years before the
   index date, added to the primary (2a) logistic models.
6. Hemoglobin trajectory after baseline by baseline grade (for an eFigure).
Sparse rule: exposure levels with <5 events are not estimated; EPV is reported.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from .paper2 import EXPO, cov_for
from .paper2_ext import GRADES, PRE_VARS
from .surgery import DATE_COLS
from .utils import RESULTS, Term, add_table, design, fit_logit, fmt_or, fmt_p, get_or, log, register_fit

LEVELS = ["None (Hb ≥12)"] + GRADES
CAP = 3 * 365
TESTS_PRE = ["0", "1–2", "3–5", "6+"]
TESTS_YR = ["0–1", "2–3", "4–11", "12+"]


def grade_of(h):
    h = np.asarray(h, float)
    return np.where(np.isnan(h), None, np.select([h >= 12, h >= 10, h >= 8], LEVELS[:3], LEVELS[3]))


def hb_long():
    long = RESULTS.get("_lab_long")
    if long is None:
        from .labs import load
        long, _ = load()
    return long[long.analyte == "Hb (g/dL)"][["mrn", "date", "v"]].sort_values(["mrn", "date"])


def base_cohort(d):
    """Eligible women without stroke before/at index or an undatable stroke; end = stroke, death or last encounter."""
    from .followup import encounters
    end_default = max(pd.to_datetime(d[c], errors="coerce").max() for c in DATE_COLS)
    x = d[~d.stroke_timing.isin([1, 2])]
    x = x[~((x.stroke_any == 1) & ~x.stroke_timing.isin([3]))].copy()
    x["idx"] = pd.to_datetime(x.index_date)
    x["sd"] = pd.to_datetime(x.stroke_date, errors="coerce")
    x["ev_any"] = ((x.stroke_any == 1) & (x.stroke_timing == 3)).astype(int)
    x["ev_isch"] = (x.ev_any.astype(bool) & (x.stroke_type == 1)).astype(int)
    lc, dth = encounters()
    lastc, dd = x.mrn.map(lc), x.mrn.map(dth)
    cens = pd.concat([lastc, dd], axis=1).min(axis=1).fillna(end_default)
    x["end"] = pd.to_datetime(np.where(x.ev_any == 1, x.sd, cens))
    x["died"] = ((x.ev_any == 0) & dd.notna() & (dd <= x.end)).astype(int)
    return x


def split_intervals(x, h, cap=CAP, lag=0):
    """One row per (woman, interval) with the Hb in force, cumulative anemic years, tests in the prior year."""
    hg = {m: (g.date.values.astype("datetime64[D]").astype("int64"), g.v.values) for m, g in h.groupby("mrn")}
    rows, n_no_hb, n_ev_gap = [], 0, 0
    for r in x.itertuples():
        if r.mrn not in hg:
            n_no_hb += 1
            continue
        dates0, vals = hg[r.mrn]
        dates = dates0 + lag
        idx = np.datetime64(r.idx, "D").astype("int64")
        end = np.datetime64(r.end, "D").astype("int64")
        if end <= idx:
            continue
        starts = np.maximum(dates, idx)
        nxt = np.append(dates[1:], np.iinfo("int64").max // 4)
        ends = np.minimum(np.minimum(nxt, dates + cap), end)
        keep = ends > starts
        if not keep.any():
            n_no_hb += 1
            continue
        cum = 0.0
        covered_end = ends[keep].max()
        for s_, e_, v_, d0 in zip(starts[keep], ends[keep], vals[keep], dates0[keep]):
            py = (e_ - s_) / 365.25
            n_prev = int(((dates0 > s_ - 365) & (dates0 <= s_)).sum())
            ev = int(r.ev_any == 1 and e_ == end)
            rows.append((r.Index, s_, e_, py, v_, cum, n_prev, ev, int(ev and r.ev_isch == 1)))
            cum += py if v_ < 12 else 0.0
        if r.ev_any == 1 and covered_end < end:
            n_ev_gap += 1
    f = pd.DataFrame(rows, columns=["i", "start", "end", "py", "hgb_t", "cum_anemic_y", "tests_prior_y", "ev_any",
                                    "ev_isch"])
    f = f[f.py > 0].copy()
    f["anemia_cat"] = grade_of(f.hgb_t)
    f["anemia"] = f.anemia_cat.map({lv: k for k, lv in enumerate(LEVELS)}).astype(float)
    f["tests_yr_cat"] = pd.cut(f.tests_prior_y, [-1, 1, 3, 11, 10 ** 6], labels=TESTS_YR).astype(str)
    return f, n_no_hb, n_ev_gap


REDUCED = [Term("age_index", "cont"), Term("bmi", "cont"), Term("htn", "bin"), Term("dm", "bin"),
           Term("dyslipidemia", "bin"), Term("migraine_any", "bin"), Term("ckd", "bin"),
           Term("haemoglobinopathy_any", "bin")]


def pois(df, oc, terms):
    df = df.dropna(subset=[t.var for t in terms])
    X, _ = design(df, terms)
    X = X.loc[:, (X != 0).any(axis=0)]
    r = sm.GLM(df[oc].values, X, family=sm.families.Poisson(), offset=np.log(df.py.values)).fit(cov_type="HC1")
    return r, df, df[oc].sum() / (X.shape[1] - 1)


def pois_stable(df, oc, expo_terms, cov, extra, label=""):
    """Pre-specified simplification: if EPV < 10 with the full covariate set, refit with the reduced set."""
    r, used, epv = pois(df, oc, expo_terms + cov + extra)
    note = ""
    if epv < 10:
        keep = [t for t in extra if t.var not in PRE_VARS]  # keep the term the model is about (e.g. test frequency)
        r, used, epv = pois(df, oc, expo_terms + REDUCED + keep)
        note = "EPV<10 with full covariates: reduced set (age, BMI, HTN, DM, dyslipidemia, migraine, CKD, haemoglobinopathy)"
        log("Paper 2 long", f"{label}: EPV<10 -> reduced covariate set (EPV {epv:.1f})")
    return r, used, epv, note


def rr(r, c):
    b, se = r.params[c], r.bse[c]
    return dict(RR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=2 * stats.norm.sf(abs(b / se)),
                txt=fmt_or(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)))


def run(d):
    """d: the Paper 2 extension frame (eligible women with comorbidity flags, y_isch, anemia_cat ...)."""
    h = hb_long()
    h = h[h.mrn.isin(d.mrn)]
    woman = d[["mrn", "index_date", "hgb", "hgb_date", "anemia_cat"]].copy()
    woman["idx"] = pd.to_datetime(woman.index_date)
    ext_a = [Term(v, "bin") for v in PRE_VARS]

    # 1. measurement pattern ---------------------------------------------------------------------------
    g = h.groupby("mrn").agg(n=("v", "size"), first=("date", "min"), last=("date", "max"))
    g["span_y"] = (g["last"] - g["first"]).dt.days / 365.25
    g["per_year"] = g.n / g.span_y.clip(lower=1 / 12)
    idx_map = woman.set_index("mrn").idx
    hh = h.join(idx_map, on="mrn")
    n_pre = hh[hh.date < hh.idx].groupby("mrn").size()
    n_post = hh[hh.date >= hh.idx].groupby("mrn").size()
    n_pre2 = hh[(hh.date < hh.idx) & (hh.date >= hh.idx - pd.Timedelta(days=730))].groupby("mrn").size()
    meas = pd.DataFrame([
        {"Item": "Eligible women", "Value": len(d)},
        {"Item": "Women with ≥1 CBC hemoglobin", "Value": len(g)},
        {"Item": "Total hemoglobin measurements", "Value": int(g.n.sum())},
        {"Item": "Measurements per woman, median (IQR)", "Value": f"{g.n.median():.0f} ({g.n.quantile(.25):.0f}–{g.n.quantile(.75):.0f})"},
        {"Item": "Women with ≥3 measurements", "Value": int((g.n >= 3).sum())},
        {"Item": "Women with ≥5 measurements", "Value": int((g.n >= 5).sum())},
        {"Item": "Span of measurements, years, median (IQR)", "Value": f"{g.span_y.median():.1f} ({g.span_y.quantile(.25):.1f}–{g.span_y.quantile(.75):.1f})"},
        {"Item": "Measurements per year (women with span ≥1 y), median (IQR)",
         "Value": f"{g[g.span_y >= 1].per_year.median():.1f} ({g[g.span_y >= 1].per_year.quantile(.25):.1f}–{g[g.span_y >= 1].per_year.quantile(.75):.1f})"},
        {"Item": "Women with a measurement before the index date", "Value": int(len(n_pre))},
        {"Item": "Women with a measurement after the index date", "Value": int(len(n_post))},
        {"Item": "Measurements in the 2 years before index, median (IQR) among all eligible",
         "Value": f"{n_pre2.reindex(d.mrn).fillna(0).median():.0f} ({n_pre2.reindex(d.mrn).fillna(0).quantile(.25):.0f}–{n_pre2.reindex(d.mrn).fillna(0).quantile(.75):.0f})"},
    ])
    add_table("P2L_measurements", meas, "Repeated CBC hemoglobin measurements (Lab_Tests_2 after cleaning), eligible women.")
    RESULTS["p2l_meas"] = meas

    # 5. healthcare-contact proxy (cross-sectional) ------------------------------------------------------
    d = d.copy()
    d["tests_2y_pre"] = n_pre2.reindex(d.mrn).fillna(0).values
    d["tests_pre_cat"] = pd.cut(d.tests_2y_pre, [-1, 0, 2, 5, 10 ** 6], labels=TESTS_PRE).astype(str)
    tp = Term("tests_pre_cat", "cat", ref="0", levels=TESTS_PRE, label="Hb tests in 2 y before index")
    hc = []
    for o, olab in [("stroke_any", "Any stroke"), ("y_isch", "Ischemic stroke"), ("y_incident", "Stroke after index")]:
        cv = cov_for(d, o, ["anemia_cat"], f"hc{o}")
        for mlab, extra in [("2a (reference)", ext_a), ("2a + Hb tests in 2 y before index", ext_a + [tp])]:
            f = register_fit(fit_logit(d, o, [EXPO["anemia_cat"]] + cv + extra, name=f"{o} {mlab}"), "P2L",
                             f"Healthcare-contact proxy: {olab}, {mlab}")
            ft = fit_logit(d, o, [Term("anemia", "cont")] + cv + extra, compute_vif=False)
            for lv in GRADES + ["Per grade (trend)"]:
                c = f"anemia_cat={lv}" if lv in GRADES else "anemia"
                q = get_or(ft if lv not in GRADES else f, c)
                hc.append({"Outcome": olab, "Model": mlab, "Level": lv, "OR (95% CI)": q["txt"], "OR": q["OR"],
                           "CI low": q["lo"], "CI high": q["hi"], "p": q["p"], "N": f.n, "Events": f.events})
            if "tests" in mlab:
                for lv in TESTS_PRE[1:]:
                    q = get_or(f, f"tests_pre_cat={lv}")
                    hc.append({"Outcome": olab, "Model": mlab, "Level": f"Hb tests before index: {lv} vs 0",
                               "OR (95% CI)": q["txt"], "OR": q["OR"], "CI low": q["lo"], "CI high": q["hi"],
                               "p": q["p"], "N": f.n, "Events": f.events})
    hc = pd.DataFrame(hc)
    hc["p (text)"] = hc.p.map(fmt_p)
    add_table("P2L_healthcare_proxy", hc, "Cross-sectional anemia-grade ORs (Paper 2 covariates + pre-Hb conditions "
              "2a) with and without adjustment for the number of CBC hemoglobin tests in the 2 years before the index "
              "date, a proxy for contact with the health system.")
    RESULTS["p2l_hc"] = hc
    RESULTS["p2l_tests_pre_dist"] = d.tests_pre_cat.value_counts().reindex(TESTS_PRE).to_dict()

    # 2/3. time-updated exposure ------------------------------------------------------------------------
    x = base_cohort(d)
    frames = {}
    for key, cap, lag in [("primary", CAP, 0), ("lag30", CAP, 30), ("cap1y", 365, 0), ("nocap", 10 ** 6, 0)]:
        f, n_no, n_gap = split_intervals(x, h, cap=cap, lag=lag)
        f = f.join(x.drop(columns=[c for c in ["ev_any", "ev_isch", "anemia_cat", "anemia", "end", "idx", "sd"] if c in x.columns]),
                   on="i")
        frames[key] = (f, n_no, n_gap)
    f0, n_no, n_gap = frames["primary"]
    RESULTS["p2l_tu_cohort"] = dict(n=int(f0.i.nunique()), py=float(f0.py.sum()), ev=int(f0.ev_any.sum()),
                                    ev_isch=int(f0.ev_isch.sum()), intervals=len(f0), n_no_hb=n_no, ev_in_gap=n_gap,
                                    excluded_prior=int(d.stroke_timing.isin([1, 2]).sum()),
                                    median_intervals=float(f0.groupby("i").size().median()))
    log("Paper 2 long", f"Time-updated cohort: {f0.i.nunique():,} women, {f0.py.sum():,.0f} person-years, "
                        f"{int(f0.ev_any.sum())} strokes/TIA ({int(f0.ev_isch.sum())} ischemic), {len(f0):,} intervals; "
                        f"{n_no} women without a usable Hb; {n_gap} strokes fell in person-time without a Hb in the "
                        f"previous 3 years and are not counted.")
    tu = []
    tests_t = Term("tests_yr_cat", "cat", ref="0–1", levels=TESTS_YR, label="Hb tests in previous year")
    specs = [("Time-updated, Paper 2 covariates", "primary", []),
             ("Time-updated + pre-Hb conditions (2a)", "primary", ext_a),
             ("Time-updated + 2a + Hb tests in previous year", "primary", ext_a + [tests_t]),
             ("Time-updated + 2a, 30-day lag", "lag30", ext_a),
             ("Time-updated + 2a, 1-year carry-forward", "cap1y", ext_a),
             ("Time-updated + 2a, no carry-forward limit", "nocap", ext_a)]
    for oc, olab in [("ev_any", "Stroke or TIA"), ("ev_isch", "Ischemic stroke")]:
        for mlab, key, extra in specs:
            f = frames[key][0]
            cv = cov_for(f, oc, ["anemia_cat"], f"tu{key}{oc}")
            r, used, epv, note = pois_stable(f, oc, [EXPO["anemia_cat"]], cv, extra, f"{olab} {mlab}")
            rt, _, _ = pois(f, oc, [Term("anemia", "cont")] + (REDUCED if note else cv + extra))
            gg = used.groupby("anemia_cat").agg(ev=(oc, "sum"), py=("py", "sum"), w=("i", "nunique"))
            for lv in LEVELS + ["Per grade (trend)"]:
                row = {"Outcome": olab, "Model": mlab, "Level": lv, "EPV": round(float(epv), 1), "Note": note,
                       "Women": int(f.i.nunique()), "Person-years (level)": round(float(gg.py.get(lv, np.nan)), 1) if lv in gg.index else np.nan,
                       "Events (level)": int(gg.ev.get(lv, 0)) if lv in gg.index else np.nan,
                       "Rate /1,000 PY": round(1000 * gg.ev.get(lv, 0) / gg.py.get(lv, np.nan), 2) if lv in gg.index else np.nan}
                if lv == LEVELS[0]:
                    row["RR (95% CI)"] = "1.00 (reference)"
                elif lv == "Per grade (trend)":
                    q = rr(rt, "anemia")
                    row.update({"RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
                elif gg.ev.get(lv, 0) < 5:
                    row["RR (95% CI)"] = "not estimated (<5 events)"
                else:
                    q = rr(r, f"anemia_cat={lv}")
                    row.update({"RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
                tu.append(row)
            if "tests in previous year" in mlab and oc == "ev_any":
                for lv in TESTS_YR[1:]:
                    q = rr(r, f"tests_yr_cat={lv}")
                    tu.append({"Outcome": olab, "Model": mlab, "Level": f"Hb tests in previous year: {lv} vs 0–1",
                               "RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
    tu = pd.DataFrame(tu)
    tu["p (text)"] = tu.p.map(lambda v: fmt_p(v) if pd.notna(v) else "")
    add_table("P2L_time_updated", tu, "Time-updated anemia exposure (most recent CBC hemoglobin, carried forward up "
              "to 3 years; person-time split at each measurement). Poisson regression, log person-time offset, HC1 "
              "SEs. Reference Hb ≥12 g/dL.")
    RESULTS["p2l_tu"] = tu
    # cumulative burden
    cum = []
    f = frames["primary"][0].copy()
    f["cum_cat"] = pd.cut(f.cum_anemic_y, [-1, 0, 1, 10 ** 6], labels=["0", ">0–1 y", ">1 y"]).astype(str)
    cc = Term("cum_cat", "cat", ref="0", levels=["0", ">0–1 y", ">1 y"], label="Years anemic so far")
    for oc, olab in [("ev_any", "Stroke or TIA"), ("ev_isch", "Ischemic stroke")]:
        cv = cov_for(f, oc, ["anemia_cat"], f"cum{oc}")
        for mlab, terms in [("Years anemic so far (per year), 2a", [Term("cum_anemic_y", "cont")] + cv + ext_a),
                            ("Years anemic so far (per year) + current grade, 2a",
                             [Term("cum_anemic_y", "cont"), EXPO["anemia_cat"]] + cv + ext_a),
                            ("Years anemic so far (categories), 2a", [cc] + cv + ext_a)]:
            r, used, epv = pois(f, oc, terms)
            note = ""
            if epv < 10:
                r, used, epv = pois(f, oc, [t for t in terms if t.var in ("cum_anemic_y", "anemia_cat", "cum_cat")] + REDUCED)
                note = "EPV<10: reduced covariate set"
            if "categories" in mlab:
                gg = used.groupby("cum_cat").agg(ev=(oc, "sum"), py=("py", "sum"))
                for lv in ["0", ">0–1 y", ">1 y"]:
                    row = {"Outcome": olab, "Model": mlab, "Term": f"Years anemic so far: {lv}", "EPV": round(float(epv), 1), "Note": note,
                           "Events (level)": int(gg.ev.get(lv, 0)), "Person-years (level)": round(float(gg.py.get(lv, 0)), 1)}
                    if lv == "0":
                        row["RR (95% CI)"] = "1.00 (reference)"
                    else:
                        q = rr(r, f"cum_cat={lv}")
                        row.update({"RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
                    cum.append(row)
            else:
                q = rr(r, "cum_anemic_y")
                cum.append({"Outcome": olab, "Model": mlab, "Term": "Per additional year anemic", "EPV": round(float(epv), 1), "Note": note,
                            "RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
                if "current grade" in mlab:
                    for lv in GRADES:
                        q = rr(r, f"anemia_cat={lv}")
                        cum.append({"Outcome": olab, "Model": mlab, "Term": f"Current grade: {lv}", "EPV": round(float(epv), 1),
                                    "RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
    cum = pd.DataFrame(cum)
    cum["p (text)"] = cum.p.map(lambda v: fmt_p(v) if pd.notna(v) else "")
    add_table("P2L_cumulative", cum, "Cumulative anemia burden: person-years with Hb <12 g/dL accrued before each "
                                     "interval (time-varying), Poisson regression with 2a adjustment.")
    RESULTS["p2l_cum"] = cum

    # 4. persistence and landmark ------------------------------------------------------------------------
    base = x[x.hgb.notna()].copy()
    base["hbd"] = pd.to_datetime(base.hgb_date, errors="coerce")
    base = base[base.hbd.notna()]
    hg = {m: (gq.date.values, gq.v.values) for m, gq in h.groupby("mrn")}
    pers, lm_rows = [], []
    for win in [182, 365, 730]:
        last_v, has = {}, {}
        for r in base.itertuples():
            if r.mrn not in hg:
                continue
            dts, vs = hg[r.mrn]
            m = (dts > np.datetime64(r.hbd)) & (dts <= np.datetime64(r.hbd) + np.timedelta64(win, "D"))
            if m.any():
                last_v[r.Index] = vs[m][-1]
        lv = pd.Series(last_v)
        base[f"last_{win}"] = lv.reindex(base.index)
        for gr in GRADES:
            s = base[base.anemia_cat == gr]
            l = s[f"last_{win}"]
            pers.append({"Baseline grade": gr, "Window (months)": round(win / 30.4), "Women": len(s),
                         "With a repeat Hb": int(l.notna().sum()), "Repeat Hb %": round(100 * l.notna().mean(), 1),
                         "Still <12 g/dL at last value": int((l < 12).sum()),
                         "Still anemic % (of those with repeat)": round(100 * (l < 12).sum() / max(l.notna().sum(), 1), 1),
                         "Median change in Hb, g/dL": round(float((l - s.hgb).median()), 2) if l.notna().any() else np.nan})
    pers = pd.DataFrame(pers)
    add_table("P2L_persistence", pers, "Among women anemic at baseline (Codebook Hb): last CBC hemoglobin within "
                                       "6, 12 and 24 months of the baseline value.")
    RESULTS["p2l_pers"] = pers
    LM_LEVELS = ["Not anemic at baseline", "Not anemic at baseline, Hb <12 within window (new anemia)",
                 "Anemic at baseline, resolved (last Hb ≥12)", "Anemic at baseline, persistent (last Hb <12)",
                 "Anemic at baseline, no repeat Hb"]
    for win in [365, 730]:
        b = base.copy()
        b["L"] = b.hbd + pd.Timedelta(days=win)
        l = b[f"last_{win}"]
        an = b.anemia_cat.isin(GRADES)
        b["lm_grp"] = np.select([~an & ~(l < 12), ~an & (l < 12), an & (l >= 12), an & (l < 12), an & l.isna()],
                                LM_LEVELS, None)
        b = b[b.lm_grp.notna() & (b.end > b.L)].copy()
        b = b[~((b.ev_any == 1) & (b.sd <= b.L))]
        b["py"] = (b.end - b.L).dt.days / 365.25
        RESULTS[f"p2l_lm{win}_n"] = dict(n=len(b), ev=int(b.ev_any.sum()), py=float(b.py.sum()))
        lt = Term("lm_grp", "cat", ref=LM_LEVELS[0], levels=LM_LEVELS, label="Anemia status at landmark")
        for oc, olab in [("ev_any", "Stroke or TIA"), ("ev_isch", "Ischemic stroke")]:
            cv = cov_for(b, oc, ["lm_grp"], f"lm{win}{oc}")
            r, used, epv, note = pois_stable(b, oc, [lt], cv, ext_a, f"Landmark {win} d {olab}")
            gg = used.groupby("lm_grp").agg(ev=(oc, "sum"), py=("py", "sum"), w=("mrn", "size"))
            # persistent vs resolved contrast
            bp, br_ = r.params[f"lm_grp={LM_LEVELS[3]}"], r.params[f"lm_grp={LM_LEVELS[2]}"]
            V = r.cov_params()
            se = np.sqrt(V.loc[f"lm_grp={LM_LEVELS[3]}", f"lm_grp={LM_LEVELS[3]}"] + V.loc[f"lm_grp={LM_LEVELS[2]}", f"lm_grp={LM_LEVELS[2]}"]
                         - 2 * V.loc[f"lm_grp={LM_LEVELS[3]}", f"lm_grp={LM_LEVELS[2]}"])
            pv = 2 * stats.norm.sf(abs((bp - br_) / se))
            for lv in LM_LEVELS:
                row = {"Landmark (months)": round(win / 30.4), "Outcome": olab, "Group": lv, "Women": int(gg.w.get(lv, 0)),
                       "Events": int(gg.ev.get(lv, 0)), "Person-years": round(float(gg.py.get(lv, 0)), 1),
                       "Rate /1,000 PY": round(1000 * gg.ev.get(lv, 0) / gg.py.get(lv, np.nan), 2), "EPV": round(float(epv), 1),
                       "Note": note}
                if lv == LM_LEVELS[0]:
                    row["Adjusted RR (95% CI)"] = "1.00 (reference)"
                elif gg.ev.get(lv, 0) < 5:
                    row["Adjusted RR (95% CI)"] = "not estimated (<5 events)"
                else:
                    q = rr(r, f"lm_grp={lv}")
                    row.update({"Adjusted RR (95% CI)": q["txt"], "RR": q["RR"], "CI low": q["lo"], "CI high": q["hi"], "p": q["p"]})
                lm_rows.append(row)
            lm_rows.append({"Landmark (months)": round(win / 30.4), "Outcome": olab, "Group": "Persistent vs resolved",
                            "Adjusted RR (95% CI)": fmt_or(np.exp(bp - br_), np.exp(bp - br_ - 1.96 * se), np.exp(bp - br_ + 1.96 * se)),
                            "RR": np.exp(bp - br_), "CI low": np.exp(bp - br_ - 1.96 * se), "CI high": np.exp(bp - br_ + 1.96 * se), "p": pv})
    lm = pd.DataFrame(lm_rows)
    lm["p (text)"] = lm.p.map(lambda v: fmt_p(v) if pd.notna(v) else "")
    add_table("P2L_landmark", lm, "Landmark analysis: anemia status defined by the last CBC hemoglobin within 12 (or "
                                  "24) months of the baseline Hb; follow-up from the landmark to stroke/TIA, death or last "
                                  "encounter. Poisson, Paper 2 covariates + pre-Hb conditions (2a).")
    RESULTS["p2l_lm"] = lm

    # 6. trajectory after baseline by baseline grade -----------------------------------------------------
    bb = base[["mrn", "hbd", "anemia_cat"]]
    hb2 = h.merge(bb, on="mrn")
    hb2["m"] = (hb2.date - hb2.hbd).dt.days / 30.44
    hb2 = hb2[(hb2.m > 0) & (hb2.m <= 36)]
    hb2["bin"] = pd.cut(hb2.m, [0, 3, 6, 12, 24, 36], labels=["0–3", "4–6", "7–12", "13–24", "25–36"])
    tr = hb2.groupby(["anemia_cat", "bin"], observed=True).agg(n=("v", "size"), women=("mrn", "nunique"),
                                                               mean=("v", "mean"), sd=("v", "std")).reset_index()
    tr["se"] = tr.sd / np.sqrt(tr.n)
    tr = tr.rename(columns={"anemia_cat": "Baseline grade", "bin": "Months after baseline", "mean": "Mean Hb", "sd": "SD"})
    base_mean = base.groupby("anemia_cat").hgb.agg(["mean", "size"]).reset_index()
    RESULTS["p2l_traj_base"] = base_mean
    add_table("P2L_trajectory", tr, "Mean CBC hemoglobin after the baseline value, by baseline anemia grade and months "
                                    "since baseline (all measurements, one woman may contribute several per bin).")
    RESULTS["p2l_traj"] = tr
    return frames["primary"][0]
