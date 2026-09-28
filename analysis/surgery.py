"""Supplementary analysis (Paper 2 supplement): fibroid procedures and
subsequent stroke among women with fibroids.

Design (pre-specified before fitting)
- Cohort: eligible women with fibroids == 1.
- Outcome: stroke after index (stroke_timing == 3), dated by stroke_date.
  Women with a stroke before/at index (timing 1/2) are excluded; strokes with
  unknown or blank timing are excluded (cannot be placed in time).
  Secondary outcome: ischaemic stroke (stroke_type == 1) after index; other
  stroke types are censored at their date.
- Follow-up: index_date to the first of stroke or END, where END is the latest
  date recorded anywhere in the dataset. No death / transfer-out dates exist,
  so everyone is assumed to be followed to END (stated as a limitation).
- Exposure is time-varying: person-time is split at each procedure date.
  Procedures: myomectomy (myomectomy_date, any myomectomy in the record) and the
  surgical_tx procedure (surgical_tx_date: hysterectomy, uterine artery
  embolisation, endometrial ablation, other gynaecological surgery). A woman's
  current category is her most recent procedure (same-day ties: hysterectomy >
  myomectomy > UAE > ablation > other). Procedures dated before index expose
  from index. This avoids immortal-time bias: time before surgery is unexposed.
- "Any fibroid procedure" = myomectomy, hysterectomy or UAE (procedures that
  remove or devascularise fibroids); ablation / other surgery are not fibroid
  procedures and are analysed as their own categories. Sensitivity: any
  gynaecological procedure including ablation / other.
- Model: Poisson regression, log person-years offset, HC1 robust SEs; rate
  ratios (RR). Adjustment: age at index, race, BMI, hypertension, diabetes,
  dyslipidaemia, smoking (3-level), migraine (yes/no), atrial fibrillation,
  heavy/abnormal bleeding, adenomyosis, endometriosis (all at baseline).
- Head-to-head: among women with myomectomy or hysterectomy, follow-up from the
  first of these procedures; myomectomy vs hysterectomy.
- Sparse rule: any exposure level with < 5 events is reported as counts and
  rates only ("not estimated"); models with < 10 events per parameter fall back
  to a reduced set (age, BMI, hypertension, diabetes, dyslipidaemia, AF), then
  to age only if still < 10.
- "Other gynaecological surgery" and endometrial ablation do not treat fibroids
  and serve as informal negative-control exposures for healthy-candidate bias.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from .utils import RESULTS, Term, add_table, design, fmt_or, fmt_p, log

DATE_COLS = ["index_date", "stroke_date", "myomectomy_date", "hgb_date", "ferritin_date", "surgical_tx_date",
             "lab_date"]
TYPES = {2: "Hysterectomy", 3: "Uterine artery embolisation", 4: "Endometrial ablation", 5: "Other gynaecological surgery"}
PRIORITY = {"Hysterectomy": 0, "Myomectomy": 1, "Uterine artery embolisation": 2, "Endometrial ablation": 3,
            "Other gynaecological surgery": 4}
LEVELS = ["No procedure", "Myomectomy", "Hysterectomy", "Uterine artery embolisation", "Endometrial ablation",
          "Other gynaecological surgery"]
FIBROID_PROC = {"Myomectomy", "Hysterectomy", "Uterine artery embolisation"}
COV = [Term("age_index", label="Age (per year)"),
       Term("race4", "cat", ref="White", levels=["White", "Black", "Asian", "Other/unknown"], label="Race"),
       Term("bmi", label="BMI"), Term("htn", "bin"), Term("dm", "bin"), Term("dyslipidemia", "bin"),
       Term("smoking3", "cat", ref="Never", levels=["Never", "Ever", "Unknown"]), Term("migraine_any", "bin"),
       Term("afib", "bin"), Term("uterine_bleeding", "bin"), Term("adenomyosis", "bin"), Term("endometriosis", "bin")]
REDUCED = [t for t in COV if t.var in ("age_index", "bmi", "htn", "dm", "dyslipidemia", "afib")]


def build_cohort(d):
    end = max(pd.to_datetime(d[c], errors="coerce").max() for c in DATE_COLS)
    f = d[(d.eligible == 1) & (d.fibroids == 1)].copy()
    n0 = len(f)
    n_prior = int(f.stroke_timing.isin([1, 2]).sum())
    n_unk = int(((f.stroke_any == 1) & ~f.stroke_timing.isin([1, 2, 3])).sum())
    f = f[~f.stroke_timing.isin([1, 2])]
    f = f[~((f.stroke_any == 1) & ~f.stroke_timing.isin([3]))]
    f["t0"] = pd.to_datetime(f.index_date)
    sd = pd.to_datetime(f.stroke_date, errors="coerce")
    f["ev_any"] = ((f.stroke_any == 1) & (f.stroke_timing == 3)).astype(int)
    f["ev_isch"] = (f.ev_any.astype(bool) & (f.stroke_type == 1)).astype(int)
    f["t1"] = pd.to_datetime(np.where(f.ev_any == 1, sd, end))
    bad = int((f.t1 <= f.t0).sum())
    f = f[f.t1 > f.t0]
    RESULTS["surg_cohort"] = dict(end=end.date().isoformat(), n_fibroid=n0, n_prior=n_prior, n_unknown=n_unk,
                                  n_zero_time=bad, n=len(f), events=int(f.ev_any.sum()),
                                  events_isch=int(f.ev_isch.sum()))
    log("Surgery", f"Fibroid cohort {n0:,}; excluded {n_prior} with stroke before/at index, {n_unk} with unknown "
                   f"timing, {bad} with stroke on the index date; analysed {len(f):,}, {int(f.ev_any.sum())} incident "
                   f"strokes. Follow-up assumed to {end.date()} (latest date in the dataset).")
    return f


def procedures(r):
    ev = []
    if r.myomectomy == 1 and pd.notna(r.myomectomy_date):
        ev.append((pd.Timestamp(r.myomectomy_date), "Myomectomy"))
    if r.surgical_tx in TYPES and pd.notna(r.surgical_tx_date):
        ev.append((pd.Timestamp(r.surgical_tx_date), TYPES[r.surgical_tx]))
    return sorted(ev, key=lambda x: (x[0], -PRIORITY[x[1]]))  # same day: higher priority applied last


def split(f):
    rows = []
    for i, r in f.iterrows():
        cur, start = "No procedure", r.t0
        anyfib, anygyn = 0, 0
        segs = []
        for t, typ in procedures(r):
            if t <= r.t0:
                cur = typ
                anyfib |= typ in FIBROID_PROC
                anygyn = 1
                continue
            if t >= r.t1:
                break
            segs.append((start, t, cur, anyfib, anygyn))
            start, cur = t, typ
            anyfib |= typ in FIBROID_PROC
            anygyn = 1
        segs.append((start, r.t1, cur, anyfib, anygyn))
        for k, (s, e, c, a, g) in enumerate(segs):
            last = k == len(segs) - 1
            rows.append((i, s, e, c, a, g, int(r.ev_any and last), int(r.ev_isch and last)))
    x = pd.DataFrame(rows, columns=["i", "start", "end", "proc", "any_fibroid_proc", "any_gyn_proc", "ev_any", "ev_isch"])
    x["py"] = (x.end - x.start).dt.days / 365.25
    x = x[x.py > 0]
    return x.join(f.drop(columns=["ev_any", "ev_isch"]), on="i")


def poisson(x, outcome, expo_term, covs):
    need = [expo_term.var] + [t.var for t in covs]
    x = x.dropna(subset=need)
    X, cmap = design(x, [expo_term] + covs)
    X = X.loc[:, (X != 0).any(axis=0)]
    r = sm.GLM(x[outcome].values, X, family=sm.families.Poisson(), offset=np.log(x.py.values)).fit(cov_type="HC1")
    epv = x[outcome].sum() / (X.shape[1] - 1)
    return r, cmap[expo_term.var], x, epv


def rr_rows(x, outcome, expo_term, ref_label, label):
    """Rates per level + crude, age-adjusted and fully adjusted RRs (sparse rule applied)."""
    lv = expo_term.levels or [0, 1]
    g = x.groupby(expo_term.var).agg(events=(outcome, "sum"), py=("py", "sum"), women=("i", "nunique"))
    fits = {}
    for nm, cv in [("crude", []), ("age", COV[:1]), ("full", COV)]:
        r, cols, xs, epv = poisson(x, outcome, expo_term, cv)
        if nm == "full" and epv < 10:
            r, cols, xs, epv2 = poisson(x, outcome, expo_term, REDUCED)
            nm = "full (fallback: reduced set - age, BMI, HTN, DM, dyslipidaemia, AF)"
            if epv2 < 10:
                r, cols, xs, epv2 = poisson(x, outcome, expo_term, COV[:1])
                nm = "full (fallback: age only)"
            log("Surgery", f"{label} [{outcome}]: full model EPV {epv:.1f} < 10 -> {nm.split(': ', 1)[1]} "
                           f"(EPV {epv2:.1f})")
            epv = epv2
        fits[nm.split(" ")[0]] = (r, nm, epv)
    out = []
    for level in lv:
        key = level
        ev = int(g.events.get(key, 0))
        py = float(g.py.get(key, 0.0))
        row = {"Analysis": label, "Outcome": {"ev_any": "Any stroke", "ev_isch": "Ischaemic stroke"}[outcome],
               "Group": ref_label if level == expo_term.ref else (level if isinstance(level, str) else
                                                                     ("Yes" if level == 1 else "No")),
               "Women contributing time": int(g.women.get(key, 0)), "Events": ev, "Person-years": round(py, 1),
               "Rate /1,000 PY": round(1000 * ev / py, 2) if py else np.nan}
        if level == expo_term.ref:
            row.update({"Crude RR": "1.00 (reference)", "Age-adjusted RR": "1.00 (reference)",
                        "Adjusted RR": "1.00 (reference)"})
        else:
            col = f"{expo_term.var}={level}" if expo_term.kind == "cat" else expo_term.var
            for nm, (r, lab, epv) in fits.items():
                colname = {"crude": "Crude RR", "age": "Age-adjusted RR", "full": "Adjusted RR"}[nm]
                if ev < 5 or col not in r.params.index:
                    row[colname] = "not estimated (<5 events)"
                    continue
                b, se = r.params[col], r.bse[col]
                row[colname] = fmt_or(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se))
                if nm == "full":
                    row.update({"RR": np.exp(b), "CI low": np.exp(b - 1.96 * se), "CI high": np.exp(b + 1.96 * se),
                                "p": 2 * stats.norm.sf(abs(b / se)), "p (text)": fmt_p(2 * stats.norm.sf(abs(b / se))),
                                "Adjustment": lab, "EPV": round(epv, 1)})
        out.append(row)
    return out


def run(d):
    f = build_cohort(d)
    x = split(f)
    x["any_fibroid_proc"] = x["any_fibroid_proc"].astype(float)
    x["any_gyn_proc"] = x["any_gyn_proc"].astype(float)
    RESULTS["surg_py"] = float(x.py.sum())
    rows = []
    for oc in ["ev_any", "ev_isch"]:
        rows += rr_rows(x, oc, Term("any_fibroid_proc", "bin", ref=0.0, levels=[0.0, 1.0]),
                        "No fibroid procedure", "A. Any fibroid procedure (myomectomy/hysterectomy/UAE) vs none")
        rows += rr_rows(x, oc, Term("proc", "cat", ref="No procedure", levels=LEVELS), "No procedure",
                        "B. By most recent procedure vs no procedure")
        rows += rr_rows(x, oc, Term("any_gyn_proc", "bin", ref=0.0, levels=[0.0, 1.0]),
                        "No procedure", "S. Sensitivity: any gynaecological procedure (incl. ablation/other) vs none")
    # C. head-to-head from date of first myomectomy/hysterectomy
    h = []
    n_switch = 0
    for i, r in f.iterrows():
        ps = [(t, typ) for t, typ in procedures(r) if typ in ("Myomectomy", "Hysterectomy")]
        if not ps:
            continue
        t, typ = ps[0]
        s = max(t, r.t0)
        if s >= r.t1:
            continue
        # a myomectomy followed by hysterectomy: censor at the hysterectomy (switch)
        e = r.t1
        ev_any, ev_isch = r.ev_any, r.ev_isch
        later = [tt for tt, ty in ps[1:] if ty != typ and tt > s]
        if later and later[0] < e:
            e, ev_any, ev_isch = later[0], 0, 0
            n_switch += 1
        h.append((i, s, e, typ, ev_any, ev_isch))
    hh = pd.DataFrame(h, columns=["i", "start", "end", "first_proc", "ev_any", "ev_isch"])
    hh["py"] = (hh.end - hh.start).dt.days / 365.25
    hh = hh[hh.py > 0].join(f.drop(columns=["ev_any", "ev_isch"]), on="i")
    RESULTS["surg_h2h_n"] = int(hh.i.nunique())
    RESULTS["surg_h2h_switch"] = n_switch
    for oc in ["ev_any", "ev_isch"]:
        rows += rr_rows(hh, oc, Term("first_proc", "cat", ref="Hysterectomy", levels=["Hysterectomy", "Myomectomy"]),
                        "Hysterectomy", "C. Head-to-head from surgery date: myomectomy vs hysterectomy")
    tab = pd.DataFrame(rows)
    add_table("SUPP_surgery_stroke", tab,
              "Fibroid procedures and subsequent stroke (women with fibroids; strokes after index). Poisson rate "
              "ratios with time-varying exposure, HC1 SEs. Adjusted = age, race, BMI, HTN, DM, dyslipidaemia, smoking, "
              "migraine, AF, bleeding, adenomyosis, endometriosis. Levels with <5 events not estimated.")
    # baseline characteristics by ever-procedure type (descriptive, explains healthy-candidate selection)
    f2 = f.copy()
    f2["ever"] = ["Hysterectomy" if r.surgical_tx == 2 else
                  ("Myomectomy" if r.myomectomy == 1 else TYPES.get(r.surgical_tx, "No procedure"))
                  for r in f2.itertuples()]
    base = f2.groupby("ever").agg(n=("age_index", "size"), age=("age_index", "median"),
                                  htn=("htn", "mean"), dm=("dm", "mean"), afib=("afib", "mean"),
                                  bmi=("bmi", "median"), bleeding=("uterine_bleeding", "mean")).reset_index()
    for c in ["htn", "dm", "afib", "bleeding"]:
        base[c] = (100 * base[c]).round(1)
    base.columns = ["Ever procedure (hysterectomy > myomectomy > other)", "n", "Median age", "Hypertension %",
                    "Diabetes %", "AF %", "Median BMI", "Bleeding %"]
    add_table("SUPP_surgery_baseline", base, "Baseline characteristics by procedure ever received (healthy-candidate "
                                             "selection: operated women are younger and healthier).")
    RESULTS["surg_tab"] = tab
    RESULTS["surg_base"] = base
    return tab
