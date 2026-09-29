"""Shared helpers: data loading/derivation, logistic models with HC1 SEs,
restricted cubic splines, diagnostics (VIF, EPV, sparse cells), MICE with
Rubin's rules, Table 1 builder, and a results registry used by the reports.

Every number that appears in the markdown reports is written into RESULTS
(or a table in TABLES) by code in this package; the report templates only
read from there.
"""
from __future__ import annotations

import math
import re
import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"
MASTER_CSV = DATA_DIR / "stroke_analysis_master.csv"
WORKBOOK = DATA_DIR / "stroke.xlsx"

PHI_COLUMNS = ["mrn", "dob"]

# ---------------------------------------------------------------------------
# Registries: tables (-> xlsx sheets), scalar results (-> markdown), log lines
# ---------------------------------------------------------------------------
TABLES: "OrderedDict[str, pd.DataFrame]" = OrderedDict()
TABLE_NOTES: dict[str, str] = {}
RESULTS: dict = {}
LOG: list[tuple[str, str]] = []


def add_table(name: str, df: pd.DataFrame, note: str | None = None) -> None:
    assert len(name) <= 31, name
    bad = [c for c in df.columns if str(c).lower() in PHI_COLUMNS]
    assert not bad, f"PHI column in output table {name}: {bad}"
    TABLES[name] = df.reset_index(drop=True)
    if note:
        TABLE_NOTES[name] = note


def log(section: str, text: str) -> None:
    LOG.append((section, text))
    print(f"[{section}] {text}")


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------
def fmt_or(o, lo, hi, digits=2) -> str:
    if any(pd.isna(v) for v in (o, lo, hi)):
        return ""
    return f"{o:.{digits}f} ({lo:.{digits}f}–{hi:.{digits}f})"


def fmt_p(p) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def fmt_n_pct(n, d, digits=1) -> str:
    if d == 0:
        return f"{n}/0"
    return f"{n:,} ({100 * n / d:.{digits}f}%)"


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return p, c - h, c + h


def fmt_prev(k, n) -> str:
    p, lo, hi = wilson(k, n)
    return f"{k:,}/{n:,} = {100 * p:.2f}% (95% CI {100 * lo:.2f}–{100 * hi:.2f})"


# ---------------------------------------------------------------------------
# Loading and derivation
# ---------------------------------------------------------------------------
def load_master() -> pd.DataFrame:
    df = pd.read_csv(MASTER_CSV, low_memory=False)
    return df


def race4(row_race, black) -> str:
    r = str(row_race).strip() if pd.notna(row_race) else ""
    if black == 1:
        return "Black"
    if r == "White":
        return "White"
    if r.lower().startswith("asian"):
        return "Asian"
    return "Other/unknown"


def derive(df: pd.DataFrame) -> pd.DataFrame:
    """Analysis variables. Every derivation here follows the Codebook and the
    coding rules in the analysis brief; none changes a Codebook definition."""
    d = df.copy()
    d["race4"] = [race4(r, b) for r, b in zip(d["race"], d["black_race"])]
    # smoking: 3-level never / ever (current or former) / unknown (code 9)
    d["smoking3"] = d["smoking"].map({0: "Never", 1: "Ever", 2: "Ever", 9: "Unknown"})
    # migraine (PI decision 2026-09-28): type and aura are not used. Any migraine
    # diagnosis (codes 1 without aura, 2 with aura, 9 aura/type not stated) = 1,
    # code 0 = 0. Used as the covariate in Papers 1-2 and the outcome in Paper 3.
    d["anemia_cat"] = d["anemia"].map({0: "None (Hb ≥12)", 1: "Mild (10–11.9)",
                                       2: "Moderate (8–9.9)", 3: "Severe (<8)"})
    d["anemia_any"] = np.where(d["anemia"].isna(), np.nan, (d["anemia"] >= 1).astype(float))
    d["mcv_cat"] = pd.Series(np.select(
        [d["mcv"] < 80, d["mcv"] <= 100, d["mcv"] > 100],
        ["Microcytic (<80)", "Normal (80–100)", "Macrocytic (>100)"], default=None),
        index=d.index).where(d["mcv"].notna())
    d["plt_cat"] = pd.Series(np.select(
        [d["platelets"] < 150, d["platelets"] <= 400, d["platelets"] > 400],
        ["Low (<150)", "Normal (150–400)", "Thrombocytosis (>400)"], default=None),
        index=d.index).where(d["platelets"].notna())
    d["microcytosis"] = np.where(d["mcv"].isna(), np.nan, (d["mcv"] < 80).astype(float))
    d["thrombocytosis"] = np.where(d["platelets"].isna(), np.nan, (d["platelets"] > 400).astype(float))
    d["dxgrp"] = d["uterine_dx_group"].map({1: "Fibroids only", 2: "Adenomyosis only",
                                            3: "Endometriosis only", 4: ">1 condition"})
    d["hormonal_type_cat"] = d["hormonal_type"].map({
        0: "None", 1: "Combined OC", 2: "Progestin-only", 3: "LNG-IUD",
        4: "GnRH agonist/antagonist", 5: "Menopausal HT", 6: "Other/multiple"})
    d["obesity"] = np.where(d["bmi"].isna(), np.nan, (d["bmi"] >= 30).astype(float))
    d["smoking_ever"] = d["smoking"].map({0: 0.0, 1: 1.0, 2: 1.0})
    d["migraine_any"] = d["migraine"].map({0: 0.0, 1: 1.0, 2: 1.0, 9: 1.0})
    for c in ["htn", "dm", "dyslipidemia", "cad", "afib", "chf", "vte_history",
              "thrombophilia", "uterine_bleeding", "fibroids", "adenomyosis",
              "endometriosis", "stroke_any"]:
        d[c] = d[c].where(d[c] != 9)
    d["ischaemic"] = np.where(d["stroke_any"] == 1, (d["stroke_type"] == 1).astype(float), 0.0)
    d["age_band"] = pd.cut(d["age_index"], [17, 24, 29, 34, 39, 44, 49, 54, 60],
                           labels=["18–24", "25–29", "30–34", "35–39", "40–44",
                                   "45–49", "50–54", "55–60"])
    d["age2"] = np.where(d["age_index"] < 40, "18–39", "40–60")
    # Collapsed versions used ONLY when a pre-specified model is too sparse
    # (EPV < 10 or a level with < 5 events). Not a redefinition of any variable.
    d["race3"] = d["race4"].replace({"Asian": "Asian/Other/unknown", "Other/unknown": "Asian/Other/unknown"})
    d["dxgrp3"] = d["dxgrp"].replace({"Adenomyosis only": "Adenomyosis only or >1",
                                      ">1 condition": "Adenomyosis only or >1"})
    d["anemia3"] = d["anemia_cat"].replace({"Moderate (8–9.9)": "Moderate/severe (<10)",
                                            "Severe (<8)": "Moderate/severe (<10)"})
    return d


# ---------------------------------------------------------------------------
# Design matrices
# ---------------------------------------------------------------------------
class Term:
    """A model term. kind: 'cont' | 'bin' | 'cat' | 'rcs'."""

    def __init__(self, var, kind="cont", ref=None, label=None, levels=None, knots=None):
        self.var, self.kind, self.ref, self.levels = var, kind, ref, levels
        self.label = label or var
        self.knots = knots

    def __repr__(self):
        return f"Term({self.var},{self.kind})"


def rcs_knots(x, k=4):
    q = {3: [10, 50, 90], 4: [5, 35, 65, 95], 5: [5, 27.5, 50, 72.5, 95]}[k]
    return np.nanpercentile(x, q)


def rcs_basis(x, knots):
    """Harrell restricted cubic spline basis (linear term + k-2 non-linear)."""
    x = np.asarray(x, float)
    t = np.asarray(knots, float)
    k = len(t)
    norm = (t[-1] - t[0]) ** 2
    cols = [x]
    for j in range(k - 2):
        c = (np.maximum(x - t[j], 0) ** 3
             - np.maximum(x - t[k - 2], 0) ** 3 * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
             + np.maximum(x - t[k - 1], 0) ** 3 * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2]))
        cols.append(c / norm)
    return np.column_stack(cols)


def design(df: pd.DataFrame, terms: list[Term]):
    """Return X (with const) and a map term.var -> list of column names."""
    cols = OrderedDict()
    colmap = OrderedDict()
    for t in terms:
        s = df[t.var]
        if t.kind in ("cont", "bin"):
            cols[t.var] = s.astype(float).values
            colmap[t.var] = [t.var]
        elif t.kind == "rcs":
            kn = t.knots if t.knots is not None else rcs_knots(s.values)
            t.knots = kn
            B = rcs_basis(s.values, kn)
            names = [f"{t.var}__rcs{i}" for i in range(B.shape[1])]
            for i, n in enumerate(names):
                cols[n] = B[:, i]
            colmap[t.var] = names
        elif t.kind == "cat":
            levels = t.levels if t.levels is not None else sorted(s.dropna().unique(), key=str)
            ref = t.ref if t.ref is not None else levels[0]
            names = []
            for lv in levels:
                if lv == ref:
                    continue
                n = f"{t.var}={lv}"
                cols[n] = (s == lv).astype(float).values
                names.append(n)
            colmap[t.var] = names
        else:
            raise ValueError(t.kind)
    X = pd.DataFrame(cols, index=df.index)
    X.insert(0, "const", 1.0)
    return X, colmap


def complete_rows(df, outcome, terms):
    vars_ = [outcome] + [t.var for t in terms]
    return df.dropna(subset=vars_)


# ---------------------------------------------------------------------------
# Logistic regression with HC1 SEs + diagnostics
# ---------------------------------------------------------------------------
class Fit:
    pass


def fit_logit(df, outcome, terms, name="", robust=True, compute_vif=True,
              min_level_events=5):
    d = complete_rows(df, outcome, terms)
    X, colmap = design(d, terms)
    # drop all-zero dummy columns (empty level) and report
    empty = [c for c in X.columns if c != "const" and X[c].abs().sum() == 0]
    y = d[outcome].astype(float).values
    f = Fit()
    f.name, f.outcome, f.terms, f.colmap = name, outcome, terms, colmap
    f.n, f.events = len(d), int(y.sum())
    f.empty_levels = empty
    f.data_index = d.index
    X = X.drop(columns=empty)
    for k in colmap:
        colmap[k] = [c for c in colmap[k] if c not in empty]
    f.X = X
    f.df_model = X.shape[1] - 1
    f.epv = f.events / f.df_model if f.df_model else np.nan
    f.epv_min = min(f.events, f.n - f.events) / f.df_model if f.df_model else np.nan
    # sparse cells: indicator columns with < min_level_events events
    sparse = []
    for c in X.columns[1:]:
        v = X[c].values
        if set(np.unique(v)) <= {0.0, 1.0}:
            ev = int(y[v == 1].sum())
            nn = int((v == 1).sum())
            nonev = nn - ev
            if ev < min_level_events or nonev < min_level_events:
                sparse.append(f"{c} (n={nn}, events={ev})")
    f.sparse = sparse
    model = sm.GLM(y, X, family=sm.families.Binomial())
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = model.fit(cov_type="HC1" if robust else "nonrobust", maxiter=100)
        f.converged = bool(res.converged)
    except Exception as exc:  # perfect separation etc.
        f.converged = False
        f.error = str(exc)
        f.res = None
        f.flags = [f"fit failed: {exc}"]
        return f
    f.res = res
    f.params = res.params
    f.cov = res.cov_params()
    f.bse = res.bse
    big = [c for c in X.columns[1:] if abs(res.params[c]) > 10 or res.bse[c] > 10]
    f.separation = big
    f.vif = {}
    if compute_vif and X.shape[1] > 2:
        Xv = X.values
        for i, c in enumerate(X.columns):
            if c == "const":
                continue
            try:
                f.vif[c] = variance_inflation_factor(Xv, i)
            except Exception:
                f.vif[c] = np.nan
    f.max_vif = max(f.vif.values()) if f.vif else np.nan
    flags = []
    if not f.converged:
        flags.append("did not converge")
    if f.epv < 10:
        flags.append(f"EPV {f.epv:.1f} < 10")
    if sparse:
        flags.append("sparse: " + "; ".join(sparse))
    if big:
        flags.append("possible separation: " + ", ".join(big))
    if f.max_vif == f.max_vif and f.max_vif > 5:
        flags.append(f"max VIF {f.max_vif:.1f} > 5")
    f.flags = flags
    f.aic = res.aic
    return f


def wald(f, cols):
    cols = [c for c in cols if c in f.params.index]
    if not cols:
        return np.nan, 0, np.nan
    b = f.params[cols].values
    V = f.cov.loc[cols, cols].values
    chi2 = float(b @ np.linalg.pinv(V) @ b)
    df_ = len(cols)
    return chi2, df_, float(stats.chi2.sf(chi2, df_))


def or_table(f, only_vars=None, model_label=None, labels=None):
    rows = []
    if f.res is None:
        return pd.DataFrame([{"Model": model_label or f.name, "Term": "MODEL FAILED",
                              "Flags": "; ".join(f.flags)}])
    z = stats.norm.ppf(0.975)
    for t in f.terms:
        if only_vars and t.var not in only_vars:
            continue
        cols = f.colmap[t.var]
        if t.kind == "rcs":
            chi2, df_, p = wald(f, cols[1:])
            rows.append({"Term": t.label, "Level": f"RCS ({len(t.knots)} knots)",
                         "OR (95% CI)": "see spline figure", "p": p,
                         "p (text)": f"non-linearity p={fmt_p(p)}"})
            continue
        if t.kind == "cat" and t.ref is not None:
            rows.append({"Term": t.label, "Level": f"{t.ref} (ref)", "OR (95% CI)": "1.00 (reference)"})
        for c in cols:
            b, se = f.params[c], f.bse[c]
            o, lo, hi = np.exp(b), np.exp(b - z * se), np.exp(b + z * se)
            p = 2 * stats.norm.sf(abs(b / se))
            lvl = c.split("=", 1)[1] if "=" in c else ("per 1 unit" if t.kind == "cont" else "yes vs no")
            if labels and c in labels:
                lvl = labels[c]
            X = f.X
            nlev = int(X[c].sum()) if set(np.unique(X[c])) <= {0.0, 1.0} else None
            ev = int(f.res.model.endog[X[c].values == 1].sum()) if nlev is not None else None
            rows.append({"Term": t.label, "Level": lvl, "OR (95% CI)": fmt_or(o, lo, hi),
                         "OR": o, "CI low": lo, "CI high": hi, "p": p, "p (text)": fmt_p(p),
                         "n in level": nlev, "events in level": ev})
        if t.kind == "cat" and len(cols) > 1:
            chi2, df_, p = wald(f, cols)
            rows.append({"Term": t.label, "Level": f"global Wald χ²({df_})",
                         "OR (95% CI)": "", "p": p, "p (text)": fmt_p(p)})
    out = pd.DataFrame(rows)
    out.insert(0, "Model", model_label or f.name)
    out["N"] = f.n
    out["Events"] = f.events
    out["EPV"] = round(f.epv, 1)
    out["Flags"] = "; ".join(f.flags)
    return out


def get_or(f, col):
    z = stats.norm.ppf(0.975)
    b, se = f.params[col], f.bse[col]
    return dict(OR=np.exp(b), lo=np.exp(b - z * se), hi=np.exp(b + z * se),
                p=2 * stats.norm.sf(abs(b / se)), txt=fmt_or(np.exp(b), np.exp(b - z * se), np.exp(b + z * se)))


def diag_row(f, paper, label):
    return {"Paper": paper, "Model": label, "Outcome": f.outcome, "N": f.n, "Events": f.events,
            "Model df": getattr(f, "df_model", np.nan), "EPV": round(f.epv, 2) if f.epv == f.epv else np.nan,
            "Converged": getattr(f, "converged", False),
            "Max VIF": round(f.max_vif, 2) if getattr(f, "max_vif", np.nan) == getattr(f, "max_vif", np.nan) else np.nan,
            "VIF > 5": ", ".join(k for k, v in getattr(f, "vif", {}).items() if v > 5),
            "Sparse cells (<5 events or non-events)": "; ".join(getattr(f, "sparse", [])),
            "Empty levels dropped": ", ".join(getattr(f, "empty_levels", [])),
            "Flags": "; ".join(f.flags)}


DIAG: list[dict] = []


def register_fit(f, paper, label):
    DIAG.append(diag_row(f, paper, label))
    return f


def linearity_rows(df, outcome, terms, cont_vars, paper, model):
    """For each continuous term: refit with that term as RCS(4) and test the
    non-linear components (robust Wald)."""
    rows = []
    for v in cont_vars:
        t2 = [Term(t.var, "rcs", label=t.label) if t.var == v else t for t in terms]
        f = fit_logit(df, outcome, t2, compute_vif=False)
        if f.res is None:
            rows.append({"Paper": paper, "Model": model, "Variable": v, "Note": "fit failed"})
            continue
        chi2, df_, p = wald(f, f.colmap[v][1:])
        kn = [t for t in t2 if t.var == v][0].knots
        rows.append({"Paper": paper, "Model": model, "Variable": v,
                     "Knots": ", ".join(f"{k:.1f}" for k in kn),
                     "Non-linearity χ²": round(chi2, 2), "df": df_, "p": p, "p (text)": fmt_p(p),
                     "Non-linear (p<0.05)": "yes" if p < 0.05 else "no"})
    return rows


LINEARITY: list[dict] = []


# ---------------------------------------------------------------------------
# Multiple imputation by chained equations (MICE), m imputations
# ---------------------------------------------------------------------------
def _dummies(frame: pd.DataFrame) -> np.ndarray:
    parts = []
    for c in frame.columns:
        s = frame[c]
        if s.dtype == object or str(s.dtype).startswith(("category", "str", "string")):
            dm = pd.get_dummies(s.astype(str), prefix=c, drop_first=True, dtype=float)
            parts.append(dm)
        else:
            parts.append(s.astype(float).to_frame())
    M = pd.concat(parts, axis=1)
    M.insert(0, "const", 1.0)
    return M.values


def mice(df: pd.DataFrame, impute: dict, predictors: list[str], m=20, iters=10, seed=20260928,
         pmm_k=5):
    """impute: {var: 'pmm' | 'logreg' | 'mlogit'}. predictors: fully observed
    variables used as predictors (must include the outcome). Returns a list of
    m completed copies of df (only `impute` columns change)."""
    rng = np.random.default_rng(seed)
    base = df.copy()
    miss = {v: base[v].isna().values for v in impute}
    out = []
    for _ in range(m):
        cur = base.copy()
        for v in impute:  # initialise by random draws from observed values
            obs = cur.loc[~miss[v], v].values
            cur.loc[miss[v], v] = rng.choice(obs, miss[v].sum())
        for _it in range(iters):
            for v, method in impute.items():
                others = [p for p in predictors if p != v] + [w for w in impute if w != v]
                Xall = _dummies(cur[others])
                mo, mm = ~miss[v], miss[v]
                if mm.sum() == 0:
                    continue
                Xo, Xm = Xall[mo], Xall[mm]
                if method == "pmm":
                    yo = cur.loc[mo, v].astype(float).values
                    XtX_inv = np.linalg.pinv(Xo.T @ Xo)
                    bhat = XtX_inv @ Xo.T @ yo
                    resid = yo - Xo @ bhat
                    dfree = len(yo) - Xo.shape[1]
                    sigma2 = resid @ resid / rng.chisquare(dfree)
                    bstar = rng.multivariate_normal(bhat, sigma2 * XtX_inv)
                    yhat_o = Xo @ bhat
                    yhat_m = Xm @ bstar
                    order = np.argsort(yhat_o)
                    sorted_hat = yhat_o[order]
                    pos = np.searchsorted(sorted_hat, yhat_m)
                    imp = np.empty(mm.sum())
                    for i, p0 in enumerate(pos):
                        lo, hi = max(0, p0 - pmm_k), min(len(sorted_hat), p0 + pmm_k)
                        cand = order[lo:hi]
                        dist = np.abs(yhat_o[cand] - yhat_m[i])
                        donors = cand[np.argsort(dist)[:pmm_k]]
                        imp[i] = yo[rng.choice(donors)]
                    cur.loc[mm, v] = imp
                elif method == "mlogit":
                    cats = sorted(cur.loc[mo, v].unique(), key=str)
                    yo = pd.Categorical(cur.loc[mo, v], categories=cats).codes
                    mod = sm.MNLogit(yo, Xo)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        r = mod.fit(method="newton", maxiter=50, disp=0)
                    P = np.asarray(r.params)  # (p, K-1), equations in columns
                    covp = np.asarray(r.cov_params())  # ordered equation by equation
                    draw = rng.multivariate_normal(P.T.ravel(), covp).reshape(P.shape[1], P.shape[0]).T
                    eta = np.column_stack([np.zeros(len(Xm)), Xm @ draw])
                    eta -= eta.max(axis=1, keepdims=True)
                    pr = np.exp(eta)
                    pr /= pr.sum(axis=1, keepdims=True)
                    u = rng.random(len(Xm))[:, None]
                    idx = (pr.cumsum(axis=1) < u).sum(axis=1)
                    cur.loc[mm, v] = np.array(cats, dtype=object)[np.minimum(idx, len(cats) - 1)]
                elif method == "logreg":
                    yo = cur.loc[mo, v].astype(float).values
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        r = sm.GLM(yo, Xo, family=sm.families.Binomial()).fit()
                    bstar = rng.multivariate_normal(r.params, r.cov_params())
                    pr = 1 / (1 + np.exp(-(Xm @ bstar)))
                    cur.loc[mm, v] = (rng.random(len(pr)) < pr).astype(float)
        out.append(cur)
    return out


def pool_fits(fits, cols):
    """Rubin's rules for the listed coefficient names."""
    m = len(fits)
    B = np.array([[f.params[c] for c in cols] for f in fits])
    W = np.array([[f.bse[c] ** 2 for c in cols] for f in fits])
    qbar = B.mean(0)
    ubar = W.mean(0)
    b = B.var(0, ddof=1)
    T = ubar + (1 + 1 / m) * b
    r = (1 + 1 / m) * b / np.where(ubar > 0, ubar, np.nan)
    dfree = (m - 1) * (1 + 1 / np.where(r > 0, r, np.nan)) ** 2
    rows = []
    for i, c in enumerate(cols):
        se = math.sqrt(T[i])
        tq = stats.t.ppf(0.975, dfree[i]) if np.isfinite(dfree[i]) else stats.norm.ppf(0.975)
        o, lo, hi = np.exp(qbar[i]), np.exp(qbar[i] - tq * se), np.exp(qbar[i] + tq * se)
        tstat = qbar[i] / se
        p = 2 * (stats.t.sf(abs(tstat), dfree[i]) if np.isfinite(dfree[i]) else stats.norm.sf(abs(tstat)))
        fmi = (r[i] + 2 / (dfree[i] + 3)) / (r[i] + 1) if np.isfinite(dfree[i]) else 0.0
        rows.append({"coef": c, "OR": o, "CI low": lo, "CI high": hi, "OR (95% CI)": fmt_or(o, lo, hi),
                     "p": p, "p (text)": fmt_p(p), "FMI": round(fmi, 3)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
def table1(df, group, groups_order, spec, overall=True):
    """spec: list of (var, label, kind) with kind 'cont' | 'cat' | 'bin'.
    Continuous: median (IQR) + mean (SD); categorical: n (col %). Missing/unknown
    shown as its own row. p: Kruskal-Wallis / chi-square (Fisher for 2x2 with
    expected <5) on observed levels."""
    cols = (["Overall"] if overall else []) + list(groups_order)
    subsets = {"Overall": df} if overall else {}
    for g in groups_order:
        subsets[g] = df[df[group] == g]
    rows = [{"Characteristic": "N", "Level": "", **{c: f"{len(subsets[c]):,}" for c in cols}, "p": ""}]
    for var, label, kind in spec:
        s_all = df[var]
        gd = [subsets[g][var] for g in groups_order]
        if kind == "cont":
            row = {"Characteristic": label, "Level": "median (IQR)"}
            row2 = {"Characteristic": "", "Level": "mean (SD)"}
            for c in cols:
                s = subsets[c][var].dropna()
                row[c] = f"{s.median():.1f} ({s.quantile(.25):.1f}–{s.quantile(.75):.1f})" if len(s) else ""
                row2[c] = f"{s.mean():.1f} ({s.std():.1f})" if len(s) else ""
            valid = [x.dropna() for x in gd if x.notna().sum() > 0]
            try:
                p = stats.kruskal(*valid).pvalue if len(valid) > 1 else np.nan
            except ValueError:
                p = np.nan
            row["p"] = fmt_p(p)
            row2["p"] = ""
            rows += [row, row2]
            nmiss = {c: int(subsets[c][var].isna().sum()) for c in cols}
            if sum(nmiss.values()):
                rows.append({"Characteristic": "", "Level": "missing",
                             **{c: fmt_n_pct(nmiss[c], len(subsets[c])) for c in cols}, "p": ""})
        else:
            levels = [lv for lv in pd.Series(s_all.dropna().unique()).sort_values(key=lambda x: x.astype(str))]
            if kind == "bin":
                levels = [1.0] if 1.0 in levels or 1 in levels else levels
            ct = pd.crosstab(df[df[group].isin(groups_order)][var], df[df[group].isin(groups_order)][group])
            p = np.nan
            if ct.shape[0] > 1 and ct.shape[1] > 1:
                try:
                    chi2, p, dof, exp = stats.chi2_contingency(ct)
                    if ct.shape == (2, 2) and (exp < 5).any():
                        p = stats.fisher_exact(ct.values)[1]
                except ValueError:
                    p = np.nan
            first = True
            for lv in levels:
                row = {"Characteristic": label if first else "",
                       "Level": ("yes" if kind == "bin" else str(lv))}
                for c in cols:
                    s = subsets[c][var]
                    row[c] = fmt_n_pct(int((s == lv).sum()), int(s.notna().sum()))
                row["p"] = fmt_p(p) if first else ""
                first = False
                rows.append(row)
            nmiss = {c: int(subsets[c][var].isna().sum()) for c in cols}
            if sum(nmiss.values()):
                rows.append({"Characteristic": "", "Level": "missing/unknown (not in %)",
                             **{c: f"{nmiss[c]:,}" for c in cols}, "p": ""})
    return pd.DataFrame(rows)


def bh(pvals):
    p = np.asarray(pvals, float)
    n = np.sum(~np.isnan(p))
    out = np.full_like(p, np.nan)
    idx = np.where(~np.isnan(p))[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * n / (np.arange(1, n + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1)
    return out


# ---------------------------------------------------------------- US spelling (manuscript and figures)
_US_RULES = [
    (r"([Aa])naem", r"\1nem"), (r"([Hh])aem", r"\1em"), (r"([Ii])schaem", r"\1schem"), (r"aemia", "emia"),
    (r"([Gg])ynaec", r"\1ynec"), (r"([Pp])aediatr", r"\1ediatr"), (r"([Oo])estrogen", r"\1strogen".replace("\\1", "e")),
    (r"([Oo])edema", "edema"), (r"([Oo])esophag", "esophag"), (r"([Dd])iarrhoea", r"\1iarrhea"),
    (r"\b([Cc])entre", r"\1enter"), (r"\b([Ff])ibre", r"\1iber"), (r"([Tt])umour", r"\1umor"),
    (r"([Cc])olour", r"\1olor"), (r"([Bb])ehaviour", r"\1ehavior"), (r"([Ff])avour", r"\1avor"),
    (r"([Mm])odell(ed|ing)", r"\1odel\2"), (r"([Ll])abell(ed|ing)", r"\1abel\2"),
    (r"([Aa])nalys(e|ed|es|ing)\b", r"\1nalyz\2"),
    (r"\b(characteri|standardi|categori|summari|recogni|minimi|maximi|organi|randomi|hospitali|utili|priori|"
     r"stabili|generali|normali|optimi|visuali|emphasi|dichotomi|categori|harmoni|finali|real|stratifi)s(e|ed|es|ing|ation|ations)\b",
     r"\1z\2"),
    (r"\bwhilst\b", "while"), (r"\bamongst\b", "among"), (r"([Pp])rogramme", r"\1rogram"),
]


def us_spelling(text):
    """Convert British to US spelling for publication text (does not touch numbers or variable names)."""
    import re as _re
    if not isinstance(text, str):
        return text
    for pat, rep in _US_RULES:
        text = _re.sub(pat, rep, text)
    return text
