"""Figures: flow diagram, Hb spline, forest plots (300-dpi PNG + SVG)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from .utils import FIG_DIR, RESULTS  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
    "ytick.color": "#52514e", "svg.fonttype": "none",
})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # validated categorical slots 1-3 (all-pairs safe)


def save(fig, name):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [f"figures/{name}.png", f"figures/{name}.svg"]


def box(ax, x, y, w, h, text, fs=8.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.01,rounding_size=0.01",
                                fc="white", ec=MUTED, lw=0.9))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=INK, fontweight="bold" if bold else "normal",
            wrap=True)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9))


def flow_diagram():
    f = RESULTS["flow"]
    n_all, n_el = RESULTS["n_total"], RESULTS["n_eligible"]
    g = RESULTS["p1_groups"]
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0.2, 1)
    ax.axis("off")
    box(ax, 0.35, 0.94, 0.5, 0.07, f"Patients with benign uterine pathology\nin analysis_master (n = {n_all:,})", bold=True)
    excl = (f"Excluded (n = {RESULTS['n_excluded']:,})\n"
            f"Active malignancy: {f.get(1, 0):,}\nAge >60 at index: {f.get(6, 0):,}\n"
            f"Age <18 at index: {f.get(3, 0):,}\nAge not derivable: {f.get(5, 0):,}")
    box(ax, 0.8, 0.82, 0.36, 0.13, excl, fs=8)
    arrow(ax, 0.35, 0.905, 0.35, 0.745)
    ax.plot([0.35, 0.62], [0.82, 0.82], color=MUTED, lw=0.9)
    box(ax, 0.35, 0.71, 0.5, 0.07, f"Eligible: age 18–60 at index, no active malignancy\n(n = {n_el:,})", bold=True)
    # three papers
    arrow(ax, 0.35, 0.675, 0.17, 0.585)
    arrow(ax, 0.35, 0.675, 0.5, 0.585)
    arrow(ax, 0.35, 0.675, 0.83, 0.585)
    box(ax, 0.17, 0.54, 0.3, 0.085, f"Paper 1: brain imaging\nreviewed (n = {RESULTS['p1_n_imaged']:,})", fs=8)
    box(ax, 0.5, 0.54, 0.3, 0.085, f"Paper 2: whole cohort\n(n = {n_el:,}; stroke {RESULTS['p2_events']:,})", fs=8)
    box(ax, 0.83, 0.54, 0.3, 0.085, f"Paper 3: whole cohort\n(n = {n_el:,})", fs=8)
    y = 0.385
    p1 = (f"Imaged, no infarct: {g.get('Imaged, no infarct', 0):,}\nCovert infarct: {g.get('Covert infarct', 0):,}\n"
          f"Clinical (coded) stroke: {g.get('Clinical stroke', 0):,}\n"
          f"Imaging-only acute/\nsymptomatic/haemorrhagic: {g.get('Imaging-only acute/symptomatic/haemorrhagic', 0):,}")
    box(ax, 0.17, y, 0.3, 0.15, p1, fs=7.5)
    arrow(ax, 0.17, 0.497, 0.17, 0.46)
    hb = RESULTS["missing"]["anemia_cat"][0]
    p2 = (f"Hb within ±3 y: {n_el - hb:,}\nMCV: {n_el - RESULTS['missing']['mcv_cat'][0]:,}\n"
          f"Platelets: {n_el - RESULTS['missing']['plt_cat'][0]:,}\nFerritin: {RESULTS['p2_ferritin_n']:,}\n"
          f"Incident stroke (after index): {RESULTS['p2_incident_events']:,}")
    box(ax, 0.5, y, 0.3, 0.15, p2, fs=7.5)
    arrow(ax, 0.5, 0.497, 0.5, 0.46)
    gn = RESULTS["p3_group_n"]
    p3 = "\n".join(f"{k}: {gn.get(k, 0):,}" for k in ["Fibroids only", "Adenomyosis only", "Endometriosis only",
                                                      ">1 condition"])
    box(ax, 0.83, y, 0.3, 0.15, p3, fs=7.5)
    arrow(ax, 0.83, 0.497, 0.83, 0.46)
    ax.text(0.0, 0.27, "Complete-case numbers per model are given in the result tables.", fontsize=7.5, color=MUTED)
    return save(fig, "fig1_flow_diagram")


def hb_spline():
    c = RESULTS["p2_curve"]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.fill_between(c["Hb (g/dL)"], c["CI low"], c["CI high"], color=SERIES[0], alpha=0.18, lw=0)
    ax.plot(c["Hb (g/dL)"], c["OR"], color=SERIES[0], lw=2)
    ax.axhline(1, color=MUTED, lw=0.8, ls="--")
    ax.axvline(13, color=GRID, lw=0.8)
    ax.set_yscale("log")
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3, 4])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Haemoglobin closest to index (g/dL)")
    ax.set_ylabel("Adjusted odds ratio for any stroke\n(reference 13 g/dL, log scale)")
    kn = RESULTS["p2_spline_knots"]
    ax.plot(kn, [0.02] * len(kn), marker="|", ls="none", color=MUTED, ms=8, transform=ax.get_xaxis_transform())
    ax.text(kn[-1], 0.05, " knots", transform=ax.get_xaxis_transform(), fontsize=7, color=MUTED)
    ax.set_title(f"Restricted cubic spline, 4 knots; n = {RESULTS['p2_spline_n']:,}, "
                 f"events = {RESULTS['p2_spline_events']:,}", fontsize=8.5, color=MUTED, loc="left")
    ax.grid(axis="y", color=GRID, lw=0.5)
    return save(fig, "fig2_hb_spline")


def forest(df, row_label_cols, group_col, groups, title, name, xlabel="Adjusted odds ratio (95% CI, log scale)",
           width=6.8):
    rows = df[row_label_cols].drop_duplicates().values.tolist()
    n = len(rows)
    fig, ax = plt.subplots(figsize=(width, 0.32 * n * max(1, len(groups) * 0.55) + 1.2))
    off = np.linspace(0.25, -0.25, len(groups)) if len(groups) > 1 else [0]
    ylabels = []
    for i, r in enumerate(rows):
        y = n - i
        ylabels.append((y, " · ".join(str(v) for v in r)))
        for j, g in enumerate(groups):
            q = df
            for col, v in zip(row_label_cols, r):
                q = q[q[col] == v]
            q = q[q[group_col] == g] if group_col else q
            if q.empty or pd.isna(q["OR"].iloc[0]):
                continue
            o, lo, hi = q["OR"].iloc[0], q["lo"].iloc[0], q["hi"].iloc[0]
            ax.plot([lo, hi], [y + off[j]] * 2, color=SERIES[j], lw=1.6, solid_capstyle="round")
            ax.plot([o], [y + off[j]], marker="o", ms=5, color=SERIES[j], mec="white", mew=0.8,
                    label=g if i == 0 or g not in [h.get_label() for h in ax.lines[:-1]] else None)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    xs = np.concatenate([np.asarray(l.get_xdata(), float) for l in ax.lines])
    xs = xs[np.isfinite(xs) & (xs > 0)]
    lo_all, hi_all = xs.min(), xs.max()
    cand = [0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10]
    ticks = [t for t in cand if lo_all * 0.9 <= t <= hi_all * 1.1] or [1]
    if 1 not in ticks:
        ticks.append(1)
    ax.set_xticks(sorted(ticks))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([y for y, _ in ylabels])
    ax.set_yticklabels([lab for _, lab in ylabels])
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    if group_col and len(groups) > 1:
        handles, labels = ax.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        ax.legend(seen.values(), seen.keys(), frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12 - 0.3 / n),
                  ncol=len(groups), fontsize=8)
    return save(fig, name)


def forest_p1():
    comp = RESULTS["p1_comp"].dropna(subset=["OR_imaged"]).copy()
    lab = {"age_index": "Age (per year)", "htn": "Hypertension", "dm": "Diabetes", "dyslipidemia": "Dyslipidaemia",
           "smoking3=Ever": "Smoking: ever", "smoking3=Unknown": "Smoking: unknown",
           "migraine_any": "Migraine (any type)",
           "race3=Black": "Race: Black", "race4=Black": "Race: Black", "race4=Asian": "Race: Asian",
           "race4=Other/unknown": "Race: other/unknown",
           "dxgrp=Adenomyosis only": "Adenomyosis only (vs fibroids only)",
           "dxgrp=Endometriosis only": "Endometriosis only (vs fibroids only)",
           "dxgrp=>1 condition": ">1 condition (vs fibroids only)",
           "anemia_cat=Mild (10–11.9)": "Anaemia mild", "anemia_cat=Moderate (8–9.9)": "Anaemia moderate",
           "anemia_cat=Severe (<8)": "Anaemia severe", "race3=Asian/Other/unknown": "Race: Asian/other/unknown",
           "bmi": "BMI (per kg/m²)", "afib": "Atrial fibrillation",
           "dxgrp3=Endometriosis only": "Endometriosis only (vs fibroids only)",
           "dxgrp3=Adenomyosis only or >1": "Adenomyosis only or >1 (vs fibroids only)",
           "anemia3=Mild (10–11.9)": "Anaemia mild", "anemia3=Moderate/severe (<10)": "Anaemia moderate/severe"}
    rows = []
    fits_m, fits_w = RESULTS["p1_main_fits"], RESULTS["p1_whole_fits"]
    for _, r in comp.iterrows():
        if r["Term"] not in lab:
            continue
        for g, col in [("Comparator: imaged, no infarct (primary)", "imaged"),
                       ("Comparator: whole eligible cohort without stroke", "whole")]:
            o, lo, hi = r[f"OR_{col}"], r[f"lo_{col}"], r[f"hi_{col}"]
            rows.append({"Factor": lab[r["Term"]], "Comparator": g, "OR": o, "lo": lo, "hi": hi})
    df = pd.DataFrame(rows)
    return forest(df, ["Factor"], "Comparator", list(dict.fromkeys(df["Comparator"])),
                  "Factors associated with covert brain infarct: effect of the comparator group", "fig3_forest_paper1")


def forest_p2():
    f = RESULTS["p2_forest"]
    f = f[f["Outcome"].isin(["Any stroke (primary)", "Ischaemic stroke vs no stroke"])].copy()
    f["Row"] = f["Exposure"] + ": " + f["Level"]
    f = f[~f["Level"].str.contains("global")]
    return forest(f, ["Row"], "Outcome", ["Any stroke (primary)", "Ischaemic stroke vs no stroke"],
                  "Haematological indices and stroke (each exposure modelled separately, fully adjusted)",
                  "fig4_forest_paper2")


def forest_p3():
    f = RESULTS["p3_forest"]
    return forest(f, ["Outcome"], "Group", ["Adenomyosis only", "Endometriosis only", ">1 condition"],
                  "Vascular risk factors by uterine condition (reference: fibroids only; age, race, BMI adjusted)",
                  "fig5_forest_paper3")


def forest_surgery():
    t = RESULTS["surg_tab"].dropna(subset=["RR"]).copy()
    t["Row"] = t["Analysis"].str.slice(0, 2).map({"A.": "Any fibroid procedure vs none", "B.": "", "S.": None,
                                                  "C.": "Myomectomy vs hysterectomy (from surgery)"})
    t = t[t["Row"].notna()]
    t.loc[t["Row"] == "", "Row"] = t["Group"] + " vs no procedure"
    t = t.rename(columns={"RR": "OR", "CI low": "lo", "CI high": "hi"})
    return forest(t, ["Row"], "Outcome", ["Any stroke", "Ischaemic stroke"],
                  "Fibroid procedures and subsequent stroke (women with fibroids; time-varying exposure)",
                  "fig6_forest_surgery", xlabel="Adjusted rate ratio (95% CI, log scale)")

# ---------------------------------------------------------------------------
# Key-finding figures for the two recommended papers
# ---------------------------------------------------------------------------
ANAEMIA_LEVELS = ["Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"]


def fig_p2_gradient():
    """Paper 2: adjusted OR by anaemia grade, for any, ischaemic and incident stroke."""
    from .utils import get_or
    fits = RESULTS["p2_fits"]
    outs = [("stroke_any", "Any stroke"), ("y_isch", "Ischaemic stroke"), ("y_incident", "Incident stroke (after diagnosis)")]
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    x0 = np.arange(len(ANAEMIA_LEVELS) + 1)
    for j, (o, lab) in enumerate(outs):
        f = fits[(o, "anemia_cat")]
        pts = [(1.0, 1.0, 1.0)] + [(get_or(f, f"anemia_cat={lv}")["OR"], get_or(f, f"anemia_cat={lv}")["lo"],
                                    get_or(f, f"anemia_cat={lv}")["hi"]) for lv in ANAEMIA_LEVELS]
        xs = x0 + (j - 1) * 0.18
        o_, lo_, hi_ = map(np.array, zip(*pts))
        ax.vlines(xs[1:], lo_[1:], hi_[1:], color=SERIES[j], lw=1.8)
        ax.plot(xs, o_, color=SERIES[j], lw=1.2, alpha=0.6)
        ax.plot(xs, o_, "o", color=SERIES[j], ms=6, mec="white", mew=0.8,
                label=f"{lab} (p-trend {_ptrend(o)})")
    ax.axhline(1, color=MUTED, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_yticks([0.75, 1, 1.5, 2, 3])
    ax.set_ylim(0.65, 3.8)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(x0)
    ax.set_xticklabels(["None\n(Hb ≥12)", "Mild\n(10–11.9)", "Moderate\n(8–9.9)", "Severe\n(<8)"])
    ax.set_xlabel("Anaemia grade (Hb g/dL)")
    ax.set_ylabel("Adjusted odds ratio (95% CI, log scale)")
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("Stroke odds rise with anaemia severity", fontsize=9, loc="left", color=INK)
    return save(fig, "fig7_p2_anaemia_gradient")


def _ptrend(o):
    from .utils import fmt_p
    t = fmt_p(RESULTS[f"p2_trend_{o}"]["p"])
    return t if t.startswith("<") else f"= {t}"


def fig_p2_robustness():
    """Paper 2: moderate and severe anaemia ORs across primary and sensitivity analyses."""
    from .utils import get_or
    fits = RESULTS["p2_fits"]
    lt = RESULTS["p2_lab_timing"]
    lt = lt[lt.Outcome == "Any stroke (primary)"].set_index("Level")
    mi = RESULTS["p2_mi"]
    mi = mi[mi.Outcome == "Any stroke (primary)"].set_index("coef")
    rows = []
    for lv in ["Moderate (8–9.9)", "Severe (<8)"]:
        c = f"anemia_cat={lv}"
        for lab, (o, lo, hi) in [
            ("Primary (Hb within ±3 y)", _t(get_or(fits[("stroke_any", "anemia_cat")], c))),
            ("Hb within ±1 y of diagnosis", (lt.loc[lv, "OR ±1y"], lt.loc[lv, "CI low ±1y"], lt.loc[lv, "CI high ±1y"])),
            ("Multiple imputation (m = 20)", (mi.loc[c, "OR"], mi.loc[c, "CI low"], mi.loc[c, "CI high"])),
            ("Joint model (+ MCV, platelets)", _t(get_or(fits[("stroke_any", "joint")], c))),
            ("Ischaemic stroke only", _t(get_or(fits[("y_isch", "anemia_cat")], c))),
            ("Incident stroke only", _t(get_or(fits[("y_incident", "anemia_cat")], c))),
        ]:
            rows.append({"Analysis": lab, "Grade": lv.split(" ")[0] + " anaemia", "OR": o, "lo": lo, "hi": hi})
    df = pd.DataFrame(rows)
    return forest(df, ["Analysis"], "Grade", ["Moderate anaemia", "Severe anaemia"],
                  "Anaemia–stroke association across sensitivity analyses (reference: Hb ≥12 g/dL)",
                  "fig8_p2_robustness", width=6.4)


def _t(g):
    return g["OR"], g["lo"], g["hi"]


def fig_p3_migraine_prev():
    """Paper 3: age-standardised migraine prevalence by uterine condition."""
    pv = RESULTS["p3_prev"]
    pv = pv[pv.Outcome == "Migraine (any)"].set_index("Group")
    groups = ["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"]
    m3 = RESULTS["p3_main"]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    y = pv.loc[groups, "Age-std %"].values
    lo = pv.loc[groups, "Age-std low"].values
    hi = pv.loc[groups, "Age-std high"].values
    cols = [MUTED] + SERIES
    ax.bar(range(4), y, color=cols, width=0.6, edgecolor="white", linewidth=2)
    ax.errorbar(range(4), y, yerr=[y - lo, hi - y], fmt="none", ecolor=INK, capsize=3, lw=0.9)
    for i, g in enumerate(groups):
        n = int(pv.loc[g, "N (outcome known)"])
        lab = f"{y[i]:.1f}%"
        if g != "Fibroids only":
            r = m3[(m3.Outcome == "Migraine (any)") & (m3["Group vs fibroids only"] == g)].iloc[0]
            lab += f"\nOR {r['OR']:.2f}"
        ax.text(i, hi[i] + 0.8, lab, ha="center", va="bottom", fontsize=8, color=INK)
        ax.text(i, 0.8, f"n={n:,}", ha="center", va="bottom", fontsize=7, color="white")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Fibroids\nonly", "Adenomyosis\nonly", "Endometriosis\nonly", ">1\ncondition"])
    ax.set_ylabel("Age-standardised migraine prevalence, % (95% CI)")
    ax.set_ylim(0, max(hi) + 7)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title("Migraine by uterine condition (OR vs fibroids only; age, race, BMI adjusted)", fontsize=9,
                 loc="left", color=INK)
    return save(fig, "fig9_p3_migraine_prevalence")


def fig_p3_migraine_age():
    """Paper 3: migraine ORs overall and by age stratum."""
    m3 = RESULTS["p3_main"]
    st = RESULTS["p3_strata"]
    rows = []
    for g in ["Adenomyosis only", "Endometriosis only", ">1 condition"]:
        r = m3[(m3.Outcome == "Migraine (any)") & (m3["Group vs fibroids only"] == g)].iloc[0]
        rows.append({"Stratum": "All ages", "Group": g, "OR": r["OR"], "lo": r["CI low"], "hi": r["CI high"]})
        for ag in ["18–39", "40–60"]:
            r = st[(st.Outcome == "Migraine (any)") & (st["Group vs fibroids only"] == g) & (st.Age == ag)].iloc[0]
            rows.append({"Stratum": f"Age {ag}", "Group": g, "OR": r["OR"], "lo": r["CI low"], "hi": r["CI high"]})
    df = pd.DataFrame(rows)
    return forest(df, ["Stratum"], "Group", ["Adenomyosis only", "Endometriosis only", ">1 condition"],
                  "Migraine vs fibroids only, overall and by age (age, race, BMI adjusted)",
                  "fig10_p3_migraine_by_age", width=6.0)


P2X_SHORT = {
    "Primary: Any stroke": "Primary model",
    "Hb ≥30 d before stroke: Any stroke": "Hb measured ≥30 d before stroke",
    "Excluding Hb within 30 d of stroke: Any stroke": "Excluding Hb within 30 d of stroke",
    "Primary: Incident stroke": "Incident strokes only",
    "2a Pre-Hb conditions: Any stroke": "+ conditions documented before Hb (2a)",
    "2b Over-adjustment check (+ undated/post-stroke): Any stroke": "+ undated / post-stroke factors (2b)",
    "Hb ≥30 d before stroke + 2a: Any stroke": "Hb before stroke + 2a",
    "Excluding haemoglobinopathies: Any stroke": "Excluding haemoglobinopathies",
    "Excluding all anaemia-causing conditions: Any stroke": "Excluding all anaemia-causing conditions",
    "Excluding deaths within 1 y of Hb (+2a): Any stroke": "Excluding deaths within 1 y of Hb (+2a)",
    "Excluding cancer, heart failure, CKD, liver disease, HIV (+2a): Any stroke": "Excluding serious chronic illness (+2a)",
    "Primary: Ischaemic stroke": "Ischaemic stroke only",
}


def fig_p2x():
    t = RESULTS["p2x_tab"].set_index("Analysis")
    rows = []
    for a, lab in P2X_SHORT.items():
        if a not in t.index:
            continue
        r = t.loc[a]
        rows.append({"Analysis": lab, "Estimate": "Moderate anaemia (8–9.9 g/dL)", "OR": r["Moderate OR"],
                     "lo": r["Moderate lo"], "hi": r["Moderate hi"]})
        rows.append({"Analysis": lab, "Estimate": "Per anaemia grade (trend)", "OR": r["trend OR"],
                     "lo": r["trend lo"], "hi": r["trend hi"]})
    tt = RESULTS["p2x_tte"]
    tt = tt[tt.Outcome == "Any stroke"].set_index("Anaemia grade")
    for key, est in [("Moderate (8–9.9)", "Moderate anaemia (8–9.9 g/dL)"), ("Per grade (trend)", "Per anaemia grade (trend)")]:
        rows.append({"Analysis": "Time-to-event from Hb date (rate ratio, 2a)", "Estimate": est,
                     "OR": tt.loc[key, "RR"], "lo": tt.loc[key, "CI low"], "hi": tt.loc[key, "CI high"]})
    df = pd.DataFrame(rows)
    return forest(df, ["Analysis"], "Estimate", ["Moderate anaemia (8–9.9 g/dL)", "Per anaemia grade (trend)"],
                  "Anaemia and any stroke: temporality, confounding and restriction analyses (OR; last row RR)",
                  "fig11_p2_extended", width=6.6)


def fig_p2_per_hb():
    """Absolute stroke rates by anaemia grade (A) and risk per 1 g/dL lower Hb (B)."""
    ab = RESULTS["p2x_abs"]
    ab = ab[ab.Outcome == "Any stroke"].reset_index(drop=True)
    per = RESULTS["p2x_per"]
    per = per[per.Adjustment == "+ pre-Hb conditions (2a)"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 3.7), gridspec_kw={"width_ratios": [1, 1.25]})
    ramp = ["#c9dcf4", "#8fb6ea", "#4f8fdd", "#1d5fb0"]  # single-hue sequential (blue), light -> dark
    x = np.arange(len(ab))
    y, lo, hi = ab["Rate /1,000 PY"].values, ab["Rate low"].values, ab["Rate high"].values
    a1.bar(x, y, color=ramp, width=0.62, edgecolor="white", linewidth=2)
    a1.errorbar(x, y, yerr=[y - lo, hi - y], fmt="none", ecolor=INK, capsize=3, lw=0.9)
    a1.axhline(y[0], color=MUTED, lw=0.8, ls="--")
    for i, r in ab.iterrows():
        a1.text(i, hi[i] + 0.25, f"{y[i]:.1f}", ha="center", va="bottom", fontsize=8, color=INK)
    a1.set_xticks(x)
    a1.set_xticklabels([f"{g}\n{int(r['Events'])} / {int(r['Women']):,}" for g, (_, r) in
                        zip(["None\n(≥12)", "Mild\n(10–11.9)", "Moderate\n(8–9.9)", "Severe\n(<8)"], ab.iterrows())],
                       fontsize=7)
    a1.set_xlabel("Anaemia grade (Hb, g/dL); strokes / women")
    a1.set_ylabel("Strokes per 1,000 person-years (95% CI)")
    a1.set_ylim(0, max(hi) * 1.12)
    a1.grid(axis="y", color=GRID, lw=0.5)
    a1.set_axisbelow(True)
    a1.set_title("A  Stroke rate after the Hb measurement", fontsize=9, loc="left", color=INK)
    rows = [("Cross-sectional (odds ratio)", "Any stroke", "Odds ratio\nany stroke"),
            ("Cross-sectional (odds ratio)", "Ischaemic stroke", "Odds ratio\nischaemic stroke"),
            ("After the Hb measurement (rate ratio)", "Any stroke", "Rate ratio after Hb\nany stroke"),
            ("After the Hb measurement (rate ratio)", "Ischaemic stroke", "Rate ratio after Hb\nischaemic stroke")]
    series = [("Below 13 g/dL (per 1 g/dL lower)", "Below 13 g/dL", SERIES[0], 0.14),
              ("Whole range", "Whole Hb range", SERIES[1], -0.14)]
    for j, (des, oc, lab) in enumerate(rows):
        yy = len(rows) - j
        for rng, slab, col, off in series:
            r = per[(per.Design == des) & (per.Outcome == oc) & (per["Hb range"] == rng)].iloc[0]
            a2.plot([r["lo"], r["hi"]], [yy + off] * 2, color=col, lw=1.8, solid_capstyle="round")
            a2.plot([r["Est"]], [yy + off], "o", color=col, ms=5.5, mec="white", mew=0.8, label=slab if j == 0 else None)
            a2.text(1.265, yy + off, f"{r['Est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})", va="center", fontsize=7,
                    color=col)
    a2.axvline(1, color=MUTED, lw=0.8, ls="--")
    a2.set_xlim(0.96, 1.40)
    a2.set_yticks([len(rows) - j for j in range(len(rows))])
    a2.set_yticklabels([r[2] for r in rows], fontsize=8)
    a2.set_xlabel("Per 1 g/dL lower Hb (95% CI)")
    a2.grid(axis="x", color=GRID, lw=0.5)
    a2.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0.0, -0.36), ncol=2)
    a2.set_title("B  Risk per 1 g/dL lower haemoglobin", fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    return save(fig, "fig13_p2_per_hb")


def fig_p2_hb_curves():
    """Stroke by Hb: cross-sectional OR (A) and rate ratio after the Hb measurement (B), both adjusted (2a)."""
    c, info = RESULTS["p2x_curves"], RESULTS["p2x_curve_info"]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.5), sharey=True)
    panels = [("Cross-sectional (odds ratio)", "A  Odds of stroke", "Adjusted odds ratio"),
              ("After the Hb measurement (rate ratio)", "B  Stroke rate after the Hb measurement", "Adjusted rate ratio")]
    for ax, (des, title, yl), col in zip(axes, panels, SERIES[:2]):
        cc = c[c.Design == des]
        ax.fill_between(cc["Hb (g/dL)"], cc["lo"], cc["hi"], color=col, alpha=0.18, lw=0)
        ax.plot(cc["Hb (g/dL)"], cc["Est"], color=col, lw=2)
        ax.axhline(1, color=MUTED, lw=0.8, ls="--")
        ax.axvline(13, color=GRID, lw=0.8)
        ax.set_yscale("log")
        ax.set_ylim(0.4, 4.5)
        ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3, 4])
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_xlabel("Haemoglobin (g/dL)")
        ax.set_title(title, fontsize=9, loc="left", color=INK)
        i = info[des]
        ax.text(0.98, 0.97, f"n = {i['n']:,}; strokes = {i['events']:,}\nnon-linearity P = "
                + (f"{i['p_nonlin']:.3f}" if i["p_nonlin"] >= 0.001 else "<0.001"),
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=MUTED)
        if ax is axes[0]:
            ax.set_ylabel("Adjusted ratio vs Hb 13 g/dL (log scale)")
    return save(fig, "fig14_p2_hb_curves")


def fig_p2_flow():
    """Participant flow for the anaemia manuscript."""
    f = RESULTS["flow"]
    n_all, n_el = RESULTS["n_total"], RESULTS["n_eligible"]
    n_nohb = RESULTS["missing"]["anemia_cat"][0]
    fa = RESULTS["p2_fits"][("stroke_any", "anemia_cat")]
    tc = RESULTS["p2x_tte_cohort"]
    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box(ax, 0.38, 0.93, 0.56, 0.08, f"Women with benign uterine pathology\n(fibroids, adenomyosis, endometriosis), n = {n_all:,}", bold=True)
    box(ax, 0.82, 0.78, 0.34, 0.14, f"Excluded, n = {RESULTS['n_excluded']:,}\nActive malignancy: {f.get(1, 0):,}\n"
                                    f"Age >60: {f.get(6, 0):,}\nAge <18: {f.get(3, 0):,}\nAge not derivable: {f.get(5, 0):,}", fs=8)
    arrow(ax, 0.38, 0.89, 0.38, 0.72)
    ax.plot([0.38, 0.65], [0.78, 0.78], color=MUTED, lw=0.9)
    box(ax, 0.38, 0.68, 0.56, 0.07, f"Eligible women aged 18–60, n = {n_el:,}", bold=True)
    box(ax, 0.82, 0.55, 0.34, 0.08, f"No Hb within ±3 years\nof index, n = {n_nohb:,}", fs=8)
    arrow(ax, 0.38, 0.645, 0.38, 0.475)
    ax.plot([0.38, 0.65], [0.55, 0.55], color=MUTED, lw=0.9)
    box(ax, 0.38, 0.44, 0.56, 0.07, f"Hb available, n = {n_el - n_nohb:,}", bold=True)
    arrow(ax, 0.25, 0.405, 0.2, 0.27)
    arrow(ax, 0.51, 0.405, 0.58, 0.27)
    box(ax, 0.2, 0.2, 0.36, 0.14, f"Primary analysis\n(complete covariates)\nn = {fa.n:,}; strokes = {fa.events:,}", fs=8)
    box(ax, 0.62, 0.2, 0.4, 0.14, f"Time-to-event cohort\n(no stroke before Hb/index)\nn = {tc['n']:,}; strokes = {tc['ev']}\n"
                                  f"{tc['py']:,.0f} person-years", fs=8)
    return save(fig, "fig12_p2_flow")


def run():
    out = {}
    out["flow"] = flow_diagram()
    out["spline"] = hb_spline()
    out["p1"] = forest_p1()
    out["p2"] = forest_p2()
    out["p3"] = forest_p3()
    out["surgery"] = forest_surgery()
    out["p2_gradient"] = fig_p2_gradient()
    out["p2_robust"] = fig_p2_robustness()
    out["p3_prev"] = fig_p3_migraine_prev()
    out["p3_age"] = fig_p3_migraine_age()
    out["p2x"] = fig_p2x()
    out["p2_flow"] = fig_p2_flow()
    out["p2_per_hb"] = fig_p2_per_hb()
    out["p2_hb_curves"] = fig_p2_hb_curves()
    RESULTS["figures"] = out
    return out
