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
           "migraine3=Without aura": "Migraine without aura", "migraine3=With aura": "Migraine with aura",
           "race3=Black": "Race: Black", "race3=Asian/Other/unknown": "Race: Asian/other/unknown",
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


def run():
    out = {}
    out["flow"] = flow_diagram()
    out["spline"] = hb_spline()
    out["p1"] = forest_p1()
    out["p2"] = forest_p2()
    out["p3"] = forest_p3()
    RESULTS["figures"] = out
    return out
