"""Export the originally reported (Python pipeline) estimates from outputs/_state.pkl to a flat CSV,
so the R reproduction compares against machine-read values rather than retyped numbers.
Usage (from reanalysis/): python3 tools/export_original_results.py ../outputs/_state.pkl output/original_reported.csv"""
import pickle
import re
import sys

import pandas as pd
sys.path.insert(0, "..")  # the pickle references classes from the original analysis package

state, out = sys.argv[1], sys.argv[2]
R = pickle.load(open(state, "rb"))["RESULTS"]
T = pickle.load(open(state, "rb"))["TABLES"]
rows = []


def add(analysis, outcome, level, txt, n=None, ev=None):
    m = re.match(r"\s*([\d.]+)\s*\(([\d.]+)[–-]([\d.]+)\)", str(txt))
    if m:
        rows.append(dict(analysis=analysis, outcome=outcome, level=level, est=float(m[1]), lo=float(m[2]),
                         hi=float(m[3]), N=n, events=ev))


mm = T["P2_main_models"]
mm = mm[(mm["Model"] == "separate") & (mm["Term"] == "Anaemia grade")]
for _, r in mm.iterrows():
    add("Cross-sectional logistic, crude", r["Outcome"], r["Level"], r["Crude OR (95% CI)"], r["N"], r["Events"])
    add("Cross-sectional logistic, adjusted", r["Outcome"], r["Level"], r["Adjusted OR (95% CI)"], r["N"], r["Events"])
ss = R["p2x_tab"] if "p2x_tab" in R else T["P2X_sensitivity_summary"]
ss = T["P2X_sensitivity_summary"]
for _, r in ss.iterrows():
    for lv in ["Mild", "Moderate", "Severe", "Per grade (trend)"]:
        add("Sensitivity: " + r["Analysis"], r["Analysis"].split(": ")[-1], lv, r[lv], r["N"], r["Events"])
for _, r in R["p2x_tte"].iterrows():
    add("Original time-to-event, Paper 2 covariates", r["Outcome"], r["Anaemia grade"], str(r["Adjusted (Paper 2 covariates) RR"]).split(";")[0], r["Women"], r["Events"])
    add("Original time-to-event, + 2a", r["Outcome"], r["Anaemia grade"], str(r["Adjusted + pre-Hb conditions (2a) RR"]).split(";")[0], r["Women"], r["Events"])
for _, r in R["p2x_per"].iterrows():
    add(f"Per 1 g/dL: {r['Design']}, {r['Adjustment']}", r["Outcome"], r["Hb range"], r["Estimate (95% CI)"], r["N"], r["Events"])
for _, r in R["p2_spline_pts"].iterrows():
    add("Spline (Paper 2 covariates), OR vs 13 g/dL", "Any stroke", f"Hb {r['Hb (g/dL)']}", r["OR vs 13 g/dL (95% CI)"])
for _, r in R["p2l_tu"].iterrows():
    add("Time-updated: " + r["Model"], r["Outcome"], r["Level"], r["RR (95% CI)"])
c = R["p2x_tte_cohort"]; t = R["p2l_tu_cohort"]
pd.DataFrame(rows).to_csv(out, index=False)
pd.DataFrame([dict(item="TTE women", value=c["n"]), dict(item="TTE person-years", value=c["py"]), dict(item="TTE events", value=c["ev"]),
              dict(item="TU women", value=t["n"]), dict(item="TU person-years", value=t["py"]), dict(item="TU events", value=t["ev"]),
              dict(item="Cross-sectional any-stroke N", value=R["p2_n"]), dict(item="Cross-sectional any-stroke events", value=R["p2_events"])]
             ).to_csv(out.replace(".csv", "_cohorts.csv"), index=False)
print(len(rows), "rows")
