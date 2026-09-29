"""Blood count results from the Lab Tests extracts (2: haemoglobin, 3: MCV, 4: MCH/MCHC, 7: RDW, 8: platelets).

Cleaning rules (every exclusion is counted in the LAB_cleaning table):
- analyte assigned from the test description (files 4, 8) or, for file 7 (no description), from
  unit / normal range / code: RDW-CV (%) vs RDW-SD (fL);
- estimates, comments, 'large platelets' and non-CBC rows are not counts and are dropped;
- blank values, sentinel values (>= 9,999,999) and sentinel/implausible dates (year < 1990) are dropped;
- values outside physiological plausibility limits (PLAUSIBLE) are dropped, not truncated;
- repeated results for the same woman, analyte and timestamp are averaged.
For each woman the plausible value nearest in time to her Hb measurement (hgb_date) is taken
("nearest plausible value"); the gap in days is kept so analyses can restrict to the same blood count.
Clinic Number (MRN) is used only in memory for linkage and never written to an output.
"""
import glob

import numpy as np
import pandas as pd

from .utils import DATA_DIR

PLAUSIBLE = {  # physiological limits (values outside are treated as errors)
    "Hb (g/dL)": (3, 20),
    "MCV (fL)": (50, 150),
    "MCH (pg)": (12, 50),
    "MCHC (g/dL)": (20, 45),
    "RDW-CV (%)": (8, 40),
    "RDW-SD (fL)": (20, 150),
    "Platelets (x10^9/L)": (1, 2000),
}
SENTINEL = 9_999_999
# files verified to be exact row-level subsets of another extract (not read twice)
DUPLICATE_FILES = {6: "not used: every row identical to a row in Lab_Tests_7"}


def _files():
    return {int(f.split("(")[-1].split(")")[0]): f for f in glob.glob(str(DATA_DIR / "labs" / "*.csv"))}


HB_EXCL = (r"\bQL\b|\bU\b|URINE|UPHB|O2|OXY|METHEM|A1C|\bA\b|\bA2\b|\bF\b|VARIANT|RETIC|PLASMA|, P$|ARTERIAL|"
           r"VENOUS|GAS|ABG|VBG|POC|TOTAL|THB|COMMENT|SOURCE|\bCBC\b|DONOR")


def _analyte_2(desc):
    """Venous whole-blood (CBC) haemoglobin only: urine, blood-gas/co-oximetry, point-of-care, HbA1c,
    haemoglobin fractions, reticulocyte and plasma haemoglobin are excluded."""
    d = desc.str.upper().str.strip().fillna("")
    hb = d.str.match(r"^(EXTI |EXTM )?(NON-GHS )?(HEMOGLOBIN|HGB|HB)\b") & ~d.str.contains(HB_EXCL)
    return np.where(hb, "Hb (g/dL)", "not a CBC haemoglobin")


def _analyte_3(desc):
    d = desc.str.upper().str.strip().fillna("")
    mcv = d.str.contains(r"\bMCV\b|MEAN CORPUSCULAR VOLUME") & ~d.str.contains(r"PUB'S")
    return np.where(mcv, "MCV (fL)", "not an MCV")


def _analyte_4(desc):
    d = desc.str.upper().fillna("")
    mchc = d.str.contains(r"\bMCHC\b|HGB CONCENTRATION|HEMOGLOBIN CONCENTRATION")
    mch = d.str.contains(r"\bMCH\b|MEAN CORPUSCULAR (?:HEMOGLOBIN|HGB)\)") & ~mchc
    return np.select([mchc, mch], ["MCHC (g/dL)", "MCH (pg)"], "not a red-cell index")


def _analyte_8(desc):
    d = desc.str.upper().fillna("")
    drop = d.str.contains(r"ESTIMATE|\bEST\b|EST\.|LARGE|COMMENT|UNIT INFO")
    plt = d.str.contains(r"PLATELET|\bPLT\b")
    return np.select([drop, plt], ["platelet estimate/comment", "Platelets (x10^9/L)"], "not a platelet count")


def _analyte_7(x):
    unit = x.LAB_UNIT_OF_MEASURE_TXT.fillna("").str.strip().str.lower()
    rng_lo = pd.to_numeric(x.LAB_NORMAL_RANGE_TXT.fillna("").str.extract(r"^\s*([\d.]+)")[0], errors="coerce")
    code_med = x.groupby("Lab Test Code").v.transform(lambda s: s[s < SENTINEL].median())
    sd = unit.str.contains("fl") | ((unit == "") & (rng_lo >= 30)) | ((unit == "") & rng_lo.isna() & (code_med > 30))
    cv = (unit.str.contains("%")) | ((unit == "") & (rng_lo < 30)) | ((unit == "") & rng_lo.isna() & (code_med <= 30))
    return np.select([sd, cv], ["RDW-SD (fL)", "RDW-CV (%)"], "unrecognised unit")


def load():
    """Return (clean long table [mrn, analyte, date, value], cleaning audit)."""
    fs = _files()
    parts, audit = [], []
    for n, fn in sorted(fs.items()):
        if n in DUPLICATE_FILES:
            k = len(pd.read_csv(fn, dtype=str, usecols=["Clinic Number"]))
            audit.append(pd.DataFrame([{"File": f"Lab_Tests_{n}", "analyte": "RDW (duplicate)",
                                        "reason": DUPLICATE_FILES[n], "rows": k}]))
            continue
        x = pd.read_csv(fn, dtype=str)
        x["v"] = pd.to_numeric(x["Value (Numeric)"], errors="coerce")
        if n == 2:
            x["analyte"] = _analyte_2(x["Lab Test Description"])
        elif n == 3:
            x["analyte"] = _analyte_3(x["Lab Test Description"])
        elif n == 4:
            x["analyte"] = _analyte_4(x["Lab Test Description"])
        elif n == 8:
            x["analyte"] = _analyte_8(x["Lab Test Description"])
        elif n == 7:
            x["analyte"] = _analyte_7(x)
        else:
            x["analyte"] = "unknown file"
        x["mrn"] = pd.to_numeric(x["Clinic Number"], errors="coerce")
        x["date"] = pd.to_datetime(x["Result Date"].str.slice(0, 10), errors="coerce")
        lo = x.analyte.map(lambda a: PLAUSIBLE.get(a, (np.nan, np.nan))[0])
        hi = x.analyte.map(lambda a: PLAUSIBLE.get(a, (np.nan, np.nan))[1])
        reason = np.select(
            [~x.analyte.isin(PLAUSIBLE.keys()), x.v.isna(), x.v >= SENTINEL, x.date.isna() | (x.date.dt.year < 1990),
             (x.v < lo) | (x.v > hi), x.mrn.isna()],
            ["not a count / other test", "blank value", "sentinel 9,999,999", "missing or sentinel date",
             "outside plausible range", "no clinic number"], "kept")
        x["reason"] = reason
        a = x.groupby(["analyte", "reason"]).size().rename("rows").reset_index()
        a.insert(0, "File", f"Lab_Tests_{n}")
        audit.append(a)
        parts.append(x.loc[x.reason == "kept", ["mrn", "analyte", "date", "v"]])
    long = pd.concat(parts, ignore_index=True)
    long = long.groupby(["mrn", "analyte", "date"], as_index=False).v.mean()
    audit = pd.concat(audit, ignore_index=True)
    return long, audit


def nearest_to_hb(df, long):
    """Per woman and analyte: nearest plausible value to her Hb date, gap (days, value date - Hb date)."""
    hb = df[["mrn", "hgb_date"]].copy()
    hb["hbd"] = pd.to_datetime(hb.hgb_date, errors="coerce")
    hb = hb.dropna(subset=["hbd", "mrn"])
    m = long.merge(hb[["mrn", "hbd"]], on="mrn", how="inner")
    m["gap"] = (m.date - m.hbd).dt.days
    m["agap"] = m.gap.abs()
    m = m.sort_values(["mrn", "analyte", "agap", "gap"])  # ties: earlier result first
    near = m.drop_duplicates(["mrn", "analyte"])
    wide_v = near.pivot(index="mrn", columns="analyte", values="v")
    wide_g = near.pivot(index="mrn", columns="analyte", values="gap")
    return wide_v, wide_g


def coverage(df, wide_v, wide_g):
    """Coverage of the nearest value relative to the Hb date (eligible women with an Hb value)."""
    e = df[(df.eligible == 1) & df.hgb.notna()]
    rows = []
    for a in PLAUSIBLE:
        if a not in wide_v:
            continue
        g = e.mrn.map(wide_g[a])
        v = e.mrn.map(wide_v[a])
        ag = g.abs()
        rows.append({"Analyte": a, "Women with Hb": len(e), "Any plausible value": int(v.notna().sum()),
                     "Same day as Hb": int((ag == 0).sum()), "Within 7 d": int((ag <= 7).sum()),
                     "Within 30 d": int((ag <= 30).sum()), "Within 365 d": int((ag <= 365).sum()),
                     "Median |gap| (d)": float(ag.median()) if ag.notna().any() else np.nan,
                     "Median value (same day)": float(v[ag == 0].median()) if (ag == 0).any() else np.nan,
                     "P5–P95 (same day)": (f"{v[ag == 0].quantile(.05):.1f}–{v[ag == 0].quantile(.95):.1f}"
                                           if (ag == 0).any() else "")})
    return pd.DataFrame(rows)
