"""Step 0: Codebook summary, flow, missingness and data-quality flags.

Nothing here modifies the data. Problems are counted and flagged in the log.
"""
import re

import numpy as np
import pandas as pd

from .utils import (RESULTS, WORKBOOK, add_table, derive, load_master, log)

# Variables used anywhere in Papers 1-3 (plus the provenance fields they need)
USED_VARS = [
    "age_index", "race", "black_race", "bmi", "fibroids", "adenomyosis", "endometriosis",
    "uterine_dx_group", "uterine_bleeding", "hormonal_tx", "hormonal_type", "htn", "dm",
    "dyslipidemia", "cad", "afib", "smoking", "migraine", "vte_history", "thrombophilia",
    "hgb", "hgb_days", "mcv", "platelets", "lab_date", "ferritin", "ferritin_days", "anemia",
    "iron_deficiency", "stroke_any", "stroke_any_incl_imaging", "stroke_type", "stroke_date",
    "stroke_timing", "stroke_confirmed_imaging", "stroke_etiology", "location",
    "vascular_territory", "num_infarct", "notes", "myomectomy", "surgical_tx",
]

CODED_WITH_9 = ["uterine_dx_group", "uterine_bleeding", "hormonal_type", "htn", "dm", "dyslipidemia",
                "cad", "afib", "smoking", "migraine", "vte_history", "thrombophilia", "stroke_any",
                "stroke_type", "stroke_timing", "stroke_confirmed_imaging", "stroke_etiology",
                "vascular_territory", "num_infarct", "hormonal_tx", "fibroids", "adenomyosis",
                "endometriosis"]


def read_codebook():
    cb = pd.read_excel(WORKBOOK, sheet_name="Codebook", header=None)
    hdr = cb.iloc[0].tolist()
    cb = cb.iloc[1:].copy()
    cb.columns = [str(h) if pd.notna(h) else f"col{i}" for i, h in enumerate(hdr)]
    # rows appended at the bottom of the Codebook are shifted one column left
    # (no '#' number) - realign them
    fixed = []
    for _, r in cb.iterrows():
        vals = r.tolist()
        if isinstance(vals[0], str) and not str(vals[0]).isdigit() and vals[1] not in (None,) \
                and isinstance(vals[1], str) and re.match(r"^[a-z_0-9]+$", str(vals[1])):
            vals = [np.nan] + vals[:-1]
        fixed.append(vals)
    cb = pd.DataFrame(fixed, columns=cb.columns)
    cb = cb[cb["Variable"].notna()]
    return cb


def run():
    raw = load_master()
    RESULTS["n_total"] = len(raw)
    assert raw["study_id"].is_unique
    RESULTS["study_id_unique"] = True
    df = derive(raw)
    elig = df[df["eligible"] == 1].copy()
    RESULTS["n_eligible"] = len(elig)

    # ---------------- Flow ----------------
    reasons = {0: "Not excluded (eligible)", 1: "Gynecologic or other active malignancy",
               2: "Uterine pathology not confirmed", 3: "Age <18 at index",
               4: "Male / no uterus documented", 5: "Insufficient records (age not derivable)",
               6: "Other (age >60 at index)"}
    vc = df["exclusion_reason"].value_counts(dropna=False)
    flow = []
    for k, lab in reasons.items():
        flow.append({"exclusion_reason": k, "Description": lab, "n": int(vc.get(k, 0))})
    if df["exclusion_reason"].isna().any():
        flow.append({"exclusion_reason": "blank", "Description": "blank", "n": int(df["exclusion_reason"].isna().sum())})
    flow = pd.DataFrame(flow)
    flow["% of all records"] = (100 * flow["n"] / len(df)).round(2)
    add_table("flow", flow, "Participant flow from all records in analysis_master.")
    RESULTS["flow"] = {int(r.exclusion_reason) if r.exclusion_reason != "blank" else r.exclusion_reason: int(r.n)
                       for r in flow.itertuples()}
    RESULTS["n_excluded"] = int((df["eligible"] != 1).sum())
    xt = pd.crosstab(df["eligible"], df["exclusion_reason"])
    consistent = bool(((df["eligible"] == 1) == (df["exclusion_reason"] == 0)).all())
    RESULTS["eligible_reason_consistent"] = consistent
    log("Flow", f"{len(df):,} records; eligible==1: {len(elig):,}; excluded {RESULTS['n_excluded']:,}. "
                f"eligible==1 exactly matches exclusion_reason==0: {consistent}.")
    for r in flow.itertuples():
        log("Flow", f"exclusion_reason {r.exclusion_reason} ({r.Description}): {r.n:,}")

    # ---------------- Codebook summary ----------------
    cb = read_codebook()
    rows = []
    for v in USED_VARS:
        c = cb[cb["Variable"] == v]
        s = elig[v] if v in elig else pd.Series(dtype=float)
        n = len(elig)
        blank = int(s.isna().sum())
        n9 = int((s == 9).sum()) if v in CODED_WITH_9 else 0
        if s.dtype == object or str(s.dtype).startswith(("str", "string")):
            typ_obs = "text"
            obs = f"{s.nunique()} distinct values"
        elif s.dropna().nunique() <= 12:
            typ_obs = "coded"
            obs = "; ".join(f"{k:g}: {v2:,}" for k, v2 in s.value_counts().sort_index().items())
        else:
            typ_obs = "numeric"
            obs = (f"median {s.median():.1f} (IQR {s.quantile(.25):.1f}–{s.quantile(.75):.1f}), "
                   f"range {s.min():g}–{s.max():g}")
        rows.append({
            "Variable": v,
            "Codebook type": c["Type"].iloc[0] if len(c) else "(not in Codebook)",
            "Codebook values": c["Permitted values / units"].iloc[0] if len(c) else "",
            "Observed type": typ_obs,
            "Observed distribution (eligible)": obs,
            "Blank n": blank, "Blank %": round(100 * blank / n, 2),
            "Code 9 n": n9 if v in CODED_WITH_9 else "", "Code 9 %": round(100 * n9 / n, 2) if v in CODED_WITH_9 else "",
        })
    cbs = pd.DataFrame(rows)
    add_table("codebook_summary", cbs, f"Variables used in Papers 1-3, eligible population (N={len(elig):,}).")

    # Codebook variables not present in analysis_master
    absent = [v for v in cb["Variable"].astype(str) if v not in raw.columns]
    RESULTS["codebook_absent"] = absent
    log("Data problem", "Codebook variables absent from analysis_master (cannot be used): " + ", ".join(absent))

    # ---------------- Missingness (eligible) ----------------
    miss_rows = []
    model_vars = {
        "age_index": "all", "race4": "all", "bmi": "all", "htn": "all", "dm": "all", "dyslipidemia": "all",
        "smoking3": "all", "migraine3": "all", "thrombophilia": "all", "afib": "all", "cad": "all",
        "vte_history": "all", "dxgrp": "all", "anemia_cat": "P1,P2", "mcv_cat": "P2", "plt_cat": "P2",
        "iron_deficiency": "P2", "hormonal_type_cat": "P2", "uterine_bleeding": "P2", "stroke_any": "all",
        "smoking_ever": "P3 outcome", "migraine_any": "P3 outcome", "migraine_aura": "P3 outcome",
        "obesity": "P3 outcome",
    }
    for v, use in model_vars.items():
        s = elig[v]
        miss = int(s.isna().sum())
        miss_rows.append({"Variable": v, "Used in": use, "Missing n": miss,
                          "Missing %": round(100 * miss / len(elig), 2),
                          "> 10% missing": "yes" if miss / len(elig) > 0.10 else "no",
                          "Source of missingness": {
                              "migraine3": "code 9 (treated as missing)",
                              "smoking3": "none: code 9 retained as 'Unknown' level (per brief)",
                              "smoking_ever": "code 9 (unknown) excluded",
                              "migraine_any": "code 9 excluded", "migraine_aura": "code 9 excluded",
                              "hormonal_type_cat": "code 9 (treated as missing)",
                              "anemia_cat": "no Hb within ±3 y", "iron_deficiency": "no ferritin within ±3 y",
                              "obesity": "BMI blank"}.get(v, "blank")})
    mt = pd.DataFrame(miss_rows)
    add_table("missingness", mt, "Missingness of analysis variables, eligible population.")
    RESULTS["missing"] = {r["Variable"]: (r["Missing n"], r["Missing %"]) for r in miss_rows}

    # ---------------- Data-quality flags (do not fix) ----------------
    flags = []

    def flag(key, desc, n_all, n_elig, handling):
        flags.append({"Check": key, "Description": desc, "n (all records)": n_all,
                      "n (eligible)": n_elig, "Handling": handling})
        RESULTS[f"flag_{key}"] = (n_all, n_elig)
        log("Data problem", f"{key}: {desc} — all records {n_all:,}, eligible {n_elig:,}. Handling: {handling}")

    e = elig
    flag("myomectomy_no_fibroids", "myomectomy == 1 but fibroids == 0",
         int(((df.myomectomy == 1) & (df.fibroids == 0)).sum()), int(((e.myomectomy == 1) & (e.fibroids == 0)).sum()),
         "Not changed. fibroids/uterine_dx_group used as coded.")
    flag("surgtx_myomectomy_no_fibroids", "surgical_tx == 1 (myomectomy) but fibroids == 0",
         int(((df.surgical_tx == 1) & (df.fibroids == 0)).sum()), int(((e.surgical_tx == 1) & (e.fibroids == 0)).sum()),
         "Not changed.")
    flag("myomectomy_before_index", "first myomectomy dated before index_date (index = first uterine dx)",
         int((df.myomectomy_timing == "before index").sum()), int((e.myomectomy_timing == "before index").sum()),
         "Not changed; implies index_date is not the first fibroid documentation for these patients.")
    flag("stroke_no_timing", "stroke_any == 1 with blank stroke_timing",
         int(((df.stroke_any == 1) & df.stroke_timing.isna()).sum()), int(((e.stroke_any == 1) & e.stroke_timing.isna()).sum()),
         "Kept as cases in primary analyses; neither case nor control in incident-only sensitivity analyses.")
    flag("stroke_timing_9", "stroke_timing == 9 (unknown)", int((df.stroke_timing == 9).sum()), int((e.stroke_timing == 9).sum()),
         "As above.")
    flag("migraine_9_is_migraine",
         "migraine == 9: per migraine_terms sheet, code 9 was assigned to 'migraine, aura status not stated' "
         "(e.g. 'Chronic Migraine'), i.e. these patients HAVE migraine; Codebook label is 'Unknown'",
         int((df.migraine == 9).sum()), int((e.migraine == 9).sum()),
         "Codebook/brief definition kept (9 = unknown -> missing). Sensitivity analyses keep 9 as its own level. "
         "QUESTION FOR PI: should migraine_any include code 9?")
    ld = pd.to_datetime(e.lab_date, errors="coerce")
    idx = pd.to_datetime(e.index_date, errors="coerce")
    lab_days = (ld - idx).dt.days
    flag("lab_date_sentinel", "lab_date == 1000-01-01 (sentinel)", int((df.lab_date == "1000-01-01").sum()),
         int((e.lab_date == "1000-01-01").sum()), "Not used; MCV/platelet timing cannot be verified for these rows.")
    flag("lab_date_outside_12m", "lab_date more than 365 days from index_date (Codebook: ±12 months)",
         -1, int((lab_days.abs() > 365).sum()),
         "Not changed. Codebook says lab_date was back-filled with the d-dimer date, so it does not date MCV/platelets.")
    flag("mcv_no_lab_date", "MCV present but lab_date blank", -1, int((e.mcv.notna() & e.lab_date.isna()).sum()),
         "Not changed; MCV/platelet timing unverifiable. ±1-year lab-timing sensitivity applies to Hb (hgb_days) only.")
    RESULTS["lab_days_median"] = float(lab_days.median())
    flag("race_zero", "race == '0' (not a valid label)", int((df.race.astype(str) == "0").sum()),
         int((e.race.astype(str) == "0").sum()), "Collapsed into Other/unknown.")
    flag("race_free_text", "race is free text (Codebook: coded 1-9)", int(df.race.nunique()), int(e.race.nunique()),
         "Collapsed to White / Black (black_race) / Asian (text starts 'Asian') / Other-unknown. Counts are distinct values.")
    flag("bmi_out_of_range", "BMI outside Codebook range 10–80 kg/m²",
         int(((df.bmi < 10) | (df.bmi > 80)).sum()), int(((e.bmi < 10) | (e.bmi > 80)).sum()), "Not changed; used as recorded.")
    flag("hgb_at_bounds", "Hb exactly at plausibility bound (3.0 or 20.0 g/dL)",
         int(df.hgb.isin([3.0, 20.0]).sum()), int(e.hgb.isin([3.0, 20.0]).sum()), "Not changed.")
    flag("hormonal_type_9", "hormonal_type == 9", int((df.hormonal_type == 9).sum()), int((e.hormonal_type == 9).sum()),
         "Treated as missing (excluded from Paper 2 complete-case models).")
    flag("etiology_code4", "stroke_etiology == 4 (Codebook: code 4 'not separately derivable')",
         int((df.stroke_etiology == 4).sum()), int((e.stroke_etiology == 4).sum()), "Descriptive only.")
    flag("etiology_nonischaemic", "stroke_etiology populated where stroke_type != 1 (Codebook: only if ischaemic)",
         int((df.stroke_etiology.notna() & (df.stroke_type != 1)).sum()),
         int((e.stroke_etiology.notna() & (e.stroke_type != 1)).sum()), "Descriptive only; reported as recorded.")
    flag("imaging_positive_no_stroke", "stroke_confirmed_imaging == 1 but stroke_any == 0",
         int(((df.stroke_confirmed_imaging == 1) & (df.stroke_any == 0)).sum()),
         int(((e.stroke_confirmed_imaging == 1) & (e.stroke_any == 0)).sum()), "None expected.")
    flag("stroke_any_9", "stroke_any == 9 or blank", int((~df.stroke_any.isin([0, 1])).sum()),
         int((~e.stroke_any.isin([0, 1])).sum()), "None.")
    n_incl_diff = int((e.stroke_any_incl_imaging == 1).sum() - (e.stroke_any == 1).sum())
    flag("incl_imaging_diff",
         "stroke_any_incl_imaging − stroke_any (Codebook describes this as +408 imaging-only patients)",
         int((df.stroke_any_incl_imaging == 1).sum() - (df.stroke_any == 1).sum()), n_incl_diff,
         "Imaging-only cases have already been merged into stroke_any; the Codebook text for "
         "stroke_any_incl_imaging is out of date.")
    # derived-variable consistency
    a = pd.cut(e.hgb, [-1, 8, 10, 12, 99], right=False, labels=[3, 2, 1, 0]).astype(float)
    RESULTS["anemia_mismatch"] = int(((a != e.anemia) & e.hgb.notna()).sum())
    RESULTS["irondef_mismatch"] = int((((e.ferritin < 30).astype(float) != e.iron_deficiency) & e.ferritin.notna()).sum())
    g = np.select([(e.fibroids + e.adenomyosis + e.endometriosis) > 1, e.fibroids == 1, e.adenomyosis == 1,
                   e.endometriosis == 1], [4, 1, 2, 3], 0)
    RESULTS["dxgroup_mismatch"] = int((g != e.uterine_dx_group).sum())
    log("Derivation check", f"anemia vs hgb mismatches: {RESULTS['anemia_mismatch']}; iron_deficiency vs ferritin "
                            f"mismatches: {RESULTS['irondef_mismatch']}; uterine_dx_group vs flags mismatches: "
                            f"{RESULTS['dxgroup_mismatch']}.")
    ldl = e.ldl.dropna()
    flag("ldl_distribution", f"LDL distribution implausible for this population: median {ldl.median():.0f}, "
                             f"P75 {ldl.quantile(.75):.0f}, >190 mg/dL in {int((ldl > 190).sum()):,}",
         -1, int((ldl > 190).sum()), "Not used in any model; possible unit/test mix-up (e.g. total cholesterol).")
    add_table("data_flags", pd.DataFrame(flags), "Data problems flagged, not fixed. -1 = not computed for all records.")

    RESULTS["eligible_race"] = elig["race4"].value_counts().to_dict()
    RESULTS["n_stroke_before_or_same"] = int(elig.stroke_timing.isin([1, 2]).sum())
    RESULTS["n_stroke_before"] = int((elig.stroke_timing == 1).sum())
    return df, elig
