"""Classify the supplementary Epic diagnosis extracts (Diagnosis_22-26) into
comorbidity flags for the anaemia paper.

Code rules are prefix-based on ICD-10-CM / ICD-9-CM codes (dots removed for
matching). Rows whose code system is Mayo 'HIC' carry no ICD code and are
classified from the description text. Rows that match no rule are listed in
the log for review; they never count as a condition.
"""
import glob
import re

import numpy as np
import pandas as pd

from .utils import DATA_DIR

# (flag, list of code prefixes without dots). ICD-9 and ICD-10 together.
RULES = [
    # --- haemoglobinopathies
    ("sickle_disease", ["D570", "D571", "D572", "D574", "D578", "2826"]),
    ("sickle_trait", ["D573", "2825"]),
    ("thal_minor", ["D563", "28246", "D564"]),          # thalassaemia minor/trait, HPFH
    ("thal_other", ["D560", "D561", "D562", "D565", "D568", "D569", "2824"]),
    ("other_haemoglobinopathy", ["D580", "D581", "D582", "D588", "D589", "D55", "2827", "2820", "2821", "2822",
                                 "2823", "2829"]),
    # --- kidney
    ("ckd_esrd", ["N186", "N185", "I120", "I1311", "I132", "Z992", "Z940", "Z4901", "Z4902", "Z4931", "Z4932",
                  "5856", "5855", "40301", "40311", "40391", "40402", "40403", "40412", "40413", "40492", "40493",
                  "V4511", "V420", "V560", "V568"]),
    ("ckd", ["N18", "I12", "I13", "E0822", "E0922", "E1022", "E1122", "E1322", "D631", "28521", "585", "403", "404"]),
    # --- liver
    ("cirrhosis", ["K703", "K704", "K717", "K72", "K74", "K766", "K767", "I85", "5712", "5715", "5716", "5722", "5723",
                   "5724", "4560", "4561", "4562", "5728"]),
    ("chronic_liver", ["K70", "K73", "K75", "K760", "B18", "5710", "5711", "5713", "5714", "5718", "5719", "0702",
                       "0703", "07044", "07054"]),
    # --- alcohol
    ("alcohol_disorder", ["F101", "F102", "F1092", "F1093", "F1094", "F1095", "F1096", "F1098", "F1099", "K292",
                          "K852", "K860", "G621", "I426", "G312", "303", "3050", "291", "5353", "4255", "V113"]),
    # --- GI bleeding / IBD / malabsorption
    ("gi_bleed", ["K920", "K921", "K922", "K625", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266",
                  "K270", "K272", "K274", "K276", "K280", "K282", "K284", "K286", "K31811", "K5521", "K5701",
                  "K5711", "K5713", "K5721", "K5731", "K5733", "K5741", "K5751", "K5753", "K5781", "K5791", "K5793",
                  "5780", "5781", "5789", "5693", "5310", "5312", "5314", "5316", "5320", "5322", "5324", "5326",
                  "5330", "5332", "5334", "5336", "5340", "5342", "5344", "5346"]),
    ("ibd", ["K50", "K51", "555", "556"]),
    ("malabsorption", ["K900", "K904", "K909", "5790", "Z9884", "V4586"]),
    # --- HIV (infection only; screening/exposure/PrEP excluded below)
    ("hiv", ["B20", "Z21", "042", "V08", "O987"]),
    # --- menopause (descriptions mentioning perimenopause / hot flash are dropped in classify_rows)
    ("menopause", ["Z780", "N951", "E2831", "E2839", "E894", "6272", "6274", "2563", "V4981"]),
    # --- pregnancy / delivery / postpartum
    ("pregnancy", ["O", "Z331", "Z34", "Z3A", "Z37", "Z39", "V22", "V23", "V24", "V27"] +
                  [str(c) for c in range(630, 680)]),
]
# Second family: coded anaemia type (a code can carry a comorbidity flag above AND an anaemia type)
ANAEMIA_RULES = [
    ("anaemia_pregnancy", ["O990", "O9081", "6482"]),
    ("anaemia_iron", ["D50", "280"]),
    ("anaemia_b12_folate", ["D51", "D52", "281"]),
    ("anaemia_nutritional", ["D53"]),
    ("anaemia_blood_loss_acute", ["D62", "2851"]),
    ("anaemia_neoplastic_chemo", ["D630", "D6481", "28522", "2853"]),
    ("anaemia_chronic_disease", ["D631", "D638", "28521", "28529"]),
    ("anaemia_haemolytic_aplastic", ["D59", "D60", "D61", "283", "284"]),
    ("anaemia_unspecified", ["D649", "D6489", "D648", "2859", "2858"]),
]
# explicit non-qualifying codes (kept out on purpose; reported in the log)
EXCLUDE = {"7948": "abnormal liver function test", "5738": "other disorders of liver (nonspecific)",
           "Z0283": "drug/alcohol screen", "Z811": "family history of alcohol abuse", "F1090": "alcohol use, no disorder",
           "F1091": "alcohol use, no disorder, in remission", "Z7141": "alcohol counselling",
           "Z114": "HIV screening", "R75": "inconclusive HIV test", "Z206": "HIV exposure", "Z2981": "HIV PrEP",
           "79571": "nonspecific HIV serology", "V016": "HIV exposure", "D75A": "G6PD deficiency without anaemia",
           "V780": "anaemia screening", "V781": "anaemia screening"}
# HIC (internal) rows: description regexes
HIC_RULES = [
    ("sickle_disease", r"sickle.?cell (disease|anemia|anaemia)|hb.?ss"),
    ("sickle_trait", r"sickle.?cell trait"),
    ("thal_minor", r"thalass?emia (minor|trait)"),
    ("thal_other", r"thalass?emia"),
    ("ckd_esrd", r"dialysis|end.?stage renal|kidney transplant|renal transplant"),
    ("ckd", r"chronic (renal|kidney)"),
    ("cirrhosis", r"cirrhosis|portal hypertension|varices"),
    ("chronic_liver", r"fatty,? liver|liver, fatty|steatohepatitis|chronic hepatitis|alcoholic liver"),
    ("alcohol_disorder", r"alcoholi(c|sm)|alcohol (abuse|dependence|addiction)|abuse, alcohol"),
    ("gi_bleed", r"melena|hematemesis|gastrointestinal (hemorrhage|bleed)|rectal (bleed|hemorrhage)"),
    ("ibd", r"crohn|ulcerative colitis"),
    ("malabsorption", r"celiac|sprue|malabsorption|bariatric"),
    ("pregnancy", r"deliver|pregnan|postpartum|puerper|cesarean|labor|gestation|abortion|miscarriage"),
    ("hiv", r"hiv test positive|human immunodeficiency virus (disease|infection|positive)|\baids\b"),
    ("menopause", r"^menopause"),
    ("anaemia_iron", r"anemia, (iron deficiency|blood loss|microcytic)|iron deficiency, secondary to bleeding"),
    ("anaemia_unspecified", r"^anemia"),
]
FLAGS = [r[0] for r in RULES] + [r[0] for r in ANAEMIA_RULES]
PERI_RE = re.compile(r"perimenopaus|hot flash|hot flush", re.I)


def _clean(code):
    return re.sub(r"[^A-Z0-9]", "", str(code).upper())


def classify_code(code):
    """Return the first matching comorbidity flag (or 'EXCLUDED' / None)."""
    c = _clean(code)
    if c in EXCLUDE or any(c.startswith(e) for e in EXCLUDE):
        return "EXCLUDED"
    for flag, prefixes in RULES:
        for p in prefixes:
            if c.startswith(p):
                return flag
    return None


def anaemia_type(code):
    c = _clean(code)
    for flag, prefixes in ANAEMIA_RULES:
        if any(c.startswith(p) for p in prefixes):
            return flag
    return None


def load_extracts():
    files = sorted(glob.glob(str(DATA_DIR / "dx" / "MDE_Workflow_Results_for_Diagnosis_*.csv")))
    parts = []
    for f in files:
        x = pd.read_csv(f, dtype=str)
        x["source_file"] = f.rsplit("_", 1)[-1].replace(".csv", "")
        parts.append(x)
    x = pd.concat(parts, ignore_index=True)
    x.columns = ["mrn", "system", "code", "desc", "date", "source_file"]
    x["mrn"] = pd.to_numeric(x["mrn"], errors="coerce")
    x["date"] = pd.to_datetime(x["date"], errors="coerce")
    x = x.drop_duplicates(subset=["mrn", "system", "code", "desc", "date"])
    return x, files


def classify_rows(x):
    rows = []
    for r in x.itertuples(index=False):
        flags = set()
        if str(r.system).upper() == "HIC":
            d = str(r.desc).lower()
            for flag, pat in HIC_RULES:
                if re.search(pat, d):
                    flags.add(flag)
                    break
            if not flags:
                flags.add("UNMATCHED")
        else:
            codes = [c.strip() for c in str(r.code).split(",")]
            ks = [classify_code(c) for c in codes]
            hits = {k for k in ks if k and k != "EXCLUDED"} | {a for a in map(anaemia_type, codes) if a}
            if "menopause" in hits and PERI_RE.search(str(r.desc)):
                hits.discard("menopause")
            flags = hits or ({"EXCLUDED"} if "EXCLUDED" in ks else {"UNMATCHED"})
        for fl in flags:
            rows.append((r.mrn, r.system, r.code, r.desc, r.date, r.source_file, fl))
    return pd.DataFrame(rows, columns=["mrn", "system", "code", "desc", "date", "source_file", "flag"])
