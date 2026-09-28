"""Follow-up dates and anaemia treatment from the Encounters and Medications
Administered extracts (linked by MRN in memory only).

- last_contact: 'Create Date' in Encounters_8 - verified to be on or after every
  stroke, surgery and Hb date recorded for each patient, i.e. the date of the
  latest encounter.
- death: 'Death Date' (Encounters_6 / _8).
- treatment: administered IV iron, oral iron (therapeutic products only;
  multivitamins, prenatal vitamins, OC packs with iron placebo and
  spironolactone - matched on the letters 'iron' - are excluded) and ESA.
"""
import glob
import re

import numpy as np
import pandas as pd

from .utils import DATA_DIR

GIVEN = {"given", "done", "new bag", "given by other", "auth (verified)", "restarted", "rate change", "rate verify",
         "stopped", "started during downtime"}
EXCLUDE_MED = re.compile(r"multivit|prenatal|pediatric|norethin|ethinyl|spironol|buspir", re.I)
IV_IRON = re.compile(r"iron sucrose|iron dextran|carboxymaltose|ferumoxytol|ferric gluconate|derisomaltose|injectafer|"
                     r"venofer|infed|feraheme", re.I)
ORAL_IRON = re.compile(r"ferrous|ferric citrate|carbonyl|polysaccharide iron|elemental iron|iron bisglycinate|"
                       r"iron fum|fe-fa|ferric maltol|iron polysaccharide", re.I)
ESA = re.compile(r"epoetin|darbepoetin|methoxy polyethylene", re.I)


def encounters():
    fs = sorted(glob.glob(str(DATA_DIR / "enc" / "MDE_Workflow_Results_for_Encounters_*.csv")))
    last, death = [], []
    for f in fs:
        x = pd.read_csv(f, dtype=str)
        x["mrn"] = pd.to_numeric(x["Clinic Number"], errors="coerce")
        if "Death Date" in x:
            death.append(pd.DataFrame({"mrn": x.mrn, "death": pd.to_datetime(x["Death Date"].str.slice(0, 10),
                                                                              errors="coerce")}))
        if "Arrive Date" in x and "Create Date" in x:  # Encounters_8 layout
            last.append(pd.DataFrame({"mrn": x.mrn, "last": pd.to_datetime(x["Create Date"].str.slice(0, 10),
                                                                            errors="coerce")}))
    d = pd.concat(death).dropna().groupby("mrn").death.min() if death else pd.Series(dtype="datetime64[ns]")
    lc = pd.concat(last).dropna().groupby("mrn")["last"].max() if last else pd.Series(dtype="datetime64[ns]")
    return lc, d


def treatments():
    fs = sorted(glob.glob(str(DATA_DIR / "meds" / "MDE_Workflow_Results_for_Medications_Administered_*.csv")))
    if not fs:
        return pd.DataFrame(columns=["mrn", "date", "cls"]), pd.DataFrame()
    x = pd.concat([pd.read_csv(f, dtype=str) for f in fs], ignore_index=True)
    x["mrn"] = pd.to_numeric(x["Clinic Number"], errors="coerce")
    x["date"] = pd.to_datetime(x["Administered Date"].str.slice(0, 10), errors="coerce")
    name = (x["Medication Generic Name"].fillna("") + " | " + x["Medication Name"].fillna(""))
    route = x["Administered Route"].fillna("").str.strip().str.lower()
    status = x["Administered Status"].fillna("").str.strip().str.lower()
    cls = np.select([name.str.contains(EXCLUDE_MED), name.str.contains(ESA),
                     name.str.contains(IV_IRON) | (name.str.contains("iron", case=False) &
                                                   route.str.contains("intraven|iv|ivpb")),
                     name.str.contains(ORAL_IRON)],
                    ["excluded", "esa", "iv_iron", "oral_iron"], default="unclassified")
    x["cls"] = cls
    x["given"] = status.isin(GIVEN)
    audit = x.groupby(["cls", "Medication Generic Name"]).agg(rows=("mrn", "size"), patients=("mrn", "nunique"),
                                                               given_rows=("given", "sum")).reset_index() \
        .sort_values(["cls", "patients"], ascending=[True, False])
    keep = x[x.given & x.cls.isin(["iv_iron", "oral_iron", "esa"])][["mrn", "date", "cls"]].dropna()
    return keep, audit
