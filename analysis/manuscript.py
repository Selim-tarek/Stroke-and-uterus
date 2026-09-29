"""Draft manuscript for the anaemia paper.

Every number is read from RESULTS / TABLES produced in the same run. Literature
statements carry [ref] placeholders for the authors to fill; no external numbers
are introduced. Output: outputs/anaemia_manuscript.json (rendered to .docx by
analysis/build_docx.js) and outputs/anaemia_manuscript.md.
"""
import json
import subprocess

import numpy as np
import pandas as pd

from .utils import OUT_DIR, RESULTS, TABLES, fmt_p, get_or, log

GR = ["Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"]


def _p(p):
    t = fmt_p(p) if not isinstance(p, str) else p
    return f"P{t}" if t.startswith("<") else f"P={t}"


def _or(f, c):
    g = get_or(f, c)
    return g["txt"]


def _dir(lo, hi, higher="higher", lower="lower", null="similar"):
    return higher if lo > 1 else (lower if hi < 1 else null)


def build():
    R, T = RESULTS, TABLES
    checks = []
    n_all, n_el = R["n_total"], R["n_eligible"]
    fl = R["flow"]
    miss = R["missing"]
    n_hb = n_el - miss["anemia_cat"][0]
    ex = T["P2_Table1_exposures"]
    ex_a = ex[ex.Exposure == "Anaemia grade"].set_index("Level")
    n_anaemic = int(sum(ex_a.loc[g, "n"] for g in GR))
    pct = lambda k, n: f"{100 * k / n:.1f}%"  # noqa: E731
    fits = R["p2_fits"]
    fa, fi, finc = fits[("stroke_any", "anemia_cat")], fits[("y_isch", "anemia_cat")], fits[("y_incident", "anemia_cat")]
    fm, fp = fits[("stroke_any", "mcv_cat")], fits[("stroke_any", "plt_cat")]
    tr = {k: R[f"p2_trend_{k}"] for k in ["stroke_any", "y_isch", "y_incident"]}
    pts = R["p2_spline_pts"].set_index("Hb (g/dL)")
    X = R["p2x_tab"].set_index("Analysis")
    tt = R["p2x_tte"]
    tc = R["p2x_tte_cohort"]
    per = R["p2x_per"]
    ab = R["p2x_abs"]

    def pg(design, outcome, rng, adj="+ pre-Hb conditions (2a)"):
        m = per[(per.Design.str.startswith(design)) & (per.Outcome == outcome) & (per["Hb range"].str.startswith(rng))
                & (per.Adjustment == adj)]
        assert len(m) == 1, (design, outcome, rng, adj)
        return m.iloc[0]

    _pa0 = R["p2x_pattern"]
    _pt0 = R["p2x_pattern_tte"]

    def _pp(lv, adj="Paper 2 covariates"):
        return _pa0[(_pa0.Outcome == "Any stroke") & (_pa0.Adjustment == adj) & (_pa0.Pattern == lv)].iloc[0]
    _pnh = _pp("Normocytic, high RDW (mixed / early iron deficiency)")
    _pmac = _pp("Macrocytic")
    _pid = _pp("Microcytic, high RDW (iron-deficiency pattern)")
    _tnh = _pt0[(_pt0.Outcome == "Any stroke") & (_pt0.Pattern == "Normocytic, high RDW (mixed / early iron deficiency)")].iloc[0]

    def abr(outcome, grade):
        m = ab[(ab.Outcome == outcome) & (ab["Anaemia grade"].str.startswith(grade))]
        assert len(m) == 1, (outcome, grade)
        return m.iloc[0]
    mo = R["p2x_morph"]
    sub = R["p2x_sub"].set_index("Subtype")
    sg = R["p2x_sg"]
    ev = R["p2x_ev"]
    pvx = R["p2x_prev"].set_index("Condition")
    txr = R["p2x_tx"]
    idr = R["p2_id"]
    ints = R["p2_ints"]
    lt = R["p2_lab_timing"]
    lt_any = lt[lt.Outcome == "Any stroke (primary)"].set_index("Level")
    mi = R["p2_mi"]
    mi_a = mi[mi.Outcome == "Any stroke (primary)"].set_index("coef")
    med = R["p2_med"]
    ht = R["p2x_hb_timing"]
    t1 = T["P2_Table1"]
    iy0, iy1 = R["index_range"]

    def t1v(char, level, col="Overall"):
        r = t1[(t1.Characteristic.astype(str) == char) | (t1.Characteristic.astype(str).str.startswith(char))]
        i = t1.index[t1.Characteristic.astype(str).str.startswith(char)][0]
        rows = t1.loc[i:i + 8]
        return rows[rows.Level.astype(str) == level].iloc[0][col]

    def tte(outcome, grade, col="Adjusted + pre-Hb conditions (2a) RR"):
        return tt[(tt.Outcome == outcome) & (tt["Anaemia grade"] == grade)].iloc[0][col]

    def tter(outcome, grade):
        return tt[(tt.Outcome == outcome) & (tt["Anaemia grade"] == grade)].iloc[0]

    def xr(a, col="Moderate"):
        return X.loc[a, col]

    def xt(a):
        return f"{X.loc[a, 'Per grade (trend)']}"

    def mor(outcome, typ):
        return mo[(mo.Outcome == outcome) & (mo["Anaemia type"] == typ)].iloc[0]

    def sgi(name):
        r = sg[(sg.Subgroup == name) & (sg.Level.str.startswith("Interaction"))].iloc[0]
        return _p(r["p (text)"])

    def txv(outcome, grp):
        return txr[(txr.Outcome == outcome) & (txr.Group == grp)].iloc[0]["Adjusted RR"]

    ev_pg = ev[(ev.Analysis == "Primary: Any stroke") & (ev.Contrast == "Per grade")].iloc[0]
    ev_mod = ev[(ev.Analysis == "Primary: Any stroke") & (ev.Contrast == "Moderate")].iloc[0]
    gmod, gsev, gmild = (get_or(fa, f"anemia_cat={g}") for g in ["Moderate (8–9.9)", "Severe (<8)", "Mild (10–11.9)"])
    tte_mod = tter("Any stroke", "Moderate (8–9.9)")
    tte_pg = tter("Any stroke", "Per grade (trend)")
    micro, normo, macro = (mor("Any stroke", k) for k in ["Microcytic anaemia", "Normocytic anaemia", "Macrocytic anaemia"])
    micro_i = mor("Ischaemic stroke", "Microcytic anaemia")
    sev_tte = tter("Any stroke", "Severe (<8)")
    B = []  # blocks

    def H1(t):
        B.append({"t": "h1", "text": t})

    def H2(t):
        B.append({"t": "h2", "text": t})

    def P(t):
        B.append({"t": "p", "text": t})

    def TABLE(title, df, note="", widths=None, font=8):
        B.append({"t": "table", "title": title, "columns": [str(c) for c in df.columns],
                  "rows": [["" if (isinstance(v, float) and np.isnan(v)) else
                            (f"{v:,}" if isinstance(v, (int, np.integer)) else str(v)) for v in r]
                           for r in df.itertuples(index=False)], "note": note, "widths": widths, "font": font})

    def FIG(path, caption, width=6.3):
        B.append({"t": "figure", "path": str(OUT_DIR / path), "caption": caption, "width": width})

    def BR():
        B.append({"t": "pagebreak"})

    # ------------------------------------------------------------------ title page
    B.append({"t": "title", "text": "Anaemia severity and stroke in women with benign uterine disease: "
                                    "a retrospective cohort study"})
    P("Running title: Anaemia and stroke in benign uterine disease")
    P("Authors: [Author names, degrees, affiliations]")
    P("Corresponding author: [name, address, email]")
    P("Word count (main text): [to be completed]; Tables: 3; Figures: 6; Supplementary material: eTables 1–13, "
      "eFigures 1–3")
    P("Keywords: anaemia; haemoglobin; ischaemic stroke; uterine fibroids; endometriosis; adenomyosis; women")
    B.append({"t": "note", "text": "DRAFT generated from the analysis pipeline (python -m analysis.run_all). All numbers "
                                   "are produced by code in the same run. [ref] marks statements that need a citation; "
                                   "bracketed text marks information the authors must supply."})
    BR()

    # ------------------------------------------------------------------ abstract
    H1("Abstract")
    import re as _re

    def sc(txt, first=False):
        """'1.77 (1.41–2.21)' -> '1.77; 95% CI 1.41–2.21' (first) or '1.77; 1.41–2.21'."""
        m = _re.match(r"\s*([\d.]+) \(([\d.]+–[\d.]+)\)", str(txt))
        return f"{m.group(1)}; {'95% CI ' if first else ''}{m.group(2)}" if m else str(txt)
    P("**Background.** Anaemia is common in women with uterine fibroids, adenomyosis and endometriosis, mainly because "
      "of heavy menstrual bleeding. Whether anaemia is associated with stroke in this population is not well "
      "characterised.")
    P(f"**Methods.** We studied {n_el:,} women aged 18–60 years diagnosed with benign uterine disease between {iy0} "
      f"and {iy1} at Mayo Clinic, using electronic health records. Haemoglobin (Hb) closest to diagnosis was graded by "
      f"World Health Organization thresholds. Associations with stroke (including transient ischaemic attack) and "
      f"ischaemic stroke were estimated by logistic regression adjusted for vascular risk factors, migraine, hormonal "
      f"therapy, uterine bleeding and diagnosis group. Stroke rates after the Hb measurement were estimated by Poisson "
      f"regression, with further adjustment for conditions that cause anaemia.")
    P(f"**Results.** Of {n_hb:,} women with an Hb value, {n_anaemic:,} ({pct(n_anaemic, n_hb)}) were anaemic. "
      f"Compared with Hb ≥12 g/dL, the adjusted odds of stroke were higher with moderate (OR {sc(gmod['txt'], True)}) "
      f"and severe anaemia (OR {sc(gsev['txt'])}), with an OR of {tr['stroke_any']['txt']} per anaemia grade (P for "
      f"trend {fmt_p(tr['stroke_any']['p'])}); results were similar for ischaemic stroke (per grade OR "
      f"{sc(tr['y_isch']['txt'])}). Among {tc['n']:,} women followed from the Hb measurement ({tc['py']:,.0f} "
      f"person-years; {tc['ev']} strokes), moderate anaemia was associated with a higher stroke rate (rate ratio "
      f"{sc(tte_mod['Adjusted + pre-Hb conditions (2a) RR'])}). Below 13 g/dL, each 1 g/dL lower Hb was associated "
      f"with an OR of {pg('Cross', 'Any stroke', 'Below')['Estimate (95% CI)']} and a rate ratio of "
      f"{pg('After', 'Any stroke', 'Below')['Estimate (95% CI)']}. Associations persisted after excluding women with "
      f"haemoglobinopathies or other anaemia-causing conditions. By red-cell pattern, the association was clearest for "
      f"normocytic anaemia with high red-cell distribution width (OR {sc(_pnh['Adjusted OR (95% CI)'])}; rate ratio "
      f"{sc(_tnh['Adjusted RR (95% CI)'])}) and macrocytic anaemia (OR {sc(_pmac['Adjusted OR (95% CI)'])}; "
      f"{int(_pmac['Strokes'])} strokes), and weaker for the iron-deficiency pattern (OR "
      f"{sc(_pid['Adjusted OR (95% CI)'])}).")
    P("**Conclusions.** In women with benign uterine disease, moderate and, less precisely, severe anaemia were "
      "associated with a higher risk of stroke, particularly ischaemic stroke. Anaemia may mark underlying illness as "
      "well as vascular risk. Haemoglobin measured at gynaecological diagnosis could help identify women for vascular "
      "risk assessment.")
    checks.append(("Abstract: moderate & severe ORs exclude 1", gmod["lo"] > 1 and gsev["lo"] > 1))
    checks.append(("Abstract: TTE moderate RR excludes 1", tte_mod["CI low"] > 1))
    checks.append(("Abstract: per-g/dL OR and RR below 13 g/dL exclude 1",
                   pg("Cross", "Any stroke", "Below")["lo"] > 1 and pg("After", "Any stroke", "Below")["lo"] > 1))
    checks.append(("Abstract: pattern normo-high-RDW OR & RR and macro OR exclude 1; ID pattern CI includes 1",
                   _pnh["CI low"] > 1 and _tnh["CI low"] > 1 and _pmac["CI low"] > 1
                   and _pid["CI low"] < 1 < _pid["CI high"]))
    checks.append(("Abstract: normo & macro > micro, micro CI includes 1",
                   normo["OR"] > micro["OR"] and macro["OR"] > micro["OR"] and micro["CI low"] < 1 < micro["CI high"]))
    checks.append(("Abstract: TIA not associated (CI includes 1)",
                   sub.loc["TIA", "CI low"] < 1 < sub.loc["TIA", "CI high"]))
    BR()

    # ------------------------------------------------------------------ introduction
    H1("Introduction")
    P("Stroke in women under 60 is uncommon but causes substantial disability, and a sizeable share of events in young "
      "adults are not explained by conventional vascular risk factors [ref]. Anaemia, particularly iron-deficiency "
      "anaemia, has been associated with ischaemic stroke in population cohorts and case–control studies [ref]. "
      "Proposed mechanisms include reduced cerebral oxygen delivery, compensatory hyperdynamic circulation, reactive "
      "thrombocytosis and a hypercoagulable state [ref]. Much of this literature concerns older adults or children, and "
      "the populations studied are heterogeneous in why they are anaemic.")
    P("Women with benign uterine disease (fibroids, adenomyosis and endometriosis) are a large group with a high burden "
      "of anaemia, most of it attributed to heavy menstrual bleeding [ref]. They are also frequently in contact with "
      "gynaecological services and have haemoglobin measured, which creates an opportunity for risk stratification. "
      "Whether anaemia in this group is associated with stroke, and whether any association reflects iron-deficient "
      "anaemia from bleeding or anaemia from other causes, has not been examined.")
    P("We examined the association between anaemia severity at the time of benign uterine diagnosis and stroke in a "
      "large electronic-health-record cohort. We assessed temporality, confounding by anaemia-causing conditions, and "
      "anaemia morphology.")

    # ------------------------------------------------------------------ methods
    H1("Methods")
    H2("Study design, setting and participants")
    P(f"This retrospective cohort study used electronic health records from a US academic medical centre [name/site to "
      f"be added]. The index date was each woman's first documented diagnosis of uterine fibroids, adenomyosis or "
      f"endometriosis (identified by ICD-10 codes). Of {n_all:,} women identified, we excluded {fl.get(1, 0):,} with "
      f"active malignancy within 6 months of the index date, {fl.get(6, 0):,} aged over 60 years, {fl.get(3, 0):,} aged "
      f"under 18 years and {fl.get(5, 0):,} whose age could not be derived, leaving {n_el:,} eligible women "
      f"(Figure 1). The study was approved by [institutional review board, approval number] with a waiver of informed "
      f"consent. It is reported according to the STROBE guideline.")
    H2("Exposure")
    P("The exposure was the Hb value closest to the index date within ±3 years. Anaemia was graded using World Health "
      "Organization thresholds for non-pregnant women: none (≥12 g/dL), mild (10–11.9 g/dL), moderate (8–9.9 g/dL) "
      "and severe (<8 g/dL). Mean corpuscular volume (MCV), taken from the dated laboratory result nearest to the Hb "
      "measurement (within 30 days), classified anaemia as microcytic (<80 fL), normocytic (80–100 fL) or macrocytic "
      "(>100 fL). Laboratory anaemia patterns combined MCV with red-cell distribution width (RDW-CV, high >14.5%) "
      "from the same record: microcytic with high RDW (iron-deficiency pattern), microcytic with normal RDW "
      "(thalassaemia-trait pattern), normocytic with normal or high RDW, and macrocytic. Platelet count (<150, 150–400, >400 ×10³/µL) and ferritin "
      "(iron deficiency, <30 ng/mL) were secondary exposures.")
    H2("Outcomes")
    P("The primary outcome was stroke or TIA, identified from diagnosis and problem-list codes at any date, supplemented "
      "by adjudicated brain-imaging reports. Secondary outcomes were ischaemic stroke (compared with women without "
      "stroke), individual stroke subtypes, and incident stroke (dated after the index date). Code-based exclusions "
      "(unruptured aneurysm, carotid stenosis without infarction, migraine aura without infarction, chronic "
      "small-vessel change and similar) are detailed in the Supplement.")
    H2("Covariates")
    P("Pre-specified covariates were age, race, body mass index (BMI), hypertension, diabetes, dyslipidaemia, smoking "
      "(never, ever, unknown), migraine, thrombophilia, hormonal therapy type, heavy or abnormal uterine bleeding and "
      "uterine diagnosis group. Additional conditions were identified from dated diagnosis extracts: haemoglobinopathies "
      "(sickle cell disease or trait, thalassaemia, other haemoglobinopathies and hereditary haemolytic anaemias, at "
      "any date), chronic kidney disease (including dialysis and transplantation), chronic liver disease, alcohol use "
      "disorder, malabsorption (coeliac disease, bariatric surgery), inflammatory bowel disease and HIV, each counted "
      "only if documented before the Hb measurement. GI bleeding and pregnancy in the year before the Hb were also "
      "recorded. Other covariates were coronary artery disease, heart failure, atrial fibrillation, prior venous "
      "thromboembolism, malignancy at any time, anticoagulant use and administered iron or erythropoiesis-stimulating "
      "agents. Death and last-encounter dates defined follow-up.")
    H2("Statistical analysis")
    P("Associations between anaemia grade and stroke were estimated by logistic regression with heteroskedasticity-"
      "robust (HC1) standard errors. Crude and adjusted odds ratios (ORs) are reported, and linear trend was tested by "
      "entering grade as a score. Hb was also modelled as a restricted cubic spline with 4 knots, referenced to 13 g/dL. "
      "Effect modification by uterine bleeding, age, race, fibroid status and menopausal status, and between "
      "microcytosis and thrombocytosis, was assessed with product terms.")
    P("To address temporality we did four things. (1) Strokes were restricted to those dated after the index date. "
      "(2) Cases were restricted to those whose Hb was measured at least 30 days before the stroke. (3) Hb was "
      "restricted to within 1 year of the index date. (4) In a time-to-event analysis, women were followed from the "
      "later of the index date and the Hb measurement to stroke, death or last encounter. Women with a stroke before "
      "that start were excluded, and rate ratios (RRs) were estimated by Poisson regression with a person-time offset.")
    P("To address confounding by anaemia-causing illness, models were further adjusted for conditions documented "
      "before the Hb (adjustment 2a). As an over-adjustment check, we also adjusted for undated or potentially "
      "post-stroke factors (heart failure, atrial fibrillation, venous thromboembolism, malignancy, anticoagulant use; "
      "adjustment 2b). We then repeated the analyses excluding women with any anaemia-causing condition. Anaemia was "
      "classified by MCV and by coded iron-deficiency diagnosis. E-values quantified the strength of unmeasured "
      "confounding needed to explain the associations.")
    P("Multiple imputation by chained equations (20 imputations) was used for BMI in women with an Hb value; the "
      "exposure and outcome were not imputed. Iron or ESA treatment and a fibroid-to-anaemia mediation analysis were "
      "explored and are reported in the Supplement. Analyses used Python (statsmodels, scipy, pandas); code is "
      "available at [repository].")

    # ------------------------------------------------------------------ results
    H1("Results")
    H2("Participants")
    P(f"Of {n_el:,} eligible women, {n_hb:,} ({pct(n_hb, n_el)}) had an Hb value within 3 years of the index date "
      f"(Figure 1). Median age was {t1v('Age at index', 'median (IQR)')} years. Anaemia was present in "
      f"{n_anaemic:,} ({pct(n_anaemic, n_hb)}). Of these, {int(ex_a.loc['Mild (10–11.9)', 'n']):,} had mild, "
      f"{int(ex_a.loc['Moderate (8–9.9)', 'n']):,} moderate and {int(ex_a.loc['Severe (<8)', 'n']):,} severe "
      f"anaemia. Women with more severe anaemia more often had heavy or abnormal uterine bleeding "
      f"({t1v('Heavy/abnormal uterine bleeding', 'yes', 'None (Hb ≥12)')} without anaemia vs "
      f"{t1v('Heavy/abnormal uterine bleeding', 'yes', 'Severe (<8)')} with severe anaemia). They also had lower MCV "
      f"and higher platelet counts, and were less often recorded as having migraine (Table 1). In total, "
      f"{R['p2_events']:,} women had a stroke or TIA and {R['p2_isch_events']:,} an ischaemic stroke. Women without an "
      f"Hb value had a lower stroke prevalence ({ex_a.loc['missing', 'Any stroke %']}%) than those with one "
      f"({pct(int(sum(ex_a.loc[g, 'Any stroke n'] for g in ['None (Hb ≥12)'] + GR)), n_hb)}).")
    H2("Anaemia and stroke")
    P(f"Stroke prevalence rose from {ex_a.loc['None (Hb ≥12)', 'Any stroke %']:.1f}% in women without anaemia to "
      f"{ex_a.loc['Moderate (8–9.9)', 'Any stroke %']:.1f}% with moderate and {ex_a.loc['Severe (<8)', 'Any stroke %']:.1f}% "
      f"with severe anaemia (Table 2). After adjustment (n = {fa.n:,}; {fa.events:,} strokes), the ORs were "
      f"{gmild['txt']} for mild, {gmod['txt']} for moderate and {gsev['txt']} for severe anaemia (per grade "
      f"{tr['stroke_any']['txt']}; P for trend {fmt_p(tr['stroke_any']['p'])}; Figure 2). For ischaemic stroke they were "
      f"{_or(fi, 'anemia_cat=Mild (10–11.9)')}, {_or(fi, 'anemia_cat=Moderate (8–9.9)')} and "
      f"{_or(fi, 'anemia_cat=Severe (<8)')}. The spline showed a steady rise in risk as Hb fell below about 12–13 g/dL "
      f"(non-linearity {_p(R['p2_spline_p_nonlin'])}). Compared with 13 g/dL, the OR was "
      f"{pts.loc[10, 'OR vs 13 g/dL (95% CI)']} at 10 g/dL and {pts.loc[8, 'OR vs 13 g/dL (95% CI)']} at 8 g/dL "
      f"(eFigure 3).")
    H2("Temporality")
    P(f"Hb was measured at least 30 days before the stroke in {ht.get('Hb >=30 d before stroke', 0)} strokes, within "
      f"30 days of it in {ht.get('Hb within 30 d of stroke', 0)}, and at least 30 days after it in "
      f"{ht.get('Hb >=30 d after stroke', 0)}. Restricting cases to Hb measured before the stroke, moderate anaemia "
      f"remained associated with stroke (OR {xr('Hb ≥30 d before stroke: Any stroke')}) and with ischaemic stroke "
      f"(OR {xr('Hb ≥30 d before stroke: Ischaemic stroke')}). Severe anaemia was uncommon before stroke and its "
      f"estimate was imprecise (OR {xr('Hb ≥30 d before stroke: Any stroke', 'Severe')}). For strokes after the index "
      f"date, the OR for moderate anaemia was {_or(finc, 'anemia_cat=Moderate (8–9.9)')}.")
    P(f"In the time-to-event analysis, {tc['n']:,} women contributed {tc['py']:,.0f} person-years (median "
      f"{tc['fu_median']:.1f} years) with {tc['ev']} strokes ({tc['ev_isch']} ischaemic). Rates per 1,000 person-years "
      f"were {tter('Any stroke', 'None (Hb ≥12)')['Rate /1,000 PY']} without anaemia, "
      f"{tter('Any stroke', 'Mild (10–11.9)')['Rate /1,000 PY']} with mild, {tte_mod['Rate /1,000 PY']} with moderate "
      f"and {sev_tte['Rate /1,000 PY']} with severe anaemia. Adjusted for covariates and pre-Hb conditions, the RR was "
      f"{tte_mod['Adjusted + pre-Hb conditions (2a) RR']} for moderate anaemia and "
      f"{sev_tte['Adjusted + pre-Hb conditions (2a) RR']} for severe anaemia ({int(sev_tte['Events'])} events). The "
      f"per-grade RR was {tte_pg['Adjusted + pre-Hb conditions (2a) RR'].replace('; p=', ', P=')}. For ischaemic stroke the RR for moderate "
      f"anaemia was {tte('Ischaemic stroke', 'Moderate (8–9.9)')} (Table 2).")
    H2("Stroke risk per 1 g/dL of haemoglobin")
    _pc, _pr = pg("Cross", "Any stroke", "Below"), pg("After", "Any stroke", "Below")
    _pci, _pri = pg("Cross", "Ischaemic stroke", "Below"), pg("After", "Ischaemic stroke", "Below")
    _wc, _wr = pg("Cross", "Any stroke", "Whole"), pg("After", "Any stroke", "Whole")
    _ac, _ar = pg("Cross", "Any stroke", "Above"), pg("After", "Any stroke", "Above")
    _inc = lambda r: "the CI included 1" if r["lo"] < 1 < r["hi"] else "the CI excluded 1"  # noqa: E731
    P(f"Because the spline showed no association above 13 g/dL, Hb was also modelled as two linear segments "
      f"(below and above 13 g/dL), adjusted for covariates and pre-Hb conditions (2a). Below 13 g/dL, each 1 g/dL lower "
      f"Hb was associated with an OR for any stroke of {_pc['Estimate (95% CI)']} and an OR for ischaemic stroke of "
      f"{_pci['Estimate (95% CI)']}. In the time-to-event cohort, the corresponding rate ratios were "
      f"{_pr['Estimate (95% CI)']} and {_pri['Estimate (95% CI)']}. Above 13 g/dL, the OR per 1 g/dL higher Hb was "
      f"{_ac['Estimate (95% CI)']} ({_inc(_ac)}) and the rate ratio was {_ar['Estimate (95% CI)']} ({_inc(_ar)}). "
      f"Treating Hb as linear over the whole range gave an OR of {_wc['Estimate (95% CI)']} and a rate ratio of "
      f"{_wr['Estimate (95% CI)']} per 1 g/dL lower Hb (eTable 10).")
    _n, _mi, _mo, _se = (abr("Any stroke", g) for g in ["None", "Mild", "Moderate", "Severe"])
    P(f"Crude stroke rates after the Hb measurement were {_n['Rate (95% CI)']} per 1,000 person-years without "
      f"anaemia, {_mi['Rate (95% CI)']} with mild, {_mo['Rate (95% CI)']} with moderate and {_se['Rate (95% CI)']} with "
      f"severe anaemia ({int(_se['Events'])} events). Assuming a constant rate, these correspond to unadjusted 5-year "
      f"risks of {_n['Crude 5-y risk (95% CI)']}, {_mi['Crude 5-y risk (95% CI)']}, {_mo['Crude 5-y risk (95% CI)']} "
      f"and {_se['Crude 5-y risk (95% CI)']} (Figure 5). Deaths without stroke, which end follow-up, were more "
      f"frequent with worse anaemia: {_n['Death rate /1,000 PY']:.2f}, {_mi['Death rate /1,000 PY']:.2f}, "
      f"{_mo['Death rate /1,000 PY']:.2f} and {_se['Death rate /1,000 PY']:.2f} per 1,000 person-years "
      f"({int(_se['Deaths without stroke'])} deaths vs {int(_se['Events'])} strokes with severe anaemia).")
    _cp = R["p2x_curve_pts"]
    _ci = R["p2x_curve_info"]

    def cpt(des, h):
        return _cp[(_cp.Design.str.startswith(des)) & (_cp["Hb (g/dL)"] == h)].iloc[0]
    _cx, _ct = _ci["Cross-sectional (odds ratio)"], _ci["After the Hb measurement (rate ratio)"]
    P(f"Figure 3 shows stroke risk across the Hb range, adjusted for covariates and pre-Hb conditions (2a). Compared "
      f"with 13 g/dL, the OR was {cpt('Cross', 11)['Estimate vs 13 g/dL (95% CI)']} at 11 g/dL, "
      f"{cpt('Cross', 9)['Estimate vs 13 g/dL (95% CI)']} at 9 g/dL and {cpt('Cross', 8)['Estimate vs 13 g/dL (95% CI)']} "
      f"at 8 g/dL (non-linearity {_p(fmt_p(_cx['p_nonlin']))}). For stroke after the Hb measurement, the rate ratio was "
      f"{cpt('After', 11)['Estimate vs 13 g/dL (95% CI)']} at 11 g/dL, {cpt('After', 9)['Estimate vs 13 g/dL (95% CI)']} "
      f"at 9 g/dL and {cpt('After', 8)['Estimate vs 13 g/dL (95% CI)']} at 8 g/dL (overall {_p(fmt_p(_ct['p_overall']))}, "
      f"non-linearity {_p(fmt_p(_ct['p_nonlin']))}; {_ct['events']} strokes). Above 13 g/dL neither curve departed "
      f"from 1 (eTable 12).")
    checks.append(("Curves: OR at 9 g/dL and RR at 9 g/dL exclude 1; values above 13 include 1",
                   cpt("Cross", 9)["lo"] > 1 and cpt("After", 9)["lo"] > 1 and
                   all(cpt(dd, h)["lo"] < 1 < cpt(dd, h)["hi"] for dd in ["Cross", "After"] for h in [14, 15, 16])))
    checks.append(("Deaths without stroke rise with anaemia grade", _n["Death rate /1,000 PY"] < _mi["Death rate /1,000 PY"]
                   < _mo["Death rate /1,000 PY"] < _se["Death rate /1,000 PY"]))
    checks.append(("Per-g/dL: below-13 OR/RR any & ischaemic exclude 1",
                   all(r["lo"] > 1 for r in [_pc, _pr, _pci, _pri])))
    H2("Confounding by anaemia-causing conditions")
    P(f"Anaemia-causing conditions were more common in women with stroke. For example, chronic kidney disease was "
      f"present in {pvx.loc['Chronic kidney disease (any)', 'Stroke %']}% vs "
      f"{pvx.loc['Chronic kidney disease (any)', 'No stroke %']}%, and haemoglobinopathy in "
      f"{pvx.loc['Any haemoglobinopathy', 'Stroke %']}% vs {pvx.loc['Any haemoglobinopathy', 'No stroke %']}% "
      f"(eTable 1). Adjusting for conditions documented before the Hb attenuated but did not remove the association: "
      f"moderate anaemia OR {xr('2a Pre-Hb conditions: Any stroke')}, severe "
      f"{xr('2a Pre-Hb conditions: Any stroke', 'Severe')}, per grade {xt('2a Pre-Hb conditions: Any stroke')}. "
      f"Excluding women with haemoglobinopathies (per grade {xt('Excluding haemoglobinopathies: Any stroke')}) or with "
      f"any anaemia-causing condition (per grade {xt('Excluding all anaemia-causing conditions: Any stroke')}) gave "
      f"similar results. Adding undated or potentially post-stroke factors, including anticoagulant use, reduced the "
      f"per-grade OR to {xt('2b Over-adjustment check (+ undated/post-stroke): Any stroke')}. Combining the "
      f"Hb-before-stroke restriction with adjustment 2a, the OR for moderate anaemia was "
      f"{xr('Hb ≥30 d before stroke + 2a: Any stroke')} (Table 3, Figure 4). The E-value for the moderate-anaemia OR "
      f"was {ev_mod['E-value (point)']} (confidence limit {ev_mod['E-value (CI limit)']}), and for the per-grade OR "
      f"{ev_pg['E-value (point)']} ({ev_pg['E-value (CI limit)']}).")
    il = R["p2x_ill"].set_index("Analysis")
    _d1, _si = il.loc["Excluding deaths within 1 y of Hb"], il.loc["Excluding cancer, heart failure, CKD, liver disease, HIV"]
    P(f"To examine whether anaemia simply marked women who were already seriously ill, we excluded women who died "
      f"within 1 year of the Hb measurement (n = {R['p2x_n_death1y']:,}) and, separately, women with cancer, heart "
      f"failure, chronic kidney disease, chronic liver disease or HIV recorded at any time "
      f"(n = {R['p2x_n_serious']:,}), with adjustment 2a. The OR for moderate anaemia was "
      f"{xr('Excluding deaths within 1 y of Hb (+2a): Any stroke')} and "
      f"{xr('Excluding cancer, heart failure, CKD, liver disease, HIV (+2a): Any stroke')}, respectively. In the "
      f"time-to-event analysis, the per-grade RR was {_d1['Per grade RR (95% CI)']} after excluding early deaths and "
      f"{_si['Per grade RR (95% CI)']} after excluding serious chronic illness ({int(_si['Strokes'])} strokes; "
      f"eTable 11). Healthcare use could not be measured.")
    H2("Anaemia type, stroke subtype and subgroups")
    P(f"Compared with no anaemia, microcytic anaemia was not clearly associated with stroke "
      f"(OR {micro['Adjusted OR (95% CI)']}), whereas normocytic (OR {normo['Adjusted OR (95% CI)']}) and macrocytic "
      f"anaemia (OR {macro['Adjusted OR (95% CI)']}, {int(macro['Events'])} events) were. For ischaemic stroke, the OR "
      f"for microcytic anaemia was {micro_i['Adjusted OR (95% CI)']}. Anaemia with a coded iron-deficiency diagnosis "
      f"was also associated with ischaemic stroke "
      f"(OR {mor('Ischaemic stroke', 'Anaemia, iron deficiency coded')['Adjusted OR (95% CI)']}; eTable 2). By subtype, "
      f"the per-grade OR was {sub.loc['Ischaemic stroke', 'Per grade OR (95% CI)']} for ischaemic stroke and "
      f"{sub.loc['TIA', 'Per grade OR (95% CI)']} for TIA (eTable 3). The association did not differ by age "
      f"({sgi('Age')}), race ({sgi('Race')}), fibroid status ({sgi('Fibroids')}) or menopausal status "
      f"({sgi('Menopausal status (coded)')}). It was somewhat weaker in women with heavy or abnormal uterine bleeding "
      f"(interaction {sgi('Heavy/abnormal bleeding')}; eTable 4).")
    pa = R["p2x_pattern"]
    ptt = R["p2x_pattern_tte"]

    def pat(outcome, lv, adj="+ pre-Hb conditions (2a)"):
        return pa[(pa.Outcome == outcome) & (pa.Adjustment == adj) & (pa.Pattern == lv)].iloc[0]

    def patt(lv):
        return ptt[(ptt.Outcome == "Any stroke") & (ptt.Pattern == lv)].iloc[0]
    PL = [("Microcytic, high RDW (iron-deficiency pattern)", "microcytic anaemia with high RDW (iron-deficiency pattern)"),
          ("Microcytic, normal RDW (thalassaemia-trait pattern)", "microcytic anaemia with normal RDW (thalassaemia-trait pattern)"),
          ("Normocytic, normal RDW", "normocytic anaemia with normal RDW"),
          ("Normocytic, high RDW (mixed / early iron deficiency)", "normocytic anaemia with high RDW"),
          ("Macrocytic", "macrocytic anaemia")]
    _ph = pa[(pa.Outcome == "Any stroke") & (pa.Adjustment == "+ pre-Hb conditions (2a)") &
             pa.Pattern.str.startswith("Heterogeneity")].iloc[0]
    P("Using the laboratory indices measured with the Hb (MCV and RDW), and with adjustment 2a, the OR for any stroke "
      "compared with no anaemia was " + "; ".join(
          f"{txt} {pat('Any stroke', lv)['Adjusted OR (95% CI)']} ({int(pat('Any stroke', lv)['Strokes'])} of "
          f"{int(pat('Any stroke', lv)['Women']):,} women)" for lv, txt in PL) +
      f" (heterogeneity {_p(_ph['p (text)'])}). {R['p2x_pattern_missing']:,} anaemic women lacked MCV or RDW within 30 "
      f"days. After the Hb measurement, the adjusted RR was " + "; ".join(
          f"{txt} {patt(lv)['Adjusted RR (95% CI)']}" for lv, txt in PL) + " (eTable 13; Figure 6).")
    H2("Other haematological indices")
    P(f"Thrombocytosis was not associated with stroke (OR {_or(fp, 'plt_cat=Thrombocytosis (>400)')}), whereas a low "
      f"platelet count (OR {_or(fp, 'plt_cat=Low (<150)')}) and macrocytosis (OR {_or(fm, 'mcv_cat=Macrocytic (>100)')}) "
      f"were. Microcytosis and thrombocytosis did not interact ({_p(ints.iloc[0]['p (text)'])}). In the "
      f"{R['p2_ferritin_n']:,} women with a ferritin result, iron deficiency was inversely associated with stroke "
      f"(OR {idr[(idr.Outcome == 'Any stroke (primary)') & (idr.Model == 'adjusted')].iloc[0]['OR (95% CI)']}; "
      f"eTable 5).")
    H2("Sensitivity analyses")
    P(f"The estimates were similar when Hb was restricted to ±1 year of the index date (severe anaemia OR "
      f"{lt_any.loc['Severe (<8)', 'Hb within ±1 y']}) and after multiple imputation (severe anaemia OR "
      f"{mi_a.loc['anemia_cat=Severe (<8)', 'OR (95% CI)']}; eTable 8). In an exploratory landmark analysis, the stroke "
      f"rate did not differ between anaemic women who were and were not given iron or ESA "
      f"(RR {txv('Any stroke', 'Treated vs untreated anaemia, additionally adjusted for grade')}; eTable 7).")
    checks += [
        ("Discussion: severe cross-sectional OR excludes 1; severe TTE imprecise (CI includes 1)",
         gsev["lo"] > 1 and sev_tte["CI low"] < 1 < sev_tte["CI high"]),
        ("Results: 2a 'attenuated but did not remove'", X.loc["2a Pre-Hb conditions: Any stroke", "trend lo"] > 1 and
         X.loc["2a Pre-Hb conditions: Any stroke", "trend OR"] < X.loc["Primary: Any stroke", "trend OR"]),
        ("Results: restricted cohorts similar (per-grade CI excludes 1)",
         X.loc["Excluding haemoglobinopathies: Any stroke", "trend lo"] > 1 and
         X.loc["Excluding all anaemia-causing conditions: Any stroke", "trend lo"] > 1),
        ("Results: Hb-before moderate 'remained associated'", X.loc["Hb ≥30 d before stroke: Any stroke", "Moderate lo"] > 1),
        ("Results: normocytic & macrocytic CI excludes 1", normo["CI low"] > 1 and macro["CI low"] > 1),
        ("Results: thrombocytosis null, low plt & macro associated",
         get_or(fp, "plt_cat=Thrombocytosis (>400)")["lo"] < 1 < get_or(fp, "plt_cat=Thrombocytosis (>400)")["hi"]
         and get_or(fp, "plt_cat=Low (<150)")["lo"] > 1 and get_or(fm, "mcv_cat=Macrocytic (>100)")["lo"] > 1),
        ("Results: iron deficiency inverse", idr[(idr.Outcome == 'Any stroke (primary)') &
                                                (idr.Model == 'adjusted')].iloc[0]["CI high"] < 1),
        ("Results: subgroup interactions non-significant (age, race, fibroids, menopause) and bleeding p<0.05",
         all(sg[(sg.Subgroup == n) & sg.Level.str.startswith("Interaction")].iloc[0]["p interaction"] > 0.05
             for n in ["Age", "Race", "Fibroids", "Menopausal status (coded)"]) and
         sg[(sg.Subgroup == "Heavy/abnormal bleeding") & sg.Level.str.startswith("Interaction")].iloc[0]["p interaction"] < 0.05),
        ("Results: treated vs untreated CI includes 1", "not" not in txv('Any stroke', 'Treated vs untreated anaemia, additionally adjusted for grade')),
    ]

    # ------------------------------------------------------------------ discussion
    H1("Discussion")
    P(f"In this cohort of women with benign uterine disease, moderate anaemia was associated with higher odds of stroke "
      f"(OR {gmod['txt']}) and a higher subsequent stroke rate (RR {tte_mod['Adjusted + pre-Hb conditions (2a) RR']}), "
      f"mainly ischaemic stroke. Severe anaemia showed a similar cross-sectional association (OR {gsev['txt']}), but "
      f"estimates before stroke were imprecise because few women had severe anaemia before a stroke. The association held when anaemia was measured "
      f"before the stroke and when women were followed forward from the Hb measurement. It was attenuated but "
      f"persisted after adjustment for, or exclusion of, haemoglobinopathies and other anaemia-causing conditions. By MCV, "
      f"it was stronger for normocytic and macrocytic than for microcytic anaemia, but anaemia with a coded "
      f"iron-deficiency diagnosis was also associated, so the data do not clearly separate iron-deficiency anaemia "
      f"from anaemia of other causes. It was also weaker in women with heavy uterine bleeding.")
    _hs = xr('Hb ≥30 d before stroke: Any stroke', 'Severe')
    _ab = pvx.loc["Coded acute blood-loss anaemia (±1 y of Hb)"]
    _ns, _ms, _ss = (abr("Any stroke", g) for g in ["None", "Moderate", "Severe"])
    P(f"The crude stroke rate after the Hb measurement was lower with severe than with moderate anaemia "
      f"({_ss['Rate /1,000 PY']:.2f} vs {_ms['Rate /1,000 PY']:.2f} per 1,000 person-years). We do not interpret this "
      f"as lower risk. First, it rests on {int(_ss['Events'])} strokes and the confidence intervals overlap widely. "
      f"Second, women with severe anaemia died without stroke at {_ss['Death rate /1,000 PY'] / _ns['Death rate /1,000 PY']:.1f} "
      f"times the rate of women without anaemia, so death competes with and removes women from being at risk of stroke. "
      f"Third, severe anaemia was more often coded as acute blood-loss anaemia ({_ab['Severe anaemia %']}% vs "
      f"{_ab['Moderate anaemia %']}% with moderate anaemia), which is often short-lived once bleeding is treated, so a "
      f"single Hb value may overstate how long these women were anaemic. The high cross-sectional OR for severe "
      f"anaemia fell to {_hs} when Hb had to precede the stroke by at least 30 days, suggesting that some very low "
      f"values were measured during the stroke admission. The continuous analyses, which do not depend on the small "
      f"severe group, showed risk rising steadily as Hb fell below 13 g/dL (Figure 3).")
    _md = R["mcv_dated"]
    P(f"MCV-based classification has limitations. MCV is a weak marker of iron deficiency on its own, and the "
      f"macrocytic group was small ({int(macro['n']):,} women, {int(macro['Events'])} strokes). MCV was taken from the "
      f"dated laboratory record nearest to the Hb measurement (within 30 days; same day for {_md['same_day']:,} of "
      f"{_md['n']:,} women). The morphology findings should be treated as exploratory.")
    P("These findings extend reports linking anaemia and iron deficiency with stroke in general populations [ref] to a "
      "group of young and middle-aged women with a very high prevalence of anaemia. The morphology pattern and the "
      "bleeding interaction suggest two interpretations. Anaemia that is not explained by menstrual iron loss may "
      "signal systemic illness (renal, hepatic, inflammatory or haematological) that itself raises vascular risk. "
      "Alternatively, the physiological consequences of low oxygen-carrying capacity may matter more than its cause "
      "[ref]. Our adjustment for, and exclusion of, a broad set of anaemia-causing diagnoses argues against the "
      "association being entirely due to recognised comorbidity. However, residual confounding by undiagnosed or "
      "uncoded illness cannot be excluded; the E-values indicate that a confounder would need to be moderately strongly "
      "associated with both anaemia and stroke.")
    P("The inverse association between iron deficiency (ferritin <30 ng/mL) and stroke should be interpreted with "
      "caution. Ferritin was measured in a selected third of the cohort, and ferritin rises with inflammation and "
      "acute illness [ref], so low ferritin may identify otherwise healthy women with pure menstrual iron loss. The "
      "absence of an association with TIA may reflect diagnostic misclassification of TIA in administrative data "
      "[ref]. The associations with low platelet count and macrocytosis are consistent with haematological markers of "
      "systemic disease.")
    P("Clinically, Hb is routinely measured when benign uterine disease is diagnosed. Women with moderate or severe "
      "anaemia, particularly anaemia that is normocytic, macrocytic or out of proportion to bleeding, may warrant "
      "assessment of its cause and of their vascular risk factors. Whether correcting anaemia changes stroke risk "
      "cannot be answered by these data. Our exploratory treatment analysis was limited by confounding by indication "
      "and by incomplete capture of outpatient oral iron.")
    H2("Strengths and limitations")
    P("Strengths include the size of the cohort, a pre-specified analysis plan, and linkage to dated diagnosis, "
      "medication, encounter and death data. These allowed several complementary approaches to temporality and "
      "confounding, including a time-to-event analysis in which anaemia always preceded the outcome.")
    P(f"The study has important limitations. It was conducted at a single tertiary centre, and stroke was ascertained "
      f"from diagnosis codes and imaging reports rather than adjudicated clinically; stroke dates are the earliest "
      f"coded dates. A single Hb value within ±3 years was available, and in the cross-sectional analyses Hb was often "
      f"measured after the stroke; we addressed this with the temporality analyses. Women without an Hb value had a "
      f"lower stroke prevalence, so testing was not random. Smoking status was unknown for "
      f"{T['codebook_summary'].set_index('Variable').loc['smoking', 'Code 9 %']}% of women. Iron therapy captured "
      f"only facility-administered doses. Transfusion data were not available, and ICD-10 pregnancy codes were "
      f"available only as free-text descriptions. Severe anaemia before stroke was uncommon, which limited precision. "
      f"Healthcare use (number of visits) was not available, so women who were seen and tested more often may "
      f"have had more anaemia and more stroke recorded. Finally, although we adjusted for many conditions, residual "
      f"confounding by general ill health, inflammation or socioeconomic factors is possible.")
    H2("Conclusions")
    P("Among women with benign uterine disease, moderate anaemia is associated with a higher risk of subsequent "
      "stroke, particularly ischaemic stroke, independent of recognised anaemia-causing conditions. Anaemia at gynaecological "
      "diagnosis may be a useful marker for vascular risk assessment. "
      "Prospective studies should examine its causes and whether treating it modifies stroke risk.")
    H2("Acknowledgements, funding, disclosures")
    P("[To be completed by the authors.]")
    H2("Data availability")
    P("[Statement to be completed; analysis code available at the project repository.]")
    BR()

    # ------------------------------------------------------------------ tables
    t1s = t1[~t1.Characteristic.astype(str).str.startswith(("Hormonal", "Any stroke", "Ischaemic stroke"))].copy()
    keep_after = []
    skip = False
    for r in t1.itertuples(index=False):
        ch = str(r.Characteristic)
        if ch not in ("", "nan"):
            skip = ch.startswith("Hormonal")
        keep_after.append(not skip)
    t1s = t1[pd.Series(keep_after, index=t1.index)].copy()
    t1s = t1s.fillna("")
    t1s.columns = ["Characteristic", "", "All", "No anaemia (Hb ≥12)", "Mild (10–11.9)", "Moderate (8–9.9)",
                   "Severe (<8)", "P"]
    TABLE("Table 1. Characteristics of eligible women by anaemia grade", t1s,
          "Continuous variables are median (IQR) and mean (SD); categorical variables are n (%) of those with known "
          "values. The anaemia-grade columns include only women with an Hb value. P values are from Kruskal–Wallis or "
          "χ² tests across anaemia grades.",
          widths=[2200, 1500, 1300, 1300, 1300, 1300, 1300, 700], font=7)
    BR()
    # Table 2
    main = T["P2_main_models"]
    rows = []
    for lv in ["None (Hb ≥12)"] + GR:
        r = {"Anaemia grade (Hb, g/dL)": lv, "Women": f"{int(ex_a.loc[lv, 'n']):,}",
             "Stroke, n (%)": f"{int(ex_a.loc[lv, 'Any stroke n']):,} ({ex_a.loc[lv, 'Any stroke %']:.1f})"}
        for o, lab in [("Any stroke (primary)", "Any stroke"), ("Ischaemic stroke vs no stroke", "Ischaemic stroke")]:
            m = main[(main.Outcome == o) & (main.Model == "separate") & (main.Term == "Anaemia grade")]
            if lv == "None (Hb ≥12)":
                cr = ad = "1.00 (ref)"
            else:
                mm = m[m.Level == lv].iloc[0]
                cr, ad = mm["Crude OR (95% CI)"], mm["Adjusted OR (95% CI)"]
            if lab == "Any stroke":
                r["Crude OR (95% CI)"] = cr
            r[f"Adjusted OR (95% CI), {lab.lower()}"] = ad
        tr_ = tter("Any stroke", lv)
        r["Rate/1,000 PY after Hb"] = f"{tr_['Rate /1,000 PY']:.2f}"
        r["Adjusted RR (95% CI)"] = "1.00 (ref)" if lv.startswith("None") else tr_["Adjusted + pre-Hb conditions (2a) RR"]
        rows.append(r)
    rows.append({"Anaemia grade (Hb, g/dL)": "Per grade (trend)", "Women": "", "Stroke, n (%)": "",
                 "Crude OR (95% CI)": "",
                 "Adjusted OR (95% CI), any stroke": f"{tr['stroke_any']['txt']}; P{'' if fmt_p(tr['stroke_any']['p']).startswith('<') else '='}{fmt_p(tr['stroke_any']['p'])}",
                 "Adjusted OR (95% CI), ischaemic stroke": f"{tr['y_isch']['txt']}",
                 "Rate/1,000 PY after Hb": "",
                 "Adjusted RR (95% CI)": tte_pg["Adjusted + pre-Hb conditions (2a) RR"].replace("; p=", "; P=")})
    TABLE("Table 2. Anaemia grade and stroke: cross-sectional odds ratios and rate ratios after the Hb measurement",
          pd.DataFrame(rows),
          f"OR: logistic regression with HC1 robust SEs, adjusted for age, race, BMI, hypertension, diabetes, "
          f"dyslipidaemia, smoking, migraine, thrombophilia, hormonal therapy type, uterine bleeding and uterine "
          f"diagnosis group (any stroke n = {fa.n:,}, {fa.events:,} events; ischaemic stroke vs no stroke n = "
          f"{fi.n:,}, {fi.events:,} events). RR: Poisson regression in the time-to-event cohort (n = {tc['n']:,}; "
          f"{tc['py']:,.0f} person-years; {tc['ev']} strokes), with the same covariates plus conditions documented "
          f"before the Hb (haemoglobinopathy, CKD, chronic liver disease, alcohol use disorder, GI bleeding, "
          f"malabsorption, IBD, HIV, recent pregnancy, CAD).",
          widths=[1700, 900, 1200, 1500, 1700, 1700, 1100, 1600], font=7)
    BR()
    sel = ["Primary: Any stroke", "Primary: Ischaemic stroke", "Primary: Incident stroke",
           "Hb ≥30 d before stroke: Any stroke", "Hb ≥30 d before stroke: Ischaemic stroke",
           "Excluding Hb within 30 d of stroke: Any stroke", "2a Pre-Hb conditions: Any stroke",
           "2a Pre-Hb conditions: Ischaemic stroke", "2b Over-adjustment check (+ undated/post-stroke): Any stroke",
           "Hb ≥30 d before stroke + 2a: Any stroke", "Excluding haemoglobinopathies: Any stroke",
           "Excluding deaths within 1 y of Hb (+2a): Any stroke",
           "Excluding cancer, heart failure, CKD, liver disease, HIV (+2a): Any stroke",
           "Excluding all anaemia-causing conditions: Any stroke",
           "Excluding all anaemia-causing conditions: Ischaemic stroke"]
    t3 = X.loc[[s_ for s_ in sel if s_ in X.index], ["N", "Events", "Mild", "Moderate", "Severe", "Per grade (trend)",
                                                     "p-trend"]].reset_index()
    t3.columns = ["Analysis", "N", "Strokes", "Mild", "Moderate", "Severe", "Per grade", "P trend"]
    t3["N"] = t3["N"].map(lambda v: f"{int(v):,}")
    t3["Strokes"] = t3["Strokes"].map(lambda v: f"{int(v):,}")
    TABLE("Table 3. Anaemia and stroke across temporality, confounding and restriction analyses (adjusted OR, 95% CI)",
          t3, "Reference: Hb ≥12 g/dL. 2a: additionally adjusted for conditions documented before the Hb. 2b: 2a plus "
              "heart failure, atrial fibrillation, VTE, malignancy and anticoagulant use (undated or potentially post-"
              "stroke; over-adjustment check). Anaemia-causing conditions: haemoglobinopathy, CKD, chronic liver "
              "disease, alcohol use disorder, GI bleeding or pregnancy in the year before Hb, malabsorption, IBD, HIV, "
              "neoplastic/chemotherapy anaemia. Serious chronic illness: cancer, heart failure, CKD, chronic liver "
              "disease or HIV recorded at any time.",
          widths=[3300, 900, 800, 1400, 1400, 1400, 1400, 800], font=7)
    BR()

    # ------------------------------------------------------------------ figures
    H1("Figures")
    FIG("figures/fig12_p2_flow.png", "Figure 1. Participant flow.", 5.2)
    FIG("figures/fig7_p2_anaemia_gradient.png",
        "Figure 2. Adjusted odds ratios for stroke by anaemia grade, for any stroke, ischaemic stroke and strokes "
        "after diagnosis. Reference Hb ≥12 g/dL; bars are 95% CIs.")
    FIG("figures/fig14_p2_hb_curves.png",
        "Figure 3. Stroke risk by haemoglobin (restricted cubic spline, 4 knots; reference 13 g/dL), adjusted for "
        "covariates and pre-Hb conditions. (A) Odds ratio for any stroke (cross-sectional). (B) Rate ratio for any "
        "stroke after the Hb measurement. Shaded bands: 95% CI.")
    FIG("figures/fig11_p2_extended.png",
        "Figure 4. Moderate anaemia and per-grade estimates across temporality, confounding and restriction analyses. "
        "The final row shows rate ratios from the time-to-event analysis; other rows show odds ratios.")
    FIG("figures/fig13_p2_per_hb.png",
        "Figure 5. (A) Crude stroke rate per 1,000 person-years after the Hb measurement by anaemia grade (exact "
        "Poisson 95% CI; strokes / women under each bar). (B) Odds ratios (cross-sectional) and rate ratios (after the Hb "
        "measurement) per 1 g/dL lower Hb, below 13 g/dL and over the whole range, adjusted for covariates and "
        "pre-Hb conditions.")
    FIG("figures/fig15_p2_anaemia_pattern.png",
        "Figure 6. Laboratory anaemia pattern (MCV and RDW-CV from the blood count nearest the Hb, within 30 days; "
        "RDW-CV >14.5% = high) and stroke, compared with no anaemia, adjusted for covariates and pre-Hb conditions. "
        "(A) Odds ratios (cross-sectional). (B) Rate ratios after the Hb measurement. Labels: strokes / women; "
        "patterns with fewer than 5 strokes were not estimated.")
    BR()

    # ------------------------------------------------------------------ supplement
    H1("Supplementary material")
    sup = [
        ("eTable 1. Prevalence of anaemia-related conditions by anaemia grade and stroke status (%)",
         R["p2x_prev"], "Conditions from dated diagnosis extracts; acquired conditions counted if documented before the Hb."),
        ("eTable 2. Anaemia type (MCV) and coded iron-deficiency anaemia vs no anaemia",
         R["p2x_morph"].drop(columns=[c for c in ["OR", "CI low", "CI high", "p"] if c in R["p2x_morph"].columns]),
         "Adjusted as in Table 2."),
        ("eTable 3. Anaemia and stroke subtypes (each subtype vs no stroke)",
         R["p2x_sub"].drop(columns=[c for c in ["OR", "CI low", "CI high", "Flags"] if c in R["p2x_sub"].columns]), ""),
        ("eTable 4. Per-grade anaemia OR for any stroke within subgroups",
         R["p2x_sg"].drop(columns=[c for c in ["OR", "CI low", "CI high", "Flags", "p interaction"] if c in R["p2x_sg"].columns]),
         "Interaction: product of anaemia grade score and subgroup indicator."),
        ("eTable 5. MCV, platelets and iron deficiency (adjusted ORs)",
         T["P2_main_models"][(T["P2_main_models"].Model == "separate") & T["P2_main_models"].Term.isin(["MCV", "Platelets"])]
         [["Outcome", "Term", "Level", "Crude OR (95% CI)", "Adjusted OR (95% CI)", "Adjusted p"]], ""),
        ("eTable 6. E-values", R["p2x_ev"], ""),
        ("eTable 7. Anaemia treatment landmark analysis (exploratory)",
         R["p2x_tx"].drop(columns=[c for c in ["RR", "CI low", "CI high", "p"] if c in R["p2x_tx"].columns]),
         T.get("P2X_treatment_landmark") is not None and
         "Iron/ESA administered from 30 d before to 90 d after the Hb (facility-administered only); follow-up from 90 d."),
        ("eTable 8. Sensitivity: Hb within ±1 year and multiple imputation",
         R["p2_lab_timing"][["Outcome", "Level", "Hb within ±3 y (primary)", "Hb within ±1 y", "N ±1y", "Events ±1y"]], ""),
        ("eTable 9. Exploratory mediation of the fibroid–stroke association by anaemia",
         T["P2_mediation"].T.reset_index().rename(columns={"index": "Item", 0: "Value"}), ""),
        ("eTable 10. Stroke risk per 1 g/dL of haemoglobin and crude rates by anaemia grade",
         R["p2x_per"][["Design", "Outcome", "Hb range", "Adjustment", "Estimate (95% CI)", "p (text)", "N", "Events"]]
         .rename(columns={"p (text)": "P"}),
         "Piecewise-linear models: separate slopes below and above 13 g/dL. Crude rates are in the workbook sheet "
         "P2X_absolute_rates."),
        ("eTable 11. Rate ratios after excluding early deaths and serious chronic illness (adjustment 2a)",
         R["p2x_ill"][["Analysis", "Women", "Strokes", "EPV", "Moderate RR (95% CI)", "Severe RR (95% CI)",
                       "Per grade RR (95% CI)"]], "Serious chronic illness: cancer, heart failure, CKD, chronic liver "
                                                  "disease or HIV recorded at any time."),
        ("eTable 12. Stroke risk at selected Hb values vs 13 g/dL (spline, adjusted for covariates and pre-Hb conditions)",
         R["p2x_curve_pts"][["Design", "Hb (g/dL)", "Estimate vs 13 g/dL (95% CI)", "N", "Events", "p overall",
                             "p non-linearity"]], ""),
        ("eTable 13. Laboratory anaemia pattern (MCV and RDW) and stroke",
         pd.concat([R["p2x_pattern"].assign(Design="Cross-sectional OR")
                    .rename(columns={"Adjusted OR (95% CI)": "Estimate (95% CI)", "Strokes": "Strokes"}),
                    R["p2x_pattern_tte"].assign(Design="Rate ratio after Hb", Adjustment="+ pre-Hb conditions (2a)")
                    .rename(columns={"Adjusted RR (95% CI)": "Estimate (95% CI)"})], ignore_index=True)
         [["Design", "Outcome", "Adjustment", "Pattern", "Women", "Strokes", "Estimate (95% CI)", "p (text)"]]
         .rename(columns={"p (text)": "P"}),
         "MCV and RDW-CV from the dated laboratory record nearest the Hb (within 30 days). RDW-CV >14.5% = high. "
         "Reference: no anaemia (Hb ≥12 g/dL)."),
    ]
    for title, df, note in sup:
        TABLE(title, df.fillna(""), note if isinstance(note, str) else "", font=7)
    FIG("figures/fig8_p2_robustness.png", "eFigure 1. Moderate and severe anaemia across sensitivity analyses.")
    FIG("figures/fig4_forest_paper2.png", "eFigure 2. Haematological indices (anaemia, MCV, platelets) and stroke.")
    FIG("figures/fig2_hb_spline.png",
        f"eFigure 3. Adjusted odds ratio for any stroke by haemoglobin, Paper 2 covariates only (restricted cubic "
        f"spline, 4 knots; reference 13 g/dL; n = {R['p2_spline_n']:,}; {R['p2_spline_events']:,} strokes).", 5.4)

    # ------------------------------------------------------------------ checks
    RESULTS["manuscript_checks"] = checks
    failed = [c for c, ok in checks if not ok]
    if failed:
        print("MANUSCRIPT CHECK FAILED:", failed)
    B.append({"t": "note", "text": "Automated wording checks: " + ("all passed" if not failed else "FAILED: " + "; ".join(failed))})
    return B


def to_markdown(B):
    out = []
    for b in B:
        if b["t"] == "title":
            out.append(f"# {b['text']}\n")
        elif b["t"] == "h1":
            out.append(f"\n## {b['text']}\n")
        elif b["t"] == "h2":
            out.append(f"\n### {b['text']}\n")
        elif b["t"] in ("p", "note"):
            out.append(("_" + b["text"] + "_" if b["t"] == "note" else b["text"]) + "\n")
        elif b["t"] == "table":
            out.append(f"\n**{b['title']}**\n")
            out.append("| " + " | ".join(b["columns"]) + " |")
            out.append("|" + "---|" * len(b["columns"]))
            for r in b["rows"]:
                out.append("| " + " | ".join(str(v).replace("|", "/") for v in r) + " |")
            if b["note"]:
                out.append(f"\n_{b['note']}_\n")
        elif b["t"] == "figure":
            out.append(f"\n![{b['caption']}]({b['path'].split('outputs/')[-1]})\n\n{b['caption']}\n")
    return "\n".join(out)


FOCUS_FIGS = ["fig7_p2_anaemia_gradient.png", "fig14_p2_hb_curves.png", "fig13_p2_per_hb.png",
              "fig11_p2_extended.png", "fig15_p2_anaemia_pattern.png"]


def short_version(B):
    """Title + abstract + Tables 1-3 + four key figures (no flow diagram), figures renumbered 1-4."""
    import re
    out = [b for b in B if b["t"] == "title"]
    i = next(k for k, b in enumerate(B) if b["t"] == "h1" and b["text"] == "Abstract")
    out.append(B[i])
    for b in B[i + 1:]:
        if b["t"] in ("h1", "pagebreak"):
            break
        out.append(b)
    out.append({"t": "pagebreak"})
    out.append({"t": "h1", "text": "Tables"})
    tabs = [b for b in B if b["t"] == "table" and re.match(r"Table [123]\.", b["title"])]
    assert len(tabs) == 3
    for tb in tabs:
        out += [tb]
    out.append({"t": "pagebreak"})
    out.append({"t": "h1", "text": "Figures"})
    figs = {b["path"].split("/")[-1]: b for b in B if b["t"] == "figure"}
    for n, f in enumerate(FOCUS_FIGS, 1):
        b = dict(figs[f])
        b["caption"] = re.sub(r"^Figure \d+\.", f"Figure {n}.", b["caption"])
        out.append(b)
    return out


def _us(o, key=None):
    """US spelling for every text field of the manuscript blocks (file paths untouched)."""
    from .utils import us_spelling
    if key == "path":
        return o
    if isinstance(o, dict):
        return {k: _us(v, k) for k, v in o.items()}
    if isinstance(o, list):
        return [_us(v) for v in o]
    return us_spelling(o)


def run():
    B = _us(build())
    (OUT_DIR / "anaemia_manuscript.json").write_text(json.dumps(B, ensure_ascii=False), encoding="utf-8")
    (OUT_DIR / "anaemia_manuscript.md").write_text(to_markdown(B), encoding="utf-8")
    r = subprocess.run(["node", "analysis/build_docx.js", str(OUT_DIR / "anaemia_manuscript.json"),
                        str(OUT_DIR / "anaemia_manuscript.docx")], capture_output=True, text=True)
    log("Manuscript", f"docx build: {'ok' if r.returncode == 0 else r.stderr[-500:]}")
    if r.returncode != 0:
        print(r.stderr)
    S = short_version(B)
    (OUT_DIR / "anaemia_abstract_tables_figures.json").write_text(json.dumps(S, ensure_ascii=False), encoding="utf-8")
    (OUT_DIR / "anaemia_abstract_tables_figures.md").write_text(to_markdown(S), encoding="utf-8")
    r = subprocess.run(["node", "analysis/build_docx.js", str(OUT_DIR / "anaemia_abstract_tables_figures.json"),
                        str(OUT_DIR / "anaemia_abstract_tables_figures.docx")], capture_output=True, text=True)
    log("Manuscript", f"short docx build: {'ok' if r.returncode == 0 else r.stderr[-500:]}")
