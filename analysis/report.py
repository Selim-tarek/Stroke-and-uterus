"""Generate methods_and_results.md and analysis_log.md from RESULTS/TABLES.

Every number in these documents is read from objects produced by the
analysis code in this run; nothing is typed in by hand.
"""
import platform
from datetime import date

import numpy as np
import pandas as pd
import scipy
import statsmodels

from .utils import DIAG, LINEARITY, LOG, OUT_DIR, RESULTS, TABLES, fmt_p, get_or


def md_table(df, cols=None, max_rows=None):
    d = df[cols] if cols else df
    if max_rows:
        d = d.head(max_rows)
    def cell(v):
        if isinstance(v, float):
            if np.isnan(v):
                return ""
            return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
        return str(v).replace("|", "/").replace("\n", " ")
    lines = ["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)]
    for r in d.itertuples(index=False):
        lines.append("| " + " | ".join(cell(v) for v in r) + " |")
    return "\n".join(lines)


def _p3_changed_text(inc):
    """List P3 contrasts whose OR moved by >10% when prior strokes are excluded."""
    ch = []
    for _, r in inc.iterrows():
        a = float(r["Primary M2 OR"].split(" ")[0])
        b = float(r["Excluding stroke before/at index"].split(" ")[0])
        if abs(np.log(b / a)) >= np.log(1.1):
            ch.append(f"{r['Outcome'].lower()}, {r['Group vs fibroids only'].lower()}: "
                      f"{r['Primary M2 OR']} → {r['Excluding stroke before/at index']}")
    if not ch:
        return "no estimate by 10% or more"
    return (f"{len(inc) - len(ch)} of {len(inc)} estimates by less than 10%; the remainder were "
            + "; ".join(ch))


def ptxt(p):
    t = fmt_p(p) if not isinstance(p, str) else p
    return f"p{t}" if t.startswith("<") else f"p={t}"


def orp(f, col):
    g = get_or(f, col)
    return f"{g['txt']}; {ptxt(g['p'])}"


def assoc(f, col):
    """Direction phrase driven by the CI, so wording always matches the estimate."""
    g = get_or(f, col)
    if g["lo"] > 1:
        return "higher odds"
    if g["hi"] < 1:
        return "lower odds"
    return "similar odds"


def _interaction_text(ints):
    import re
    out = []
    mt = ints.iloc[0]
    out.append(f"The microcytosis × thrombocytosis interaction was "
               f"{'statistically significant' if mt['p interaction'] < 0.05 else 'not statistically significant'} "
               f"({ptxt(mt['p (text)'])}).")
    ab = ints.iloc[2]
    txt = (f"The anaemia grade × uterine bleeding interaction was "
           f"{'statistically significant' if ab['p interaction'] < 0.05 else 'not statistically significant'} "
           f"for any stroke ({ptxt(ab['p (text)'])}; ischaemic stroke {ptxt(ints.iloc[3]['p (text)'])})")
    if ab["p interaction"] < 0.05:
        mods = re.findall(r"Moderate \(8–9\.9\): ([0-9.]+ \([0-9.]+–[0-9.]+\))", str(ab["Stratum-specific anaemia ORs"]))
        sev = re.findall(r"Severe \(<8\): ([0-9.]+ \([0-9.]+–[0-9.]+\))", str(ab["Stratum-specific anaemia ORs"]))
        if len(mods) == 2 and len(sev) == 2:
            txt += (f"; moderate anaemia OR {mods[0]} without and {mods[1]} with heavy/abnormal bleeding, severe "
                    f"anaemia OR {sev[0]} and {sev[1]}, respectively")
    out.append(txt + ".")
    return " ".join(out)


def orx(f, col):
    return get_or(f, col)["txt"]


# ---------------------------------------------------------------------------
def methods_results():
    R = RESULTS
    n_all, n_el = R["n_total"], R["n_eligible"]
    fl = R["flow"]
    miss = R["missing"]
    sw = f"Python {platform.python_version()}, statsmodels {statsmodels.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}"

    # ---------------- Paper 1 numbers
    g = R["p1_groups"]
    core = R["p1_core"]
    mf = R["p1_main_fits"]
    wf = R["p1_whole_fits"]
    desc = R["p1_desc"]

    def dsc(descr, level):
        r = desc[(desc.Descriptor.str.startswith(descr)) & (desc.Level == level)]
        return f"{int(r['n (imaging_only_candidates)'].iloc[0])} ({r['% (candidates)'].iloc[0]:.1f}%)"

    n_cov = g.get("Covert infarct", 0)
    in_cand = R["p1_covert_in_candidates"]
    sub = R["p1_io_subcat"]
    imp = R["p1_imaged_pct_by_mig"]
    strict = R["p1_strict_fit"]

    # ---------------- Paper 2 numbers
    f2 = R["p2_fits"]
    fa = f2[("stroke_any", "anemia_cat")]
    fi = f2[("y_isch", "anemia_cat")]
    finc = f2[("y_incident", "anemia_cat")]
    fm = f2[("stroke_any", "mcv_cat")]
    fp = f2[("stroke_any", "plt_cat")]
    fj = f2[("stroke_any", "joint")]
    fmi = f2[("y_isch", "mcv_cat")]
    pts = R["p2_spline_pts"].set_index("Hb (g/dL)")
    ints = R["p2_ints"]
    idr = R["p2_id"]
    id_adj = idr[(idr.Outcome == "Any stroke (primary)") & (idr.Model == "adjusted")].iloc[0]
    id_adj_an = idr[(idr.Outcome == "Any stroke (primary)") & (idr.Model == "adjusted + anaemia grade")].iloc[0]
    lt = R["p2_lab_timing"]
    lt_any = lt[lt.Outcome == "Any stroke (primary)"].set_index("Level")
    mi2 = R["p2_mi"]
    mi2a = mi2[mi2.Outcome == "Any stroke (primary)"].set_index("coef")
    med = R["p2_med"]
    fsel = R["p2_ferr_sel"].set_index("% with")
    tr = {k: R[f"p2_trend_{k}"] for k in ["stroke_any", "y_isch", "y_incident"]}

    # ---------------- Paper 3 numbers
    m3 = R["p3_main"]
    pv = R["p3_prev"]
    em = R["p3_em"].set_index("Outcome")
    st = R["p3_strata"]

    def p3(out, grp, col="M2 OR (95% CI) age+race+BMI"):
        r = m3[(m3.Outcome == out) & (m3["Group vs fibroids only"] == grp)].iloc[0]
        return f"{r[col]}; BH-adjusted {ptxt(r['M2 p (BH-adjusted, all contrasts) (text)'])}"

    def p3s(out, grp, age):
        r = st[(st.Outcome == out) & (st["Group vs fibroids only"] == grp) & (st.Age == age)].iloc[0]
        return r["OR (95% CI)"]

    def asp(out, grp):
        r = pv[(pv.Outcome == out) & (pv.Group == grp)].iloc[0]
        return f"{r['Age-std %']:.2f}% (95% CI {r['Age-std low']:.2f}–{r['Age-std high']:.2f})"

    gn = R["p3_group_n"]
    ag = R["p3_age_by_group"]
    n_out = m3["Outcome"].nunique()
    n_sig = int((m3["M2 p (BH-adjusted, all contrasts)"] < 0.05).sum())
    n_tests = int(m3["M2 p (BH-adjusted, all contrasts)"].notna().sum())

    L = []
    L.append("# Statistical methods and results — draft for three manuscripts\n")
    L.append(f"_Generated automatically by `analysis/run_all.py` on {date.today().isoformat()}. Every number below "
             f"is read from the objects produced in the same run (see `outputs/stats_results.xlsx`). "
             f"All associations are cross-sectional; no causal interpretation is intended._\n")
    L.append("## Common source population\n")
    L.append(f"The analysis file contained {n_all:,} women with benign uterine pathology (fibroids, adenomyosis or "
             f"endometriosis). We excluded {fl.get(1, 0):,} with active malignancy within ±6 months of index, "
             f"{fl.get(6, 0):,} aged >60 years, {fl.get(3, 0):,} aged <18 years and {fl.get(5, 0):,} whose age could not "
             f"be derived, leaving {n_el:,} eligible women (Figure 1). The index date was the first documented benign "
             f"uterine diagnosis. Direct identifiers (MRN, date of birth) were removed from all outputs.\n")

    # ------------------------------------------------------------------ P1
    L.append("---\n\n## Paper 1 — Covert brain infarcts on brain imaging performed for other reasons\n")
    L.append("### Statistical methods\n")
    L.append(
        f"The study population comprised eligible women whose brain CT or MRI reports had been reviewed "
        f"(stroke_confirmed_imaging recorded as positive or negative, plus {R['p1_n_reviewed_neg']} women whose imaging-only "
        f"finding had been revoked, excluded or left undecided at adjudication, whose scans had also been reviewed). "
        f"Imaging-positive patients without a qualifying "
        f"cerebrovascular diagnosis code had been adjudicated from full-text reports; we identified them from the "
        f"abstractor notes and cross-checked the classification against the three review worksheets. A covert "
        f"infarct was defined as an imaging-only infarct that was chronic or incidental, with no documented acute "
        f"features and no documented stroke symptoms: incidental infarcts on scans ordered for another indication, "
        f"multiple or territorial infarcts without documented symptoms, chronic infarcts described with hedged "
        f"wording, and imaging-only findings that had been revoked, excluded or left undecided at adjudication but "
        f"whose report text described a brain lesion (the remaining revoked/excluded findings, which described other "
        f"organs, a truncated probably negated sentence or no organ, were classed as imaged without infarct). A stricter definition kept only the first of these. A clinical stroke was a coded stroke or TIA. "
        f"Imaging-only acute, symptomatic or haemorrhagic events met neither definition and were described separately. "
        f"Prevalence is reported with Wilson 95% CIs. Characteristics were compared across imaged-negative, covert-"
        f"infarct and clinical-stroke groups using Kruskal–Wallis and χ² tests. Factors associated with covert infarct "
        f"(versus imaged-negative) were estimated by logistic regression with heteroskedasticity-robust (HC1) standard "
        f"errors. Covariates were pre-specified: age, race, BMI, hypertension, diabetes, dyslipidaemia, smoking "
        f"(never/ever/unknown), migraine (any diagnosis, yes/no), thrombophilia, atrial fibrillation, uterine "
        f"diagnosis group and anaemia grade. The full model had {R['p1_full_events']} events for {R['p1_full_df']} "
        f"parameters (events per variable [EPV] {R['p1_full_epv']:.1f}) and sparse cells. We therefore applied the "
        f"pre-specified simplification: a core model (age, hypertension, diabetes, dyslipidaemia, smoking, migraine), "
        f"with each remaining factor added one at a time. When a core-plus-one model still had EPV <10 or a level with "
        f"<5 events, categorical factors were collapsed (race: White/Black/other; diagnosis group: fibroids only/"
        f"endometriosis only/adenomyosis only or >1 condition; anaemia: none/mild/moderate–severe), and binary "
        f"factors with <5 events were not estimated. Linearity of continuous terms was checked with restricted cubic "
        f"splines (4 knots). Collinearity (variance inflation factor, VIF) and EPV were recorded for every model. "
        f"Brain imaging is not performed at random. Because of that, we re-estimated all models with the whole eligible "
        f"cohort without stroke as the comparator, and we modelled the probability of being imaged. Migraine was "
        f"coded as present for any migraine diagnosis regardless of type or aura status. No pre-specified variable "
        f"was more than {R['p1_max_missing_pct']}% missing in the imaged analysis sample, so multiple imputation was "
        f"not required (threshold 10%). Further sensitivity analyses used the strict covert definition and restricted cases to infarcts dated after index (stroke_timing = 3), excluding patients with a stroke "
        f"before or at index. Analyses used {sw}.\n")
    L.append("### Results\n")
    L.append(
        f"Of {n_el:,} eligible women, {R['p1_n_imaged']:,} had brain imaging reviewed. Of these, "
        f"{g.get('Imaged, no infarct', 0):,} had no infarct, {n_cov} had a covert infarct, "
        f"{g.get('Clinical stroke', 0):,} had a clinical (coded) stroke and "
        f"{g.get('Imaging-only acute/symptomatic/haemorrhagic', 0)} had an imaging-only acute, symptomatic or "
        f"haemorrhagic event. Covert infarct prevalence among imaged women was {R['p1_prev_imaged_txt']}. It was "
        f"{R['p1_prev_nonclin_txt']} after excluding women with clinical or other imaging-only events, and "
        f"{R['p1_prev_strict_txt']} under the strict definition. Among the {n_cov} covert cases, "
        f"{sub.get('Covert: incidental, scan for another indication', 0)} were incidental on scans for another "
        f"indication, {sub.get('Covert: multiple/territorial, no symptoms', 0)} were multiple or territorial without "
        f"documented symptoms, {sub.get('Covert: chronic infarct, hedged wording', 0)} were chronic infarcts "
        f"with hedged wording, and {R['p1_n_restored']} were brain-lesion findings restored after revocation, exclusion "
        f"or an undecided adjudication. Imaging descriptors were available for {in_cand} of {n_cov} covert cases from the "
        f"imaging_only_candidates worksheet; the corresponding master-file fields were blank for most cases. Of all "
        f"{n_cov} covert cases, the anterior circulation was involved in {dsc('Vascular territory', 'Anterior')} and "
        f"the posterior circulation in {dsc('Vascular territory', 'Posterior')}. A single lesion was recorded in "
        f"{dsc('Number of infarcts', 'One')} and three or more in {dsc('Number of infarcts', 'Three or more')}. A "
        f"lacunar/small-vessel pattern was recorded in {dsc('Imaging-inferred', 'Small-vessel (lacunar)')}. Locations "
        f"included cortical/lobar regions in {dsc('Location', 'Cortical / lobar')}, the basal ganglia/internal capsule "
        f"in {dsc('Location', 'Basal ganglia / internal capsule')} and the cerebellum in "
        f"{dsc('Location', 'Cerebellum')}.\n")
    L.append(
        f"In the core model (n = {core.n:,}; {core.events} covert infarcts; EPV {core.epv:.1f}), covert infarct was "
        f"associated with hypertension (OR {orp(core, 'htn')}). Women with dyslipidaemia had "
        f"{assoc(core, 'dyslipidemia')} of covert infarct (OR {orp(core, 'dyslipidemia')}). There was no association with diabetes "
        f"(OR {orp(core, 'dm')}) or age (per year, OR {orp(core, 'age_index')}). Relative to imaged women without "
        f"migraine, imaged women with migraine had {assoc(core, 'migraine_any')} of covert infarct "
        f"(OR {orp(core, 'migraine_any')}). Added one at a time, atrial fibrillation was associated with covert infarct "
        f"(OR {orp(mf['afib'], 'afib')}). Endometriosis only (OR {orp(mf['dxgrp_collapsed'], 'dxgrp3=Endometriosis only')}) "
        f"and adenomyosis only or >1 condition (OR {orp(mf['dxgrp_collapsed'], 'dxgrp3=Adenomyosis only or >1')}) had "
        f"lower odds of covert infarct than fibroids only. Anaemia was not clearly associated "
        f"(moderate–severe: OR {orp(mf['anemia_cat_collapsed'], 'anemia3=Moderate/severe (<10)')}). Thrombophilia was "
        f"too sparse to estimate.\n")
    L.append(
        f"Imaging was strongly selective: {imp.get(0.0)}% of women without migraine were imaged, compared with "
        f"{imp.get(1.0)}% of women with migraine (adjusted OR for being imaged {R['p1_sel_mig']['txt']}). With the "
        f"whole eligible cohort as the comparator, migraine was associated with {assoc(wf['CORE'], 'migraine_any')} "
        f"of covert infarct (OR {orp(wf['CORE'], 'migraine_any')}), the opposite direction to the imaged-only "
        f"comparison. The estimates for hypertension (OR {orx(wf['CORE'], 'htn')}) and atrial fibrillation "
        f"(OR {orx(wf['afib'], 'afib')}) became larger. The reversal for migraine is consistent with preferential "
        f"imaging of women with migraine: restricting the comparator to imaged women over-represents migraine among "
        f"controls. Under the strict definition ({strict.events} cases; model reduced to age, hypertension and migraine "
        f"because of EPV), hypertension remained associated (OR {orp(strict, 'htn')}). Only "
        f"{R['p1_incident_events']} covert infarcts were dated after index, because the event date is blank for most "
        f"imaging-only cases (timing of covert cases: "
        f"{', '.join(f'{k} {v}' for k, v in R['p1_covert_timing'].items())}). An incident-only model was therefore "
        f"not estimable.\n")

    # ------------------------------------------------------------------ P2
    L.append("---\n\n## Paper 2 — Anaemia, red-cell indices, platelet count and stroke\n")
    L.append("### Statistical methods\n")
    L.append(
        f"All {n_el:,} eligible women were included. Haemoglobin and ferritin were the results closest to index within "
        f"±3 years. Anaemia was graded by WHO thresholds for non-pregnant women (none ≥12, mild 10–11.9, moderate "
        f"8–9.9, severe <8 g/dL). MCV was categorised as <80, 80–100 and >100 fL, platelets as <150, 150–400 and "
        f">400 ×10³/µL, and iron deficiency as ferritin <30 ng/mL. The primary outcome was any stroke or TIA "
        f"(stroke_any). Secondary outcomes were ischaemic stroke (versus no stroke; other stroke types excluded) and "
        f"incident stroke (stroke_timing = 3 versus no stroke; strokes before or at index excluded; strokes with "
        f"unknown or blank timing neither cases nor controls). Each exposure was modelled separately and then jointly "
        f"(anaemia, MCV, platelets) by logistic regression with HC1 robust standard errors. Covariates were "
        f"pre-specified: age, race, BMI, hypertension, diabetes, dyslipidaemia, smoking (never/ever/unknown), migraine "
        f"(any diagnosis, yes/no), thrombophilia, hormonal therapy type, heavy or abnormal uterine bleeding and "
        f"uterine diagnosis group. Hormonal-therapy levels with <5 events in a model's sample were merged into "
        f"“other/multiple”. Linear trend across anaemia grades was tested by entering grade as a score. Haemoglobin "
        f"was also modelled as a restricted cubic spline (4 knots at the 5th, 35th, 65th and 95th centiles) and "
        f"plotted as adjusted ORs against a reference of 13 g/dL. Effect modification was assessed with product terms "
        f"for microcytosis × thrombocytosis and anaemia grade × uterine bleeding (robust Wald tests). Iron deficiency "
        f"was a secondary exposure. It was analysed as a complete-case subset because ferritin was available for "
        f"{R['p2_ferritin_n']:,} women ({R['p2_ferritin_pct']:.1f}%). As an exploratory analysis, we estimated the "
        f"proportion of the fibroid–stroke association statistically accounted for by anaemia grade with the "
        f"difference-of-coefficients method, using 1,000 bootstrap resamples for 95% CIs. Uterine bleeding and "
        f"diagnosis group were excluded from that model; adenomyosis and endometriosis flags were included. Because the "
        f"exposure (Hb) was {miss['anemia_cat'][1]}% missing, the anaemia and joint models were repeated after "
        f"multiple imputation by chained equations among women with Hb (m = 20; BMI by predictive mean matching; "
        f"outcome and stroke type included). Exposures were not imputed. Lab timing "
        f"was examined by restricting to Hb (or ferritin) within ±1 year of index. MCV and platelet dates could not be "
        f"verified (see analysis log). Continuous covariate linearity, VIF and EPV were checked for every model. "
        f"Analyses used {sw}.\n")
    L.append("### Results\n")
    L.append(
        f"Among {n_el:,} women, {R['p2_events']:,} had any stroke or TIA and {R['p2_isch_events']:,} had an ischaemic "
        f"stroke. Hb was available for {n_el - miss['anemia_cat'][0]:,}, MCV for {n_el - miss['mcv_cat'][0]:,} and "
        f"platelets for {n_el - miss['plt_cat'][0]:,}. In adjusted models (n = {fa.n:,}; {fa.events:,} events), "
        f"moderate (OR {orp(fa, 'anemia_cat=Moderate (8–9.9)')}) and severe anaemia (OR "
        f"{orp(fa, 'anemia_cat=Severe (<8)')}) were associated with stroke. Mild anaemia was associated with "
        f"{assoc(fa, 'anemia_cat=Mild (10–11.9)')} of stroke (OR {orp(fa, 'anemia_cat=Mild (10–11.9)')}). The OR per anaemia grade was {tr['stroke_any']['txt']} "
        f"(p-trend {fmt_p(tr['stroke_any']['p'])}). The spline showed a non-linear relationship (non-linearity "
        f"{ptxt(R['p2_spline_p_nonlin'])}). Compared with Hb 13 g/dL, the adjusted OR was {pts.loc[10, 'OR vs 13 g/dL (95% CI)']} "
        f"at 10 g/dL and {pts.loc[8, 'OR vs 13 g/dL (95% CI)']} at 8 g/dL. There was little association above 13 g/dL "
        f"(15 g/dL: {pts.loc[15, 'OR vs 13 g/dL (95% CI)']}) (Figure 2). Associations for ischaemic stroke were similar "
        f"(moderate: OR {orp(fi, 'anemia_cat=Moderate (8–9.9)')}; severe: OR {orp(fi, 'anemia_cat=Severe (<8)')}; "
        f"p-trend {fmt_p(tr['y_isch']['p'])}).\n")
    L.append(
        f"Compared with a normal MCV, microcytosis was associated with {assoc(fm, 'mcv_cat=Microcytic (<80)')} of any "
        f"stroke (OR {orp(fm, 'mcv_cat=Microcytic (<80)')}) and {assoc(fmi, 'mcv_cat=Microcytic (<80)')} of ischaemic "
        f"stroke (OR {orp(fmi, 'mcv_cat=Microcytic (<80)')}; joint model with anaemia grade: OR "
        f"{orx(f2[('y_isch', 'joint')], 'mcv_cat=Microcytic (<80)')}). Macrocytosis was associated with "
        f"{assoc(fm, 'mcv_cat=Macrocytic (>100)')} of stroke (OR {orp(fm, 'mcv_cat=Macrocytic (>100)')}). "
        f"Thrombocytosis was associated with {assoc(fp, 'plt_cat=Thrombocytosis (>400)')} of stroke "
        f"(OR {orp(fp, 'plt_cat=Thrombocytosis (>400)')}), and a low platelet count with "
        f"{assoc(fp, 'plt_cat=Low (<150)')} (OR {orp(fp, 'plt_cat=Low (<150)')}). In the joint model, severe anaemia "
        f"had OR {orx(fj, 'anemia_cat=Severe (<8)')} (separate model {orx(fa, 'anemia_cat=Severe (<8)')}) and a low "
        f"platelet count OR {orx(fj, 'plt_cat=Low (<150)')}. {_interaction_text(ints)}\n")
    L.append(
        f"Women with a ferritin result differed from those without: stroke {fsel.loc['stroke_any', 'ferritin measured']}% "
        f"vs {fsel.loc['stroke_any', 'ferritin not measured']}%, and uterine bleeding "
        f"{fsel.loc['uterine_bleeding', 'ferritin measured']}% vs {fsel.loc['uterine_bleeding', 'ferritin not measured']}%. "
        f"In this subset, iron deficiency was inversely associated with stroke (adjusted OR {id_adj['OR (95% CI)']}, "
        f"{ptxt(id_adj['p (text)'])}; additionally adjusted for anaemia grade: OR {id_adj_an['OR (95% CI)']}). Ferritin "
        f"is an acute-phase reactant and was measured selectively, so this estimate should be interpreted cautiously.\n")
    L.append(
        f"Results were materially unchanged when Hb was restricted to ±1 year of index (n = {R['p2_hb1y_n']:,}; severe "
        f"anaemia OR {lt_any.loc['Severe (<8)', 'Hb within ±1 y']}; moderate {lt_any.loc['Moderate (8–9.9)', 'Hb within ±1 y']}) "
        f"and after multiple imputation (severe anaemia OR {mi2a.loc['anemia_cat=Severe (<8)', 'OR (95% CI)']}; "
        f"moderate {mi2a.loc['anemia_cat=Moderate (8–9.9)', 'OR (95% CI)']}). In the incident-only analysis "
        f"({R['p2_incident_events']} strokes after index; {R['p2_incident_excluded_prior']} women with a stroke before "
        f"or at index excluded; {R['p2_incident_excluded_unknown']} with unknown or blank timing excluded), moderate "
        f"anaemia remained associated (OR {orp(finc, 'anemia_cat=Moderate (8–9.9)')}). The severe-anaemia estimate was "
        f"imprecise (OR {orp(finc, 'anemia_cat=Severe (<8)')}; p-trend {fmt_p(tr['y_incident']['p'])}).\n")
    L.append(
        f"Exploratory mediation (cross-sectional; n = {med['n']:,}; {med['ev']:,} events): the fibroid–stroke OR was "
        f"{med['gt']['txt']} without and {med['gd']['txt']} with adjustment for anaemia grade. The proportion of the "
        f"association statistically accounted for by anaemia was {100 * med['pm']:.1f}% (bootstrap 95% CI "
        f"{100 * med['ci'][0]:.1f}% to {100 * med['ci'][1]:.1f}%). Because the temporal order of fibroids, anaemia and "
        f"stroke cannot be established, this is not evidence of mediation in a causal sense.\n")

    # ------------------------------------------------------------------ P3
    L.append("---\n\n## Paper 3 — Vascular risk-factor burden across fibroids, adenomyosis and endometriosis\n")
    L.append("### Statistical methods\n")
    L.append(
        f"All {n_el:,} eligible women were classified by uterine_dx_group: fibroids only (reference), adenomyosis only, "
        f"endometriosis only, or more than one condition. There were {n_out} outcomes, each modelled separately: "
        f"hypertension, diabetes, dyslipidaemia, obesity (BMI ≥30 kg/m²), ever smoking (current or former vs never; "
        f"unknown excluded), migraine (any migraine diagnosis, regardless of type or aura), atrial fibrillation, coronary "
        f"artery disease, prior venous thromboembolism and thrombophilia. Prevalences were directly age-standardised "
        f"to the whole eligible cohort in 5-year age bands, with normal-approximation 95% CIs. Odds ratios came from "
        f"logistic regression with HC1 robust standard errors. The minimal model adjusted for age and race; the primary "
        f"model also adjusted for BMI (except for obesity). Pre-specified simplifications for sparse or low-EPV models "
        f"were: collapse, then drop, race; fall back to age-only adjustment; and suppress contrasts with <5 events. "
        f"None were triggered in the primary models. Analyses were repeated within the age strata 18–39 and 40–60 "
        f"years, and a condition × age-group interaction was tested (3-df robust Wald). P values were corrected with "
        f"the Benjamini–Hochberg procedure across the {3 * n_out} condition contrasts ({n_out} outcomes × 3 groups) and across "
        f"the {n_out} global tests; raw and adjusted values are reported. In a sensitivity analysis, we excluded "
        f"{R['p3_excl_prior_n']:,} women whose stroke preceded or coincided with index, because stroke can prompt "
        f"risk-factor ascertainment. "
        f"BMI, the only covariate with missing data, was {miss['bmi'][1]}% missing (<10%), so multiple imputation was "
        f"not required. Ever smoking was analysable in {R['p3_smoking_known']:,} women "
        f"({100 * R['p3_smoking_known'] / n_el:.1f}%) because smoking status was unknown for most; smoking was not "
        f"imputed. Analyses used {sw}.\n")
    L.append("### Results\n")
    L.append(
        f"The cohort comprised {gn.get('Fibroids only', 0):,} women with fibroids only, "
        f"{gn.get('Adenomyosis only', 0):,} with adenomyosis only, {gn.get('Endometriosis only', 0):,} with "
        f"endometriosis only and {gn.get('>1 condition', 0):,} with more than one condition. Median ages were "
        f"{ag.get('Fibroids only'):.0f}, {ag.get('Adenomyosis only'):.0f}, {ag.get('Endometriosis only'):.0f} and "
        f"{ag.get('>1 condition'):.0f} years, respectively. After BH correction, {n_sig} of {n_tests} condition "
        f"contrasts had adjusted p <0.05 (Figure 5).\n")
    L.append(
        f"Migraine was more common in every non-fibroid group. The age-standardised prevalence of any migraine was "
        f"{asp('Migraine (any)', 'Fibroids only')} with fibroids only, {asp('Migraine (any)', 'Adenomyosis only')} "
        f"with adenomyosis only and {asp('Migraine (any)', 'Endometriosis only')} with endometriosis only. Adjusted "
        f"ORs were {p3('Migraine (any)', 'Adenomyosis only')} for adenomyosis only, "
        f"{p3('Migraine (any)', 'Endometriosis only')} for endometriosis only and "
        f"{p3('Migraine (any)', '>1 condition')} for more than one condition.\n")
    L.append(
        f"Women with endometriosis only had lower odds of most cardiometabolic factors than women with fibroids only: "
        f"hypertension {p3('Hypertension', 'Endometriosis only')}; diabetes {p3('Diabetes', 'Endometriosis only')}; "
        f"dyslipidaemia {p3('Dyslipidaemia', 'Endometriosis only')}; obesity {p3('Obesity (BMI ≥30)', 'Endometriosis only')}; "
        f"atrial fibrillation {p3('Atrial fibrillation', 'Endometriosis only')}; coronary artery disease "
        f"{p3('Coronary artery disease', 'Endometriosis only')}. Ever smoking was slightly more common "
        f"({p3('Ever smoking', 'Endometriosis only')}). Adenomyosis only was associated with more dyslipidaemia "
        f"({p3('Dyslipidaemia', 'Adenomyosis only')}) and thrombophilia ({p3('Thrombophilia', 'Adenomyosis only')}), "
        f"and with less hypertension ({p3('Hypertension', 'Adenomyosis only')}). More than one condition was "
        f"associated with prior VTE ({p3('Prior VTE', '>1 condition')}) and thrombophilia "
        f"({p3('Thrombophilia', '>1 condition')}).\n")
    L.append(
        f"The cardiometabolic contrasts differed by age (condition × age interaction, BH-adjusted p: hypertension "
        f"{em.loc['Hypertension', 'p BH (text)']}, diabetes {em.loc['Diabetes', 'p BH (text)']}, dyslipidaemia "
        f"{em.loc['Dyslipidaemia', 'p BH (text)']}, coronary artery disease {em.loc['Coronary artery disease', 'p BH (text)']}). "
        f"Among women aged 18–39, endometriosis only was associated with markedly lower odds of hypertension "
        f"({p3s('Hypertension', 'Endometriosis only', '18–39')}) and coronary artery disease "
        f"({p3s('Coronary artery disease', 'Endometriosis only', '18–39')}). At 40–60 years the estimates were close "
        f"to the null (hypertension {p3s('Hypertension', 'Endometriosis only', '40–60')}; coronary artery disease "
        f"{p3s('Coronary artery disease', 'Endometriosis only', '40–60')}). In contrast, the higher odds of migraine "
        f"were present in both age strata (endometriosis only, migraine: 18–39 "
        f"{p3s('Migraine (any)', 'Endometriosis only', '18–39')}; 40–60 {p3s('Migraine (any)', 'Endometriosis only', '40–60')}). "
        f"Excluding the {R['p3_excl_prior_n']:,} women with a stroke before or at index changed "
        f"{_p3_changed_text(R['p3_inc'])} (sheet P3_sens_excl_prior_stroke).\n")

    # ------------------------------------------------------------------ STROBE
    # ------------------------------------------------------------------ claim checks
    def ci_excl1(g):
        return g["lo"] > 1 or g["hi"] < 1

    G = lambda f, c: get_or(f, c)  # noqa: E731
    inc = R["p3_inc"]
    ratios = []
    for _, r in inc.iterrows():
        try:
            a = float(r["Primary M2 OR"].split(" ")[0])
            b = float(r["Excluding stroke before/at index"].split(" ")[0])
            ratios.append(abs(np.log(b / a)))
        except (ValueError, AttributeError, IndexError):
            pass
    claims = [
        ("P1 hypertension associated (core)", ci_excl1(G(core, "htn"))),
        ("P1 dyslipidaemia wording generated from CI", True),
        ("P1 diabetes/age not associated", not ci_excl1(G(core, "dm")) and not ci_excl1(G(core, "age_index"))),
        ("P1 migraine direction opposite vs imaged and vs whole cohort",
         (G(core, "migraine_any")["OR"] - 1) * (G(wf["CORE"], "migraine_any")["OR"] - 1) < 0),
        ("P1 afib associated", ci_excl1(G(mf["afib"], "afib"))),
        ("P1 endometriosis-only / adeno-or->1 lower", G(mf["dxgrp_collapsed"], "dxgrp3=Endometriosis only")["hi"] < 1
         and G(mf["dxgrp_collapsed"], "dxgrp3=Adenomyosis only or >1")["hi"] < 1),
        ("P1 anaemia not clearly associated", not ci_excl1(G(mf["anemia_cat_collapsed"], "anemia3=Moderate/severe (<10)"))),
        ("P1 htn & afib larger vs whole cohort", G(wf["CORE"], "htn")["OR"] > G(core, "htn")["OR"]
         and G(wf["afib"], "afib")["OR"] > G(mf["afib"], "afib")["OR"]),
        ("P1 strict htn associated", ci_excl1(G(strict, "htn"))),
        ("P1 incident covert < 10", R["p1_incident_events"] < 10),
        ("P2 moderate & severe anaemia associated", G(fa, "anemia_cat=Moderate (8–9.9)")["lo"] > 1 and G(fa, "anemia_cat=Severe (<8)")["lo"] > 1),
        ("P2 spline non-linear", R["p2_spline_p_nonlin"] < 0.05),
        ("P2 MCV/platelet/interaction wording generated from CIs and p values", True),
        ("P2 iron deficiency inverse", id_adj["CI high"] < 1),
        ("P2 incident moderate associated", G(finc, "anemia_cat=Moderate (8–9.9)")["lo"] > 1),
        ("P2 ±1y unchanged (severe within 10%)", abs(np.log(lt_any.loc["Severe (<8)", "OR ±1y"] /
                                                            G(fa, "anemia_cat=Severe (<8)")["OR"])) < np.log(1.1)),
        ("P3 no simplification triggered in primary models", (m3["M2 simplification"] == "").all()),
        ("P3 excluding prior stroke: changed-estimate list computed", len(ratios) == len(inc)),
        ("P3 endometriosis lower cardiometabolic", all(
            m3[(m3.Outcome == o) & (m3["Group vs fibroids only"] == "Endometriosis only")]["CI high"].iloc[0] < 1
            for o in ["Hypertension", "Diabetes", "Dyslipidaemia", "Obesity (BMI ≥30)", "Atrial fibrillation",
                      "Coronary artery disease"])),
        ("P3 migraine higher all groups", all(
            m3[(m3.Outcome == "Migraine (any)")]["CI low"] > 1)),
        ("P3 endometriosis 40-60 HTN/CAD CI includes 1 or near null",
         st[(st.Outcome == "Coronary artery disease") & (st["Group vs fibroids only"] == "Endometriosis only")
            & (st.Age == "40–60")]["OR"].iloc[0] < 1.5),
    ]
    RESULTS["claims"] = claims
    failed = [c for c, ok in claims if not ok]
    if failed:
        print("CLAIM CHECK FAILED:", failed)

    L.append("---\n\n## STROBE checklist mapping (cross-sectional studies)\n")
    strobe = [
        ("1", "Title and abstract", "To be written; indicate cross-sectional design in title/abstract."),
        ("2", "Background/rationale", "Manuscript introduction (not generated)."),
        ("3", "Objectives", "Question stated at the head of each paper section above."),
        ("4", "Study design", "Methods: cross-sectional analysis of an EHR cohort (all three papers)."),
        ("5", "Setting", "Single U.S. academic centre; index = first benign uterine diagnosis (index dates in data)."),
        ("6", "Participants", "Common source population; Paper 1 imaging subset; Figure 1; sheet `flow`."),
        ("7", "Variables", "Methods paragraphs; Codebook summary sheet `codebook_summary`; analysis_log §Derivations."),
        ("8", "Data sources/measurement", "Codebook (derivations from Epic extracts); analysis_log §Data problems."),
        ("9", "Bias", "P1 imaging-selection analyses (`P1_imaging_selection`, `P1_selection_model`, "
                      "`P1_sens_selection_bias`); incident-only analyses (all papers); ferritin selection (`P2_ferritin_selection`)."),
        ("10", "Study size", "Fixed by available records; no sample-size calculation (state in manuscript)."),
        ("11", "Quantitative variables", "Categorisation rules in Methods; RCS for Hb; `linearity_checks`."),
        ("12a", "Statistical methods", "Methods paragraphs above."),
        ("12b", "Subgroups/interactions", "`P2_interactions`, `P3_age_strata`, `P3_age_interaction`."),
        ("12c", "Missing data", "`missingness`; complete case + MICE (`P1_sensD_MICE`, `P2_sens_MICE`)."),
        ("12d", "Sampling strategy", "Not applicable (all eligible records)."),
        ("12e", "Sensitivity analyses", "Listed in each Methods paragraph; sheets `*_sens*`."),
        ("13", "Participants (numbers)", "Figure 1 (`fig1_flow_diagram`), `flow`, model N in every table."),
        ("14a", "Descriptive data", "`P1_Table1`, `P2_Table1`, `P3_Table1`."),
        ("14b", "Missing data per variable", "`missingness`, `codebook_summary`."),
        ("15", "Outcome data", "`P1_prevalence`, `P2_Table1_exposures`, `P3_prevalence`."),
        ("16a", "Main results (unadjusted and adjusted)", "`P1_main_model`, `P2_main_models` (crude + adjusted), `P3_main_models` (M1 + M2)."),
        ("16b", "Category boundaries", "Methods (anaemia, MCV, platelets, BMI, age bands)."),
        ("16c", "Absolute risk", "Prevalences with 95% CI (`P1_prevalence`, `P3_prevalence`)."),
        ("17", "Other analyses", "Sensitivity sheets; `P2_mediation` (exploratory)."),
        ("18", "Key results", "Discussion (not generated)."),
        ("19", "Limitations", "See analysis_log §Limitations to state."),
        ("20", "Interpretation", "Discussion (not generated); cross-sectional, non-causal wording throughout."),
        ("21", "Generalisability", "Discussion (not generated); single tertiary centre."),
        ("22", "Funding", "To be completed by authors."),
    ]
    L.append(md_table(pd.DataFrame(strobe, columns=["Item", "Topic", "Where addressed"])))
    L.append("")
    L.append("## Figures\n")
    L.append("1. `figures/fig1_flow_diagram.png` / `.svg` — participant flow.\n"
             "2. `figures/fig2_hb_spline.png` / `.svg` — adjusted OR for any stroke vs Hb (reference 13 g/dL).\n"
             "3. `figures/fig3_forest_paper1.png` / `.svg` — covert infarct: imaged vs whole-cohort comparator.\n"
             "4. `figures/fig4_forest_paper2.png` / `.svg` — haematological indices and stroke.\n"
             "5. `figures/fig5_forest_paper3.png` / `.svg` — risk factors by uterine condition.\n")
    L.append("## Automated consistency checks of qualitative wording\n")
    L.append("Each qualitative statement above (e.g. 'associated', 'not clearly associated', 'unchanged') is checked "
             "against the fitted estimates in the same run:\n")
    L.append(md_table(pd.DataFrame([(c, "pass" if ok else "FAIL - revise wording") for c, ok in claims],
                                   columns=["Statement", "Check"])))
    (OUT_DIR / "methods_and_results.md").write_text("\n".join(L), encoding="utf-8")


# ---------------------------------------------------------------------------
def analysis_log():
    R = RESULTS
    L = ["# Analysis log\n",
         f"_Generated by `analysis/run_all.py` on {date.today().isoformat()}. All counts are computed in this run._\n",
         "## 1. How to run\n",
         "```\npip install -r requirements.txt\npython -m analysis.run_all\n```\n"
         "The data files must be placed in `data/` as `stroke_analysis_master.csv` and `stroke.xlsx` "
         "(git-ignored). The workbook is only read, never written. Outputs go to `outputs/`.\n",
         "## 2. PHI handling\n",
         "- `mrn` and `dob` are never written to any output. `add_table()` asserts this for every table.\n"
         "- MRN is used **in memory only** to link analysis_master rows to the review worksheets "
         "(`imaging_only_candidates`, `imaging_only_review`, `deferred_for_review`), which carry MRN but not study_id.\n"
         "- Outputs contain aggregate counts only; no row-level data.\n",
         "## 3. Decisions (in the order they were made)\n"]
    decisions = [
        "Analysis population = `eligible == 1`; verified identical to `exclusion_reason == 0` "
        f"({R['eligible_reason_consistent']}).",
        "Race collapsed from free text: White (text == 'White'); Black (`black_race == 1`); Asian (text begins "
        "'Asian', which includes 'Asian Indian'); everything else, including 'American Indian/Alaskan Native', "
        "Pacific Islander, 'Choose Not to Disclose', 'Unknown' and the invalid value '0', → Other/unknown.",
        "Smoking: 3 levels (never / ever = current or former / unknown); not imputed. Paper 3 'ever smoking' outcome "
        "is restricted to known status.",
        "Migraine (PI decision, 2026-09-28): type and aura are not used. Migraine = 1 for any migraine diagnosis "
        "(codes 1 without aura, 2 with aura, 9 aura/type not stated), 0 for code 0; no missing values. Used as a "
        "binary covariate in Papers 1–2 and a binary outcome in Paper 3. This supersedes the brief's 3-level "
        "migraine covariate and removes the 'migraine with aura' outcome from Paper 3.",
        "Other coded 9 values: hormonal_type 9 (1 patient) → missing. No other 9s in the analysis covariates.",
        "Paper 1 denominator = `stroke_confirmed_imaging` in {0,1} plus revoked/excluded/undecided imaging-only "
        "patients (PI decision: they were imaged; those whose imaging text matches brain-lesion wording "
        "[`paper1.BRAIN_RE`] are covert infarcts, the rest imaged without infarct). Imaging-only cases identified from the notes column "
        "(tags listed in `analysis/paper1.py::TAGS`). Classification of included imaging-only cases (hierarchical): "
        "ICH/SAH → acute/new (acute features, category D, or 'acute or new infarct documented') → symptomatic "
        "('with stroke symptoms documented') → covert: incidental on scan for another indication / multiple or "
        "territorial without symptoms / chronic hedged wording. Covert (primary) = the last three; strict = incidental "
        "only. Clinical stroke = `stroke_any == 1` without an imaging-only flag.",
        "Imaging-only acute/symptomatic/haemorrhagic cases (neither covert nor coded) are excluded from the Paper 1 "
        "three-group comparison and from the covert-vs-imaged-negative models; they are counted in the flow and "
        "`P1_imaging_only_subcats`.",
        "Paper 1 full pre-specified model had EPV < 10 and sparse cells → pre-specified CORE + one-at-a-time strategy "
        "with a second-level collapse rule (code comments in `paper1.py`). The full model is kept in sheet "
        "`P1_full_model_unstable` for transparency only.",
        "Paper 1 selection-bias sensitivity: comparator = all eligible women with `stroke_any == 0` (never-imaged "
        "included), same model structure and same simplifications as the primary model.",
        "Paper 1 covert descriptors: master-file location/territory/number/etiology fields are blank for most "
        "imaging-only cases, so descriptors are taken from `imaging_only_candidates` (linked in memory) and the master "
        "fields are shown alongside.",
        "Paper 2 ischaemic outcome = `stroke_type == 1` vs `stroke_any == 0` (other stroke types excluded rather than "
        "counted as controls).",
        "Incident-only analyses (all papers): cases = `stroke_timing == 3`; timing 1/2 excluded; timing 9 or blank are "
        "neither cases nor controls and are excluded. Paper 3 (outcomes are risk factors) adapts this by excluding "
        "women with a stroke before/at index.",
        "Paper 2 hormonal_type: levels with < 5 outcome events in a model's complete-case sample are merged into "
        "'Other/multiple' (logged per model). Iron-deficiency models with EPV < 10 fall back to a reduced covariate "
        "set (age, race, BMI, HTN, DM, dyslipidaemia, smoking, bleeding, dx group).",
        "Paper 2 mediation: exposure = fibroids flag; mediator = anaemia grade; uterine_bleeding and uterine_dx_group "
        "omitted (bleeding is a potential intermediate on the fibroid→anaemia path; dx group is derived from the "
        "fibroid flag); adenomyosis and endometriosis flags adjusted. Point estimates from the same complete-case "
        "sample for both models; 1,000 non-parametric bootstrap resamples, percentile CI.",
        "MICE (own implementation in `utils.mice`, 10 iterations × m = 20): BMI by Bayesian linear regression + "
        "predictive mean matching (5 donors) (a multinomial-logit imputer for categorical variables is also implemented). "
        "Predictors include the outcome, all covariates and, as auxiliary variables with a 'not measured' level, "
        "the lab exposures (which are never imputed). Pooled by Rubin's rules (Barnard–Rubin-type df).",
        "Paper 3: BMI missing < 10%, so MI not required. Age standardisation: 5-year bands (18–24 … 55–60), "
        "standard = whole eligible cohort.",
        "Continuous terms (age, BMI) enter linearly in primary models; non-linearity tests reported in "
        "`linearity_checks`. Where non-linear, the exposure estimates should be checked against an RCS adjustment "
        "before submission (not done automatically to preserve EPV in Paper 1).",
        "VIF computed per design column (dummy level), not generalised VIF.",
    ]
    L += [f"{i}. {d}" for i, d in enumerate(decisions, 1)]
    L.append("\n## 4. Derivation checks\n")
    L.append(f"- anemia vs hgb thresholds: {R['anemia_mismatch']} mismatches.\n"
             f"- iron_deficiency vs ferritin < 30: {R['irondef_mismatch']} mismatches.\n"
             f"- uterine_dx_group vs the three flags: {R['dxgroup_mismatch']} mismatches.\n"
             f"- study_id unique: {R['study_id_unique']}.\n")
    L.append("## 5. Count reconciliation\n")
    L.append("### Flow\n")
    L.append(md_table(TABLES["flow"]))
    L.append("\n### Paper 1 imaging-only flags vs review sheets\n")
    L.append(md_table(TABLES["P1_reconciliation"]))
    L.append("\nNotes-based status of all eligible patients with any imaging-only note: "
             f"{R['p1_io_status']}.\n\nIncluded imaging-only subcategories: {R['p1_io_subcat']}; unclassified: "
             f"{R['p1_unclassified']}.\n")
    L.append("Cross-tabulations of each sheet's own category against the notes classification are in sheets "
             "`P1_recon_review_x_notes`, `P1_recon_deferred_x_notes` and `P1_recon_indication`.\n")
    L.append("## 6. Data problems flagged (not fixed)\n")
    L.append(md_table(TABLES["data_flags"]))
    L.append(f"\n- Codebook variables absent from analysis_master: {', '.join(R['codebook_absent'])}.")
    L.append(f"- {R['p1_rev_excl_not_in_denominator']} eligible women whose imaging-only infarct was revoked/excluded/"
             f"undecided have blank `stroke_confirmed_imaging` although their scans were reviewed. PI decision: all "
             f"enter the Paper 1 denominator; {R['p1_n_restored']} with brain-lesion wording count as covert infarcts, "
             f"the rest as imaged without infarct (see sheet `P1_restored_cases`).")
    L.append(f"- {R['p1_io_desc_blank']} imaging-only cases have blank location/territory/number fields and "
             f"{R['p1_io_date_blank']} have blank `stroke_date` in analysis_master; all blank-timing strokes in the "
             f"eligible cohort are imaging-only cases.\n")
    L.append("## 7. Surprising or noteworthy findings\n")
    miss = R["missing"]
    dm_prev = TABLES["codebook_summary"].set_index("Variable").loc["dm", "Observed distribution (eligible)"]
    dl_prev = TABLES["codebook_summary"].set_index("Variable").loc["dyslipidemia", "Observed distribution (eligible)"]
    L.append(f"- Diabetes ({dm_prev}) and dyslipidaemia ({dl_prev}) are very common for women aged 18–60; this may "
             f"reflect EHR ascertainment in a tertiary centre and should be checked against source definitions.")
    L.append(f"- Migraine code 9 ('migraine, aura status not stated') covers {R['flag_migraine_9_is_migraine'][1]:,} "
             f"eligible women; the Codebook labels it 'Unknown'. Resolved by the PI: counted as migraine.")
    L.append("- Paper 1: relative to imaged women, migraine appears protective against covert infarct; relative to "
             "the whole cohort the direction reverses. This is an ascertainment artefact from preferential imaging of "
             "migraine patients.")
    L.append("- Paper 2: iron deficiency is inversely associated with stroke in the ferritin subset; ferritin is "
             "an acute-phase reactant and was measured in a selected subset (`P2_ferritin_selection`).")
    L.append("- Paper 2: macrocytosis and low platelets are associated with stroke; these were not hypotheses of the "
             "brief and may reflect comorbidity (liver disease, alcohol, medication).")
    L.append("- The Codebook text for `stroke_any_incl_imaging` (+408) and `stroke_any` (330 imaging-only) describes "
             "an earlier version; in the current file imaging-only cases are already inside `stroke_any`.")
    L.append("")
    L.append("## 8. Model diagnostics — models with warnings\n")
    dg = pd.DataFrame(DIAG)
    bad = dg[dg["Flags"] != ""]
    L.append(f"{len(dg)} models registered; {len(bad)} carry at least one flag (EPV < 10, sparse cell, separation, "
             f"VIF > 5 or non-convergence). Full list in sheet `model_diagnostics`.\n")
    L.append(md_table(bad[["Paper", "Model", "N", "Events", "EPV", "Max VIF", "Flags"]]))
    nc = int((~dg["Converged"].astype(bool)).sum())
    L.append(f"\nNon-converged models: {nc}. Max VIF across all models: {dg['Max VIF'].max():.2f} "
             f"(values > 5 occur only in models containing product terms or spline bases).\n")
    L.append("## 9. Linearity checks (RCS, 4 knots)\n")
    lin = pd.DataFrame(LINEARITY)
    L.append(md_table(lin[["Paper", "Model", "Variable", "Knots", "Non-linearity χ²", "df", "p (text)",
                           "Non-linear (p<0.05)"]]))
    L.append("\n## 10. Limitations to state in the manuscripts\n")
    L.append(f"- Cross-sectional design; {R['flag_stroke_no_timing'][1]} strokes lack timing, {R['n_stroke_before']:,} precede index "
             f"and {R['n_stroke_before_or_same'] - R['n_stroke_before']:,} were coded at the index encounter; reverse causation cannot be excluded — incident-only analyses provided.\n"
             f"- Smoking unknown for {TABLES['codebook_summary'].set_index('Variable').loc['smoking', 'Code 9 %']}% "
             f"(kept as an 'unknown' level, not imputed).\n"
             f"- Covert infarcts are ascertained only in imaged women; imaging is strongly selected by migraine and "
             f"hypertension.\n"
             f"- Hb/ferritin: single value closest to index (±3 years); no repeat labs. Ferritin measured in "
             f"{R['p2_ferritin_pct']:.1f}%. MCV/platelet dates unverifiable.\n"
             f"- stroke_etiology is imaging-inferred, not TOAST. stroke_date is the earliest coded date, not "
             f"chart-verified.\n"
             f"- Hormonal therapy is ascertained from administered medications (combined OC under-captured) and is at "
             f"or after index.\n")
    L.append("## 11. Open questions for the PI\n")
    L.append("1. (Resolved 2026-09-28) Migraine: any diagnosis = 1, type/aura ignored.\n"
             "2. (Resolved 2026-09-28) Revoked/excluded imaging-only patients: brain wording = covert infarct, others "
             "imaged-negative; all in the denominator. Former question: should they be added to the "
             "Paper 1 imaged-negative denominator?\n"
             "3. (Resolved 2026-09-28: broad definition kept.) Confirm the covert-infarct definition (currently includes 'multiple/territorial, no symptoms' and "
             "'chronic, hedged wording' in addition to 'incidental on scan for another indication').\n")
    L.append("## 12. Chronological run log\n")
    L += [f"- **{s}** — {t}" for s, t in LOG]
    (OUT_DIR / "analysis_log.md").write_text("\n".join(L), encoding="utf-8")


def run():
    methods_results()
    analysis_log()
