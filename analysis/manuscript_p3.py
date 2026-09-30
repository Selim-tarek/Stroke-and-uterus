"""Draft manuscript for Paper 3: migraine and vascular risk-factor profiles across fibroids, adenomyosis and
endometriosis. Every number is read from RESULTS / TABLES produced in the same run; literature statements carry
[ref] placeholders. Outputs: outputs/migraine_vascular_manuscript.{json,md,docx} and
outputs/migraine_vascular_abstract_tables_figures.{json,md,docx}. US spelling.
"""
import json
import re
import subprocess

import pandas as pd

from .manuscript import _us, to_markdown
from .utils import OUT_DIR, RESULTS, TABLES, fmt_p, log

GROUPS = ["Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition"]
CARDIO = ["Hypertension", "Diabetes", "Dyslipidaemia", "Obesity (BMI ≥30)", "Atrial fibrillation",
          "Coronary artery disease"]
FIGS = [("fig9_p3_migraine_prevalence.png", 5.6), ("fig5_forest_paper3.png", 6.3), ("fig16_p3_age_strata.png", 6.5),
        ("fig10_p3_migraine_by_age.png", 5.8), ("fig17_p4_incidence.png", 6.5)]


def build():
    R, T = RESULTS, TABLES
    checks = []
    main, prev, st, em = R["p3_main"], R["p3_prev"], R["p3_strata"], R["p3_em"].set_index("Outcome")
    ms, inc = R["p3_msens"], R["p3_inc"]
    gn = R["p3_group_n"]
    age = R["p3_age_by_group"]
    iy0, iy1 = R["index_range"]
    n = R["p3_n"]

    def m(o, g):
        return main[(main.Outcome == o) & (main["Group vs fibroids only"] == g)].iloc[0]

    def orx(o, g):
        return m(o, g)["M2 OR (95% CI) age+race+BMI"]

    def pv(o, g):
        return prev[(prev.Outcome == o) & (prev.Group == g)].iloc[0]

    def sx(o, a, g):
        return st[(st.Outcome == o) & (st.Age == a) & (st["Group vs fibroids only"] == g)].iloc[0]

    def msx(model, g):
        return ms[(ms.Model == model) & (ms["Group vs fibroids only"] == g)].iloc[0]

    def P_(t):
        t = str(t)
        return f"P{t}" if t.startswith("<") else f"P={t}"

    pc, pm = R["p4_cohort"], R["p4_models"]
    pmig = R["p4_mig_event"].set_index("Event type")

    def q4(oc, model, g):
        return pm[(pm.Outcome == oc) & (pm.Model == model) & (pm["Group vs fibroids only"] == g)].iloc[0]

    def sc(txt, first=False):
        mm = re.match(r"\s*([\d.]+) \(([\d.]+–[\d.]+)\)", str(txt))
        return f"{mm.group(1)}; {'95% CI ' if first else ''}{mm.group(2)}" if mm else str(txt)

    def pstd(o, g):
        return f"{float(pv(o, g)['Age-standardised % (95% CI)'].split(' ')[0]):.1f}%"

    n_mig = int(T["P3_Table1"].set_index("Characteristic").loc["Migraine (any)", "Overall"].split(" ")[0]
                .replace(",", ""))
    mig = "Migraine (any)"
    B = []
    H1 = lambda t: B.append({"t": "h1", "text": t})  # noqa: E731
    H2 = lambda t: B.append({"t": "h2", "text": t})  # noqa: E731
    P = lambda t: B.append({"t": "p", "text": t})  # noqa: E731

    def TABLE(title, df, note="", widths=None, font=8):
        B.append({"t": "table", "title": title, "columns": [str(c) for c in df.columns],
                  "rows": [["" if (isinstance(v, float) and pd.isna(v)) else
                            (f"{v:,}" if isinstance(v, int) else str(v)) for v in r] for r in df.itertuples(index=False)],
                  "note": note, "widths": widths, "font": font})

    def FIG(path, caption, width=6.3):
        B.append({"t": "figure", "path": str(OUT_DIR / "figures" / path), "caption": caption, "width": width})

    BR = lambda: B.append({"t": "pagebreak"})  # noqa: E731

    # ------------------------------------------------------------------ title page
    B.append({"t": "title", "text": "Migraine, cardiometabolic risk profiles and stroke in women with uterine "
                                    "fibroids, adenomyosis and endometriosis"})
    P("Running title: Migraine and vascular risk in benign uterine disease")
    P("Authors: [Author names, degrees, affiliations]")
    P("Corresponding author: [name, address, email]")
    P("Word count (main text): [to be completed]; Tables: 3; Figures: 5; Supplementary material: eTables 1–7")
    P("Keywords: migraine; endometriosis; adenomyosis; uterine fibroids; hypertension; cardiovascular risk; women")
    B.append({"t": "note", "text": "DRAFT generated from the analysis pipeline (python -m analysis.run_all). All numbers "
                                   "are produced by code in the same run. [ref] marks statements that need a citation."})
    BR()

    # ------------------------------------------------------------------ abstract
    endo, aden, multi = "Endometriosis only", "Adenomyosis only", ">1 condition"
    H1("Abstract")
    P("**Background.** Fibroids, adenomyosis and endometriosis are common, hormonally driven conditions, and "
      "endometriosis has been linked to migraine and to cardiovascular disease. Whether vascular risk profiles differ "
      "between these conditions within one population is not well characterized.")
    P(f"**Methods.** We studied {n:,} women aged 18–60 years diagnosed with benign uterine disease between {iy0} and "
      f"{iy1} at Mayo Clinic, grouped as fibroids only, adenomyosis only, endometriosis only or more than one "
      f"condition. Ten vascular risk factors, including migraine, were compared with fibroids only by logistic "
      f"regression adjusted for age, race and body mass index, with Benjamini–Hochberg correction for 30 comparisons. "
      f"Prevalences were directly age-standardized. Stroke or transient ischemic attack (TIA) after diagnosis was "
      f"compared by Poisson regression.")
    P(f"**Results.** The cohort included {gn['Fibroids only']:,} women with fibroids only, {gn[aden]:,} with "
      f"adenomyosis only, {gn[endo]:,} with endometriosis only and {gn[multi]:,} with more than one condition. "
      f"Age-standardized migraine prevalence was {pstd(mig, 'Fibroids only')} with fibroids only, "
      f"{pstd(mig, aden)} with adenomyosis only and {pstd(mig, endo)} with endometriosis only. Compared with fibroids "
      f"only, migraine was more common with adenomyosis (OR {sc(orx(mig, aden), True)}), endometriosis (OR "
      f"{sc(orx(mig, endo))}) and more than one condition (OR {sc(orx(mig, multi))}), in both age groups and after "
      f"adjustment for hormonal therapy and uterine bleeding. Endometriosis only was associated with lower odds of "
      f"hypertension (OR {sc(orx('Hypertension', endo))}), diabetes (OR {sc(orx('Diabetes', endo))}), dyslipidemia "
      f"(OR {sc(orx('Dyslipidaemia', endo))}) and coronary artery disease (OR {sc(orx('Coronary artery disease', endo))}). "
      f"These differences were concentrated in women aged 18–39 years (hypertension OR "
      f"{sc(sx('Hypertension', '18–39', endo)['OR (95% CI)'])}) and were close to null at 40–60 years (OR "
      f"{sc(sx('Hypertension', '40–60', endo)['OR (95% CI)'])}). Adenomyosis only was associated with more "
      f"dyslipidemia (OR {sc(orx('Dyslipidaemia', aden))}) and thrombophilia (OR {sc(orx('Thrombophilia', aden))}). "
      f"During {pc['py']:,.0f} person-years of follow-up ({pc['ev']} strokes or TIAs), adjusted stroke rates did not "
      f"differ between endometriosis and fibroids (rate ratio {sc(q4('Stroke or TIA', 'M2 + cardiometabolic', endo)['RR (95% CI)'])}). "
      f"The higher rate with adenomyosis (rate ratio {sc(q4('Stroke or TIA', 'M2 + cardiometabolic', aden)['RR (95% CI)'])}; "
      f"{int(q4('Stroke or TIA', 'M2 + cardiometabolic', aden)['Events in group'])} events) was not significant for "
      f"ischemic stroke, and women with TIA more often had migraine than women with ischemic stroke "
      f"({pmig.loc['TIA', 'Migraine %']:.0f}% vs {pmig.loc['Ischemic stroke', 'Migraine %']:.0f}%).")
    P("**Conclusions.** Women with adenomyosis or endometriosis had a higher migraine burden than women with fibroids, "
      "whereas cardiometabolic risk factors were less common in young women with endometriosis. Vascular risk "
      "assessment in benign uterine disease may need to account for migraine as well as conventional risk factors, "
      "particularly in adenomyosis and endometriosis. Stroke rates after diagnosis were similar across conditions; "
      "the excess of TIA diagnoses in adenomyosis may partly reflect migraine.")
    checks.append(("Abstract: endometriosis stroke RR CI includes 1; adenomyosis RR > 1; adenomyosis ischemic CI includes 1; "
                   "TIA migraine % > ischemic, P<0.05",
                   q4("Stroke or TIA", "M2 + cardiometabolic", endo)["CI low"] < 1 < q4("Stroke or TIA", "M2 + cardiometabolic", endo)["CI high"]
                   and q4("Stroke or TIA", "M2 + cardiometabolic", aden)["CI low"] > 1
                   and q4("Ischemic stroke", "M2 + cardiometabolic", aden)["CI low"] < 1
                   and pmig.loc["TIA", "Migraine %"] > pmig.loc["Ischemic stroke", "Migraine %"]
                   and R["p4_mig_event_p"] < 0.05))
    checks.append(("Abstract: migraine ORs exclude 1 in all three groups",
                   all(m(mig, g)["CI low"] > 1 for g in [aden, endo, multi])))
    checks.append(("Abstract: endometriosis lower HTN/DM/dyslipidemia/CAD",
                   all(m(o, endo)["CI high"] < 1 for o in ["Hypertension", "Diabetes", "Dyslipidaemia",
                                                           "Coronary artery disease"])))
    checks.append(("Abstract: endometriosis HTN 18–39 excludes 1; 40–60 attenuated (OR closer to 1)",
                   sx("Hypertension", "18–39", endo)["CI high"] < 1 and
                   abs(sx("Hypertension", "40–60", endo)["OR"] - 1) < abs(sx("Hypertension", "18–39", endo)["OR"] - 1)))
    checks.append(("Abstract: migraine higher in both age groups and after HT/bleeding adjustment",
                   all(sx(mig, a, g)["CI low"] > 1 for a in ["18–39", "40–60"] for g in [aden, endo, multi]) and
                   all(msx("M2 + hormonal therapy + uterine bleeding", g)["CI low"] > 1 for g in [aden, endo, multi])))
    checks.append(("Abstract: adenomyosis dyslipidemia and thrombophilia > 1",
                   m("Dyslipidaemia", aden)["CI low"] > 1 and m("Thrombophilia", aden)["CI low"] > 1))
    # AAN-format abstract (five headings, ≤300 words in the body, no institution name) for the submission file
    _q_endo = q4("Stroke or TIA", "M2 + cardiometabolic", endo)
    _q_aden = q4("Stroke or TIA", "M2 + cardiometabolic", aden)
    aan = [
        ("OBJECTIVE", "To compare migraine, vascular risk factors and stroke across uterine fibroids, adenomyosis and "
                      "endometriosis."),
        ("BACKGROUND", "Endometriosis has been linked to migraine and cardiovascular disease, but vascular risk "
                       "profiles have rarely been compared across benign uterine conditions within one population."),
        ("DESIGN/METHODS", f"We studied {n:,} women aged 18–60 years diagnosed with benign uterine disease between "
                           f"{iy0} and {iy1} at one academic center, grouped as fibroids only (reference), adenomyosis "
                           f"only, endometriosis only or more than one condition. Ten vascular risk factors, including "
                           f"migraine, were compared by logistic regression adjusted for age, race and BMI, with "
                           f"Benjamini–Hochberg correction; prevalences were age-standardized. Stroke (ischemic, "
                           f"hemorrhagic or cerebral venous thrombosis) or TIA after diagnosis was compared by Poisson "
                           f"regression."),
        ("RESULTS", f"Age-standardized migraine prevalence was {pstd(mig, 'Fibroids only')} with fibroids only "
                    f"(n={gn['Fibroids only']:,}), {pstd(mig, aden)} with adenomyosis only (n={gn[aden]:,}) and "
                    f"{pstd(mig, endo)} with endometriosis only (n={gn[endo]:,}); adjusted ORs versus fibroids were "
                    f"{orx(mig, aden).replace(' (', ' (95% CI ', 1)} and {orx(mig, endo)}, consistent across age groups and after "
                    f"adjustment for hormonal therapy and uterine bleeding. Endometriosis only was associated with lower "
                    f"odds of hypertension (OR {sc(orx('Hypertension', endo))}), diabetes "
                    f"({sc(orx('Diabetes', endo))}) and coronary artery disease "
                    f"({sc(orx('Coronary artery disease', endo))}), mainly before age 40 (hypertension OR "
                    f"{sc(sx('Hypertension', '18–39', endo)['OR (95% CI)'])} at 18–39 vs "
                    f"{sc(sx('Hypertension', '40–60', endo)['OR (95% CI)'])} at 40–60 years). Adenomyosis only was "
                    f"associated with more dyslipidemia ({sc(orx('Dyslipidaemia', aden))}) and thrombophilia "
                    f"({sc(orx('Thrombophilia', aden))}). During {pc['py']:,.0f} person-years ({pc['ev']} strokes or TIAs), "
                    f"stroke rates did not differ between endometriosis and fibroids (rate ratio "
                    f"{sc(_q_endo['RR (95% CI)'])}). The higher rate with adenomyosis ({sc(_q_aden['RR (95% CI)'])}; "
                    f"{int(_q_aden['Events in group'])} events) was not significant for ischemic stroke, and women with "
                    f"TIA more often had migraine than those with ischemic stroke "
                    f"({pmig.loc['TIA', 'Migraine %']:.0f}% vs {pmig.loc['Ischemic stroke', 'Migraine %']:.0f}%)."),
        ("CONCLUSIONS", "Adenomyosis and endometriosis carried a higher migraine burden than fibroids, whereas young "
                        "women with endometriosis had fewer cardiometabolic risk factors. Stroke rates after diagnosis "
                        "were similar across conditions; excess TIA diagnoses in adenomyosis may partly reflect "
                        "migraine."),
    ]
    _w = sum(len(t.split()) for _, t in aan)
    RESULTS["aan_p3"], RESULTS["aan_p3_words"] = aan, _w
    checks.append((f"AAN abstract body ≤300 words (now {_w})", _w <= 300))
    BR()

    # ------------------------------------------------------------------ introduction
    H1("Introduction")
    P("Uterine fibroids, adenomyosis and endometriosis affect a large proportion of women of reproductive age and "
      "share hormonal drivers, heavy or painful bleeding and frequent hormonal treatment [ref]. Endometriosis has been "
      "associated with migraine [ref] and, in prospective cohorts, with coronary heart disease and stroke [ref]. "
      "Fibroids have been linked to hypertension and metabolic risk factors [ref]. Adenomyosis has been studied far "
      "less.")
    P("Most studies compare women with one condition with women from the general population, which mixes the effect "
      "of the condition with differences in age, health-care contact and diagnostic pathway. Comparing the "
      "conditions with one another, within one health system, gives a more direct view of how their vascular risk "
      "profiles differ. We compared the prevalence of migraine and nine other vascular risk factors across fibroids "
      "only, adenomyosis only, endometriosis only and more than one condition, and examined whether differences "
      "varied with age.")

    # ------------------------------------------------------------------ methods
    H1("Methods")
    H2("Design and population")
    P(f"This cross-sectional study used electronic health records of women aged 18–60 years with a first diagnosis of "
      f"uterine fibroids, adenomyosis or endometriosis between {iy0} and {iy1} at Mayo Clinic (index date = first "
      f"diagnosis). Of {R['n_total']:,} records, {n:,} met the eligibility criteria. The study is reported following "
      f"STROBE guidance for cross-sectional studies.")
    H2("Exposure")
    P("Women were grouped by recorded diagnoses as fibroids only (reference), adenomyosis only, endometriosis only, or "
      "more than one condition.")
    H2("Outcomes")
    P("Ten vascular risk factors were examined: hypertension, diabetes, dyslipidemia, obesity (body mass index "
      "[BMI] ≥30 kg/m²), ever smoking (current or former vs never; unknown status excluded), migraine (any migraine "
      "diagnosis, regardless of type or aura), atrial fibrillation, coronary artery disease, prior venous "
      "thromboembolism (VTE) and thrombophilia. All were recorded as present or absent from diagnosis codes and "
      "clinical records at any time.")
    H2("Statistical analysis")
    P("Prevalences were directly age-standardized to the whole cohort in 5-year age bands. Each risk factor was "
      "modeled by logistic regression with robust (HC1) standard errors. Model 1 adjusted for age (linear) and race; "
      "the primary model (Model 2) additionally adjusted for BMI (except for obesity). Pre-specified rules simplified "
      "models with sparse covariate levels or fewer than 10 events per variable; none were needed for the primary "
      "models. P values were corrected with the Benjamini–Hochberg method across the 30 condition contrasts (10 risk "
      "factors × 3 groups). Effect modification by age (18–39 vs 40–60 years) was tested with condition × age-group "
      "interaction terms, and models were repeated within each age group. Linearity of age and BMI was checked with "
      "restricted cubic splines.")
    P(f"Sensitivity analyses excluded women whose stroke occurred before or at the index date "
      f"(n = {R['p3_excl_prior_n']:,}), because a stroke may prompt risk-factor testing. For migraine, models were "
      f"further adjusted for hormonal therapy and heavy or abnormal uterine bleeding, and repeated in women without "
      f"hormonal therapy and in women without heavy bleeding, because both are linked to the uterine conditions and "
      f"to headache.")
    H2("Stroke after diagnosis")
    P(f"Women without a stroke or TIA before or at the index date were followed from the index date to the first of "
      f"stroke or TIA, death or last recorded encounter. Women with an undatable stroke (n = {pc['n_unknown']:,}) or "
      f"no follow-up time (n = {pc['n_zero']:,}) were excluded. Crude rates were calculated with exact Poisson "
      f"confidence intervals. Rate ratios versus fibroids only were estimated by Poisson regression with a log "
      f"person-time offset and robust standard errors, adjusted for age and race, then additionally for BMI, "
      f"hypertension, diabetes, dyslipidemia, smoking, atrial fibrillation and coronary artery disease (primary), then "
      f"for migraine, and finally for hormonal therapy and uterine bleeding. Outcomes were stroke or TIA (primary), "
      f"ischemic stroke, and stroke excluding TIA. Because migraine can mimic TIA, migraine prevalence was compared "
      f"between women with TIA and with ischemic stroke. Analyses used Python (pandas, statsmodels).")

    # ------------------------------------------------------------------ results
    H1("Results")
    H2("Participants")
    t1 = T["P3_Table1"].set_index("Characteristic")

    def t1v(ch, col):
        return t1.loc[ch, col] if not isinstance(t1.loc[ch, col], pd.Series) else t1.loc[ch, col].iloc[0]
    P(f"Of {n:,} women, {gn['Fibroids only']:,} had fibroids only, {gn[aden]:,} adenomyosis only, {gn[endo]:,} "
      f"endometriosis only and {gn[multi]:,} more than one condition (Table 1). Women with endometriosis only were "
      f"younger (median age {age[endo]:.0f} vs {age['Fibroids only']:.0f} years with fibroids only), had lower BMI "
      f"and less often reported heavy or abnormal bleeding ({t1v('Heavy/abnormal uterine bleeding', endo)} vs "
      f"{t1v('Heavy/abnormal uterine bleeding', 'Fibroids only')}). Migraine was recorded in {n_mig:,} women overall.")
    H2("Migraine")
    P(f"Age-standardized migraine prevalence was {pv(mig, 'Fibroids only')['Age-standardised % (95% CI)']}% with "
      f"fibroids only, {pv(mig, aden)['Age-standardised % (95% CI)']}% with adenomyosis only, "
      f"{pv(mig, endo)['Age-standardised % (95% CI)']}% with endometriosis only and "
      f"{pv(mig, multi)['Age-standardised % (95% CI)']}% with more than one condition (Figure 1). Compared with "
      f"fibroids only, the adjusted OR for migraine was {orx(mig, aden)} for adenomyosis only, {orx(mig, endo)} for "
      f"endometriosis only and {orx(mig, multi)} for more than one condition (all Benjamini–Hochberg-adjusted "
      f"{P_(m(mig, endo)['M2 p (BH-adjusted, all contrasts) (text)'])}; Table 2).")
    hb = "M2 + hormonal therapy + uterine bleeding"
    P(f"The associations were present in women aged 18–39 years (adenomyosis only "
      f"{sx(mig, '18–39', aden)['OR (95% CI)']}; endometriosis only {sx(mig, '18–39', endo)['OR (95% CI)']}) and "
      f"40–60 years ({sx(mig, '40–60', aden)['OR (95% CI)']}; {sx(mig, '40–60', endo)['OR (95% CI)']}), although "
      f"they were somewhat weaker at older ages (interaction {P_(em.loc[mig, 'p BH (text)'])}; Figure 4). Further "
      f"adjustment for hormonal therapy and uterine bleeding gave ORs of {msx(hb, aden)['OR (95% CI)']}, "
      f"{msx(hb, endo)['OR (95% CI)']} and {msx(hb, multi)['OR (95% CI)']}. Among women without hormonal therapy "
      f"the ORs were {msx('M2, women without hormonal therapy', aden)['OR (95% CI)']}, "
      f"{msx('M2, women without hormonal therapy', endo)['OR (95% CI)']} and "
      f"{msx('M2, women without hormonal therapy', multi)['OR (95% CI)']}, and among women without heavy bleeding "
      f"{msx('M2, women without heavy/abnormal bleeding', aden)['OR (95% CI)']}, "
      f"{msx('M2, women without heavy/abnormal bleeding', endo)['OR (95% CI)']} and "
      f"{msx('M2, women without heavy/abnormal bleeding', multi)['OR (95% CI)']} (eTable 4).")
    checks.append(("Results: all migraine BH P <0.001",
                   all(m(mig, g)["M2 p (BH-adjusted, all contrasts) (text)"] == "<0.001" for g in [aden, endo, multi])))
    checks.append(("Results: migraine ORs > 1 in every sensitivity model",
                   bool((ms.dropna(subset=["CI low"])["CI low"] > 1).all())))
    checks.append(("Results: migraine weaker at older ages (40–60 OR < 18–39 OR for adenomyosis and endometriosis)",
                   all(sx(mig, "40–60", g)["OR"] < sx(mig, "18–39", g)["OR"] for g in [aden, endo])))
    H2("Cardiometabolic risk factors")
    P(f"Compared with fibroids only, endometriosis only was associated with lower odds of hypertension "
      f"({orx('Hypertension', endo)}), diabetes ({orx('Diabetes', endo)}), dyslipidemia ({orx('Dyslipidaemia', endo)}), "
      f"obesity ({orx('Obesity (BMI ≥30)', endo)}), atrial fibrillation ({orx('Atrial fibrillation', endo)}) and "
      f"coronary artery disease ({orx('Coronary artery disease', endo)}) (Figure 2). Adenomyosis only was associated "
      f"with lower odds of hypertension ({orx('Hypertension', aden)}) and atrial fibrillation "
      f"({orx('Atrial fibrillation', aden)}), but higher odds of dyslipidemia ({orx('Dyslipidaemia', aden)}). Diabetes "
      f"was slightly more common ({orx('Diabetes', aden)}; Benjamini–Hochberg-adjusted "
      f"P={m('Diabetes', aden)['M2 p (BH-adjusted, all contrasts) (text)']}, not significant after correction) and "
      f"obesity was similar ({orx('Obesity (BMI ≥30)', aden)}). More than one condition "
      f"showed smaller differences in the same direction as endometriosis.")
    checks.append(("Results: >1 condition cardiometabolic ORs between endometriosis OR and 1",
                   all(m(o, endo)["OR"] < m(o, multi)["OR"] < 1 for o in CARDIO)))
    checks.append(("Results: endometriosis all six cardiometabolic ORs < 1",
                   all(m(o, endo)["CI high"] < 1 for o in CARDIO)))
    checks.append(("Results: adenomyosis HTN & AF < 1; dyslipidemia > 1; diabetes & obesity CI include 1",
                   m("Hypertension", aden)["CI high"] < 1 and m("Atrial fibrillation", aden)["CI high"] < 1 and
                   m("Dyslipidaemia", aden)["CI low"] > 1 and
                   m("Diabetes", aden)["OR"] > 1 and m("Diabetes", aden)["M2 p (BH-adjusted, all contrasts)"] >= 0.05 and
                   m("Obesity (BMI ≥30)", aden)["CI low"] < 1 < m("Obesity (BMI ≥30)", aden)["CI high"]))
    H2("Differences by age")
    P(f"Differences in cardiometabolic risk factors varied with age (Benjamini–Hochberg-adjusted interaction P: "
      f"hypertension {em.loc['Hypertension', 'p BH (text)']}, diabetes {em.loc['Diabetes', 'p BH (text)']}, "
      f"dyslipidemia {em.loc['Dyslipidaemia', 'p BH (text)']}, coronary artery disease "
      f"{em.loc['Coronary artery disease', 'p BH (text)']}; Table 3, Figure 3). Among women aged 18–39 years, "
      f"endometriosis only was associated with markedly lower odds of hypertension "
      f"({sx('Hypertension', '18–39', endo)['OR (95% CI)']}), diabetes ({sx('Diabetes', '18–39', endo)['OR (95% CI)']}) "
      f"and coronary artery disease ({sx('Coronary artery disease', '18–39', endo)['OR (95% CI)']}). At 40–60 years "
      f"the corresponding ORs were {sx('Hypertension', '40–60', endo)['OR (95% CI)']}, "
      f"{sx('Diabetes', '40–60', endo)['OR (95% CI)']} and {sx('Coronary artery disease', '40–60', endo)['OR (95% CI)']}.")
    checks.append(("Results: endometriosis CAD 18–39 < 1",
                   sx("Coronary artery disease", "18–39", endo)["CI high"] < 1))
    H2("Thrombotic risk factors and smoking")
    P(f"Prior VTE was less common with endometriosis only ({orx('Prior VTE', endo)}) and more common with more than "
      f"one condition ({orx('Prior VTE', multi)}). Thrombophilia was more common with adenomyosis only "
      f"({orx('Thrombophilia', aden)}) and more than one condition ({orx('Thrombophilia', multi)}). Smoking status was "
      f"known for {R['p3_smoking_known']:,} women; among them, ever smoking was slightly more common with "
      f"endometriosis only ({orx('Ever smoking', endo)}).")
    changed = [r for _, r in inc.iterrows()]
    P(f"Excluding the {R['p3_excl_prior_n']:,} women with a stroke before or at the index date gave similar "
      f"estimates (eTable 5).")
    H2("Stroke after diagnosis")
    rt = R["p4_rates"]

    def rate(oc, g):
        return rt[(rt.Outcome == oc) & (rt.Group == g)].iloc[0]
    M2n, M3n = "M2 + cardiometabolic", "M3 + migraine"
    P(f"{pc['n']:,} women were followed for {pc['py']:,.0f} person-years (median {pc['fu_median']:.1f} years), with "
      f"{pc['ev']} strokes or TIAs ({pc['ev_isch']} ischemic strokes). Crude rates per 1,000 person-years were "
      f"{rate('Stroke or TIA', 'Fibroids only')['Rate /1,000 PY (95% CI)']} with fibroids only, "
      f"{rate('Stroke or TIA', aden)['Rate /1,000 PY (95% CI)']} with adenomyosis only, "
      f"{rate('Stroke or TIA', endo)['Rate /1,000 PY (95% CI)']} with endometriosis only and "
      f"{rate('Stroke or TIA', multi)['Rate /1,000 PY (95% CI)']} with more than one condition (Figure 5). After "
      f"adjustment for age, race and cardiometabolic factors, the rate ratio versus fibroids only was "
      f"{q4('Stroke or TIA', M2n, endo)['RR (95% CI)']} for endometriosis only, "
      f"{q4('Stroke or TIA', M2n, aden)['RR (95% CI)']} for adenomyosis only "
      f"({int(q4('Stroke or TIA', M2n, aden)['Events in group'])} events) and "
      f"{q4('Stroke or TIA', M2n, multi)['RR (95% CI)']} for more than one condition. Further adjustment for migraine "
      f"gave {q4('Stroke or TIA', M3n, endo)['RR (95% CI)']}, {q4('Stroke or TIA', M3n, aden)['RR (95% CI)']} and "
      f"{q4('Stroke or TIA', M3n, multi)['RR (95% CI)']}. For ischemic stroke the adjusted rate ratios were "
      f"{q4('Ischemic stroke', M2n, endo)['RR (95% CI)']}, {q4('Ischemic stroke', M2n, aden)['RR (95% CI)']} and "
      f"{q4('Ischemic stroke', M2n, multi)['RR (95% CI)']}, and for stroke excluding TIA "
      f"{q4('Stroke excluding TIA', M2n, endo)['RR (95% CI)']}, {q4('Stroke excluding TIA', M2n, aden)['RR (95% CI)']} "
      f"and {q4('Stroke excluding TIA', M2n, multi)['RR (95% CI)']} (eTable 6).")
    ty = R["p4_type"].set_index("Group")
    P(f"TIA accounted for {ty.loc[aden, 'TIA']} of all strokes or TIAs in women with adenomyosis only and "
      f"{ty.loc[endo, 'TIA']} with endometriosis only, compared with {ty.loc['Fibroids only', 'TIA']} with fibroids "
      f"only. Migraine was recorded in {pmig.loc['TIA', 'Migraine %']:.1f}% of women with TIA and "
      f"{pmig.loc['Ischemic stroke', 'Migraine %']:.1f}% of women with ischemic stroke "
      f"({P_(fmt_p(R['p4_mig_event_p']))}; eTable 7).")
    checks.append(("Results: >1 condition stroke/TIA RR > 1 at M2 and CI includes 1 after migraine",
                   q4("Stroke or TIA", M2n, multi)["CI low"] > 1 and q4("Stroke or TIA", M3n, multi)["CI low"] < 1))
    checks.append(("Results: adenomyosis stroke/TIA CI includes 1 after migraine",
                   q4("Stroke or TIA", M3n, aden)["CI low"] < 1))

    # ------------------------------------------------------------------ discussion
    H1("Discussion")
    P(f"In this cohort of {n:,} women with benign uterine disease, adenomyosis and endometriosis were associated with a "
      f"higher prevalence of migraine than fibroids, a difference that was present in younger and older women and "
      f"persisted after accounting for hormonal therapy and uterine bleeding. By contrast, endometriosis was "
      f"associated with a markedly lower prevalence of cardiometabolic risk factors, but mainly in women under 40 "
      f"years; at 40–60 years the profiles were broadly similar. Adenomyosis combined a higher migraine burden with "
      f"more dyslipidemia and thrombophilia.")
    P("The migraine finding agrees with reports linking endometriosis and migraine [ref] and extends it to "
      "adenomyosis, which has rarely been studied. Shared mechanisms have been proposed, including estrogen "
      "fluctuation, prostaglandin-mediated inflammation and central sensitization to pain [ref]. Because women with "
      "endometriosis and adenomyosis are often assessed for pelvic pain, headache may also be more often recorded; "
      "our design cannot separate these explanations.")
    P("The lower cardiometabolic burden in young women with endometriosis probably reflects, at least in part, "
      "differences in who is diagnosed: endometriosis is often diagnosed at a younger age and at lower BMI, and "
      "fibroids are associated with obesity and hypertension [ref]. The attenuation after age 40 suggests that "
      "these differences may not persist. Prospective studies linking endometriosis with later cardiovascular "
      "disease [ref] should therefore be interpreted alongside the migraine burden, which is itself associated with "
      "stroke in women [ref].")
    P("Despite these different risk-factor profiles, stroke rates after diagnosis did not differ between "
      "endometriosis and fibroids once age and cardiometabolic factors were accounted for. The lower crude rate with "
      "endometriosis reflected younger age. The higher rate of stroke or TIA with adenomyosis was based on few events, "
      "was not significant for ischemic stroke, and was attenuated after adjustment for migraine. TIA made up a larger "
      "share of events in adenomyosis and endometriosis, and women with TIA more often had migraine than women with "
      "ischemic stroke. Migraine with aura is a recognized TIA mimic [ref], so some TIA diagnoses in these women may "
      "represent migraine rather than cerebral ischemia. Studies of stroke in endometriosis and adenomyosis that "
      "include TIA should consider this source of misclassification.")
    H2("Strengths and limitations")
    P("Strengths include the large cohort, comparison between conditions within one health system, age "
      "standardization, pre-specified models and correction for multiple comparisons.")
    P(f"This is a cross-sectional analysis of records; the timing of risk factors relative to the uterine diagnosis "
      f"was not available, so associations cannot establish order or cause. Diagnoses were taken from records and "
      f"not adjudicated. Smoking status was unknown for most women (known for {R['p3_smoking_known']:,}), and BMI was "
      f"missing for some. Migraine may be under-recorded, and ascertainment may differ between conditions because of "
      f"different care pathways. The groups were defined by recorded diagnoses, and some women with fibroids may have "
      f"undiagnosed endometriosis or adenomyosis. Follow-up for stroke was short (median {pc['fu_median']:.1f} "
      f"years), events in the adenomyosis group were few, and TIA diagnoses were not adjudicated. Migraine was recorded at "
      f"any time, so its order relative to stroke is not established. The results come from a single tertiary center.")
    H2("Conclusions")
    P("Among women with benign uterine disease, adenomyosis and endometriosis were associated with more migraine and, "
      "in younger women with endometriosis, with fewer cardiometabolic risk factors than fibroids. Migraine should be "
      "considered alongside conventional risk factors when assessing vascular risk in these women. Stroke rates after "
      "diagnosis were similar across conditions, and TIA diagnoses in women with migraine should be interpreted "
      "with care.")
    H2("Acknowledgements, funding, disclosures")
    P("[To be completed by the authors.]")
    BR()

    # ------------------------------------------------------------------ tables
    H1("Tables")
    t1s = T["P3_Table1"].copy()
    t1s = t1s[t1s.Characteristic.fillna("") != "Any stroke"].fillna("")
    TABLE("Table 1. Characteristics by uterine diagnosis group", t1s,
          "Values are n (%) or median (IQR). P: chi-square or Kruskal–Wallis test. Smoking: unknown status excluded "
          "from percentages.", widths=[2300, 1500, 1200, 1100, 1100, 1100, 1100, 700], font=7)
    rows = []
    outs = list(dict.fromkeys(main.Outcome))
    for o in outs:
        r = {"Risk factor": o, "Fibroids only, %": pv(o, "Fibroids only")["Age-standardised % (95% CI)"]}
        for g, s in [(aden, "Adenomyosis only"), (endo, "Endometriosis only"), (multi, ">1 condition")]:
            r[f"{s}, %"] = pv(o, g)["Age-standardised % (95% CI)"]
            r[f"{s}, OR"] = orx(o, g)
        rows.append(r)
    TABLE("Table 2. Age-standardized prevalence (%) and adjusted odds ratios (vs fibroids only) of vascular risk "
          "factors", pd.DataFrame(rows),
          "Prevalence directly age-standardized (5-year bands, whole cohort). OR (95% CI) from logistic regression "
          "adjusted for age, race and BMI (BMI not adjusted for obesity); robust standard errors. Benjamini–Hochberg-"
          "adjusted P values in eTable 1.", widths=[1700, 1150, 1150, 1100, 1150, 1100, 1150, 1100], font=7)
    rows = []
    for o in CARDIO + [mig]:
        r = {"Risk factor": o}
        for g in [endo, aden]:
            for a in ["18–39", "40–60"]:
                r[f"{g.split()[0]}, {a} y"] = sx(o, a, g)["OR (95% CI)"]
        r["Interaction P (BH)"] = em.loc[o, "p BH (text)"]
        rows.append(r)
    TABLE("Table 3. Adjusted odds ratios by age group, endometriosis only and adenomyosis only vs fibroids only",
          pd.DataFrame(rows), "Models adjusted for age, race and BMI within each age group. Interaction P: condition × "
                              "age group (3 df), Benjamini–Hochberg-adjusted across 10 risk factors.",
          widths=[2000, 1500, 1500, 1500, 1500, 1360], font=7)
    BR()

    # ------------------------------------------------------------------ figures
    H1("Figures")
    caps = [
        "Figure 1. Age-standardized prevalence of migraine by uterine diagnosis group (95% CI), with adjusted odds "
        "ratios vs fibroids only.",
        "Figure 2. Adjusted odds ratios (95% CI) of vascular risk factors by uterine diagnosis group vs fibroids only "
        "(adjusted for age, race and BMI).",
        "Figure 3. Adjusted odds ratios by age group (18–39 vs 40–60 years). (A) Endometriosis only and (B) adenomyosis "
        "only, each vs fibroids only.",
        "Figure 4. Adjusted odds ratios for migraine vs fibroids only, overall and by age group.",
        "Figure 5. Stroke or TIA and ischemic stroke after diagnosis by uterine diagnosis group. (A) Crude rates per "
        "1,000 person-years (exact 95% CI). (B) Rate ratios vs fibroids only, adjusted for age, race and cardiometabolic "
        "factors, and additionally for migraine.",
    ]
    for (fn, w), cap in zip(FIGS, caps):
        FIG(fn, cap, w)
    BR()

    # ------------------------------------------------------------------ supplement
    H1("Supplementary material")
    sup = [
        ("eTable 1. Adjusted odds ratios, Model 1 and Model 2, with Benjamini–Hochberg-adjusted P values",
         main[["Outcome", "Group vs fibroids only", "M1 OR (95% CI) age+race", "M2 OR (95% CI) age+race+BMI",
               "M2 p (BH-adjusted, all contrasts) (text)", "N M2", "Events M2"]]
         .rename(columns={"M1 OR (95% CI) age+race": "Model 1 OR", "M2 OR (95% CI) age+race+BMI": "Model 2 OR",
                          "M2 p (BH-adjusted, all contrasts) (text)": "P (BH)", "N M2": "N", "Events M2": "Events"}), ""),
        ("eTable 2. Crude and age-standardized prevalence by group",
         prev[["Outcome", "Group", "n with outcome", "N (outcome known)", "Crude %", "Age-standardised % (95% CI)"]], ""),
        ("eTable 3. Age-stratified models (all risk factors)",
         st[["Outcome", "Age", "Group vs fibroids only", "n in group", "events in group", "OR (95% CI)"]], ""),
        ("eTable 4. Migraine: sensitivity to hormonal therapy and uterine bleeding",
         ms[["Model", "Group vs fibroids only", "N", "Migraine cases", "OR (95% CI)", "p"]], ""),
        ("eTable 5. Excluding women with stroke before or at the index date",
         inc[["Outcome", "Group vs fibroids only", "Primary M2 OR", "Excluding stroke before/at index", "N", "Events"]],
         ""),
        ("eTable 6. Stroke after diagnosis: rate ratios vs fibroids only by model and outcome",
         pm[["Outcome", "Model", "Group vs fibroids only", "Events in group", "RR (95% CI)", "p (text)", "N", "Events"]]
         .rename(columns={"p (text)": "P"}),
         "M1: age, race. M2: + BMI, hypertension, diabetes, dyslipidemia, smoking, atrial fibrillation, coronary artery "
         "disease. M3: + migraine. M4: + hormonal therapy, heavy/abnormal bleeding."),
        ("eTable 7. Stroke type by group, and migraine among women with TIA vs ischemic stroke",
         pd.concat([R["p4_type"], R["p4_mig_event"].rename(columns={"Event type": "Group"})], ignore_index=True), ""),
    ]
    for title, df, note in sup:
        TABLE(title, df.fillna(""), note, font=7)
    failed = [c for c, ok in checks if not ok]
    RESULTS["manuscript_p3_checks"] = checks
    B.append({"t": "note", "text": "Automated wording checks: " + ("all passed" if not failed else "FAILED: " +
                                                                   "; ".join(failed))})
    if failed:
        log("Manuscript P3", "WORDING CHECK FAILED: " + "; ".join(failed))
    return B


def short_version(B):
    out = [b for b in B if b["t"] == "title"]
    out.append({"t": "note", "text": f"AAN 2027 abstract format: body {RESULTS['aan_p3_words']} words (limit 300, "
                                     "headings excluded). Author names and affiliations are entered in the "
                                     "submission system."})
    out.append({"t": "h1", "text": "Abstract"})
    out += [{"t": "p", "text": f"**{h}:** {t}"} for h, t in RESULTS["aan_p3"]]
    out += [{"t": "pagebreak"}, {"t": "h1", "text": "Tables"}]
    out += [b for b in B if b["t"] == "table" and re.match(r"Table [123]\.", b["title"])]
    out += [{"t": "pagebreak"}, {"t": "h1", "text": "Figures"}]
    out += [b for b in B if b["t"] == "figure"]
    return out


def _write(B, stem):
    (OUT_DIR / f"{stem}.json").write_text(json.dumps(B, ensure_ascii=False), encoding="utf-8")
    (OUT_DIR / f"{stem}.md").write_text(to_markdown(B), encoding="utf-8")
    r = subprocess.run(["node", "analysis/build_docx.js", str(OUT_DIR / f"{stem}.json"), str(OUT_DIR / f"{stem}.docx")],
                       capture_output=True, text=True)
    log("Manuscript P3", f"{stem}.docx: {'ok' if r.returncode == 0 else r.stderr[-500:]}")


def run():
    B = _us(build())
    _write(B, "migraine_vascular_manuscript")
    _write(_us(short_version(B)), "migraine_vascular_abstract_tables_figures")
