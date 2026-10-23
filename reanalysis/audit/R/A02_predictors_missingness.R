# Audit Q3.4–3.5 and Q4: timing of predictor ascertainment relative to time zero and to the stroke,
# missingness, selection into the development sample, and smoking.
source("audit/R/00_audit_setup.R")
md_open("A02_predictors_missingness.R", "A02 Predictor timing, missing data and selection")
co <- add_vars(load_cohort()); d <- load_dev()
cm <- readRDS(file.path(DERIVED, "comorbid_long.rds"))

# ---- 1. documented ascertainment window of each MED-HACK predictor ------------------------------------------
win <- data.table(
  predictor = c("Migraine", "Estrogen-containing therapy", "Diabetes", "Hypertension", "Hemoglobinopathy", "Hemoglobin", "Age",
                "Coronary artery disease", "Chronic kidney disease", "BMI (14-predictor model)", "Smoking (14-predictor model)", "Dyslipidemia (14-predictor model)"),
  source = c("Codebook: Problem List + Migraine Diagnosis extracts; no date restriction stated", "Codebook: Administered Medications, window = index_date − 365 days ONWARD (no end date)",
             "Codebook: no timing stated", "Codebook: documented on or before index_date", "R/13_anchor_design.R:65 `any_date = TRUE` (diagnosis extracts 22–31, any date)",
             "R/13_anchor_design.R:70 latest CBC in the 3 years up to and including the anchor", "R/13_anchor_design.R:29 age at anchor",
             "Codebook: qualifying diagnosis on or before index_date", "R/13_anchor_design.R:66 dated on or before the anchor", "Codebook: closest value within ±12 months of index_date",
             "Codebook: latest statement on/before index_date, else EARLIEST AFTER", "dyslipidemia_update sheet: code or abnormal lipid on or before index_date"),
  relative_to_anchor = c("unknown — may postdate anchor and stroke", "postdates anchor by design", "unknown", "≤ index (index = anchor for 98.6%)", "may postdate anchor and stroke",
                         "≤ anchor", "at anchor", "≤ index", "≤ anchor", "up to 12 months after index", "may postdate anchor and stroke", "≤ index"),
  raw_dated_source_in_project = c("no", "no (meds/ file holds iron and spironolactone only)", "no", "no", "yes", "yes", "yes", "no", "yes", "no", "yes (smoking_source sheet)", "partly"))
aud_tab(win, "A02a_ascertainment_windows", "When each predictor is measured, from the Codebook notes and the scripts")

# ---- 2. hemoglobinopathy: dated codes ---------------------------------------------------------------------------
HBP <- c("sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy")
fh <- cm[flag %in% HBP, .(first_hbp = min(date)), by = study_id]
dh <- merge(d[, .(study_id, anchor, y, stroke_date, haemoglobinopathy, t_end)], fh, by = "study_id", all.x = TRUE)
hb_t <- dh[haemoglobinopathy == 1, .(women = .N, strokes = sum(y),
  first_code_on_or_before_anchor = sum(first_hbp <= anchor), first_code_after_anchor = sum(first_hbp > anchor),
  strokes_first_code_after_anchor = sum(y == 1 & first_hbp > anchor), strokes_first_code_on_or_after_stroke_minus_30d = sum(y == 1 & first_hbp >= stroke_date - 30, na.rm = TRUE))]
aud_tab(hb_t, "A02b_hemoglobinopathy_timing", "Hemoglobinopathy-positive women in the development sample: date of the first qualifying code")
d[, hbp_pre := as.integer(study_id %in% dh[first_hbp <= anchor, study_id])]
f_orig <- coxph(fml_of(PRED9), data = d); f_pre <- coxph(fml_of(sub("haemoglobinopathy", "hbp_pre", PRED9)), data = d)
hr <- function(f, v) sprintf("%.2f (%.2f–%.2f)", exp(coef(f)[v]), exp(confint(f)[v, 1]), exp(confint(f)[v, 2]))
aud_tab(data.table(definition = c("any date (as in the calculator)", "first code on or before the anchor"),
                   women_positive = c(sum(d$haemoglobinopathy), sum(d$hbp_pre)), strokes_positive = c(d[haemoglobinopathy == 1, sum(y)], d[hbp_pre == 1, sum(y)]),
                   hazard_ratio = c(hr(f_orig, "haemoglobinopathy"), hr(f_pre, "hbp_pre")),
                   harrell_c = c(f2(concordance(f_orig)$concordance, 3), f2(concordance(f_pre)$concordance, 3))),
        "A02c_hemoglobinopathy_redefined", "Nine-item Cox model (unshrunk) with hemoglobinopathy defined from any-date vs pre-anchor codes")

# ---- 3. indirect look-ahead check for predictors without dated sources ------------------------------------------
# among women without a stroke/TIA, does a predictor's prevalence rise with observed follow-up within the same anchor year?
nz <- d[ev_any == 0 & t > 0]; nz[, `:=`(log2fu = log2(t), ay = factor(year(anchor)))]
la <- rbindlist(lapply(c("migraine_any", "estrogen", "dm", "htn", "cad", "dyslipidemia", "haemoglobinopathy", "ckd_pre"), function(v) {
  g <- glm(as.formula(paste(v, "~ log2fu + ay + age")), data = nz, family = binomial)
  ci <- exp(confint.default(g)["log2fu", ])
  data.table(predictor = v, dated_on_or_before_baseline = v %in% c("htn", "cad", "dyslipidemia", "ckd_pre"),
             prevalence_fu_lt_1y = pct(nz[t < 1, mean(get(v))]), prevalence_fu_ge_3y = pct(nz[t >= 3, mean(get(v))]),
             or_per_doubling_of_followup = sprintf("%.2f (%.2f–%.2f)", exp(coef(g)["log2fu"]), ci[1], ci[2]))
}))
aud_tab(la, "A02d_lookahead_check", sprintf("Women with no stroke/TIA (n = %s). Logistic model: predictor ~ log2(follow-up) + anchor year + age. Predictors dated on or before baseline act as negative controls; a much larger OR for an undated predictor suggests it was partly recorded after baseline. This is an indirect check, not proof.", format(nrow(nz), big.mark = ",")))

# ---- 4. smoking ---------------------------------------------------------------------------------------------
ss <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "smoking_source"))[, .(mrn = as.numeric(MRN), smk_code = smoking, statement_date = as.IDate(substr(as.character(statement_date), 1, 10)))]
ds <- merge(d[, .(mrn, study_id, smoking3, smoking, anchor, y, stroke_date)], ss, by = "mrn", all.x = TRUE)
smk <- rbind(
  ds[, .(item = c("development sample", "status unknown (code 9)", "never", "current", "former"),
         n = c(.N, sum(smoking == 9), sum(smoking == 0), sum(smoking == 1), sum(smoking == 2)))][, pct := pct(n / n[1])],
  data.table(item = c("documented status with a statement date", "  statement after the anchor", "  statement more than 1 year after the anchor",
                      "women with stroke and documented status", "  statement on or after the stroke date"),
             n = c(ds[smoking != 9 & !is.na(statement_date), .N], ds[smoking != 9 & statement_date > anchor, .N], ds[smoking != 9 & statement_date > anchor + 365, .N],
                   ds[y == 1 & smoking != 9, .N], ds[y == 1 & smoking != 9 & statement_date >= stroke_date, .N]), pct = ""))
aud_tab(smk, "A02e_smoking", "Smoking status in the development sample and timing of the source statement (workbook sheet smoking_source)")
d[, smk4 := relevel(factor(c("0" = "never", "1" = "current", "2" = "former", "9" = "unknown")[as.character(smoking)]), "never")]
d[, smk_known := smoking != 9]
sm <- rbindlist(list(
  { f <- coxph(fml_of(c(PRED9, "smoke")), data = d); data.table(coding = "ever vs never-or-unknown (SAP coding)", level = "ever", hr = hr(f, "smoke"), women = nrow(d), strokes = sum(d$y)) },
  { f <- coxph(fml_of(c(PRED9, "smk4")), data = d); rbindlist(lapply(c("current", "former", "unknown"), function(l) data.table(coding = "never / current / former / unknown", level = l, hr = hr(f, paste0("smk4", l)), women = nrow(d), strokes = sum(d$y)))) },
  { dk <- d[smk_known == TRUE]; f <- coxph(fml_of(c(PRED9, "smoke")), data = dk); data.table(coding = "known status only: ever vs never", level = "ever", hr = hr(f, "smoke"), women = nrow(dk), strokes = sum(dk$y)) }))
aud_tab(sm, "A02f_smoking_models", "Smoking added to the nine-item Cox model (unshrunk), three codings")

# ---- 5. missingness of every candidate predictor (anchored cohort) --------------------------------------------
mis <- rbindlist(lapply(list(c("hb", "CBC hemoglobin in 3 y up to anchor"), c("bmi", "BMI"), c("smoking3", "smoking"), c("hormonal_type_cat", "hormonal type"),
                             c("migraine", "migraine (Codebook code 9 = aura not stated, counted as migraine)")), function(z) {
  x <- co[[z[1]]]; miss <- if (z[1] == "smoking3") x == "Unknown" else if (z[1] == "migraine") co$migraine == 9 else is.na(x)
  data.table(variable = z[2], missing_or_unknown = sum(miss), pct = pct(mean(miss)))}))
mis <- rbind(mis, data.table(variable = c("hypertension, diabetes, dyslipidemia, CAD, CKD, hemoglobinopathy, endometriosis, adenomyosis"), missing_or_unknown = NA, pct = "absence of a code is coded 0; true missingness cannot be distinguished"))
aud_tab(mis, "A02g_missingness", sprintf("Missing or unknown values among the %s women of the anchored cohort", format(nrow(co), big.mark = ",")))

# ---- 6. included vs excluded --------------------------------------------------------------------------------
co[, grp := fifelse(is.na(hb), "no Hb (excluded)", fifelse(is.na(bmi), "Hb, no BMI (excluded)", "development sample"))]
co[, `:=`(black = as.integer(race4 == "Black"), smk_doc = as.integer(smoking3 != "Unknown"), any_hormonal = as.integer(hormonal_type_cat != "None"),
          fibroid = as.integer(fibroids == 1), endo_any = as.integer(endometriosis == 1), adeno_any = as.integer(adenomyosis == 1))]
vars <- c(age = "age at anchor, mean", black = "Black race", fibroid = "fibroids", endo_any = "endometriosis", adeno_any = "adenomyosis", htn = "hypertension", dm = "diabetes",
          dyslipidemia = "dyslipidemia", cad = "coronary artery disease", migraine_any = "migraine", ckd_pre = "CKD", haemoglobinopathy = "hemoglobinopathy",
          estrogen = "estrogen-containing therapy", any_hormonal = "any hormonal therapy", smk_doc = "smoking status documented", t = "follow-up years, mean")
inc <- co[grp == "development sample"]; exc <- co[grp == "no Hb (excluded)"]
ie <- rbindlist(lapply(names(vars), function(v) data.table(characteristic = vars[v],
  development_sample = if (v %in% c("age", "t")) f2(mean(inc[[v]], na.rm = TRUE), 1) else pct(mean(inc[[v]], na.rm = TRUE)),
  excluded_no_hb = if (v %in% c("age", "t")) f2(mean(exc[[v]], na.rm = TRUE), 1) else pct(mean(exc[[v]], na.rm = TRUE)),
  standardized_difference = f2(smd(inc[[v]], exc[[v]]), 2))))
rt <- co[, .(women = .N, strokes = sum(ev_stroke), py = sum(py)), by = grp]
ie <- rbind(ie, data.table(characteristic = "stroke rate per 1,000 person-years (95% CI)", development_sample = rt[grp == "development sample", pois_ci(strokes, py)],
                           excluded_no_hb = rt[grp == "no Hb (excluded)", pois_ci(strokes, py)], standardized_difference = ""),
            data.table(characteristic = "women / strokes", development_sample = rt[grp == "development sample", paste(women, "/", strokes)],
                       excluded_no_hb = rt[grp == "no Hb (excluded)", paste(women, "/", strokes)], standardized_difference = ""))
aud_tab(ie, "A02h_included_vs_excluded", "Development sample vs anchored-cohort women excluded for having no CBC Hb (percentages unless stated). |SMD| > 0.10 is conventionally a meaningful imbalance.")
# selection model
g <- glm(I(!is.na(hb)) ~ age + black + fibroid + endo_any + adeno_any + htn + dm + dyslipidemia + cad + migraine_any + ckd_pre + haemoglobinopathy + estrogen + smk_doc, data = co, family = binomial)
sel <- data.table(term = names(coef(g))[-1], or_for_having_hb = sprintf("%.2f (%.2f–%.2f)", exp(coef(g))[-1], exp(confint.default(g))[-1, 1], exp(confint.default(g))[-1, 2]))
aud_tab(sel, "A02i_selection_model", "Logistic model for having a CBC Hb (inclusion), anchored cohort, mutually adjusted")
# outcome-associated selection: stroke hazard in excluded women relative to included, adjusted
co[, excl := as.integer(is.na(hb))]
fx <- coxph(Surv(py, ev_stroke) ~ excl + age + htn + dm + cad + migraine_any + ckd_pre + haemoglobinopathy + estrogen, data = co)
md_note("Adjusted stroke hazard ratio, women without vs with a CBC Hb (anchored cohort, Cox with the MED-HACK binary items and age): ", hr(fx, "excl"), ".")
# apply the published model to excluded women under the assumption Hb ≥13 (deficit 0) — exploratory transportability check
ex <- copy(exc); ex[, hb_def := 0]; ex[, `:=`(t = py, y = ev_stroke)]
ex <- ex[age >= 18]
ex[, p5 := risk_spec(ex)]; ex[, p3 := risk_spec(ex, h_s0 = s0_at(coxph(fml_of(PRED9), data = d), 3))]
o5 <- km_risk(ex, 5); o3 <- km_risk(ex, 3)
aud_tab(data.table(horizon = c("3 years", "5 years"), women = nrow(ex), strokes = sum(ex$y), mean_predicted_pct = pct(c(mean(ex$p3), mean(ex$p5)), 2),
                   observed_km_pct = c(sprintf("%s (%s–%s)", pct(o3[1], 2), pct(o3[2], 2), pct(o3[3], 2)), sprintf("%s (%s–%s)", pct(o5[1], 2), pct(o5[2], 2), pct(o5[3], 2))),
                   observed_over_expected = f2(c(o3[1] / mean(ex$p3), o5[1] / mean(ex$p5)), 2)),
        "A02j_model_in_excluded_women", "Exploratory: the published nine-item model applied to women excluded for no Hb, with Hb deficit set to 0 (an assumption, not data). 3-year baseline survival from the unshrunk nine-item refit.")
res_put("A02", list(hbp = hb_t, lookahead = la))
