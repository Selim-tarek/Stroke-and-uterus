# Part 3: predictor provenance, strict baseline reconstruction, outcome flags, and an exploratory detection/utilization-bias analysis.
# Writes the patient-level analysis file derived/post_audit/baseline_cohort.rds (git-ignored).
source("post_audit/R/P00_setup.R")
md_open("P02_baseline_reconstruction.R", "P02 Predictor provenance and strict baseline reconstruction")
co <- add_vars(load_cohort())
require_cols(co, c("mrn", "anchor", "index_date", "t_end", "hb", "age_anchor", "stroke_date", "stroke_type", "source_kinds", "notes", "smoking", "ev_stroke", "ev_any",
                   "migraine_any", "estrogen", "dm", "htn", "cad", "ckd_pre", "dyslipidemia", "haemoglobinopathy", "fibroids", "dxgrp"), "anchored cohort")
cm <- readRDS(file.path(DERIVED, "comorbid_long.rds")); hb <- readRDS(file.path(DERIVED, "labs_long.rds"))$hb[src == "CBC", .(study_id, date)]
co[, `:=`(t = py, y = ev_stroke)]

# ---- 1. cohort: CBC-tested, anchored at the first coded diagnosis (so index-dated covariates are at the anchor), age >= 18 -------------
# BMI is not a candidate predictor, so the BMI requirement of the published model is dropped (removes 880 needless exclusions).
b <- co[!is.na(hb) & age_anchor >= 18]
fl <- data.table(step = c("anchored cohort with CBC Hb, age >= 18", "excluded: anchor moved earlier by a procedure (index-dated covariates postdate the anchor)", "baseline cohort"),
                 n = c(nrow(b), b[anchor < index_date, .N], b[anchor == index_date, .N]),
                 strokes = c(sum(b$y), b[anchor < index_date, sum(y)], b[anchor == index_date, sum(y)]))
b <- b[anchor == index_date]
fl <- rbind(fl, data.table(step = "  of which in the published development sample (BMI recorded)", n = b[!is.na(bmi), .N], strokes = b[!is.na(bmi), sum(y)]))
aud_tab(fl, "P02a_baseline_cohort", "Baseline cohort for the post-audit candidate models")

# ---- 2. baseline-valid predictors ---------------------------------------------------------------------------------------------
HBP <- c("sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy")
fh <- cm[flag %in% HBP, .(first_hbp = min(date)), by = study_id]
b <- merge(b, fh, by = "study_id", all.x = TRUE)
b[, hbp_pre := as.integer(!is.na(first_hbp) & first_hbp <= anchor)]
ss <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "smoking_source"))
require_cols(ss, c("MRN", "smoking", "statement_date"), "smoking_source sheet")
ss <- ss[, .(mrn = as.numeric(MRN), smk_code = as.integer(smoking), statement_date = as.IDate(substr(as.character(statement_date), 1, 10)))]
b <- merge(b, ss, by = "mrn", all.x = TRUE)
b[, smoke_status_b := factor(fifelse(!is.na(statement_date) & statement_date <= anchor & smk_code %in% 0:2, c("never", "current", "former")[smk_code + 1], "unknown at baseline"),
                              levels = c("never", "current", "former", "unknown at baseline"))]
b[, smoke_b := as.integer(smoke_status_b %in% c("current", "former"))]

# ---- 3. utilization proxies (baseline: before the anchor; post-anchor rates are FUTURE information, exploratory only) -------------
cb <- merge(unique(hb), b[, .(study_id, anchor, t_end)], by = "study_id")
u_pre <- cb[date < anchor & date >= anchor - 730, .(cbc_days_2y_pre = uniqueN(date)), by = study_id]
u_post <- cb[date > anchor & date <= t_end, .(cbc_days_post = uniqueN(date)), by = study_id]
dxu <- merge(unique(cm[, .(study_id, date)]), b[, .(study_id, anchor)], by = "study_id")[date < anchor & date >= anchor - 730, .(dx_days_2y_pre = uniqueN(date)), by = study_id]
b <- Reduce(function(a, z) merge(a, z, by = "study_id", all.x = TRUE), list(b, u_pre, u_post, dxu))
for (v in c("cbc_days_2y_pre", "cbc_days_post", "dx_days_2y_pre")) b[is.na(get(v)), (v) := 0L]
b[, cbc_rate_post := cbc_days_post / pmax(t, 1 / 12)]

# ---- 4. outcome flags (see P03 for definitions) ---------------------------------------------------------------------------------
rs <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "radiology_stroke"))[, .(mrn = as.numeric(MRN), img_date = as.IDate(substr(as.character(earliest_positive_date), 1, 10)))]
nf <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "nihss_flowsheet"))[, .(mrn = as.numeric(MRN), nih_date = as.IDate(substr(as.character(first_recorded), 1, 10)))]
b <- merge(merge(b, rs, by = "mrn", all.x = TRUE), nf, by = "mrn", all.x = TRUE)
b[, `:=`(gap = as.numeric(stroke_date - anchor), acute = grepl("acute", source_kinds), imaging_only = grepl("CASE FOUND ON IMAGING ONLY", notes),
         corroborated = (!is.na(img_date) & abs(as.numeric(img_date - stroke_date)) <= 30) | (!is.na(nih_date) & abs(as.numeric(nih_date - stroke_date)) <= 30))]
b[, `:=`(y_all = y, y_acute = as.integer(y == 1 & acute & !imaging_only))]
b[, `:=`(y_acute30 = as.integer(y_acute == 1 & gap > 30), y_isch = as.integer(y_acute == 1 & stroke_type == 1),
         y_hem = as.integer(y_acute == 1 & stroke_type %in% 2:3), y_corr = as.integer(y_acute == 1 & corroborated))]
b[, y := y_acute]                                                     # pre-specified primary post-audit outcome
saveRDS(b, BASE_RDS); log_n(b, "baseline cohort, primary outcome (acute-coded stroke)")

# ---- 5. provenance quantification -----------------------------------------------------------------------------------------------
pv <- data.table(predictor = c("hemoglobinopathy", "smoking (ever)"), as_coded = c(sum(b$haemoglobinopathy), b[smoking %in% 1:2, .N]),
                 baseline_valid = c(sum(b$hbp_pre), sum(b$smoke_b)),
                 reclassified_to_0 = c(b[haemoglobinopathy == 1 & hbp_pre == 0, .N], b[smoking %in% 1:2 & smoke_b == 0, .N]),
                 strokes_among_reclassified = c(b[haemoglobinopathy == 1 & hbp_pre == 0, sum(y)], b[smoking %in% 1:2 & smoke_b == 0, sum(y)]))
aud_tab(pv, "P02b_reconstructed_predictors", "Predictors that can be re-dated with data in the project: as coded vs restricted to information dated on or before the anchor (baseline cohort; strokes = primary outcome)")
aud_tab(b[, .N, by = smoke_status_b][order(smoke_status_b)][, pct := pct(N / sum(N))], "P02c_smoking_baseline", "Smoking status known at the anchor (statement dated on or before the anchor)")

# ---- 6. exploratory detection / utilization bias (NOT a correction) ---------------------------------------------------------------
nz <- b[ev_any == 0 & t > 0.25]; nz[, ay := factor(year(anchor))]
db <- rbindlist(lapply(c("migraine_any", "dm", "estrogen", "haemoglobinopathy", "hbp_pre", "htn", "cad", "ckd_pre", "dyslipidemia"), function(v) {
  g1 <- glm(as.formula(paste(v, "~ log2(1 + cbc_days_2y_pre) + age + ay")), data = nz, family = binomial)
  g2 <- glm(as.formula(paste(v, "~ log2(1 + cbc_rate_post) + log2(1 + cbc_days_2y_pre) + age + ay")), data = nz, family = binomial)
  e <- function(g, k) sprintf("%.2f (%.2f–%.2f)", exp(coef(g)[k]), exp(confint.default(g)[k, 1]), exp(confint.default(g)[k, 2]))
  data.table(predictor = v, dated_at_or_before_anchor = v %in% c("hbp_pre", "htn", "cad", "ckd_pre", "dyslipidemia"),
             or_baseline_cbc_doubling = e(g1, 2), or_post_anchor_cbc_rate_doubling_adj_baseline = e(g2, 2)) }))
aud_tab(db, "P02d_utilization_association", sprintf("Exploratory. Women without stroke/TIA followed > 3 months (n = %s). Association of each predictor with baseline CBC testing (2 years before the anchor) and with FUTURE testing intensity (CBC days per year after the anchor, adjusted for baseline testing, age and anchor year). A predictor recorded at baseline should not depend on future testing; a strong future association suggests ascertainment after the anchor.", format(nrow(nz), big.mark = ",")))
f0 <- coxph(fml_of(CAND$C2_current_nine_as_coded), data = b); f1 <- coxph(update(fml_of(CAND$C2_current_nine_as_coded), . ~ . + log2(1 + cbc_days_2y_pre) + log2(1 + dx_days_2y_pre)), data = b)
hr <- function(f, v) sprintf("%.2f (%.2f–%.2f)", exp(coef(f)[v]), exp(confint(f)[v, 1]), exp(confint(f)[v, 2]))
aud_tab(data.table(predictor = CAND$C2_current_nine_as_coded, hr_current_nine = sapply(CAND$C2_current_nine_as_coded, hr, f = f0),
                   hr_adjusted_for_baseline_utilization = sapply(CAND$C2_current_nine_as_coded, hr, f = f1)),
        "P02e_hr_utilization_adjusted", "Exploratory. Current nine predictors (as coded) in the baseline cohort, primary outcome, with and without adjustment for baseline utilization (CBC days and coded-diagnosis days in the 2 years before the anchor)")
res_put("P02", list(flow = fl, reconstructed = pv, utilization = db))
