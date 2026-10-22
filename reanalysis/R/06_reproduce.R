# Phase 4: reproduce the ORIGINAL analyses (Python pipeline) in R, with the original definitions,
# and compare every estimate with the value exported from the original results (tools/export_original_results.py).
source("R/00_setup.R")
e  <- readRDS(file.path(DERIVED, "eligible.rds"))
cm <- readRDS(file.path(DERIVED, "comorbid_long.rds"))
labs <- readRDS(file.path(DERIVED, "labs_long.rds"))
orig <- fread(file.path(OUT, "original_reported.csv"))
LV <- c(Mild = "Mild (10–11.9)", Moderate = "Moderate (8–9.9)", Severe = "Severe (<8)")

# ---- original 2a flags (reference date = Codebook Hb date, else index; definitions exactly as in the Python code) ----
e[, ref := fifelse(is.na(hgb_date), index_date, hgb_date)]
k <- merge(cm, e[, .(study_id, ref, hgb_date)], by = "study_id")
has <- function(sel) as.integer(e$study_id %in% unique(k[sel]$study_id))
e[, `:=`(
  o_haemoglobinopathy = has(k$flag %in% c("sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy")),
  o_ckd = has(k$flag == "ckd" & k$date <= k$ref),                         # 'ckd' flag only (ESRD-only codes not included)
  o_liver = has(k$flag %in% c("cirrhosis", "chronic_liver") & k$date <= k$ref),
  o_alcohol = has(k$flag == "alcohol_disorder" & k$date <= k$ref),
  o_gibleed = has(k$flag == "gi_bleed" & !is.na(k$hgb_date) & k$date >= k$hgb_date - 365 & k$date <= k$hgb_date + 30),
  o_malabs = has(k$flag == "malabsorption" & k$date <= k$ref),
  o_ibd = has(k$flag == "ibd" & k$date <= k$ref),
  o_hiv = has(k$flag == "hiv"),                                            # any date
  o_preg = has((k$flag == "pregnancy" & !is.na(k$hgb_date) & k$date >= k$hgb_date - 365 & k$date <= k$hgb_date + 30) |
               (k$flag == "anaemia_pregnancy" & !is.na(k$hgb_date) & abs(k$date - k$hgb_date) <= 365)))]
A2 <- c("cad", "o_haemoglobinopathy", "o_ckd", "o_liver", "o_alcohol", "o_gibleed", "o_malabs", "o_ibd", "o_hiv", "o_preg")
COV0 <- c("age_index", "race4", "bmi", "htn", "dm", "dyslipidemia", "smoking3", "migraine_any", "thrombophilia")

e[, score := as.numeric(grade_cb) - 1]
e[, `:=`(y_any = as.integer(stroke_any == 1),
         y_isch = fifelse(stroke_any == 0, 0L, fifelse(stroke_type == 1, 1L, NA_integer_)),
         y_inc = fifelse(stroke_any == 0, 0L, fifelse(stroke_timing %in% 3, 1L, NA_integer_)))]

# hormonal-type sparse rule: in the complete-case sample, levels (other than None / Other) with <5 events merge into Other/multiple
ht_merge <- function(d, y) {
  ev <- d[, .(ev = sum(get(y))), by = hormonal_type_cat]
  small <- setdiff(as.character(ev[ev < 5, hormonal_type_cat]), c("None", "Other/multiple"))
  small <- union(small, setdiff(levels(d$hormonal_type_cat)[-c(1, 7)], as.character(ev$hormonal_type_cat)))
  h <- as.character(d$hormonal_type_cat); h[h %in% small] <- "Other/multiple"
  factor(h, levels = setdiff(levels(d$hormonal_type_cat), small))
}
REP <- list()
fit_xs <- function(d, y, covs, label, outcome_key, crude = FALSE) {
  need <- unique(c(y, "grade_cb", if (!crude) c(covs, "hormonal_type_cat")))
  dd <- d[complete.cases(d[, ..need])]
  rhs_cov <- if (crude) character() else c(covs, "ht")
  if (!crude) dd[, ht := ht_merge(dd, y)]
  f <- with_warnings(glm(reformulate(c("grade_cb", rhs_cov), y), binomial, data = dd), label)
  ft <- glm(reformulate(c("score", rhs_cov), y), binomial, data = dd)
  t <- robust_tab(f, vcovHC(f, type = "HC1"), paste0("grade_cb", names(LV)))
  tt <- robust_tab(ft, vcovHC(ft, type = "HC1"), "score")
  out <- rbind(t[, .(level = LV[sub("grade_cb", "", term)], est, lo, hi, p)], tt[, .(level = "Per grade (trend)", est, lo, hi, p)])
  out[, `:=`(analysis = label, outcome = outcome_key, N = nrow(dd), events = sum(dd[[y]]),
             epv = sum(dd[[y]]) / (length(coef(f)) - 1), converged = f$converged)]
  REP[[length(REP) + 1]] <<- out
  invisible(out)
}
OUTC <- list(c("y_any", "Any stroke (primary)", "Any stroke"), c("y_isch", "Ischaemic stroke vs no stroke", "Ischaemic stroke"),
             c("y_inc", "Incident stroke (after index) vs no stroke", "Incident stroke"))
for (o in OUTC) {
  fit_xs(e, o[1], COV0, "Cross-sectional logistic, crude", o[2], crude = TRUE)
  fit_xs(e, o[1], COV0, "Cross-sectional logistic, adjusted", o[2])
  fit_xs(e, o[1], c(COV0, A2), paste0("Sensitivity: 2a Pre-Hb conditions: ", o[3]), o[3])
}
fit_xs(e, "y_any", c(COV0, "uterine_bleeding", "dxgrp"), "Sensitivity: + uterine bleeding and diagnosis group: Any stroke", "Any stroke")
# cases restricted to Hb >= 30 d before the stroke (controls: all women without stroke)
e[, hb_before30 := fifelse(stroke_any == 0, 1L, as.integer(!is.na(stroke_date) & !is.na(hgb_date) & as.integer(stroke_date - hgb_date) >= 30))]
fit_xs(e[hb_before30 == 1], "y_any", COV0, "Sensitivity: Hb ≥30 d before stroke: Any stroke", "Any stroke")
fit_xs(e[hb_before30 == 1], "y_any", c(COV0, A2), "Sensitivity: Hb ≥30 d before stroke + 2a: Any stroke", "Any stroke")
fit_xs(e[o_haemoglobinopathy == 0], "y_any", COV0, "Sensitivity: Excluding haemoglobinopathies: Any stroke", "Any stroke")

# ---- per 1 g/dL (cross-sectional, 2a) and spline (Paper 2 covariates) ----------------------------------
e[, `:=`(hb_below = pmax(13 - hgb, 0), hb_above = pmax(hgb - 13, 0))]
dd <- e[complete.cases(e[, c("y_any", "hgb", COV0, A2, "hormonal_type_cat"), with = FALSE])]; dd[, ht := ht_merge(dd, "y_any")]
f <- glm(reformulate(c("hb_below", "hb_above", COV0, A2, "ht"), "y_any"), binomial, data = dd)
t <- robust_tab(f, vcovHC(f, type = "HC1"), c("hb_below", "hb_above"))
REP[[length(REP) + 1]] <- data.table(analysis = "Per 1 g/dL: Cross-sectional (odds ratio), + pre-Hb conditions (2a)", outcome = "Any stroke",
  level = c("Below 13 g/dL (per 1 g/dL lower)", "Above 13 g/dL (per 1 g/dL higher)"), t[, .(est, lo, hi, p)], N = nrow(dd), events = sum(dd$y_any),
  epv = NA, converged = f$converged)
KN <- c(9.7, 12.6, 13.5, 14.9)                       # knots as reported by the original analysis
rcs4 <- function(x) { b <- Hmisc::rcspline.eval(x, knots = KN, inclx = TRUE); colnames(b) <- paste0("s", 1:ncol(b)); b }
dd <- e[complete.cases(e[, c("y_any", "hgb", COV0, "hormonal_type_cat"), with = FALSE])]; dd[, ht := ht_merge(dd, "y_any")]
S <- rcs4(dd$hgb); dd <- cbind(dd, S)
f <- glm(reformulate(c(colnames(S), COV0, "ht"), "y_any"), binomial, data = dd); V <- vcovHC(f, type = "HC1")
grid <- c(7:12, 14:16); r13 <- rcs4(13)
sp <- rbindlist(lapply(grid, function(h) { dx <- rcs4(h) - r13; b <- sum(coef(f)[colnames(S)] * dx)
  se <- sqrt(drop(dx %*% V[colnames(S), colnames(S)] %*% t(dx))); data.table(level = paste("Hb", h), est = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = NA) }))
REP[[length(REP) + 1]] <- sp[, `:=`(analysis = "Spline (Paper 2 covariates), OR vs 13 g/dL", outcome = "Any stroke", N = nrow(dd), events = sum(dd$y_any), epv = NA, converged = f$converged)]

# ---- original time-to-event: start = later of index and Hb date ----------------------------------------------
t0 <- e[!is.na(grade_cb) & !(stroke_timing %in% c(1, 2)) & !(stroke_any == 1 & !(stroke_timing %in% 3))]
t0[, start := pmax(index_date, hgb_date)]
n_before <- t0[stroke_any == 1 & !is.na(stroke_date) & stroke_date <= start, .N]
t0 <- t0[!(stroke_any == 1 & !is.na(stroke_date) & stroke_date <= start)]
t0[, ev := as.integer(stroke_any == 1)]
t0[, end := fifelse(ev == 1, stroke_date, pmin(last, death, na.rm = TRUE))]
t0 <- t0[!is.na(end) & end > start]
t0[, py := as.numeric(end - start) / 365.25]
fit_p <- function(d, covs, label, outcome_key, y = "ev", vc = "HC1", cluster = NULL) {
  need <- unique(c(y, "grade_cb", covs, "hormonal_type_cat", "py"))
  dd <- d[complete.cases(d[, ..need])]; dd[, ht := ht_merge(dd, y)]
  f <- with_warnings(glm(reformulate(c("grade_cb", covs, "ht", "offset(log(py))"), y), poisson, data = dd), label)
  ft <- glm(reformulate(c("score", covs, "ht", "offset(log(py))"), y), poisson, data = dd)
  vf <- function(x) if (is.null(cluster)) vcovHC(x, type = vc) else vcovCL(x, cluster = dd[[cluster]], type = "HC1")
  out <- rbind(robust_tab(f, vf(f), paste0("grade_cb", names(LV)))[, .(level = names(LV)[match(term, paste0("grade_cb", names(LV)))], est, lo, hi, p)],
               robust_tab(ft, vf(ft), "score")[, .(level = "Per grade (trend)", est, lo, hi, p)])
  out[, `:=`(analysis = label, outcome = outcome_key, N = uniqueN(dd$study_id), events = sum(dd[[y]]),
             epv = sum(dd[[y]]) / (length(coef(f)) - 1), converged = f$converged)]
  REP[[length(REP) + 1]] <<- out
  invisible(list(fit = f, data = dd))
}
r_tte <- fit_p(t0, COV0, "Original time-to-event, Paper 2 covariates", "Any stroke")
fit_p(t0, c(COV0, A2), "Original time-to-event, + 2a", "Any stroke")
res_put("repro_tte_cohort", list(women = nrow(t0), py = sum(t0$py), events = sum(t0$ev), excluded_stroke_before_start = n_before))

# ---- original time-updated analysis (3-y carry-forward, lag 0 / 30; tests in previous year) ------------------
x <- e[!(stroke_timing %in% c(1, 2)) & !(stroke_any == 1 & !(stroke_timing %in% 3))]
x[, ev := as.integer(stroke_any == 1 & stroke_timing %in% 3)]
end_default <- max(c(x$last, x$death), na.rm = TRUE)
x[, end := fifelse(ev == 1, stroke_date, pmin(last, death, na.rm = TRUE))]
x[is.na(end), end := end_default]
hbl <- labs$hb[src == "CBC", .(study_id, date, val)]
split_tu <- function(cap, lag) {
  h <- merge(hbl, x[, .(study_id, index_date, end, ev)], by = "study_id")
  setorder(h, study_id, date)
  h[, `:=`(d_eff = date + lag, nxt = shift(date, -1) + lag), by = study_id]
  h[, start := pmax(d_eff, index_date)]
  h[, stop := pmin(fifelse(is.na(nxt), as.IDate("2200-01-01"), nxt), d_eff + cap, end)]
  iv <- h[stop > start]
  iv[, ev_iv := as.integer(ev == 1 & stop == end)]
  # Hb tests in the year before the interval start (counted on measurement dates)
  cnt <- hbl[iv[, .(study_id, start, lo = start - 365L)], on = .(study_id, date > lo, date <= start), .N, by = .EACHI]
  iv[, tests_prior_y := cnt$N]
  iv[, py := as.numeric(stop - start) / 365.25]
  iv[, grade_cb := factor(fifelse(val >= 12, "None", fifelse(val >= 10, "Mild", fifelse(val >= 8, "Moderate", "Severe"))), levels = GRADES)]
  iv[, score := as.numeric(grade_cb) - 1]
  iv[, tests_cat := cut(tests_prior_y, c(-1, 1, 3, 11, Inf), labels = c("0-1", "2-3", "4-11", "12+"))]
  cols <- c("study_id", setdiff(c(COV0, A2, "hormonal_type_cat"), "study_id"))
  merge(iv[, .(study_id, start, stop, py, ev = ev_iv, grade_cb, score, tests_cat)], x[, ..cols], by = "study_id")
}
tu0 <- split_tu(1095L, 0L)
res_put("repro_tu_cohort", list(women = uniqueN(tu0$study_id), py = sum(tu0$py), events = sum(tu0$ev), intervals = nrow(tu0)))
fit_p(tu0, COV0, "Time-updated: Time-updated, Paper 2 covariates", "Stroke or TIA")
fit_p(tu0, c(COV0, A2), "Time-updated: Time-updated + pre-Hb conditions (2a)", "Stroke or TIA")
fit_p(tu0, c(COV0, A2), "Time-updated + 2a, SE clustered by woman (check)", "Stroke or TIA", cluster = "study_id")
fit_p(tu0, c(COV0, A2, "tests_cat"), "Time-updated: Time-updated + 2a + Hb tests in previous year", "Stroke or TIA")
tu30 <- split_tu(1095L, 30L)
fit_p(tu30, c(COV0, A2), "Time-updated: Time-updated + 2a, 30-day lag", "Stroke or TIA")

# ---- comparison table ---------------------------------------------------------------------------------------
rep <- rbindlist(REP, fill = TRUE)
norm_level <- function(l) { l <- sub(" \\(.*", "", l); sub("Per grade.*", "Per grade (trend)", l) }
rep[, lv := norm_level(level)]
orig[, lv := norm_level(level)]
orig[, okey := fifelse(analysis %like% "^Cross-sectional", paste(analysis, outcome), analysis)]
orig[analysis %like% "^Time-updated", okey := paste(analysis, outcome)]
orig[analysis %like% "^Original time-to-event", okey := paste(analysis, outcome)]
orig[analysis %like% "^Per 1 g/dL", okey := paste(analysis, outcome)]
rep[, okey := fifelse(analysis %like% "^Cross-sectional|^Time-updated|^Original time-to-event|^Per 1 g/dL", paste(analysis, outcome), analysis)]
rep[analysis %like% "^Time-updated: ", okey := paste(analysis, fifelse(outcome == "Stroke or TIA", "Stroke or TIA", outcome))]
cmp <- merge(rep, orig[, .(okey, lv, rep_est = est, rep_lo = lo, rep_hi = hi, rep_N = N, rep_events = events)], by = c("okey", "lv"), all.x = TRUE)
cmp[, `:=`(abs_diff = round(est - rep_est, 2), rel_diff_pct = round(100 * (est / rep_est - 1), 1))]
cmp[, status := fifelse(is.na(rep_est), "no reported value (R check only)",
                fifelse(abs(log(est / rep_est)) <= 0.02 & sign(log(lo)) == sign(log(rep_lo)) & sign(log(hi)) == sign(log(rep_hi)), "reproduced (≤2%)",
                fifelse(abs(log(est / rep_est)) <= 0.10, "close (≤10%)", "NOT reproduced (>10%)")))]
setorder(cmp, analysis, outcome, lv)
out <- cmp[, .(Analysis = analysis, Outcome = outcome, Level = lv,
               `Originally reported` = fci(rep_est, rep_lo, rep_hi), `Reproduced in R` = fci(est, lo, hi),
               `Absolute difference` = abs_diff, `Relative difference %` = rel_diff_pct,
               `N (orig / R)` = fifelse(is.na(rep_N), as.character(N), paste(format(rep_N, big.mark = ","), "/", format(N, big.mark = ","))),
               `Events (orig / R)` = fifelse(is.na(rep_events), as.character(events), paste(rep_events, "/", events)),
               EPV = round(epv, 1), Converged = converged, Status = status)]
save_tab(out, "S2a_reproduction_comparison", "Original (Python) vs independent R reproduction of the original analyses")
saveRDS(list(e = e, A2 = A2, COV0 = COV0, t0 = t0), file.path(DERIVED, "repro_frames.rds"))
res_put("reproduction_status", as.list(table(cmp$status)))
log_msg("Reproduction", paste(names(table(cmp$status)), table(cmp$status), sep = ": ", collapse = "; "))

# ---- diagnostic for the time-updated discrepancy (computed on the raw rows with the original Python Hb filter) ------
h <- fread(dpath("labs", "MDE Workflow Results for Lab Tests (2).csv"), select = c(1, 3, 4, 5), colClasses = list(character = 1))
setnames(h, c("cn", "desc", "val", "date")); d <- toupper(trimws(h$desc))
h <- h[grepl("^(EXTI |EXTM )?(NON-GHS )?(HEMOGLOBIN|HGB|HB)\\b", d) & !grepl("\\bQL\\b|\\bU\\b|URINE|UPHB|O2|OXY|METHEM|A1C|\\bA\\b|\\bA2\\b|\\bF\\b|VARIANT|RETIC|PLASMA|, P$|ARTERIAL|VENOUS|GAS|ABG|VBG|POC|TOTAL|THB|COMMENT|SOURCE|\\bCBC\\b|DONOR", d) & !is.na(val) & val >= 3 & val <= 20]
h[, `:=`(day = as.IDate(substr(date, 1, 10)), gr = as.integer(cut(val, c(-Inf, 8, 10, 12, Inf), right = FALSE)))]
x <- h[, .(n = .N, ng = uniqueN(gr), sev = any(gr == 1)), by = .(cn, day)][n > 1]
res_put("tu_sameday_diag", list(woman_days_multiple_cbc = nrow(x), grade_differs_within_day = x[ng > 1, .N],
                                grade_differs_involving_severe = x[ng > 1 & sev == TRUE, .N]))
deviation("Time-updated reproduction: same-day multiple Hb results",
          sprintf("R uses the daily median (Codebook rule); the original kept every row, so the value in force differs on %s woman-days whose same-day results span grades (%s involve Hb <8). This explains the severe-grade discrepancy (6 events).",
                  format(x[ng > 1, .N], big.mark = ","), format(x[ng > 1 & sev == TRUE, .N], big.mark = ",")))
rm(h, x); gc()
