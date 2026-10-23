# PI design (2026-10-23): anchor each woman at the EARLIEST evidence of uterine pathology, follow her for stroke after
# that date, and for women who had a stroke describe their anemia before the stroke (nested case-control).
# Earliest evidence = min(first coded diagnosis from 2018 [index_date; IRB window], dated myomectomy, uterine artery
# embolisation, hysterectomy, endometrial ablation). 'Other gynaecological surgery' (tubal, ovarian, D&C, IUD removal)
# is not evidence of fibroids/adenomyosis/endometriosis and is not used.
source("R/00_setup.R"); source("R/model_fns.R")
e  <- readRDS(file.path(DERIVED, "eligible.rds"))
m  <- readRDS(file.path(DERIVED, "master_derived.rds"))[eligible == 1, .(study_id, notes, surgical_tx, surgical_tx_date, myomectomy_date)]
hb <- readRDS(file.path(DERIVED, "labs_long.rds"))$hb[src == "CBC", .(study_id, date, val)]
cm <- readRDS(file.path(DERIVED, "comorbid_long.rds"))
e  <- merge(e[, !c("notes", "surgical_tx", "surgical_tx_date", "myomectomy_date")], m, by = "study_id")

# ---- 1. anchor -----------------------------------------------------------------------------------------------
SPECIFIC <- c("myomectomy", "uterine artery embolisation", "hysterectomy", "endometrial ablation")
st_lab <- c("1" = "myomectomy", "2" = "hysterectomy", "3" = "uterine artery embolisation", "4" = "endometrial ablation")
ev <- rbind(
  e[!is.na(surgical_tx_date) & surgical_tx %in% 1:4, .(study_id, src = st_lab[as.character(surgical_tx)], date = surgical_tx_date)],
  e[!is.na(myomectomy_date), .(study_id, src = "myomectomy", date = myomectomy_date)],
  e[!is.na(notes), {
    mt <- regmatches(notes, gregexpr("(myomectomy|uterine artery embolisation|hysterectomy|endometrial ablation) \\d{4}-\\d{2}-\\d{2}", notes))[[1]]
    if (length(mt)) .(src = sub(" \\d{4}.*", "", mt), date = as.IDate(sub(".* ", "", mt))) else .(src = character(), date = as.IDate(character()))
  }, by = study_id])
ev <- ev[!is.na(date) & year(date) >= 1980]
first_ev <- ev[order(date), .SD[1], by = study_id]
e <- merge(e, first_ev[, .(study_id, proc_src = src, proc_date = date)], by = "study_id", all.x = TRUE)
e[, anchor := fifelse(!is.na(proc_date) & proc_date < index_date &
                       age_index - as.numeric(index_date - proc_date) / 365.25 >= 18, proc_date, index_date)]   # implausible procedure dates (age <18) ignored
e[, anchor_src := fifelse(anchor < index_date, proc_src, "first coded diagnosis (2018+)")]
e[, age_anchor := age_index - as.numeric(index_date - anchor) / 365.25]
shift <- e[anchor < index_date]
anc <- rbind(data.table(item = "Eligible women", value = nb <- format(nrow(e), big.mark = ",")),
  data.table(item = "Anchor = first coded diagnosis (2018 onward)", value = format(e[anchor == index_date, .N], big.mark = ",")),
  data.table(item = "Anchor moved earlier by a dated procedure", value = format(nrow(shift), big.mark = ",")),
  shift[, .(item = paste0("  of which first evidence = ", anchor_src), value = format(.N, big.mark = ",")), by = anchor_src][, .(item, value)],
  data.table(item = "  years moved earlier, median (IQR)", value = sprintf("%.1f (%.1f–%.1f)", median(as.numeric(shift$index_date - shift$anchor)) / 365.25,
    quantile(as.numeric(shift$index_date - shift$anchor), .25) / 365.25, quantile(as.numeric(shift$index_date - shift$anchor), .75) / 365.25)),
  data.table(item = "  anchor before 2018", value = format(shift[anchor < as.IDate("2018-01-01"), .N], big.mark = ",")))
save_tab(anc, "A1_anchor_definition", "Earliest evidence of uterine pathology (anchor)")

# ---- 2. stroke relative to the anchor --------------------------------------------------------------------------
s <- e[stroke_any == 1]
s[, rel := fifelse(is.na(stroke_date), "Undated", fifelse(stroke_date <= anchor, "On or before uterine diagnosis", "After uterine diagnosis"))]
TYPE <- c("1" = "Ischemic stroke", "2" = "Intracerebral hemorrhage", "3" = "Subarachnoid hemorrhage", "4" = "TIA", "5" = "Cerebral venous thrombosis", "9" = "Stroke, unknown type")
s[, type := TYPE[as.character(stroke_type)]]
rt <- dcast(s[, .N, by = .(type, rel)], type ~ rel, value.var = "N", fill = 0)
rt <- rbind(rt, cbind(data.table(type = "All"), rt[, lapply(.SD, sum), .SDcols = -1]))
save_tab(rt, "A2_stroke_vs_anchor", "Stroke/TIA by timing relative to the earliest evidence of uterine pathology")

# ---- 3. anchored cohort: women stroke-free at the anchor, followed for stroke after it ---------------------------
co <- e[!(stroke_any == 1 & (is.na(stroke_date) | stroke_date <= anchor)) & !is.na(hormonal_type_cat)]
co[, t_end := pmin(last, death, ADMIN_END, na.rm = TRUE)]
co[stroke_any == 1, t_end := pmin(t_end, stroke_date)]
co <- co[!is.na(t_end) & t_end > anchor]
co[, `:=`(py = as.numeric(t_end - anchor) / 365.25,
          ev_stroke = as.integer(stroke_any == 1 & stroke_type != 4 & stroke_date <= t_end),
          ev_any = as.integer(stroke_any == 1 & stroke_date <= t_end),
          ev_isch = as.integer(stroke_any == 1 & stroke_type == 1 & stroke_date <= t_end),
          ev_tia = as.integer(stroke_any == 1 & stroke_type == 4 & stroke_date <= t_end))]
co[, died := as.integer(ev_any == 0 & !is.na(death) & death <= t_end)]
# covariates re-dated to the anchor
k <- merge(cm, co[, .(study_id, anchor)], by = "study_id")
has <- function(fl, recent = FALSE, any_date = FALSE) {
  x <- k[flag %in% fl]; if (!any_date) x <- x[date <= anchor]; if (recent) x <- x[date >= anchor - 365]
  as.integer(co$study_id %in% x$study_id) }
co[, `:=`(haemoglobinopathy = has(c("sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy"), any_date = TRUE),
          ckd_pre = has(c("ckd", "ckd_esrd")), liver_pre = has(c("cirrhosis", "chronic_liver")), alcohol_pre = has("alcohol_disorder"),
          ibd_pre = has("ibd"), hiv_pre = has("hiv"), malabs_pre = has("malabsorption"),
          gibleed_1y = has("gi_bleed", recent = TRUE), preg_1y = has("pregnancy", recent = TRUE), age_index = age_anchor)]
# anemia at the anchor: most recent CBC Hb in the 3 years up to and including the anchor
b <- merge(hb, co[, .(study_id, anchor)], by = "study_id")[date <= anchor & date >= anchor - 1095]
setorder(b, study_id, -date); b <- b[, .SD[1], by = study_id][, .(study_id, hb_a = val)]
co <- merge(co[, !c("hb", "grade")], b, by = "study_id", all.x = TRUE)
co[, grade := factor(fifelse(hb_a >= 12, "None", fifelse(hb_a >= 10, "Mild", fifelse(hb_a >= 8, "Moderate", "Severe"))), levels = GRADES)]
co[, hb := hb_a]
saveRDS(co, file.path(DERIVED, "anchored_cohort.rds"))
cw <- prep_covs(co[!is.na(grade)], "ev_stroke")
flow <- data.table(step = c("Eligible women", "Stroke/TIA on or before the anchor, or undated (excluded)", "No follow-up after the anchor (excluded)",
                            "Anchored cohort (stroke-free at earliest uterine evidence)", "  with a CBC Hb in the 3 y up to the anchor (analysis A)",
                            "  incident strokes (excl. TIA) after the anchor", "  incident stroke or TIA after the anchor"),
                   n = c(nrow(e), e[stroke_any == 1 & (is.na(stroke_date) | stroke_date <= anchor), .N],
                         e[!(stroke_any == 1 & (is.na(stroke_date) | stroke_date <= anchor)) & !is.na(hormonal_type_cat), .N] - nrow(co),
                         nrow(co), nrow(cw), sum(co$ev_stroke), sum(co$ev_any)))
save_tab(flow, "A3_anchored_flow", "Anchored cohort")
imps <- impute_bmi(cw, "ev_stroke")
ra <- rbind(grade_models(imps, "ev_stroke", PRIMARY_COVS, "A: anchored cohort, stroke (excl. TIA)", cox = FALSE),
            grade_models(imps, "ev_isch", REDUCED_COVS, "A: anchored cohort, ischemic stroke (reduced set)", cox = FALSE),
            grade_models(imps, "ev_any", PRIMARY_COVS, "A: anchored cohort, stroke or TIA", cox = FALSE))
rates <- cw[, .(women = .N, strokes = sum(ev_stroke), py = sum(py)), by = grade][order(grade)]
rates[, rate := 1000 * strokes / py]
save_tab(ra[, .(analysis, level, estimate = fci(est, lo, hi), p = fp(p), N, events, epv = round(epv, 1))], "A4_anchored_models",
         "Anemia at the earliest uterine evidence and stroke afterwards (Poisson, robust SE, MI for BMI)")
save_tab(rates[, .(grade = GRADE_LABELS[as.character(grade)], women, strokes, person_years = round(py), rate_per_1000py = f2(rate))], "A4a_anchored_rates")

# ---- 4. nested case-control: for each woman with a stroke after the anchor, what was her anemia before the stroke? ----------
# time scale = time since the anchor; for each case, up to 10 controls still under follow-up and stroke-free at the same
# time since anchor, with age at that date within 2 years. Exposure = last CBC Hb in the year before the reference date,
# excluding the last 30 days (Hb drawn during the stroke presentation cannot define exposure); no-lag version as sensitivity.
set.seed(SEED)
co[, t_case := fifelse(ev_any == 1, as.integer(stroke_date - anchor), NA_integer_)]
co[, fu_days := as.integer(t_end - anchor)]
cases <- co[ev_any == 1 & t_case > 30]
pool <- co[, .(study_id, anchor, fu_days, t_case, age_anchor)]
sets <- rbindlist(lapply(seq_len(nrow(cases)), function(i) {
  cs <- cases[i]; t <- cs$t_case; age_ref <- cs$age_anchor + t / 365.25
  r <- pool[study_id != cs$study_id & fu_days > t & (is.na(t_case) | t_case > t) & abs(age_anchor + t / 365.25 - age_ref) <= 2]
  ctl <- if (nrow(r) > 10) r[sample(.N, 10)] else r
  rbind(data.table(set = i, study_id = cs$study_id, case = 1L, ref = cs$anchor + t),
        data.table(set = rep(i, nrow(ctl)), study_id = ctl$study_id, case = rep(0L, nrow(ctl)), ref = ctl$anchor + t))
}))
expo <- function(lag) {
  x <- merge(sets[, .(set, study_id, ref)], hb, by = "study_id", allow.cartesian = TRUE)[date <= ref - lag & date >= ref - 365]
  setorder(x, set, study_id, -date)
  last <- x[, .(hb_pre = val[1], n_tests = .N), by = .(set, study_id)]
  d <- merge(sets, last, by = c("set", "study_id"), all.x = TRUE)
  d[is.na(n_tests), n_tests := 0L]
  d[, ex := factor(fifelse(is.na(hb_pre), "No Hb in window", as.character(fifelse(hb_pre >= 12, "None", fifelse(hb_pre >= 10, "Mild", fifelse(hb_pre >= 8, "Moderate", "Severe"))))),
                   levels = c(GRADES, "No Hb in window"))]
  d[, anemic := fifelse(is.na(hb_pre), NA_integer_, as.integer(hb_pre < 12))]
  d
}
cov_ncc <- c("race4", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "dxgrp", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre")
ncc <- function(lag, label) {
  d <- merge(expo(lag), co[, c("study_id", "stroke_type", cov_ncc), with = FALSE], by = "study_id")
  d[, case_type := fifelse(case == 1, fifelse(stroke_type == 4, "TIA", "Stroke"), "Control")]
  desc <- d[, .(n = .N, with_hb = sum(!is.na(hb_pre)), anemic = sum(anemic, na.rm = TRUE), mild = sum(ex == "Mild"), moderate = sum(ex == "Moderate"),
                severe = sum(ex == "Severe"), median_hb = as.double(median(hb_pre, na.rm = TRUE)), median_tests = as.double(median(n_tests))), by = case_type]
  desc[, pct_anemic_of_tested := round(100 * anemic / with_hb, 1)]
  fits <- list()
  for (sub in c("Stroke", "TIA")) {
    dd <- d[set %in% d[case == 1 & case_type == sub, set]]
    for (spec in list(list("crude (matched on age and time since diagnosis)", character()), list("adjusted", cov_ncc),
                      list("adjusted + number of Hb tests in the window", c(cov_ncc, "n_tests")))) {
      f <- with_warnings(clogit(as.formula(paste("case ~ ex +", paste(c("strata(set)", spec[[2]]), collapse = " + "))), data = dd, method = "exact"),
                         paste("NCC", label, sub, spec[[1]]))
      ct <- summary(f)$conf.int; pv <- summary(f)$coefficients[, "Pr(>|z|)"]
      rn <- grep("^ex", rownames(ct))
      fits[[length(fits) + 1]] <- data.table(window = label, cases = sub, model = spec[[1]], level = sub("^ex", "", rownames(ct)[rn]),
                                             or = fci(ct[rn, 1], ct[rn, 3], ct[rn, 4]), p = fp(pv[rn]),
                                             n_case_sets = uniqueN(dd$set), n_cases = sum(dd$case))
    }
  }
  list(desc = desc[, window := label], fits = rbindlist(fits))
}
r1 <- ncc(30L, "Hb 31-365 days before (30-day lag)")
r0 <- ncc(0L, "Hb 1-365 days before (no lag)")
save_tab(rbind(r1$desc, r0$desc), "A5_ncc_anemia_before_stroke",
         "Women with stroke/TIA after uterine diagnosis vs matched controls: Hb in the year before the stroke date")
save_tab(rbind(r1$fits, r0$fits), "A6_ncc_models", "Conditional logistic regression (matched sets): anemia grade before the reference date vs Hb ≥12")
saveRDS(sets, file.path(DERIVED, "ncc_sets.rds"))
res_put("anchor_design", list(flow = flow, models = ra, ncc_desc = rbind(r1$desc, r0$desc), ncc = rbind(r1$fits, r0$fits)))
