# Phase 7: pre-specified sensitivity analyses S2-S11 (S1 complete case and S12 E-values are in 07_primary.R).
source("R/00_setup.R"); source("R/model_fns.R")
e <- readRDS(file.path(DERIVED, "eligible.rds"))
ic <- prep_covs(readRDS(file.path(DERIVED, "incident_cohort.rds"))[!is.na(hormonal_type_cat)], "ev_stroke")
imps <- readRDS(file.path(DERIVED, "imputed_primary.rds"))
deviation("Sensitivity analyses S2-S5, S10, S11", "Complete case for BMI (3.5% missing) instead of MI, to limit computation; S1 shows MI and complete-case estimates agree. S6-S9 reuse the primary imputations.")

# cohort builder for an alternative baseline-Hb definition (same rules as 04_cohort.R)
make_cohort <- function(hbvar, start_offset = 0L) {
  d <- e[prevalent == 0 & timing_unknown == 0 & !is.na(get(hbvar)) & !is.na(hormonal_type_cat)]
  d[, t0 := index_date + start_offset]
  d <- d[!is.na(fu_end0) & fu_end0 > t0]
  d <- d[!(stroke_any == 1 & stroke_timing %in% 3 & !is.na(ev_date) & ev_date <= t0)]
  d[, py := as.numeric(fu_end0 - t0) / 365.25]
  d[, ev_stroke := as.integer(stroke_any == 1 & stroke_timing %in% 3 & stroke_type != 4 & !is.na(ev_date) & ev_date <= fu_end0)]
  d[, grade := factor(fifelse(get(hbvar) >= 12, "None", fifelse(get(hbvar) >= 10, "Mild", fifelse(get(hbvar) >= 8, "Moderate", "Severe"))), levels = GRADES)]
  d[, hb := get(hbvar)]
  prep_covs(d, "ev_stroke")[!is.na(bmi)]
}
S <- list()
run <- function(id, label, bias, dl, covs = PRIMARY_COVS, y = "ev_stroke") {
  r <- grade_models(dl, y, covs, label, cox = FALSE)
  S[[id]] <<- r[, .(id = id, analysis = label, bias_addressed = bias, level, est, lo, hi, p, N, events, epv)]
}
cc <- ic[!is.na(bmi)]
run("S2", "Follow-up starts at day 31 (early events excluded)", "Reverse causation / detection around the index date", list(make_cohort("hb", 30L)))
run("S3", "Baseline Hb within 365 d before index", "Stale exposure (Hb up to 3 y old)", list(make_cohort("hb_1y")))
run("S4", "POC Hb included in baseline Hb", "Exposure definition (point-of-care values)", list(make_cohort("hb_poc")))
run("S5", "Events limited to women with ≥1 acute stroke code", "Outcome misdating (history/sequela codes)", list(cc), y = "ev_stroke_acute")
run("S6", "+ heavy/abnormal uterine bleeding", "Confounding vs adjustment for an exposure determinant", imps, c(PRIMARY_COVS, "uterine_bleeding"))
run("S7", "+ thrombophilia, AF, heart failure, VTE (undated)", "Over-adjustment check (factors possibly recorded after stroke)", imps,
    c(PRIMARY_COVS, "thrombophilia", "afib", "chf", "vte_history"))
imps8 <- lapply(imps, function(z) { z <- copy(z); z[, tests_cat := cut(hb_tests_2y, c(-1, 0, 2, 5, Inf), labels = c("0", "1-2", "3-5", "6+"))]; z })
run("S8", "+ number of Hb tests in the 2 y before index", "Healthcare contact / surveillance", imps8, c(PRIMARY_COVS, "tests_cat"))
run("S9", "Excluding women with any pre-index anemia-causing condition", "Residual confounding by illness causing anemia",
    lapply(imps, function(z) z[anaemia_cause_any == 0]),
    setdiff(PRIMARY_COVS, c("haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y")))

# S10 Fine-Gray (death competing), complete case
cc[, status := factor(fifelse(ev_stroke == 1, "stroke", fifelse(died == 1, "death", "censor")), levels = c("censor", "stroke", "death"))]
fg_dat <- finegray(as.formula(paste("Surv(py, status) ~ grade +", paste(c("age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad",
  "migraine_any", "ht", "dxgrp", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y"), collapse = " + "))),
  data = cc, etype = "stroke")
fg <- with_warnings(coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~ grade +", paste(PRIMARY_COVS, collapse = " + "))), weights = fgwt, data = fg_dat), "S10 Fine-Gray")
cc[, score := as.numeric(grade) - 1]
fg_dat2 <- finegray(as.formula(paste("Surv(py, status) ~ score +", paste(c("age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad",
  "migraine_any", "ht", "dxgrp", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y"), collapse = " + "))),
  data = cc, etype = "stroke")
fg2 <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~ score +", paste(PRIMARY_COVS, collapse = " + "))), weights = fgwt, data = fg_dat2)
t <- robust_tab(fg, vcov(fg), GTERMS); t2 <- robust_tab(fg2, vcov(fg2), "score")
S[["S10"]] <- rbind(t[, .(level = sub("grade", "", term), est, lo, hi, p)], t2[, .(level = "Per grade (trend)", est, lo, hi, p)])[
  , `:=`(id = "S10", analysis = "Fine-Gray subdistribution HR (death competing)", bias_addressed = "Competing risk of death", N = nrow(cc), events = sum(cc$ev_stroke), epv = NA)]

# S11 inverse-probability weights for having a baseline Hb (selection into the analyzed cohort)
c1 <- e[prevalent == 0 & timing_unknown == 0 & !is.na(fu_end0) & fu_end0 > index_date & !is.na(hormonal_type_cat)]
c1[, `:=`(tested = as.integer(!is.na(hb)), bmi_miss = as.integer(is.na(bmi)), bmi0 = fifelse(is.na(bmi), median(bmi, na.rm = TRUE), bmi))]
ps <- glm(tested ~ ns(age_index, 3) + race4 + bmi0 + bmi_miss + smoking3 + htn + dm + dyslipidemia + cad + migraine_any + hormonal_type_cat + dxgrp +
            uterine_bleeding + haemoglobinopathy + ckd_pre + liver_pre + alcohol_pre + ibd_pre + malabs_pre + gibleed_1y + preg_1y, binomial, data = c1)
c1[, pr := fitted(ps)]
c1[, sw := fifelse(tested == 1, mean(tested) / pr, NA_real_)]
q <- quantile(c1$sw, c(.01, .99), na.rm = TRUE); c1[, sw := pmin(pmax(sw, q[1]), q[2])]
ccw <- merge(cc, c1[, .(study_id, sw)], by = "study_id")
f <- with_warnings(glm(reformulate(c("grade", PRIMARY_COVS, "offset(log(py))"), "ev_stroke"), quasipoisson, data = ccw, weights = sw), "S11 IPW")
ccw[, score := as.numeric(grade) - 1]
ft <- glm(reformulate(c("score", PRIMARY_COVS, "offset(log(py))"), "ev_stroke"), quasipoisson, data = ccw, weights = sw)
t <- robust_tab(f, vcovHC(f, type = "HC0"), GTERMS); t2 <- robust_tab(ft, vcovHC(ft, type = "HC0"), "score")
S[["S11"]] <- rbind(t[, .(level = sub("grade", "", term), est, lo, hi, p)], t2[, .(level = "Per grade (trend)", est, lo, hi, p)])[
  , `:=`(id = "S11", analysis = "Inverse-probability weighted for having a baseline Hb", bias_addressed = "Selection into the tested cohort",
         N = nrow(ccw), events = sum(ccw$ev_stroke), epv = NA)]
res_put("ipw_weights", list(min = min(ccw$sw), max = max(ccw$sw), mean = mean(ccw$sw), ps_c = NA))

sens <- rbindlist(S, fill = TRUE)
prim <- fread(file.path(OUT, "tables", "Table2_primary_models.csv"))[analysis %like% "^Primary" & measure %like% "Rate"]
wide <- dcast(sens[level != "Any anemia (Hb <12)", .(id, analysis, bias_addressed, N, events, level, txt = fci(est, lo, hi))],
              id + analysis + bias_addressed + N + events ~ level, value.var = "txt")
prim_row <- data.table(id = "P", analysis = "Primary (MI, full adjustment)", bias_addressed = "—", N = prim$N[1], events = prim$events[1],
                       Mild = prim[level == "Mild", estimate], Moderate = prim[level == "Moderate", estimate], Severe = prim[level == "Severe", estimate],
                       `Per grade (trend)` = prim[level == "Per grade (trend)", estimate])
wide <- wide[order(as.integer(sub("S", "", id)))]
wide <- rbind(prim_row, wide, fill = TRUE)
setcolorder(wide, c("id", "analysis", "bias_addressed", "N", "events", "Mild", "Moderate", "Severe", "Per grade (trend)"))
wide[, `Any anemia` := c(prim[level == "Any anemia (Hb <12)", estimate], sens[level == "Any anemia (Hb <12)"][match(wide$id[-1], id), fci(est, lo, hi)])]
wide[is.na(`Any anemia`) | `Any anemia` == "", `Any anemia` := "not estimated"]
save_tab(wide, "Table5_sensitivity", "Sensitivity analyses for the primary estimand (rate ratios unless stated; S10 = subdistribution hazard ratios)")
res_put("sensitivity", sens)
