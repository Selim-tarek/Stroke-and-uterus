# Phase 3: missingness and selection into the analyzed (Hb-tested) cohort.
source("R/00_setup.R")
e <- readRDS(file.path(DERIVED, "eligible.rds"))
ic <- readRDS(file.path(DERIVED, "incident_cohort.rds"))

# women free of stroke at index, with or without a baseline Hb; follow-up computed the same way for both
c1 <- e[prevalent == 0 & timing_unknown == 0 & !is.na(fu_end0) & fu_end0 > index_date]
c1[, py := as.numeric(fu_end0 - index_date) / 365.25]
c1[, ev_stroke := as.integer(stroke_any == 1 & stroke_timing %in% 3 & stroke_type != 4 & !is.na(ev_date) & ev_date <= fu_end0)]
c1[, ev_any := as.integer(stroke_any == 1 & stroke_timing %in% 3 & !is.na(ev_date) & ev_date <= fu_end0)]
c1[, died := as.integer(ev_any == 0 & !is.na(death) & death <= fu_end0)]
c1[, tested := factor(fifelse(is.na(hb), "No baseline Hb", "Baseline Hb"), levels = c("Baseline Hb", "No baseline Hb"))]
vars <- c("age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "hormonal_type_cat",
          "dxgrp", "uterine_bleeding", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "gibleed_1y", "preg_1y",
          "anaemia_cause_any", "hb_tests_2y", "py", "ev_stroke", "ev_any", "died")
t1 <- describe(c1, "tested", vars)
rates <- c1[, .(women = .N, strokes = sum(ev_stroke), py = sum(py), rate = 1000 * sum(ev_stroke) / sum(py),
                lo = 1000 * qchisq(.025, 2 * sum(ev_stroke)) / 2 / sum(py), hi = 1000 * qchisq(.975, 2 * (sum(ev_stroke) + 1)) / 2 / sum(py)), by = tested]
save_tab(t1, "S1a_tested_vs_untested", "Women free of stroke at index: with vs without a baseline laboratory Hb")
save_tab(rates[, .(tested, women, strokes, person_years = round(py), rate_per_1000py = f2(rate), ci = sprintf("%s–%s", f2(lo), f2(hi)))],
         "S1b_rates_tested_vs_untested", "Crude incident stroke rate by baseline-Hb availability")
# adjusted rate ratio for being untested (age, race, htn, dm) - describes selection, not an exposure effect
fit <- glm(ev_stroke ~ tested + ns(age_index, 3) + race4 + htn + dm + dyslipidemia, family = poisson, offset = log(py), data = c1)
sel <- robust_tab(fit, vcovHC(fit, type = "HC0"), "testedNo baseline Hb")
res_put("selection", list(rates = rates, adj_rr_untested = sel[, .(est, lo, hi, p)]))
log_msg("Selection", sprintf("stroke rate per 1,000 PY: tested %s, untested %s; adjusted RR untested vs tested %s",
                             f2(rates[tested == "Baseline Hb", rate]), f2(rates[tested == "No baseline Hb", rate]), fci(sel$est, sel$lo, sel$hi)))

# missingness of every model variable in the incident cohort
mv <- c("hb", "age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "hormonal_type_cat",
        "dxgrp", "uterine_bleeding", "thrombophilia", "afib", "chf", "vte_history", "mcv_lab", "rdw_lab", "plt_lab", "ferritin")
miss <- rbindlist(lapply(mv, function(v) data.table(variable = v, label = if (!is.na(VAR_LABELS[v])) VAR_LABELS[v] else v,
  missing = sum(is.na(ic[[v]])), pct = round(100 * mean(is.na(ic[[v]])), 2),
  note = if (v == "smoking3") sprintf("'Unknown' (code 9) is a category, not blank: %s (%.1f%%)", format(sum(ic$smoking3 == "Unknown"), big.mark = ","),
                                       100 * mean(ic$smoking3 == "Unknown")) else "")))
save_tab(miss, "S1c_missingness_incident_cohort", "Missing values in the incident cohort")
pat <- ic[, .(bmi_missing = is.na(bmi), hormonal_missing = is.na(hormonal_type_cat), smoking_unknown = smoking3 == "Unknown")][, .N, by = .(bmi_missing, hormonal_missing, smoking_unknown)][order(-N)]
save_tab(pat, "S1d_missing_patterns", "Joint missingness pattern of covariates with missing/unknown values (incident cohort)")
res_put("missing_incident", setNames(as.list(miss$missing), miss$variable))
