# Phase 5: revised primary analysis - baseline anemia (Hb on/before index) and incident stroke after index.
source("R/00_setup.R"); source("R/model_fns.R")
ic0 <- readRDS(file.path(DERIVED, "incident_cohort.rds"))
n_ht_missing <- ic0[is.na(hormonal_type_cat), .N]
ic <- ic0[!is.na(hormonal_type_cat)]
if (n_ht_missing) deviation("Primary cohort", sprintf("%d woman with hormonal_type = 9 (missing) excluded from models (complete case for this variable).", n_ht_missing))
deviation("Sparse binary covariates", "HIV before index had 0 incident strokes and IBD 4; they are pooled into one 'IBD or HIV' indicator (analogous to the pre-specified hormonal-type rule) instead of being entered separately.")
ic <- prep_covs(ic, "ev_stroke")
if (length(attr(ic, "ht_merged"))) log_msg("Primary", paste("hormonal levels merged into Other/multiple:", paste(attr(ic, "ht_merged"), collapse = ", ")))

# ---- descriptive: Table 1 by baseline grade --------------------------------------------------------------
t1vars <- c("age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "hormonal_type_cat", "dxgrp",
            "uterine_bleeding", "thrombophilia", "afib", "chf", "vte_history", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre",
            "ibd_pre", "hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y", "anaemia_cause_any", "hb", "mcv_lab", "rdw_lab", "plt_lab",
            "hb_tests_2y", "py", "ev_stroke", "ev_any", "died")
t1 <- describe(ic0, "grade", t1vars)
setnames(t1, GRADES, GRADE_LABELS[GRADES])
save_tab(t1, "Table1_baseline_by_grade", "Baseline characteristics of the incident-stroke cohort by baseline anemia grade")

# ---- rates and cumulative incidence (death as competing event) ---------------------------------------------
ic[, status := factor(fifelse(ev_stroke == 1, "stroke", fifelse(died == 1, "death", "censor")), levels = c("censor", "stroke", "death"))]
aj <- survfit(Surv(py, status) ~ grade, data = ic)
ci5 <- summary(aj, times = 5, extend = TRUE)
cif <- data.table(grade = sub("grade=", "", ci5$strata), cif5 = ci5$pstate[, which(aj$states == "stroke")],
                  lo = ci5$lower[, which(aj$states == "stroke")], hi = ci5$upper[, which(aj$states == "stroke")])
rates <- ic[, .(women = .N, strokes = sum(ev_stroke), py = sum(py), deaths = sum(died)), by = grade][order(grade)]
rates[, `:=`(rate = 1000 * strokes / py, rate_lo = 1000 * qchisq(.025, 2 * strokes) / 2 / py,
             rate_hi = 1000 * qchisq(.975, 2 * (strokes + 1)) / 2 / py, death_rate = 1000 * deaths / py)]
rates <- merge(rates, cif, by = "grade")
rates[, grade := factor(grade, levels = GRADES)]; setorder(rates, grade)
saveRDS(aj, file.path(DERIVED, "aj_primary.rds"))

# ---- models: MI for BMI (primary), complete case (S1) -------------------------------------------------------
imps <- impute_bmi(ic, "ev_stroke", m = 20)
saveRDS(imps, file.path(DERIVED, "imputed_primary.rds"))
prim <- grade_models(imps, "ev_stroke", PRIMARY_COVS, "Primary: incident stroke, full adjustment (MI)")
crude <- grade_models(list(ic), "ev_stroke", character(), "Crude", cox = FALSE)
agerace <- grade_models(list(ic), "ev_stroke", c("ns(age_index, 3)", "race4"), "Age and race adjusted", cox = FALSE)
cc <- grade_models(list(ic[!is.na(bmi)]), "ev_stroke", PRIMARY_COVS, "S1 Complete case (BMI)", cox = FALSE)
red <- grade_models(imps, "ev_stroke", REDUCED_COVS, "Reduced pre-declared set (MI)", cox = FALSE)
allm <- rbind(prim, crude, agerace, cc, red)

# proportional hazards check (first imputation)
cz <- cox.zph(coxph(as.formula(paste("Surv(py, ev_stroke) ~ grade +", paste(PRIMARY_COVS, collapse = " + "))), data = imps[[1]]))
ph <- data.table(term = rownames(cz$table), chisq = cz$table[, 1], df = cz$table[, 2], p = cz$table[, 3])
save_tab(ph, "S2b_ph_test_primary", "Schoenfeld residual test, Cox model for the primary outcome (imputation 1)")

# ---- continuous Hb: piecewise per g/dL and restricted cubic spline ------------------------------------------------
add_pw <- function(z) { z <- copy(z); z[, `:=`(hb_below13 = pmax(13 - hb, 0), hb_above13 = pmax(hb - 13, 0))]; z }
pw <- pooled(pois_fit(lapply(imps, add_pw), "ev_stroke", c("hb_below13", "hb_above13"), PRIMARY_COVS, "Piecewise Hb"), c("hb_below13", "hb_above13"))
KN <- as.numeric(quantile(ic$hb, c(.05, .35, .65, .95)))
rcs_basis <- function(x) { b <- Hmisc::rcspline.eval(x, knots = KN, inclx = TRUE); colnames(b) <- paste0("hbs", seq_len(ncol(b))); b }
add_rcs <- function(z) { z <- copy(z); z <- cbind(z, rcs_basis(z$hb)); z }
spf <- pois_fit(lapply(imps, add_rcs), "ev_stroke", paste0("hbs", 1:3), PRIMARY_COVS, "Spline Hb")
pf <- pool_full(spf$fits, vc_rob)
sn <- paste0("hbs", 1:3)
rng <- as.numeric(quantile(ic$hb, c(.01, .99)))
grid <- seq(floor(rng[1] * 10) / 10, ceiling(rng[2] * 10) / 10, by = 0.1)
r13 <- rcs_basis(13)
curve <- rbindlist(lapply(grid, function(h) { dx <- rcs_basis(h) - r13; b <- sum(pf$b[sn] * dx)
  se <- sqrt(drop(dx %*% pf$V[sn, sn] %*% t(dx))); data.table(hb = h, rr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se)) }))
wald <- function(idx) { b <- pf$b[idx]; V <- pf$V[idx, idx, drop = FALSE]; s <- drop(t(b) %*% solve(V) %*% b); c(chisq = s, df = length(idx), p = pchisq(s, length(idx), lower.tail = FALSE)) }
w_all <- wald(sn); w_nl <- wald(sn[-1])
saveRDS(list(curve = curve, knots = KN, hb = ic$hb, range = rng), file.path(DERIVED, "spline_primary.rds"))
fwrite(curve, file.path(OUT, "tables", "Fig2_spline_curve_data.csv"))
pts <- curve[round(hb, 1) %in% c(7:12, 14:16)]

# ---- E-values ---------------------------------------------------------------------------------------------
ev_mod <- evalue(prim[measure %like% "Rate" & level == "Moderate", est], prim[measure %like% "Rate" & level == "Moderate", lo], prim[measure %like% "Rate" & level == "Moderate", hi])
ev_tr <- evalue(prim[measure %like% "Rate" & level == "Per grade (trend)", est], prim[measure %like% "Rate" & level == "Per grade (trend)", lo], prim[measure %like% "Rate" & level == "Per grade (trend)", hi])

# ---- outputs ---------------------------------------------------------------------------------------------------
save_tab(allm[, .(analysis, outcome, measure, level, estimate = fci(est, lo, hi), p = fp(p), N, events, df, epv = round(epv, 1), converged, imputations)],
         "Table2_primary_models", "Revised primary analysis: baseline anemia grade and incident stroke (excluding TIA)")
save_tab(rates[, .(grade = GRADE_LABELS[as.character(grade)], women, strokes, person_years = round(py), rate_per_1000py = fci(rate, rate_lo, rate_hi),
                   deaths_without_stroke = deaths, death_rate_per_1000py = f2(death_rate),
                   cumulative_incidence_5y_pct = fci(100 * cif5, 100 * lo, 100 * hi, 1))],
         "Table2a_rates_by_grade", "Incident stroke rates and 5-year cumulative incidence (Aalen-Johansen, death as competing event)")
save_tab(rbind(pw[, .(analysis = "Per 1 g/dL", term = c("Per 1 g/dL lower, below 13", "Per 1 g/dL higher, above 13"), estimate = fci(est, lo, hi), p = fp(p))],
               pts[, .(analysis = "Spline, RR vs 13 g/dL", term = paste0("Hb ", hb), estimate = fci(rr, lo, hi), p = "")],
               data.table(analysis = "Spline tests", term = c("overall association", "non-linearity"),
                          estimate = c(sprintf("chi2 %.2f, df %d", w_all[1], w_all[2]), sprintf("chi2 %.2f, df %d", w_nl[1], w_nl[2])),
                          p = fp(c(w_all[3], w_nl[3])))),
         "Table2b_continuous_hb", sprintf("Continuous Hb, incident stroke, full adjustment (MI). Spline knots %s g/dL; curve restricted to 1st-99th percentile (%s-%s).",
                                          paste(f2(KN, 1), collapse = ", "), f2(rng[1], 1), f2(rng[2], 1)))
res_put("primary", list(N = nrow(ic), events = sum(ic$ev_stroke), py = sum(ic$py),
  estimates = allm[analysis %like% "^Primary", .(measure, level, est, lo, hi, p)],
  rates = rates, per_g = pw, spline_knots = KN, spline_p_overall = w_all[[3]], spline_p_nonlinear = w_nl[[3]],
  evalue_moderate = ev_mod, evalue_trend = ev_tr, ph_global_p = ph[term == "GLOBAL", p]))
log_msg("Primary", sprintf("N %s, strokes %d; moderate RR %s; severe %s; per grade %s",
  format(nrow(ic), big.mark = ","), sum(ic$ev_stroke),
  prim[measure %like% "Rate" & level == "Moderate", fci(est, lo, hi)], prim[measure %like% "Rate" & level == "Severe", fci(est, lo, hi)],
  prim[measure %like% "Rate" & level == "Per grade (trend)", fci(est, lo, hi)]))
deviation("Imputation model", "mice removed the Nelson-Aalen cumulative hazard as collinear (logged event); BMI imputation used the event indicator and all covariates. BMI missing 3.5%; complete-case results (S1) are reported alongside.")
deviation("Uterine diagnosis group in the primary model", "Included as a confounder per the DAG (output/dag.md): it affects anemia through fibroid bleeding, and endometriosis may itself be associated with stroke. The original revision dropped it at the co-author's request; adding it changed the original estimates by <2% (reproduction table). Heavy bleeding stays out of the model (exposure determinant, S6).")
deviation("Primary estimand", "Changed from prevalent-or-incident stroke/TIA (cross-sectional, Hb +-3 y in either direction) to incident stroke excluding TIA after index, with Hb on/before index (PI decision 2026-10-22).")
deviation("Hb source", "Primary exposure uses laboratory CBC Hb only; the Codebook hgb also included point-of-care Hb (S4).")
