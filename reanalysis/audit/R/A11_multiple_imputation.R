# Audit Q4: multiple-imputation sensitivity analysis for missing Hb and BMI in the full anchored cohort.
# Coefficient-level comparison only. Validation with imputation repeated inside each bootstrap resample is designed below but NOT executed.
source("audit/R/00_audit_setup.R")
md_open("A11_multiple_imputation.R", "A11 Multiple-imputation sensitivity analysis")
co <- add_vars(load_cohort()[age_anchor >= 18])
co[, `:=`(na_haz = nelsonaalen(as.data.frame(co), "t", "y"), anchor_year = year(anchor), black = as.integer(race4 == "Black"),
          smk_doc = as.integer(smoking3 != "Unknown"), dx = as.integer(dxgrp))]
vars <- c("hb", "bmi", "y", "na_haz", "age", "htn", "dm", "dyslipidemia", "smoke", "smk_doc", "cad", "migraine_any", "ckd_pre", "haemoglobinopathy", "endo", "adeno",
          "estrogen", "black", "dx", "anchor_year", "hb_tests_2y")
X <- as.data.frame(co[, ..vars])
M <- 10
set.seed(SEED)
imp <- mice(X, m = M, method = ifelse(colSums(is.na(X)) > 0, "pmm", ""), printFlag = FALSE, maxit = 10)
aud_log("mice logged events: ", if (is.null(imp$loggedEvents)) "none" else paste(unique(imp$loggedEvents$out), collapse = ", "))
fits <- lapply(1:M, function(k) { z <- as.data.table(complete(imp, k)); z[, `:=`(t = co$t, hb_def = pmax(13 - hb, 0))]; coxph(fml_of(PRED9), data = z) })
Q <- sapply(fits, coef); U <- lapply(fits, vcov)
qbar <- rowMeans(Q); ubar <- Reduce(`+`, U) / M; bv <- cov(t(Q)); tv <- ubar + (1 + 1 / M) * bv
se <- sqrt(diag(tv)); fmi <- (1 + 1 / M) * diag(bv) / diag(tv)
dcc <- load_dev(); fcc <- coxph(fml_of(PRED9), data = dcc)
res <- data.table(predictor = PRED9, complete_case_hr = sprintf("%.2f (%.2f–%.2f)", exp(coef(fcc)), exp(confint(fcc)[, 1]), exp(confint(fcc)[, 2])),
                  mi_hr = sprintf("%.2f (%.2f–%.2f)", exp(qbar), exp(qbar - 1.96 * se), exp(qbar + 1.96 * se)), fraction_missing_information = round(fmi, 2))
aud_tab(res, "A11a_mi_vs_complete_case", sprintf("Nine-item Cox model (unshrunk): complete case (n = %s, %d strokes) vs multiple imputation of Hb and BMI by predictive mean matching (m = %d; n = %s, %d strokes). Imputation model includes the event indicator, Nelson-Aalen hazard, all candidate predictors, race, diagnosis group, anchor year and number of Hb tests.",
        format(nrow(dcc), big.mark = ","), sum(dcc$y), M, format(nrow(co), big.mark = ","), sum(co$y)))
imp_hb <- sapply(1:M, function(k) mean(complete(imp, k)$hb[is.na(X$hb)]))
md_note("Mean imputed Hb among women with no CBC: ", f2(mean(imp_hb), 2), " g/dL (observed mean among tested women ", f2(mean(X$hb, na.rm = TRUE), 2), "). ",
        "MAR is doubtful: having no CBC in 3 years is itself informative (A02h, A02i), and predictive mean matching can only borrow values from tested women. The MI result is therefore a sensitivity analysis, not a correction.")
md_note("Designed, not executed: imputation inside resampling. For b = 1..200: draw a bootstrap sample of the 38,332 anchored women; run mice (m = 5) on it; fit the nine-item model in each imputed set and pool; estimate shrinkage; evaluate C, calibration slope and O/E on the original data imputed once with an imputation model fitted to the bootstrap sample only (or on the observed-Hb subset). Optimism = mean (apparent − test). Expected runtime is several hours.")
res_put("A11", list(mi = res))
