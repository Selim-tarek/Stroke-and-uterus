# Part 3: does the from-scratch P05 bootstrap reproduce the cache, and what exactly were the near-infinite hbp_pre fits?
source("final_verification/R/V00_setup.R")
md_open("V06_compare_bootstrap.R", "V06 Bootstrap reproduction and sparse-fit diagnosis")
f_new <- file.path(DERIVED, "final_verification", "P05_bootstrap_rerun.rds")
if (!file.exists(f_new)) stop("Run final_verification/R/V02_rerun_p05_instrumented.R first")
old <- readRDS(file.path(DERIVED, "post_audit", "P05_bootstrap.rds")); new <- readRDS(f_new)
k <- c("b", "model", "coef_name")
m <- merge(old, new, by = k, suffixes = c("_cache", "_rerun"))
rep <- data.table(check = c("fits (resample × model) in cache", "fits in re-run", "coefficient rows matched", "max |coef cache − coef re-run|", "max |C_test cache − re-run|"),
                  value = c(uniqueN(old[, .(b, model)]), uniqueN(new[, .(b, model)]), nrow(m), signif(max(abs(m$coef_cache - m$coef_rerun)), 3), signif(max(abs(m$c_test_cache - m$c_test_rerun)), 3)))
aud_tab(rep, "V06a_reproduction", "From-scratch re-run of the P05 bootstrap (seed 20261022; 200 patient-level resamples × 7 models = 1,400 fits) vs the cached result used in the post-audit reports")
inf <- unique(new[abs(coef) > 5, .(b, model, coef_name, coef = round(coef, 1), hbp_pos_women, hbp_pos_events)])
aud_tab(inf, "V06b_near_infinite_fits", "Every fit with |β| > 5. hbp_pos_events = events among hbp_pre-positive women in that resample (0 means monotone likelihood / quasi-separation: the MLE does not exist and coxph stops at a large negative value).")
ev0 <- unique(new[, .(b, outcome = fifelse(model == "C6_baseline_only_ischemic", "ischemic (C6)", "acute-coded (all other models)"), hbp_pos_events)])
z <- ev0[hbp_pos_events == 0, .(resamples_with_zero_events = uniqueN(b), which = paste(sort(unique(b)), collapse = ", ")), by = outcome]
md_note("Resamples with zero events among hbp_pre-positive women, by outcome: ", paste(sprintf("%s: %d (resample %s)", z$outcome, z$resamples_with_zero_events, z$which), collapse = "; "),
        ". The four divergent fits are the three unpenalized models containing hbp_pre in that resample for their outcome (C3, C3s, C6 in resample 150; C6 in resample 2). The ridge model (C5) stayed finite (|β| < 5).")
bs <- unique(new[, .(b, model, c_app, c_test, slope_test)])
sens <- bs[model %in% c("C3_baseline_only", "C3s_baseline_only_age_spline", "C6_baseline_only_ischemic"),
           .(optimism_all = f2(mean(c_app - c_test), 4), optimism_excluding_affected = f2(mean((c_app - c_test)[!b %in% inf$b]), 4),
             slope_all = f2(mean(slope_test), 3), slope_excluding_affected = f2(mean(slope_test[!b %in% inf$b]), 3)), by = model]
aud_tab(sens, "V06c_sensitivity", "Effect of the affected resamples on the optimism and calibration-slope estimates")
res_put("V06", list(reproduction = rep, infinite = inf))
