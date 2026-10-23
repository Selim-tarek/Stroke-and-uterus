# Audit Q11, Q12, Q14, Q15, Q16 (development data): full-process bootstrap, model comparison with uncertainty,
# calibration (O/E by expected counts, flexible calibration, risk groups), Brier/IPA, and decision curves for the actual calculator.
source("audit/R/00_audit_setup.R"); source("audit/R/val_fns.R")
md_open("A08_internal_validation.R", "A08 Internal validation, calibration, model comparison and decision curves")
d <- load_dev()
MODELS <- list(nine_item = PRED9, fourteen = PRED14, conventional = CONV, conventional_plus_medhack = c(CONV, "migraine_any", "estrogen", "haemoglobinopathy", "hb_def"))
fits <- lapply(MODELS, fitc, dd = d)

# ---- Q12: bootstrap — fixed nine predictors (as published) vs the full process (selection repeated) ---------------
set.seed(SEED)
B <- 200; sel_count <- setNames(numeric(length(PRED14)), PRED14)
cache <- file.path(DERIVED, "audit", "A08_bootstrap.rds")       # cached so that re-running the reporting part is cheap
if (file.exists(cache)) { z <- readRDS(cache); bt <- z$bt; sel_count <- z$sel } else {
bt <- rbindlist(lapply(1:B, function(b) {
  i <- sample(nrow(d), replace = TRUE); db <- d[i]
  out <- list(b = b)
  for (m in names(MODELS)) { fb <- fitc(db, MODELS[[m]]); lo <- lpof(fb, d)
    out[[paste0("opt_", m)]] <- cstat(lpof(fb, db), db)[["c"]] - cstat(lo, d)[["c"]]
    out[[paste0("slope_", m)]] <- unname(coef(coxph(Surv(t, y) ~ lo, data = d))) }
  sp <- pipeline_select(db); sel_count[sp] <<- sel_count[sp] + 1
  fp <- fitc(db, sp); lo <- lpof(fp, d)
  out$opt_pipeline <- cstat(lpof(fp, db), db)[["c"]] - cstat(lo, d)[["c"]]
  out$slope_pipeline <- unname(coef(coxph(Surv(t, y) ~ lo, data = d)))
  out$n_selected <- length(sp); out$same_as_published <- setequal(sp, PRED9)
  as.data.table(out)
}))
saveRDS(list(bt = bt, sel = sel_count), cache) }
sel_orig <- pipeline_select(d)
app <- sapply(fits, function(f) cstat(lpof(f, d), d))
iv <- rbindlist(lapply(c(names(MODELS), "pipeline"), function(m) { a <- if (m == "pipeline") app[, "nine_item"] else app[, m]; o <- bt[[paste0("opt_", m)]]; s <- bt[[paste0("slope_", m)]]
  data.table(model = m, apparent_c = sprintf("%.3f (%.3f–%.3f)", a[1], a[1] - 1.96 * a[2], a[1] + 1.96 * a[2]), mean_optimism = sprintf("%.4f (MC SE %.4f)", mean(o), sd(o) / sqrt(B)),
             optimism_corrected_c = f2(a[1] - mean(o), 3), mean_bootstrap_slope = sprintf("%.3f (MC SE %.3f)", mean(s), sd(s) / sqrt(B))) }))
iv[model == "pipeline", model := "full MED-HACK process (14-predictor fit → points → drop zero-point items → refit), apparent C of the published nine-item model"]
aud_tab(iv, "A08a_bootstrap_validation", sprintf("Harrell optimism bootstrap, %d resamples, patient-level resampling, seed %d. Apparent C with normal-approximation 95%% CI. The process row repeats predictor selection in every resample using heuristic shrinkage (an approximation of the published bootstrap-slope shrinkage).", B, SEED))
aud_tab(data.table(predictor = PRED14, selected_pct = pct(sel_count / B), in_published_set = PRED14 %in% PRED9),
        "A08b_selection_stability", sprintf("How often each candidate received non-zero points across %d bootstrap resamples. Same set as published in %s%% of resamples; the pipeline on the original data selects: %s.", B, pct(mean(bt$same_as_published)), paste(sel_orig, collapse = ", ")))

# ---- Q11/Q16: paired comparisons in the development data -----------------------------------------------------------
set.seed(SEED + 1)
pb <- rbindlist(lapply(1:500, function(b) { i <- sample(nrow(d), replace = TRUE); x <- d[i]
  cs <- sapply(fits, function(f) cstat(lpof(f, x), x)[["c"]]); as.data.table(as.list(cs)) }))
cmpC <- function(a, b2) { dd <- pb[[a]] - pb[[b2]]; sprintf("%.4f (%.4f to %.4f)", mean(dd), quantile(dd, .025), quantile(dd, .975)) }
lr <- function(f0, f1) { a <- anova(f0, f1); sprintf("χ² %.1f, df %d, P %s", a$Chisq[2], a$Df[2], fp(a$`Pr(>|Chi|)`[2])) }
dcmp <- data.table(comparison = c("nine-item − fourteen", "nine-item − conventional", "conventional + MED-HACK items − conventional"),
                   delta_apparent_c_95ci = c(cmpC("nine_item", "fourteen"), cmpC("nine_item", "conventional"), cmpC("conventional_plus_medhack", "conventional")),
                   delta_optimism_corrected_c = f2(c(app[1, "nine_item"] - mean(bt$opt_nine_item) - (app[1, "fourteen"] - mean(bt$opt_fourteen)),
                                                     app[1, "nine_item"] - mean(bt$opt_nine_item) - (app[1, "conventional"] - mean(bt$opt_conventional)),
                                                     app[1, "conventional_plus_medhack"] - mean(bt$opt_conventional_plus_medhack) - (app[1, "conventional"] - mean(bt$opt_conventional))), 4),
                   likelihood_ratio = c(lr(fits$nine_item, fitc(d, union(PRED9, PRED14))), "not nested", lr(fits$conventional, fits$conventional_plus_medhack)))
aud_tab(dcmp, "A08c_model_comparison_dev", "Paired bootstrap (500 resamples) differences in apparent C; difference in optimism-corrected C from A08a. LR for nine vs fourteen tests the five dropped predictors jointly.")

# ---- Q14: calibration of the published calculator (development data, apparent) -----------------------------------
d[, p5 := risk_spec(d)]
f9 <- fits$nine_item; S03 <- exp(-H0at(f9, 3)); d[, p3 := risk_spec(d, h_s0 = S03)]
oe5 <- oe_counts(f9, d, 5, SPEC9$shrink); oe3 <- oe_counts(f9, d, 3, SPEC9$shrink)
km5 <- oe_km(d$p5, d, 5); km3 <- oe_km(d$p3, d, 3)
fc5 <- flex_cal(d$p5, d, 5); fc3 <- flex_cal(d$p3, d, 3)
cs5 <- cal_slope(log(-log(1 - d$p5)), d, 5)
calm <- data.table(measure = c("O/E, expected counts, 5 y (95% CI)", "O/E, expected counts, 3 y", "O/E, KM risk / mean predicted, 5 y", "O/E, KM risk / mean predicted, 3 y",
                               "calibration slope on cloglog(predicted 5-y risk)", "ICI, 5 y (mean |observed − predicted|, percentage points)", "E50 / E90, 5 y (percentage points)", "ICI, 3 y"),
                   value = c(fmt_ci(oe5[c("oe", "lo", "hi")], 2), fmt_ci(oe3[c("oe", "lo", "hi")], 2), fmt_ci(km5, 2), fmt_ci(km3, 2), fmt_ci(cs5, 2),
                             f2(100 * fc5$ici, 3), sprintf("%.3f / %.3f", 100 * fc5$e50, 100 * fc5$e90), f2(100 * fc3$ici, 3)))
aud_tab(calm, "A08d_calibration_dev", "Calibration of the published nine-item calculator in the development data (apparent; shrinkage applied). The 3-year risk uses the nine-item baseline hazard at 3 years.")
grid <- seq(0.003, 0.15, length.out = 60)
set.seed(SEED + 2)
cb <- sapply(1:200, function(b) { i <- sample(nrow(d), replace = TRUE); x <- d[i]; flex_cal(x$p5, x, 5, grid)$curve })
cc <- data.table(pred = grid, obs = flex_cal(d$p5, d, 5, grid)$curve, lo = apply(cb, 1, quantile, .025), hi = apply(cb, 1, quantile, .975))
cc <- cc[pred <= quantile(d$p5, .995)]
dec <- d[, .(p = mean(p5), o = km_risk(.SD, 5)[["est"]], n = .N), by = .(g = cut(p5, quantile(p5, 0:10 / 10), include.lowest = TRUE))]
gp <- ggplot(cc, aes(100 * pred, 100 * obs)) + geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi), alpha = .2) + geom_line() + geom_abline(lty = 3) +
  geom_point(data = dec, aes(100 * p, 100 * o), shape = 1) + geom_rug(data = d[sample(.N, 3000)], aes(x = 100 * p5), inherit.aes = FALSE, alpha = .1, sides = "b") +
  labs(x = "Predicted 5-year risk, %", y = "Observed 5-year risk, %", title = "Flexible calibration, nine-item calculator (development, apparent)",
       subtitle = "Line: RCS-cloglog Cox fit, 200-bootstrap band. Circles: deciles (KM)") + theme_minimal()
ggsave(file.path(AUD, "figures", "A08_calibration_dev.png"), gp, width = 7, height = 5.5, dpi = 130)
rg <- d[, .(women = .N, strokes = sum(y), mean_predicted_pct = pct(mean(p5), 2),
            observed_km_pct = { r <- km_risk(.SD, 5); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) }),
        by = .(calculator_band = cut(p5, c(0, .01, .025, .05, .1, 1), c("<1% (below average)", "1–2.5% (around average)", "2.5–5% (elevated)", "5–10% (high)", "≥10% (high, beyond p99 when >9.8%)"), right = FALSE))][order(calculator_band)]
aud_tab(rg, "A08e_calibration_risk_bands", "Calibration within the calculator's interpretation bands (5-year, development, apparent)")

# ---- Brier / IPA at 3 and 5 years ------------------------------------------------------------------------------------
d[, ptab := 1 - SPEC14$score_baseline_survival_5y^exp(SPEC14$score_coef * (fifelse(age >= 45, -1L, 0L) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre))]
pr5 <- list(nine_item_calculator = d$p5, points_table = d$ptab, fourteen = riskh(fits$fourteen, d, 5, SPEC14$shrinkage), conventional = riskh(fits$conventional, d, 5))
bri <- rbindlist(lapply(names(pr5), function(m) { b5 <- brier(pr5[[m]], d, 5); data.table(model = m, brier_5y = signif(b5[["brier"]], 4), null_brier_5y = signif(b5[["brier_null"]], 4), ipa_5y = sprintf("%.2f%%", 100 * b5[["ipa"]])) }))
aud_tab(bri, "A08f_brier", "IPCW Brier score at 5 years (apparent) and index of prediction accuracy (IPA = 1 − Brier/Brier of the Kaplan-Meier null model)")

# ---- Q15: decision curves for the actual calculator and the points ------------------------------------------------------
th <- c(.01, .015, .02, .03, .05, .075, .10)
dca <- rbindlist(lapply(th, function(p) data.table(threshold_pct = 100 * p, treat_all = round(nb_all(d, p, 5), 5), nine_item = round(nb_at(d$p5, d, p, 5), 5),
                                                    points_table = round(nb_at(d$ptab, d, p, 5), 5), fourteen = round(nb_at(pr5$fourteen, d, p, 5), 5),
                                                    pct_flagged_nine_item = pct(mean(d$p5 >= p)))))
set.seed(SEED + 3)
nbb <- sapply(1:200, function(b) { i <- sample(nrow(d), replace = TRUE); x <- d[i]; sapply(th, function(p) nb_at(x$p5, x, p, 5) - max(0, nb_all(x, p, 5))) })
dca[, nine_item_minus_best_default_95ci := sprintf("%.5f (%.5f to %.5f)", rowMeans(nbb), apply(nbb, 1, quantile, .025), apply(nbb, 1, quantile, .975))]
aud_tab(dca, "A08g_decision_curve_dev", "Net benefit at 5 years (development, apparent; deaths censored). Last column: nine-item net benefit minus the better of treat-all and treat-none, with 200-bootstrap 95% CI. Net benefit 0.001 = 1 additional true-positive per 1,000 women at no extra false positives.")
res_put("A08", list(bootstrap = iv, comparison = dcmp, calibration = calm, dca = dca))
