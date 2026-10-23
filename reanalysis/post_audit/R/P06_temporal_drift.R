# Part 6B–C: clean chronological validation of the pre-specified candidates and decomposition of calibration drift.
# Primary: develop 2018–2022, evaluate 2023+ at 3 years. Secondary (pre-specified 3-way): develop 2018–2022, baseline-only update on 2023–2024,
# evaluate 2025+ at 1 year (the only horizon that 2025+ follow-up supports). No cut-off, horizon or predictor set was chosen after seeing these results.
source("post_audit/R/P00_setup.R")
md_open("P06_temporal_drift.R", "P06 Temporal validation and calibration drift")
b <- load_base(); b[, ay := year(anchor)]
sz <- b[, .(women = .N, events = sum(y), events_1y = sum(y == 1 & t <= 1), events_3y = sum(y == 1 & t <= 3), max_followup = round(max(t), 2)),
        by = .(period = fifelse(ay <= 2022, "2018–2022", fifelse(ay <= 2024, "2023–2024", "2025–2026")))][order(period)]
aud_tab(sz, "P06a_period_sizes", "Baseline cohort by anchor period (primary outcome)")
dv <- b[ay <= 2022]; vv <- b[ay >= 2023]
fitc2 <- function(dd, v) coxph(fml_of(v), data = dd, x = TRUE, model = TRUE)
res <- rbindlist(lapply(names(CAND), function(n) { f <- fitc2(dv, CAND[[n]]); lp <- lpof(f, vv); cc <- cstat(lp, vv, H); sl <- cal_slope(lp, vv, H)
  oe <- oe_counts(f, vv, H); r <- riskh(f, vv, H); fc <- flex_cal(r, vv, H)
  data.table(model = n, c3 = sprintf("%.3f (%.3f–%.3f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]), slope = fmt_ci(sl, 2),
             oe3_expected_counts = fmt_ci(oe[c("oe", "lo", "hi")], 2), ici3_pp = f2(100 * fc$ici, 3), ipa3 = sprintf("%.2f%%", 100 * brier(r, vv, H)[["ipa"]])) }))
aud_tab(res, "P06b_temporal_primary", sprintf("Develop on anchors 2018–2022 (%s women, %d events), evaluate on 2023+ (%s women, %d events within 3 years). Unshrunk Cox refits.",
        format(nrow(dv), big.mark = ","), sum(dv$y), format(nrow(vv), big.mark = ","), vv[y == 1 & t <= H, .N]))
set.seed(SEED)
lps <- lapply(CAND, function(v) lpof(fitc2(dv, v), vv))
pb <- rbindlist(lapply(1:500, function(i) { ix <- sample(nrow(vv), replace = TRUE); x <- vv[ix]; as.data.table(lapply(lps, function(l) cstat(l[ix], x, H)[["c"]])) }))
dc <- function(a, r) { z <- pb[[a]] - pb[[r]]; sprintf("%.4f (%.4f to %.4f)", mean(z), quantile(z, .025), quantile(z, .975)) }
aud_tab(data.table(comparison = c("C2 current nine − C1 conventional", "C3 baseline-only − C1 conventional", "C4 reduced − C1 conventional", "C2 current nine − C3 baseline-only"),
                   delta_c3 = c(dc("C2_current_nine_as_coded", "C1_conventional_baseline"), dc("C3_baseline_only", "C1_conventional_baseline"),
                                dc("C4_reduced_baseline", "C1_conventional_baseline"), dc("C2_current_nine_as_coded", "C3_baseline_only"))),
        "P06c_temporal_delta_c", "Paired bootstrap (500 resamples of the 2023+ set) differences in 3-year C")

# ---- secondary: 3-way split at 1 year ----------------------------------------------------------------------------------------
up <- b[ay %in% 2023:2024]; ev <- b[ay >= 2025]
tw <- rbindlist(lapply(c("C3_baseline_only", "C2_current_nine_as_coded"), function(n) { f <- fitc2(dv, CAND[[n]])
  ou <- oe_counts(f, up, 1); o0 <- oe_counts(f, ev, 1)
  oe_upd <- c(o0[["O"]] / (o0[["E"]] * ou[["oe"]]), o0[c("lo", "hi")] / ou[["oe"]])
  cc <- cstat(lpof(f, ev), ev, 1)
  data.table(model = n, oe_1y_update_set_2023_24 = fmt_ci(ou[c("oe", "lo", "hi")], 2), evaluation_2025_events_1y = o0[["O"]],
             evaluation_c_1y = sprintf("%.3f (%.3f–%.3f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]),
             evaluation_oe_1y_original = fmt_ci(o0[c("oe", "lo", "hi")], 2), evaluation_oe_1y_after_baseline_update = fmt_ci(oe_upd, 2)) }))
aud_tab(tw, "P06d_three_way_1y", "Pre-specified 3-way chronological design at 1 year: develop 2018–2022; baseline-only recalibration factor (O/E) estimated on 2023–2024; both evaluated on 2025+. The update set is never used for evaluation.")

# ---- drift decomposition ------------------------------------------------------------------------------------------------------
co <- add_vars(load_cohort()); co[, ay := year(anchor)]
mix <- co[ay >= 2018 & anchor == index_date, .(anchored_women = .N, pct_with_cbc = pct(mean(!is.na(hb)))), by = ay]
mix <- merge(mix, b[, .(baseline_women = .N, mean_age = f2(mean(age), 1), pct_fibroids = pct(mean(fibroids == 1)), pct_migraine = pct(mean(migraine_any)), pct_dm = pct(mean(dm)),
                        pct_htn = pct(mean(htn)), pct_estrogen = pct(mean(estrogen), 2), pct_hbp_pre = pct(mean(hbp_pre), 2), mean_hb = f2(mean(hb), 2),
                        median_cbc_days_2y_pre = as.double(median(cbc_days_2y_pre)), rate_1y_per_1000py = { x <- trunc_h(.SD, 1); f2(1000 * sum(x$y) / sum(x$t), 2) }), by = ay], by = "ay")
aud_tab(mix[order(ay)], "P06e_drift_case_mix", "Case mix, measurement and first-year event rate by anchor year (baseline cohort; pct_with_cbc among all index-anchored women)")
oc <- b[y_all == 1, .(strokes_o1 = .N, pct_acute_code = pct(mean(acute & !imaging_only)), pct_within_30d = pct(mean(gap <= 30)), pct_ischemic = pct(mean(stroke_type == 1)),
                      pct_corroborated = pct(mean(corroborated))), by = .(period = fifelse(ay <= 2022, "2018–2022", "2023–2026"))]
aud_tab(oc, "P06f_drift_outcome_coding", "Outcome composition by anchor period (all O1 events)")
# effects drift: closed test for the baseline-only and current models
ct <- rbindlist(lapply(c("C3_baseline_only", "C2_current_nine_as_coded"), function(n) { f <- fitc2(dv, CAND[[n]]); x <- trunc_h(vv, H); x[, lp := lpof(f, vv)]
  ll0 <- coxph(Surv(t, y) ~ offset(lp), data = x)$loglik; m1 <- coxph(Surv(t, y) ~ lp, data = x)
  m2 <- coxph(as.formula(paste("Surv(t, y) ~ offset(lp) +", paste(CAND[[n]], collapse = " + "))), data = x)
  data.table(model = n, slope = f2(coef(m1), 2), lr_slope = sprintf("χ² %.2f, P %s", 2 * (m1$loglik[2] - ll0), fp(pchisq(2 * (m1$loglik[2] - ll0), 1, lower.tail = FALSE))),
             lr_revision = sprintf("χ² %.2f on %d df, P %s", 2 * (m2$loglik[2] - ll0), length(CAND[[n]]), fp(pchisq(2 * (m2$loglik[2] - ll0), length(CAND[[n]]), lower.tail = FALSE))),
             largest_hr_ratios = paste(head(sprintf("%s %.2f", names(coef(m2)), exp(coef(m2)))[order(-abs(coef(m2)))], 3), collapse = "; ")) }))
aud_tab(ct, "P06g_effect_drift", "Closed testing in 2023+ (3-year): slope update and full revision vs a baseline-only update. HR ratio = 2023+ effect / 2018–2022 effect.")
# is the drift concentrated in events within 30 days of the anchor? (outcome O3 = acute-coded stroke > 30 days after the anchor; models redeveloped on O3)
o3 <- rbindlist(lapply(c("C3_baseline_only", "C2_current_nine_as_coded"), function(n) { d3 <- copy(dv)[, y := y_acute30]; v3 <- copy(vv)[, y := y_acute30]
  f <- fitc2(d3, CAND[[n]]); oe <- oe_counts(f, v3, H); cc <- cstat(lpof(f, v3), v3, H)
  data.table(model = n, events_dev = sum(d3$y), events_val_3y = v3[y == 1 & t <= H, .N], c3 = sprintf("%.3f (%.3f–%.3f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]),
             oe3_expected_counts = fmt_ci(oe[c("oe", "lo", "hi")], 2)) }))
aud_tab(o3, "P06h_temporal_excluding_30d", "Temporal validation (develop 2018–2022, evaluate 2023+, 3 years) with the outcome restricted to acute-coded strokes more than 30 days after the anchor (O3)")
e30 <- b[, .(women = .N, events_0_30d = sum(y == 1 & gap <= 30), rate_0_30d_per_1000py = f2(1000 * sum(y == 1 & gap <= 30) / sum(pmin(t, 30 / 365.25)), 1)),
         by = .(period = fifelse(ay <= 2022, "2018–2022", fifelse(ay <= 2024, "2023–2024", "2025–2026")))][order(period)]
aud_tab(e30, "P06i_early_events", "Acute-coded strokes in the first 30 days after the anchor, by anchor period")
res_put("P06", list(primary = res, three_way = tw, drift = ct, o3 = o3))
