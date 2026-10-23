# Audit Q13 (and temporal parts of Q11, Q14–Q16): temporal validation with uncertainty, independence of the
# validation data, full-process temporal validation, explanation of O/E, and an updating analysis (closed testing; nothing is changed in the calculator).
source("audit/R/00_audit_setup.R"); source("audit/R/val_fns.R")
md_open("A09_temporal_validation.R", "A09 Temporal validation")
d <- load_dev(); CUT <- as.IDate("2023-01-01")
dv <- d[anchor < CUT]; vv <- d[anchor >= CUT]; H <- 3

# ---- 1. reproduce R/16 temporal results, with confidence intervals -----------------------------------------------
fv <- fitc(dv, PRED9); lpv <- lpof(fv, vv)
c3 <- cstat(lpv, vv, H); sl <- cal_slope(lpv, vv, H)
pv <- riskh(fv, vv, H); okm <- oe_km(pv, vv, H); oec <- oe_counts(fv, vv, H)
s0v3 <- exp(-approx(basehaz(fv, centered = TRUE)$time, basehaz(fv, centered = TRUE)$hazard, xout = 3, rule = 2)$y)
pv_lin <- 1 - s0v3^exp(lpv)                                     # R/16 interpolates the baseline linearly
tv <- data.table(measure = c("development (anchor < 2023): women / strokes", "validation (anchor ≥ 2023): women / strokes", "validation: maximum follow-up, years",
                             "validation: women followed ≥ 3 years", "validation: strokes within 3 years",
                             "C at 3 years (95% CI)", "calibration slope at 3 years (95% CI)", "O/E 3-year, KM / mean predicted (R/16 method) (95% CI)",
                             "O/E 3-year, expected counts (95% CI)", "mean predicted 3-year risk, %", "observed KM 3-year risk, % (95% CI)"),
                 value = c(sprintf("%s / %d", format(nrow(dv), big.mark = ","), sum(dv$y)), sprintf("%s / %d", format(nrow(vv), big.mark = ","), sum(vv$y)), f2(max(vv$t), 2),
                           format(vv[t >= 3, .N], big.mark = ","), vv[y == 1 & t <= 3, .N],
                           sprintf("%.3f (%.3f–%.3f)", c3[1], c3[1] - 1.96 * c3[2], c3[1] + 1.96 * c3[2]), fmt_ci(sl, 2), fmt_ci(oe_km(pv_lin, vv, H), 2),
                           fmt_ci(oec[c("oe", "lo", "hi")], 2), pct(mean(pv_lin), 3), { r <- km_risk(vv, H); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) }))
aud_tab(tv, "A09a_temporal_reproduction", "Temporal validation of the nine-item model refitted on anchors before 2023 (unshrunk, as in R/16_calculator_model.R:22–28). Reported: C 0.695, slope 0.90, O/E 1.56.")
md_note("Five-year predictions were never evaluated temporally: maximum follow-up in the validation set is ", f2(max(vv$t), 2), " years and ", vv[t >= 3, .N], " women reached 3 years.")

# ---- 2. independence of the validation data --------------------------------------------------------------------
sel_dv <- pipeline_select(dv)
fproc <- fitc(dv, sel_dv); lpp <- lpof(fproc, vv); cp <- cstat(lpp, vv, H); slp <- cal_slope(lpp, vv, H); oep <- oe_counts(fproc, vv, H)
ind <- data.table(item = c("predictor set chosen using data that include the validation years", "points (and therefore which predictors carry points) derived on the full sample",
                           "shrinkage factor estimated on the full sample", "temporal model in R/16 is a refit on <2023 data, not the published coefficients",
                           "full-process temporal validation: predictors selected on <2023 data only", "  C at 3 years", "  calibration slope", "  O/E 3-year, expected counts"),
                  finding = c("yes — the nine predictors were fixed after fitting the 14-predictor model on all 24,006 women (R/15, R/16)", "yes (R/15_risk_model.R:237–252)", "yes; R/16 temporal predictions are unshrunk",
                              "yes — it validates the procedure, not the exact published equation", paste(sel_dv, collapse = ", "),
                              sprintf("%.3f (%.3f–%.3f)", cp[1], cp[1] - 1.96 * cp[2], cp[1] + 1.96 * cp[2]), fmt_ci(slp, 2), fmt_ci(oep[c("oe", "lo", "hi")], 2)))
aud_tab(ind, "A09b_independence", "Was the 2023+ validation sample independent of development decisions?")
# the published equation (fit on all data, includes the validation women) applied to the validation set — not independent, shown for completeness
pp <- risk_spec(vv, h_s0 = exp(-H0at(fitc(d, PRED9), 3)))
md_note("Published equation (all-data fit, shrinkage applied, 3-year baseline from the all-data fit) applied to the 2023+ women — NOT independent: C ", f2(cstat(lpof(fitc(d, PRED9), vv), vv, H)[1], 3),
        ", O/E (KM/mean) ", fmt_ci(oe_km(pp, vv, H), 2), ".")

# ---- 3. why O/E = 1.56? event rates by period and time since anchor ----------------------------------------------
sp <- function(dd, lab) { x <- as.data.table(survSplit(Surv(t, y) ~ study_id, data = dd[, .(study_id, t, y)], cut = c(30 / 365.25, 1, 3), episode = "w"))
  x[, py := t - tstart][, .(period = lab, strokes = sum(y), person_years = round(sum(py)), rate = pois_ci(sum(y), sum(py))), by = .(window = c("0–30 days", "31 days–1 year", "1–3 years", ">3 years")[w])] }
rt <- rbind(sp(dv, "anchor < 2023"), sp(vv, "anchor ≥ 2023"))
aud_tab(rt[order(window, period)], "A09c_rates_by_time_since_anchor", "Stroke rate per 1,000 person-years by time since anchor and anchor period")
oe_y <- rbindlist(lapply(list(c(2023, 2023), c(2024, 2024), c(2025, 2026)), function(yr) { x <- vv[year(anchor) >= yr[1] & year(anchor) <= yr[2]]; o <- oe_counts(fv, x, H)
  data.table(anchor_years = paste(unique(yr), collapse = "–"), women = nrow(x), strokes_3y = o[["O"]], expected = round(o[["E"]], 1), oe = fmt_ci(o[c("oe", "lo", "hi")], 2)) }))
aud_tab(oe_y, "A09d_oe_by_anchor_year", "O/E (expected counts, 3-year) by anchor year within the validation set")
dvp <- oe_counts(fv, dv, H)
md_note("Apparent O/E in the development period itself (expected counts, 3-year): ", fmt_ci(dvp[c("oe", "lo", "hi")], 2), ".")

# ---- 4. model comparison in the temporal validation set ------------------------------------------------------------
MODELS <- list(nine_item = PRED9, fourteen = PRED14, conventional = CONV, conventional_plus_medhack = c(CONV, "migraine_any", "estrogen", "haemoglobinopathy", "hb_def"))
fvs <- lapply(MODELS, fitc, dd = dv); lps <- lapply(fvs, lpof, nd = vv)
vv[, pts := as.integer(age < 45) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
lps$points_published <- vv$pts
mc <- rbindlist(lapply(names(lps), function(m) { cc <- cstat(lps[[m]], vv, H)
  pr <- if (m == "points_published") NULL else riskh(fvs[[m]], vv, H)
  data.table(model = m, c_3y = sprintf("%.3f (%.3f–%.3f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]),
             slope = if (is.null(pr)) "" else fmt_ci(cal_slope(lps[[m]], vv, H), 2),
             oe_expected_counts = if (is.null(pr)) "" else fmt_ci(oe_counts(fvs[[m]], vv, H)[c("oe", "lo", "hi")], 2),
             ipa_3y = if (is.null(pr)) "" else sprintf("%.2f%%", 100 * brier(pr, vv, H)[["ipa"]])) }))
aud_tab(mc, "A09e_temporal_model_comparison", "All models refitted on anchors < 2023 and evaluated on ≥ 2023 at 3 years. The published points were derived on all data and are therefore not independent.")
set.seed(SEED)
pbt <- rbindlist(lapply(1:500, function(b) { i <- sample(nrow(vv), replace = TRUE); x <- vv[i]; as.data.table(lapply(lps, function(l) cstat(l[i], x, H)[["c"]])) }))
dc <- function(a, b2) { z <- pbt[[a]] - pbt[[b2]]; sprintf("%.4f (%.4f to %.4f)", mean(z), quantile(z, .025), quantile(z, .975)) }
aud_tab(data.table(comparison = c("nine-item − fourteen", "nine-item − conventional", "conventional + MED-HACK items − conventional", "nine-item − points"),
                   delta_c_3y_95ci = c(dc("nine_item", "fourteen"), dc("nine_item", "conventional"), dc("conventional_plus_medhack", "conventional"), dc("nine_item", "points_published"))),
        "A09f_temporal_delta_c", "Paired bootstrap (500 resamples of the validation set) differences in 3-year C")

# ---- 5. calibration in the validation set --------------------------------------------------------------------------
fc <- flex_cal(pv, vv, H, seq(0.002, 0.08, length.out = 40))
md_note("Flexible calibration in the validation set (3-year, nine-item temporal model): ICI ", f2(100 * fc$ici, 3), " percentage points, E50 ", f2(100 * fc$e50, 3), ", E90 ", f2(100 * fc$e90, 3), ".")
bands <- vv[, .(women = .N, strokes = sum(y & t <= H), mean_predicted_pct = pct(mean(pv[.I]), 2)), by = .(band = cut(pv, c(0, .005, .01, .02, .05, 1), c("<0.5%", "0.5–1%", "1–2%", "2–5%", "≥5%"), right = FALSE))]
vv[, pv := pv]
bands <- vv[, .(women = .N, strokes_3y = sum(y == 1 & t <= H), mean_predicted_pct = pct(mean(pv), 2),
                observed_km_pct = { r <- km_risk(.SD, H); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) }),
            by = .(predicted_3y_band = cut(pv, c(0, .005, .01, .02, .05, 1), c("<0.5%", "0.5–1%", "1–2%", "2–5%", "≥5%"), right = FALSE))][order(predicted_3y_band)]
aud_tab(bands, "A09g_temporal_calibration_bands", "Temporal validation: predicted vs observed 3-year risk by predicted-risk band")
gp <- ggplot(data.table(p = seq(0.002, 0.08, length.out = 40), o = fc$curve), aes(100 * p, 100 * o)) + geom_line() + geom_abline(lty = 3) +
  geom_point(data = vv[, .(p = mean(pv), o = km_risk(.SD, H)[["est"]]), by = .(g = cut(pv, unique(quantile(pv, 0:5 / 5)), include.lowest = TRUE))], aes(100 * p, 100 * o), shape = 1) +
  labs(x = "Predicted 3-year risk, %", y = "Observed 3-year risk, %", title = "Temporal validation (anchor ≥ 2023), nine-item model refitted on < 2023",
       subtitle = "Line: flexible calibration; circles: fifths (Kaplan-Meier)") + coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) + theme_minimal()
ggsave(file.path(AUD, "figures", "A09_calibration_temporal.png"), gp, width = 7, height = 5.5, dpi = 130)

# ---- 6. updating analysis (closed testing; separate from the calculator) ---------------------------------------------
x <- trunc_h(vv, H); x[, lp := lpv]
m_slope <- coxph(Surv(t, y) ~ lp, data = x)
m_rev <- coxph(as.formula(paste("Surv(t, y) ~ offset(lp) +", paste(PRED9, collapse = " + "))), data = x)
ll0 <- coxph(Surv(t, y) ~ offset(lp), data = x)$loglik
upd <- data.table(step = c("1. baseline only (calibration-in-the-large): O/E by expected counts", "2. slope vs baseline-only: LR test (1 df)", "3. revision of all nine coefficients vs baseline-only: LR test (9 df)",
                           "3-year risk of the reference woman (age 41, Hb deficit 0.67, no binary risk factor): baseline-only update vs development"),
                  result = c(fmt_ci(oec[c("oe", "lo", "hi")], 2), sprintf("χ² %.2f, P %s; slope %.2f", 2 * (m_slope$loglik[2] - ll0), fp(pchisq(2 * (m_slope$loglik[2] - ll0), 1, lower.tail = FALSE)), coef(m_slope)),
                             sprintf("χ² %.2f, P %s", 2 * (m_rev$loglik[2] - ll0), fp(pchisq(2 * (m_rev$loglik[2] - ll0), 9, lower.tail = FALSE))),
                             sprintf("%s%% vs %s%% (baseline-only update = development cumulative baseline hazard × O/E)", pct(1 - exp(-H0at(fv, H) * oec[["oe"]]), 3), pct(1 - exp(-H0at(fv, H)), 3))))
nf <- vv[htn + dm + cad + ckd_pre + migraine_any + estrogen + haemoglobinopathy == 0]; onf <- oe_counts(fv, nf, H); owf <- oe_counts(fv, vv[!study_id %in% nf$study_id], H)
upd <- rbind(upd, data.table(step = c("   O/E in validation women with no binary risk factor", "   O/E in validation women with ≥1 binary risk factor"),
                             result = c(sprintf("%s (n = %s, %d strokes)", fmt_ci(onf[c("oe", "lo", "hi")], 2), format(nrow(nf), big.mark = ","), onf[["O"]]),
                                        sprintf("%s (n = %s, %d strokes)", fmt_ci(owf[c("oe", "lo", "hi")], 2), format(nrow(vv) - nrow(nf), big.mark = ","), owf[["O"]]))))
upd <- rbind(upd, data.table(step = paste0("   revision: validation/development HR ratio, ", PRED9), result = sprintf("%.2f (%.2f–%.2f)", exp(coef(m_rev)), exp(confint(m_rev)[, 1]), exp(confint(m_rev)[, 2]))))
aud_tab(upd, "A09h_updating_closed_test", "Closed testing for model updating in the 2023+ data (Vergouwe et al. 2017 logic, Cox version). Any updated model must be evaluated in data not used for the update; this audit does not change the calculator.")

# ---- 7. temporal decision curve (3 years) -------------------------------------------------------------------------
th <- c(.005, .01, .015, .02, .03, .05)
aud_tab(rbindlist(lapply(th, function(p) data.table(threshold_pct = 100 * p, treat_all = round(nb_all(vv, p, H), 5), nine_item = round(nb_at(pv, vv, p, H), 5), pct_flagged = pct(mean(pv >= p))))),
        "A09i_temporal_decision_curve", "Net benefit at 3 years in the temporal validation set (nine-item model refitted on < 2023)")
res_put("A09", list(temporal = tv, independence = ind, updating = upd))
