# Audit Q6: follow-up, extrapolation to 5 years, censoring, and death as a competing event.
source("audit/R/00_audit_setup.R")
md_open("A04_followup_competing.R", "A04 Follow-up and competing risks")
d <- load_dev()
d[, death_ev := as.integer(died == 1)]
rkm <- survfit(Surv(t, 1 - y) ~ 1, data = d)          # reverse Kaplan-Meier (potential follow-up)
fu <- data.table(measure = c("median observed follow-up, years (IQR)", "median potential follow-up (reverse KM), years", "maximum follow-up, years", "total person-years"),
                 value = c(sprintf("%.2f (%.2f–%.2f)", median(d$t), quantile(d$t, .25), quantile(d$t, .75)), f2(summary(rkm)$table["median"], 2), f2(max(d$t), 2),
                           format(round(sum(d$t)), big.mark = ",")))
aud_tab(fu, "A04a_followup", "Follow-up in the development sample (n = 24,006)")
hz <- rbindlist(lapply(c(1, 3, 5), function(h) data.table(horizon_years = h, still_under_observation = d[t >= h, .N], pct_still_observed = pct(d[, mean(t >= h)]),
  strokes_by_horizon = d[y == 1 & t <= h, .N], deaths_by_horizon = d[death_ev == 1 & t <= h, .N], censored_before_horizon_without_event = d[y == 0 & death_ev == 0 & t < h, .N],
  km_stroke_risk_pct = { r <- km_risk(d, h); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) })))
aud_tab(hz, "A04b_horizons", "Numbers at 1, 3 and 5 years. Kaplan-Meier risk treats deaths as censored.")
# Aalen-Johansen cumulative incidence with death as a competing event
d[, st := factor(fifelse(y == 1, "stroke", fifelse(death_ev == 1, "death", "censor")), levels = c("censor", "stroke", "death"))]
aj <- survfit(Surv(t, st) ~ 1, data = d)
ajs <- summary(aj, times = c(1, 3, 5), extend = TRUE)
ci_tab <- data.table(horizon_years = c(1, 3, 5), km_stroke_pct = sapply(c(1, 3, 5), function(h) pct(km_risk(d, h)[1], 3)),
                     aj_stroke_pct = pct(ajs$pstate[, which(aj$states == "stroke")], 3), aj_death_pct = pct(ajs$pstate[, which(aj$states == "death")], 3))
aud_tab(ci_tab, "A04c_competing_cumulative_incidence", sprintf("Stroke risk: Kaplan-Meier (death censored) vs Aalen-Johansen cumulative incidence (death competing). Deaths without prior stroke: %d.", sum(d$death_ev)))
# Fine-Gray subdistribution model with the nine predictors vs cause-specific Cox
fg <- finegray(Surv(t, st) ~ ., data = d[, c("t", "st", PRED9), with = FALSE], etype = "stroke")
ffg <- coxph(fml_of(PRED9, "Surv(fgstart, fgstop, fgstatus)"), weights = fgwt, data = fg)
fcs <- coxph(fml_of(PRED9), data = d)
cmp <- data.table(predictor = PRED9, cause_specific_hr = sprintf("%.2f (%.2f–%.2f)", exp(coef(fcs)), exp(confint(fcs)[, 1]), exp(confint(fcs)[, 2])),
                  fine_gray_subdistribution_hr = sprintf("%.2f (%.2f–%.2f)", exp(coef(ffg)), exp(confint(ffg)[, 1]), exp(confint(ffg)[, 2])))
aud_tab(cmp, "A04d_cause_specific_vs_fine_gray", "Nine-item model: cause-specific Cox (as published, unshrunk) vs Fine-Gray")
# individual 5-year risk: published (cause-specific, 1-KM type) vs Fine-Gray CIF
bfg <- basehaz(ffg, centered = FALSE); H5 <- approx(bfg$time, bfg$hazard, xout = 5, rule = 2)$y
cif5 <- 1 - exp(-H5 * exp(drop(as.matrix(d[, ..PRED9]) %*% coef(ffg))))
p5 <- risk_spec(d)
aud_tab(data.table(measure = c("mean predicted 5-year risk, published model (%)", "mean 5-year cumulative incidence, Fine-Gray (%)", "median ratio published / Fine-Gray",
                               "women with absolute difference > 0.5 percentage points", "99th percentile of absolute difference (percentage points)"),
                   value = c(pct(mean(p5), 3), pct(mean(cif5), 3), f2(median(p5 / cif5), 3), format(sum(abs(p5 - cif5) > .005)), f2(100 * quantile(abs(p5 - cif5), .99), 3))),
        "A04e_individual_risk_competing", "Effect of treating death as competing on individual predicted 5-year risk")
# informative censoring: who is lost early (last encounter > 1 year before the administrative end without event)?
d[, lost := as.integer(y == 0 & death_ev == 0 & ev_any == 0 & t_end < ADMIN_END - 365)]
gl <- glm(fml_of(PRED9, "lost"), data = d, family = binomial)
lc <- data.table(predictor = PRED9, or_for_loss_to_follow_up = sprintf("%.2f (%.2f–%.2f)", exp(coef(gl))[PRED9], exp(confint.default(gl))[PRED9, 1], exp(confint.default(gl))[PRED9, 2]))
aud_tab(lc, "A04f_loss_to_followup", sprintf("Logistic model for loss to follow-up (last encounter > 1 year before 2026-09-25, no event; n = %d of %s). Strong associations indicate censoring that depends on the predictors; this is compatible with non-informative censoring conditional on them but cannot test it.", sum(d$lost), format(nrow(d), big.mark = ",")))
res_put("A04", list(horizons = hz, cif = ci_tab))
