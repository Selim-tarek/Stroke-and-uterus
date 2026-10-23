# Part 6D: which prediction horizon is supported, and how much does absolute risk depend on the censoring assumption?
source("post_audit/R/P00_setup.R")
md_open("P04_horizon_censoring.R", "P04 Prediction horizon and censoring")
b <- load_base(); CUT <- as.IDate("2023-01-01")
b[, lost := as.integer(y == 0 & ev_any == 0 & died == 0 & t_end < ADMIN_END - 365)]
hz <- rbindlist(lapply(list(c("all", "all"), c("anchor < 2023", "pre"), c("anchor ≥ 2023", "post")), function(g) {
  x <- if (g[2] == "all") b else if (g[2] == "pre") b[anchor < CUT] else b[anchor >= CUT]
  rbindlist(lapply(1:5, function(h) data.table(group = g[1], horizon_years = h, at_risk = x[t >= h, .N], pct_at_risk = pct(x[, mean(t >= h)]),
    events_by_h = x[y == 1 & t <= h, .N], censored_before_h = x[y == 0 & t < h, .N], deaths_before_h = x[died == 1 & t < h, .N],
    km_risk_pct = { r <- km_risk(x, h); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) })))
}))
aud_tab(hz, "P04a_risk_sets", "Numbers at risk, events (primary outcome), censoring and deaths by horizon; baseline cohort")
# pattern-mixture sensitivity: women lost to follow-up (last encounter > 1 y before the administrative end) have hazard delta × the KM hazard after censoring
km <- survfit(Surv(t, y) ~ 1, data = b); Sf <- stepfun(km$time, c(1, km$surv))
pm <- rbindlist(lapply(c(1, 3, 5), function(h) rbindlist(lapply(c(1, 1.5, 2, 3), function(delta) {
  obs <- b[y == 1 & t <= h, .N]; cen <- b[y == 0 & t < h]
  e <- cen[, ifelse(lost == 1, 1 - (Sf(h) / Sf(t))^delta, 1 - Sf(h) / Sf(t))]
  data.table(horizon_years = h, hazard_multiplier_after_loss = delta, risk_pct = pct((obs + sum(e)) / nrow(b), 3)) }))))
pm <- dcast(pm, horizon_years ~ hazard_multiplier_after_loss, value.var = "risk_pct"); setnames(pm, -1, paste0("delta_", names(pm)[-1]))
aud_tab(pm, "P04b_censoring_sensitivity", sprintf("Pattern-mixture sensitivity of absolute risk (%%). delta = 1 corresponds to non-informative censoring (≈ Kaplan-Meier). %d women were lost to follow-up (no encounter in the last year, no event or death). Untestable assumption; shown to bound the dependence on censoring.", sum(b$lost)))
# competing death at each horizon
b[, st := factor(fifelse(y == 1, "stroke", fifelse(died == 1, "death", "censor")), levels = c("censor", "stroke", "death"))]
aj <- summary(survfit(Surv(t, st) ~ 1, data = b), times = c(1, 3, 5), extend = TRUE)
aud_tab(data.table(horizon_years = c(1, 3, 5), aj_stroke_pct = pct(aj$pstate[, 2], 3), km_stroke_pct = sapply(c(1, 3, 5), function(h) pct(km_risk(b, h)[1], 3)), aj_death_pct = pct(aj$pstate[, 3], 3)),
        "P04c_competing_death", "Aalen-Johansen cumulative incidence (death competing) vs Kaplan-Meier, primary outcome")
res_put("P04", list(risk_sets = hz, pattern_mixture = pm))
