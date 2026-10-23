# Part 4: outcome definitions and their effect on coefficients, discrimination and absolute risk (baseline cohort).
source("post_audit/R/P00_setup.R")
md_open("P03_outcome_definitions.R", "P03 Outcome definitions")
b <- load_base()
DEF <- c(y_all = "O1 any coded stroke excl. TIA (published definition)", y_acute = "O2 acute-coded stroke (no imaging-only, no history/sequela-only) — PRIMARY",
         y_acute30 = "O3 O2 occurring > 30 days after the anchor", y_isch = "O4 O2, ischemic", y_hem = "O5 O2, ICH or SAH",
         y_corr = "O6 O2 with dated imaging or NIHSS flowsheet within ±30 days (corroborated)")
cnt <- rbindlist(lapply(names(DEF), function(v) { x <- copy(b)[, y := get(v)]; r <- km_risk(x, H); r1 <- km_risk(x, 1)
  data.table(definition = DEF[v], events = sum(x$y), events_within_3y = x[y == 1 & t <= H, .N], km_1y_pct = pct(r1[1], 2),
             km_3y_pct = sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2))) }))
aud_tab(cnt, "P03a_outcome_counts", sprintf("Event counts and Kaplan-Meier risk under each definition (baseline cohort, n = %s). Non-qualifying events end follow-up (censored), because only the first coded event per woman is recorded.", format(nrow(b), big.mark = ",")))
cmp <- b[y_all == 1, .N, by = .(acute_code = acute, imaging_only, within_30d = gap <= 30, corroborated)][order(-N)]
aud_tab(cmp, "P03b_outcome_composition", "Composition of the O1 events")
fits <- rbindlist(lapply(c("y_all", "y_acute", "y_acute30", "y_isch"), function(v) rbindlist(lapply(c("C2_current_nine_as_coded", "C3_baseline_only"), function(m) {
  x <- copy(b)[, y := get(v)]; f <- coxph(fml_of(CAND[[m]]), data = x); cc <- cstat(predict(f, type = "lp"), x, H)
  data.table(definition = sub(" —.*| \\(.*", "", DEF[v]), model = m, events = sum(x$y), predictor = names(coef(f)),
             hr = sprintf("%.2f (%.2f–%.2f)", exp(coef(f)), exp(confint(f)[, 1]), exp(confint(f)[, 2])), c3_apparent = sprintf("%.3f", cc[1])) }))))
w <- dcast(fits, model + predictor ~ definition, value.var = "hr")
w2 <- dcast(unique(fits[, .(model, predictor = "Harrell C at 3 y (apparent)", definition, c3_apparent)]), model + predictor ~ definition, value.var = "c3_apparent")
aud_tab(rbind(w, w2), "P03c_models_by_outcome", "Unshrunk hazard ratios (95% CI) and apparent 3-year C by outcome definition. O5 (hemorrhagic) and O6 (corroborated) have too few events for 8–9-parameter models and are not fitted.")
md_note("Chart validation is not possible with the project data (no notes or charts). The prepared 100-event sample (derived/audit/chart_review_sample_events.csv) and the adjudication protocol are described in ISSUE_RESOLUTION_MATRIX.md and MISSING_DATA_REQUEST.md.")
res_put("P03", list(counts = cnt))
