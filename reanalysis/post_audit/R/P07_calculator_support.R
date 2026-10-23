# Part 7: empirical support for calculator outputs (combinations, risk range, age point). The production calculator is not changed.
source("post_audit/R/P00_setup.R")
md_open("P07_calculator_support.R", "P07 Empirical support for calculator outputs")
d <- load_dev()                       # the published calculator's development sample
BIN7 <- c("migraine_any", "estrogen", "dm", "htn", "haemoglobinopathy", "cad", "ckd_pre")
cmb <- d[, .(women = .N, strokes = sum(y)), by = BIN7]
sup <- data.table(measure = c("possible combinations of the 7 yes/no items", "combinations observed at least once", "combinations with ≥ 20 women", "combinations with ≥ 5 strokes",
                              "women in combinations with < 20 women"),
                  value = c(128, nrow(cmb), cmb[women >= 20, .N], cmb[strokes >= 5, .N], sprintf("%d (%s%%)", cmb[women < 20, sum(women)], pct(cmb[women < 20, sum(women)] / nrow(d)))))
d[, pts := as.integer(age < 45) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
d[, exact5 := risk_spec(d)]
sup <- rbind(sup, data.table(measure = c("published points range in the calculator", "highest points observed", "women with ≥ 10 points (strokes)", "exact 5-y risk, 1st / 50th / 99th percentile (%)",
                                         "women with exact 5-y risk > 10%"),
                             value = c("0–19", max(d$pts), sprintf("%d (%d)", d[pts >= 10, .N], d[pts >= 10, sum(y)]),
                                       paste(pct(quantile(d$exact5, c(.01, .5, .99)), 2), collapse = " / "), sprintf("%d (%d strokes)", d[exact5 > .1, .N], d[exact5 > .1, sum(y)]))))
aud_tab(sup, "P07a_support", "Empirical support for the current calculator's inputs and outputs (development sample, 24,006 women, 345 strokes)")
# the age point: observed 3-year risk by age band vs what the points imply
ab <- d[, .(women = .N, strokes_3y = sum(y == 1 & t <= 3), km_3y_pct = { r <- km_risk(.SD, 3); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) },
            age_point = fifelse(age[1] < 45, "+1", "0")), by = .(age_band = cut(age, c(18, 30, 35, 45, 55, 61), c("18–29", "30–34", "35–44", "45–54", "55–60"), right = FALSE))][order(age_band)]
aud_tab(ab, "P07b_age_point", "Observed 3-year risk (Kaplan-Meier) by age band vs the calculator's age point (+1 under 45)")
# 'unknown' in the source data: how often is a calculator item genuinely undocumented?
unk <- data.table(item = c("migraine: code 9 (aura status not stated, counted as yes)", "smoking (not in the calculator): unknown", "hormonal therapy: unknown (code 9)",
                           "diagnosis-based items (diabetes, hypertension, CAD, CKD, hemoglobinopathy)"),
                  women = c(d[migraine == 9, .N], d[smoking == 9, .N], d[hormonal_type == 9, .N], NA),
                  note = c("", "", "", "absence of a code is coded 'no'; the data cannot distinguish 'no' from 'not documented'"))
aud_tab(unk, "P07c_unknown_status", "How 'unknown' was represented in the development data")
res_put("P07", list(support = sup, age = ab))
