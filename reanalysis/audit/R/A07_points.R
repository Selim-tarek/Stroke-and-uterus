# Audit Q10: points derivation (Sullivan), the +1 shift, score-to-risk table, internal consistency and information loss.
source("audit/R/00_audit_setup.R")
md_open("A07_points.R", "A07 Points score")
d <- load_dev()
b <- unlist(SPEC14$coef) * SPEC14$shrinkage; B <- b["htn"] / 2
# ---- 1. reproduce the points from the 14-predictor shrunken coefficients ---------------------------------------
band <- function(lab, mids, beta, fx) data.table(item = lab, midpoint = mids, raw = beta * (fx(mids) - fx(mids[1])) / B)[, points := round(raw)]
pt <- rbind(band(paste("Age", c("18–34", "35–44", "45–54", "55–60")), c(30, 40, 50, 57.5), b["age"], identity),
            band(paste("Hb", c("≥12", "10–11.9", "8–9.9", "<8")), c(13.2, 11, 9, 7), b["hb_def"], function(x) pmax(13 - x, 0)),
            band(paste("BMI", c("<25", "25–29.9", "30–34.9", "≥35")), c(22, 27.5, 32.5, 39), b["bmi5"], function(x) x / 5),
            data.table(item = names(b)[!names(b) %in% c("age", "hb_def", "bmi5")], midpoint = NA, raw = b[!names(b) %in% c("age", "hb_def", "bmi5")] / B)[, points := round(raw)])
pt[, raw := round(raw, 3)]
pt[, half_point_boundary := fifelse(abs(raw - round(raw)) > 0.4, "within 0.1 of a rounding boundary", "")]
aud_tab(pt, "A07a_points_reproduced", sprintf("Points = round(shrunken β × distance from reference / B), B = β_HTN/2 = %.4f. Reproduced from risk_model_spec.json; published calculator points: migraine 2, estrogen 2, diabetes 1, HTN 2, hemoglobinopathy 3, CAD 4, CKD 2, Hb 10–11.9 and 8–9.9 = 1, <8 = 2, age ≥45 = −1 (shown as age <45 = +1).", B))

# ---- 2. the shift, and the score-to-risk equation ---------------------------------------------------------------
d[, pts_dev := fifelse(age >= 45, -1L, 0L) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
d[, pts_pub := as.integer(age < 45) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
fs <- coxph(Surv(t, y) ~ pts_dev, data = d); S0s <- exp(-approx(basehaz(fs, centered = FALSE)$time, basehaz(fs, centered = FALSE)$hazard, xout = 5, rule = 2)$y)
r6 <- fread(file.path("output", "tables", "R6_score_to_risk.csv"))
st <- d[, .(women = .N, strokes = sum(y)), by = .(published_points = pts_pub)][order(published_points)]
st[, `:=`(risk_reproduced = 100 * (1 - S0s^exp(coef(fs) * (published_points - 1))),
          risk_calculator = 100 * (1 - SPEC14$score_baseline_survival_5y^exp(SPEC14$score_coef * (published_points - 1))),
          km_observed = sapply(published_points, function(p) { r <- km_risk(d[pts_pub == p], 5); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) }))]
st <- merge(st, r6[, .(published_points = score + 1, r6_women = women, r6_risk = predicted_5y_risk_pct)], by = "published_points", all = TRUE)
st[, `:=`(risk_reproduced = round(risk_reproduced, 3), risk_calculator = round(risk_calculator, 3))]
aud_tab(st, "A07b_score_to_risk", sprintf("Score-to-risk: reproduced Cox model on the score (β = %.5f, S0(5) = %.6f) vs calculator constants and R6. Observed = Kaplan-Meier 5-year risk in each score group (95%% CI; NA if no one reached 5 years).", coef(fs), S0s))
md_note("Development score = published score − 1 for every woman: ", all(d$pts_dev == d$pts_pub - 1), ". Harrell C of the points: ", f2(harrell(d$pts_pub, d)$concordance, 4),
        " (identical under the shift, since a constant does not change ranks). Score HR per point ", sprintf("%.3f (%.3f–%.3f)", exp(coef(fs)), exp(confint(fs)[1]), exp(confint(fs)[2])), ".")

# ---- 3. what a point means: implied vs model hazard ratio ---------------------------------------------------------
b9 <- unlist(SPEC9$coef) * SPEC9$shrink
imp <- data.table(item = c("migraine", "estrogen", "diabetes", "hypertension", "hemoglobinopathy", "CAD", "CKD", "Hb 10–11.9 (vs ≥13: deficit 2 at midpoint)", "Hb <8 (deficit 6 at midpoint)", "age 30 vs 50 (20 years)"),
                  points = c(2, 2, 1, 2, 3, 4, 2, 1, 2, 1),
                  hr_nine_item_model_shrunk = round(exp(c(b9[c("migraine_any", "estrogen", "dm", "htn", "haemoglobinopathy", "cad", "ckd_pre")], 2 * b9["hb_def"], 6 * b9["hb_def"], -20 * b9["age"])), 2))
imp[, hr_implied_by_points := round(exp(coef(fs) * points), 2)]
aud_tab(imp, "A07c_point_meaning", "Relative hazard implied by the points (score HR^points) vs the nine-item model's shrunken HR. A point is a fixed relative-hazard step only in the score model; absolute risk per point grows with the starting risk.")

# ---- 4. every combination of the nine items -----------------------------------------------------------------
g <- CJ(migraine_any = 0:1, estrogen = 0:1, dm = 0:1, htn = 0:1, haemoglobinopathy = 0:1, cad = 0:1, ckd_pre = 0:1, hb = c(13, 11, 9, 7), age = c(30, 50))
g[, hb_def := pmax(13 - hb, 0)]
g[, pts := as.integer(age < 45) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
g[, `:=`(exact = 100 * risk_spec(g), table = 100 * (1 - SPEC14$score_baseline_survival_5y^exp(SPEC14$score_coef * (pts - 1))))]
cmb <- g[, .(combinations = .N, exact_min = round(min(exact), 2), exact_max = round(max(exact), 2), table_risk = round(table[1], 2)), by = pts][order(pts)]
aud_tab(cmb, "A07d_all_combinations", sprintf("All %d combinations (7 binary items × Hb 13/11/9/7 × age 30/50): range of exact-model risk within each points total. Spearman(points, exact risk) = %.3f.", nrow(g), cor(g$pts, g$exact, method = "spearman")))
md_note("Pairs of combinations where the points rank and the exact-risk rank disagree: ", { o <- outer(g$pts, g$pts, "-"); e <- outer(g$exact, g$exact, "-"); sprintf("%.1f%%", 100 * mean((o > 0 & e < 0)[upper.tri(o)] | (o < 0 & e > 0)[upper.tri(o)])) },
        ". Note: a 30-year-old woman with no risk factor scores 1 point and a 50-year-old 0 points, while the exact model gives ", f2(100 * risk_spec(data.table(migraine_any = 0, estrogen = 0, dm = 0, htn = 0, haemoglobinopathy = 0, hb_def = 0, age = 30, cad = 0, ckd_pre = 0)), 2),
        "% and ", f2(100 * risk_spec(data.table(migraine_any = 0, estrogen = 0, dm = 0, htn = 0, haemoglobinopathy = 0, hb_def = 0, age = 50, cad = 0, ckd_pre = 0)), 2), "%.")

# ---- 5. two risks per patient, and information loss ------------------------------------------------------------
d[, `:=`(exact = risk_spec(d), table = 1 - SPEC14$score_baseline_survival_5y^exp(SPEC14$score_coef * (pts_pub - 1)))]
d[, diffpp := 100 * (table - exact)]
aud_tab(data.table(measure = c("median (IQR) table − exact, percentage points", "women with |difference| > 1 percentage point", "women with table risk < half or > double the exact risk",
                               "Harrell C, exact nine-item model (apparent)", "Harrell C, points (apparent)"),
                   value = c(sprintf("%.2f (%.2f to %.2f)", median(d$diffpp), quantile(d$diffpp, .25), quantile(d$diffpp, .75)), sprintf("%d (%s%%)", sum(abs(d$diffpp) > 1), pct(mean(abs(d$diffpp) > 1))),
                             sprintf("%d (%s%%)", sum(d$table < d$exact / 2 | d$table > 2 * d$exact), pct(mean(d$table < d$exact / 2 | d$table > 2 * d$exact))),
                             f2(harrell(d$exact, d)$concordance, 4), f2(harrell(d$pts_pub, d)$concordance, 4))),
        "A07e_two_risks", "Agreement between the points-table risk (headline number in the calculator) and the exact-model risk, development sample")
set.seed(SEED)
dc <- replicate(500, { i <- sample(nrow(d), replace = TRUE); x <- d[i]; harrell(x$exact, x)$concordance - harrell(x$pts_pub, x)$concordance })
md_note("Paired bootstrap (500 resamples) difference in apparent C, exact − points: ", f2(mean(dc), 4), " (95% CI ", f2(quantile(dc, .025), 4), " to ", f2(quantile(dc, .975), 4), ").")
res_put("A07", list(points = pt, score = st))
