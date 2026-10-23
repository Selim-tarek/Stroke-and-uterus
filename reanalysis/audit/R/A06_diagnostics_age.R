# Audit Q8 (regression diagnostics) and Q9 (the inverse age association).
source("audit/R/00_audit_setup.R")
md_open("A06_diagnostics_age.R", "A06 Regression diagnostics and the age association")
d <- load_dev(); co <- add_vars(load_cohort())
f <- coxph(fml_of(PRED9), data = d, x = TRUE)
hrci <- function(f, v, k = 1) sprintf("%.2f (%.2f–%.2f)", exp(k * coef(f)[v]), exp(k * confint(f)[v, 1]), exp(k * confint(f)[v, 2]))

# ---- Q8.1 proportional hazards, per predictor -----------------------------------------------------------------
z <- cox.zph(f, transform = "km")
aud_tab(data.table(term = rownames(z$table), chisq = round(z$table[, 1], 2), df = z$table[, 2], p = fp(z$table[, 3])), "A06a_schoenfeld", "Schoenfeld residual tests (KM time transform), nine-item model")
png(file.path(AUD, "figures", "A06a_schoenfeld.png"), width = 1500, height = 1000, res = 130); par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
for (i in seq_along(PRED9)) { plot(z[i], resid = FALSE, main = PRED9[i]); abline(h = coef(f)[i], lty = 3, col = "red") }; dev.off()

# ---- Q8.2 functional forms --------------------------------------------------------------------------------------
ff <- function(base, alt, lab) { a <- anova(base, alt); data.table(test = lab, lr_chisq = round(a$Chisq[2], 2), df = a$Df[2], p = fp(a$`Pr(>|Chi|)`[2])) }
d[, hb_c := hb]
fa  <- coxph(update(fml_of(PRED9), . ~ . - age + rcs(age, 4)), data = d)
fh  <- coxph(update(fml_of(PRED9), . ~ . - hb_def + rcs(hb_c, 4)), data = d)
fl  <- coxph(update(fml_of(PRED9), . ~ . - hb_def + hb_c), data = d)
f14 <- coxph(fml_of(PRED14), data = d)
fb  <- coxph(update(fml_of(PRED14), . ~ . - bmi5 + rcs(bmi, 4)), data = d)
fx  <- coxph(update(fml_of(PRED9), . ~ . + migraine_any:estrogen), data = d)
func <- rbind(ff(f, fa, "age: linear vs restricted cubic spline (4 knots)"),
              data.table(test = "Hb: hinge at 13 (published) vs RCS(Hb, 4 knots) — AIC difference (hinge − spline)", lr_chisq = round(AIC(f) - AIC(fh), 2), df = NA, p = ""),
              data.table(test = "Hb: hinge at 13 vs linear Hb — AIC difference (hinge − linear)", lr_chisq = round(AIC(f) - AIC(fl), 2), df = NA, p = ""),
              ff(f14, fb, "BMI (14-predictor model): linear per 5 vs RCS(BMI, 4 knots)"),
              ff(f, fx, "migraine × estrogen interaction (pre-specified clinically: Codebook migraine note)"))
aud_tab(func, "A06b_functional_form", "Functional-form and interaction checks (likelihood-ratio tests; AIC difference > 0 favours the alternative)")
md_note("Migraine × estrogen: women with both = ", d[migraine_any == 1 & estrogen == 1, .N], ", strokes = ", d[migraine_any == 1 & estrogen == 1, sum(y)],
        "; interaction HR ", hrci(fx, "migraine_any:estrogen"), ".")
# adjusted curves for age and Hb
curve <- function(fit, var, grid, ref, other) { nd <- d[rep(1, length(grid))][, (PRED9) := lapply(PRED9, function(v) if (v %in% names(other)) other[[v]] else 0)]
  nd[, (var) := grid]; nd[, hb_c := if (var == "hb_c") grid else 13]; nd[, age := if (var == "age") grid else 41]
  r <- d[1][, (PRED9) := lapply(PRED9, function(v) 0)][, `:=`(hb_c = 13, age = 41)]; r[, (var) := ref]
  p <- predict(fit, newdata = nd, type = "lp", se.fit = TRUE); p0 <- predict(fit, newdata = r, type = "lp")
  data.table(x = grid, hr = exp(p$fit - p0), lo = exp(p$fit - p0 - 1.96 * p$se.fit), hi = exp(p$fit - p0 + 1.96 * p$se.fit)) }
ca <- curve(fa, "age", 18:60, 41, list()); ch <- curve(fh, "hb_c", seq(6, 16, .25), 13, list())
g1 <- ggplot(ca, aes(x, hr)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) + geom_line() + geom_hline(yintercept = 1, lty = 3) + scale_y_log10() +
  labs(x = "Age at anchor, years", y = "Adjusted HR vs age 41 (RCS, 4 knots)", title = "Age") + theme_minimal()
g2 <- ggplot(ch, aes(x, hr)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) + geom_line() + geom_hline(yintercept = 1, lty = 3) + scale_y_log10() +
  labs(x = "Hemoglobin, g/dL", y = "Adjusted HR vs Hb 13 (RCS, 4 knots)", title = "Hemoglobin") + theme_minimal()
ggsave(file.path(AUD, "figures", "A06b_age_hb_curves.png"), patchwork::wrap_plots(g1, g2), width = 10, height = 4, dpi = 130)
aud_tab(rbind(ca[x %in% c(20, 25, 30, 35, 45, 50, 55, 60), .(variable = "age", value = x, hr = sprintf("%.2f (%.2f–%.2f)", hr, lo, hi))],
              ch[x %in% c(7, 8, 9, 10, 11, 12, 14, 15), .(variable = "Hb g/dL", value = x, hr = sprintf("%.2f (%.2f–%.2f)", hr, lo, hi))]),
        "A06c_spline_hr", "Adjusted HR from spline models (reference age 41, Hb 13); figure audit/output/figures/A06b_age_hb_curves.png")

# ---- Q8.3 influence, collinearity, sparse cells ------------------------------------------------------------------
dfb <- residuals(f, type = "dfbetas")
top <- order(-apply(abs(dfb), 1, max))[1:20]
f_drop <- coxph(fml_of(PRED9), data = d[-top])
infl <- data.table(predictor = PRED9, max_abs_dfbetas = round(apply(abs(dfb), 2, max), 3), hr_all = sapply(PRED9, hrci, f = f), hr_without_20_most_influential = sapply(PRED9, hrci, f = f_drop))
aud_tab(infl, "A06d_influence", "Influence: largest scaled dfbeta per predictor, and the model refitted without the 20 most influential women")
X <- as.matrix(d[, ..PRED9]); vif <- sapply(seq_along(PRED9), function(j) 1 / (1 - summary(lm(X[, j] ~ X[, -j]))$r.squared))
sp <- rbindlist(lapply(c("migraine_any", "estrogen", "dm", "htn", "haemoglobinopathy", "cad", "ckd_pre"), function(v) data.table(predictor = v, women_yes = sum(d[[v]]), strokes_yes = d[get(v) == 1, sum(y)])))
sp <- merge(sp, data.table(predictor = PRED9, vif = round(vif, 2)), by = "predictor", all = TRUE)
aud_tab(sp, "A06e_sparsity_collinearity", "Events among predictor-positive women and variance inflation factors (VIF > 5 would indicate collinearity)")
cmb <- d[, .(women = .N, strokes = sum(y)), by = .(cad, ckd_pre, haemoglobinopathy, estrogen)][cad + ckd_pre + haemoglobinopathy + estrogen >= 2][order(-women)]
aud_tab(cmb, "A06f_rare_combinations", "Women with two or more of the rare predictors (CAD, CKD, hemoglobinopathy, estrogen)")

# ---- Q9 age -------------------------------------------------------------------------------------------------------
br <- c(18, 30, 35, 40, 45, 50, 55, 61); lb <- c("18–29", "30–34", "35–39", "40–44", "45–49", "50–54", "55–60")
d[, ab := cut(age, br, lb, right = FALSE)]
crude <- d[, .(women = .N, strokes = sum(y), person_years = round(sum(t)), rate = pois_ci(sum(y), sum(t))), by = ab][order(ab)]
# attained age: split person-time
sp_d <- survSplit(Surv(age, age + t, y) ~ study_id, data = d[, .(study_id, age, t, y)], cut = c(30, 35, 40, 45, 50, 55, 60, 65, 70), episode = "band")
sp_d <- as.data.table(sp_d)[, py := tstop - age]
att <- sp_d[, .(strokes = sum(y), person_years = round(sum(py)), rate = pois_ci(sum(y), sum(py))), by = .(band = c("<30", "30–34", "35–39", "40–44", "45–49", "50–54", "55–59", "60–64", "65–69", "70+")[band])]
aud_tab(crude[, .(age_at_anchor = ab, women, strokes, person_years, rate_per_1000py = rate)], "A06g_age_crude", "Crude stroke rate by age at anchor (exact Poisson 95% CI)")
aud_tab(att[, .(attained_age = band, strokes, person_years, rate_per_1000py = rate)], "A06h_age_attained", "Stroke rate by attained age (person-time split by age during follow-up)")
young <- d[y == 1 & age < 30]
md_note("Strokes in women under 30 at the anchor (n = ", nrow(young), "): ischemic ", young[stroke_type == 1, .N], ", ICH ", young[stroke_type == 2, .N], ", SAH ", young[stroke_type == 3, .N],
        ", CVT ", young[stroke_type == 5, .N], ", unspecified ", young[stroke_type == 9, .N], "; within 30 days of the anchor ", young[stroke_date - anchor <= 30, .N],
        "; no acute code ", young[!grepl("acute", source_kinds), .N], "; pregnancy code in the year before the anchor ", young[preg_1y == 1, .N],
        "; with migraine ", young[migraine_any == 1, .N], "; with hemoglobinopathy ", young[haemoglobinopathy == 1, .N], ".")
# sensitivity of the age coefficient (HR per 10 years)
asens <- function(dd, form, lab) { g <- coxph(form, data = dd); data.table(analysis = lab, women = nrow(dd), strokes = sum(dd$y), hr_per_10_years = hrci(g, "age", 10)) }
dd30 <- copy(d)[, y := as.integer(y == 1 & stroke_date - anchor > 30)]
dda <- copy(d)[, y := as.integer(y == 1 & grepl("acute", source_kinds))]
ddi <- copy(d)[, y := as.integer(y == 1 & stroke_type == 1)]
coall <- co[age_anchor >= 18]; coall[, `:=`(y = ev_stroke, t = py)]
sens <- rbind(
  data.table(analysis = "crude (age only)", women = nrow(d), strokes = sum(d$y), hr_per_10_years = hrci(coxph(Surv(t, y) ~ age, data = d), "age", 10)),
  asens(d, fml_of(PRED9), "nine-item model (published)"),
  asens(d, update(fml_of(PRED9), . ~ . + hb_tests_2y), "+ number of CBC Hb tests in 2 y before index (utilization proxy)"),
  asens(d, update(fml_of(PRED9), . ~ . + dyslipidemia + smoke + bmi5 + endo + adeno), "14-predictor model"),
  asens(d[endometriosis == 0], fml_of(PRED9), "women without endometriosis"),
  asens(d[fibroids == 0], fml_of(PRED9), "women without fibroids"),
  asens(dd30, fml_of(PRED9), "events > 30 days after anchor only"),
  asens(dda, fml_of(PRED9), "events with an acute code only"),
  asens(ddi, fml_of(PRED9), "ischemic stroke only"),
  asens(coall, fml_of(setdiff(PRED9, "hb_def")), "all anchored women, no Hb requirement (Hb term dropped)"),
  NULL)
aud_tab(sens, "A06i_age_sensitivity", "Hazard ratio per 10 years of age at anchor under alternative specifications (unshrunk)")
# attained-age time scale (age absorbed into the baseline hazard; left-truncated entry at age at anchor), other HRs compared
d[, `:=`(age_in = age, age_out = age + t)]
gat <- coxph(Surv(age_in, age_out, y) ~ migraine_any + estrogen + dm + htn + haemoglobinopathy + hb_def + cad + ckd_pre, data = d)
P8 <- setdiff(PRED9, "age")
aud_tab(data.table(predictor = P8, hr_time_since_anchor_scale = sapply(P8, hrci, f = f), hr_attained_age_scale = sapply(P8, hrci, f = gat)),
        "A06j_attained_age_model", "Nine-item model on the published time scale (time since anchor, age as covariate) vs attained-age time scale (no age covariate)")
# fragility of the sparse binary predictors: remove the most influential predictor-positive strokes one at a time until the 95% CI includes 1
frag <- rbindlist(lapply(c("estrogen", "haemoglobinopathy", "cad", "ckd_pre"), function(v) {
  j <- which(PRED9 == v); cand <- which(d[[v]] == 1 & d$y == 1); ordr <- cand[order(-dfb[cand, j])]
  k <- 0; dd <- d
  repeat { g <- coxph(fml_of(PRED9), data = dd); if (exp(confint(g)[v, 1]) <= 1 || k >= length(ordr)) break; k <- k + 1; dd <- d[-ordr[1:k]] }
  data.table(predictor = v, positive_women = sum(d[[v]]), positive_strokes = length(cand), hr = hrci(f, v),
             strokes_removed_until_ci_includes_1 = k, hr_at_that_point = hrci(g, v))
}))
aud_tab(frag, "A06k_fragility", "Fragility: number of the most influential predictor-positive strokes (largest dfbeta) that must be removed before the 95% CI includes 1")
res_put("A06", list(zph = z$table, func = func, age = sens))
