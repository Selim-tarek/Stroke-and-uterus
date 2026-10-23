# SAP addendum 2: 5-year stroke risk model and points score for women with benign uterine disease
# (development + internal and temporal validation). Writes the model specification used by the calculator.
source("R/00_setup.R")
co <- readRDS(file.path(DERIVED, "anchored_cohort.rds"))
set.seed(SEED)
d <- co[!is.na(grade) & !is.na(bmi)]
d[, `:=`(age = age_anchor, bmi5 = bmi / 5, smoke = as.integer(smoking3 == "Ever"), endo = as.integer(endometriosis %in% 1),
         adeno = as.integer(adenomyosis %in% 1), estrogen = as.integer(hormonal_type_cat %in% c("Combined OC", "Menopausal HT")),
         hb_def = pmax(13 - hb, 0), y = ev_stroke, t = py)]
PRED <- c("age", "bmi5", "htn", "dm", "dyslipidemia", "smoke", "cad", "migraine_any", "ckd_pre", "haemoglobinopathy", "endo", "adeno", "estrogen", "hb_def")
LAB <- c(age = "Age (per year)", bmi5 = "BMI (per 5 kg/m²)", htn = "Hypertension", dm = "Diabetes", dyslipidemia = "Dyslipidemia",
         smoke = "Ever smoker", cad = "Coronary artery disease", migraine_any = "Migraine", ckd_pre = "Chronic kidney disease",
         haemoglobinopathy = "Hemoglobinopathy (sickle cell / thalassemia)", endo = "Endometriosis", adeno = "Adenomyosis",
         estrogen = "Estrogen-containing hormonal therapy", hb_def = "Hemoglobin, per 1 g/dL below 13")
fml <- as.formula(paste("Surv(t, y) ~", paste(PRED, collapse = " + ")))
H <- 5

fit_model <- function(dd) coxph(fml, data = dd, x = TRUE)
risk_at <- function(f, newd, h, shrink = 1, s0 = NULL, means = NULL) {
  if (is.null(means)) means <- f$means          # must match the centring used by basehaz(centered = TRUE)
  if (is.null(s0)) { bh <- basehaz(f, centered = TRUE); s0 <- exp(-approx(bh$time, bh$hazard, xout = h, rule = 2)$y) }
  X <- as.matrix(newd[, ..PRED])
  lp <- drop(sweep(X, 2, means) %*% coef(f)) * shrink
  1 - s0^exp(lp)
}
cidx <- function(lp, dd, h = Inf) { yy <- dd$y * (dd$t <= h); tt <- pmin(dd$t, h)
  concordance(Surv(tt, yy) ~ lp, reverse = TRUE)$concordance }

# ---- development on the full cohort ---------------------------------------------------------------------------
f <- fit_model(d)
ph <- cox.zph(f)$table["GLOBAL", 3]
nl <- anova(f, coxph(update(fml, . ~ . + pspline(age, df = 3) - age), data = d))   # age non-linearity check
lp <- predict(f, type = "lp")
app_c <- cidx(lp, d)
boot <- rbindlist(lapply(1:200, function(b) {
  i <- sample(nrow(d), replace = TRUE); db <- d[i]; fb <- fit_model(db)
  lpb <- predict(fb, newdata = db, type = "lp"); lpo <- predict(fb, newdata = d, type = "lp")
  data.table(opt = cidx(lpb, db) - cidx(lpo, d), slope = coef(coxph(Surv(t, y) ~ lpo, data = d)))
}))
shrink <- mean(boot$slope)
bh <- basehaz(f, centered = TRUE); S0 <- exp(-approx(bh$time, bh$hazard, xout = H, rule = 2)$y)
means <- f$means                                     # coxph centring (binary predictors at 0), matches basehaz()
d[, risk5 := risk_at(f, d, H, shrink, S0, means)]
race_check <- anova(f, coxph(update(fml, . ~ . + race4), data = d))
coef_tab <- data.table(predictor = LAB[names(coef(f))], hazard_ratio = fci(exp(coef(f)), exp(confint(f)[, 1]), exp(confint(f)[, 2])),
                       p = fp(summary(f)$coefficients[, 5]), beta = round(coef(f), 4), beta_shrunk = round(coef(f) * shrink, 4),
                       mean_in_cohort = round(means, 4))
save_tab(coef_tab, "R1_model_coefficients", sprintf("Cox model for 5-year stroke risk (n = %s, %d strokes). Uniform shrinkage %.3f.", format(nrow(d), big.mark = ","), sum(d$y), shrink))

# ---- calibration (development, bootstrap-shrunk) -----------------------------------------------------------------
calib <- function(dd, pr, h) {
  dd <- copy(dd); dd[, pr := pr]; dd[, dec := cut(pr, unique(quantile(pr, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
  dd[, { s <- summary(survfit(Surv(t, y) ~ 1), times = h, extend = TRUE)
    .(women = .N, strokes = sum(y), predicted = 100 * mean(pr), observed = 100 * (1 - s$surv), lo = 100 * (1 - s$upper), hi = 100 * (1 - s$lower)) }, by = dec][order(dec)]
}
cal_dev <- calib(d, d$risk5, H)

# ---- temporal validation: develop on anchors before 2023, validate on 2023 onward (3-year horizon) ------------------
cut <- as.IDate("2023-01-01")
dv <- d[anchor < cut]; vv <- d[anchor >= cut]
fv <- fit_model(dv)
bt <- rbindlist(lapply(1:100, function(b) { i <- sample(nrow(dv), replace = TRUE); fb <- fit_model(dv[i])
  data.table(slope = coef(coxph(Surv(t, y) ~ predict(fb, newdata = dv, type = "lp"), data = dv))) }))
shv <- mean(bt$slope)
bhv <- basehaz(fv, centered = TRUE); S0v3 <- exp(-approx(bhv$time, bhv$hazard, xout = 3, rule = 2)$y)
pv <- risk_at(fv, vv, 3, shv, S0v3, fv$means)
lpv <- drop(sweep(as.matrix(vv[, ..PRED]), 2, fv$means) %*% coef(fv))
tv_c <- cidx(lpv, vv, 3)
tv_slope <- coef(coxph(Surv(pmin(t, 3), y * (t <= 3)) ~ lpv, data = vv))
cal_tv <- calib(vv, pv, 3)
oe <- {s <- summary(survfit(Surv(t, y) ~ 1, data = vv), times = 3, extend = TRUE); (1 - s$surv) / mean(pv)}

# ---- decision curve (5-year) ---------------------------------------------------------------------------------
nb <- function(pr, dd, pt, h) {
  hi <- pr >= pt; n <- nrow(dd)
  if (!any(hi)) return(0)
  s <- summary(survfit(Surv(t, y) ~ 1, data = dd[hi]), times = h, extend = TRUE)$surv
  tp <- (1 - s) * mean(hi); fp_ <- s * mean(hi)
  tp - fp_ * pt / (1 - pt)
}
s_all <- summary(survfit(Surv(t, y) ~ 1, data = d), times = H, extend = TRUE)$surv
dca <- rbindlist(lapply(seq(0.01, 0.10, 0.01), function(pt) data.table(
  threshold_pct = 100 * pt, model = nb(d$risk5, d, pt, H), treat_all = (1 - s_all) - s_all * pt / (1 - pt), treat_none = 0,
  pct_above_threshold = 100 * mean(d$risk5 >= pt))))

# ---- points score (Sullivan) ------------------------------------------------------------------------------------
b <- coef(f) * shrink
B <- b["htn"] / 2      # SAP specified 1 point = 5 years of age, but age had no independent (slightly inverse) association,
                       # so the constant is set so that hypertension = 2 points, as in CHA2DS2-VASc (logged deviation)
deviation("Points-score constant", sprintf("Pre-specified constant (5 years of age) unusable: adjusted age coefficient %.4f per year (HR %.2f). Constant set to half the hypertension coefficient (hypertension = 2 points).", b["age"], exp(b["age"])))
age_bands <- data.table(lab = c("18–34", "35–44", "45–54", "55–60"), lo = c(18, 35, 45, 55), hi = c(35, 45, 55, 61), mid = c(30, 40, 50, 57.5))
bmi_bands <- data.table(lab = c("<25", "25–29.9", "30–34.9", "≥35"), lo = c(0, 25, 30, 35), hi = c(25, 30, 35, 200), mid = c(22, 27.5, 32.5, 39))
hb_bands  <- data.table(lab = c("≥12", "10–11.9", "8–9.9", "<8"), lo = c(12, 10, 8, 0), hi = c(99, 12, 10, 8), mid = c(13.2, 11, 9, 7))
pts_band <- function(tab, beta, fx) { ref <- fx(tab$mid[1]); tab[, points := round(beta * (fx(mid) - ref) / B)]; tab }
age_bands <- pts_band(age_bands, b["age"], identity)
bmi_bands <- pts_band(bmi_bands, b["bmi5"], function(x) x / 5)
hb_bands  <- pts_band(hb_bands, b["hb_def"], function(x) pmax(13 - x, 0))
bin <- c("htn", "dm", "dyslipidemia", "smoke", "cad", "migraine_any", "ckd_pre", "haemoglobinopathy", "endo", "adeno", "estrogen")
bin_pts <- data.table(var = bin, label = LAB[bin], points = round(b[bin] / B))
band_pts <- function(x, tab) tab$points[pmax(1L, findInterval(x, tab$lo))]
d[, score := band_pts(age, age_bands) + band_pts(bmi, bmi_bands) + hb_bands$points[fifelse(hb >= 12, 1L, fifelse(hb >= 10, 2L, fifelse(hb >= 8, 3L, 4L)))] +
             rowSums(sweep(as.matrix(.SD), 2, bin_pts$points, `*`)), .SDcols = bin]
fs <- coxph(Surv(t, y) ~ score, data = d)
bhs <- basehaz(fs, centered = FALSE); S0s <- exp(-approx(bhs$time, bhs$hazard, xout = H, rule = 2)$y)
rng <- range(d$score)
score_risk <- data.table(score = rng[1]:rng[2])[, risk5_pct := 100 * (1 - S0s^exp(coef(fs) * score))]
obs <- d[, .(women = .N, strokes = sum(y)), by = score]
score_risk <- merge(score_risk, obs, by = "score", all.x = TRUE)
score_c <- cidx(d$score, d)
cal_score <- calib(d, 1 - S0s^exp(coef(fs) * d$score), H)

# ---- outputs --------------------------------------------------------------------------------------------------
val <- data.table(
  measure = c("Women (development)", "Strokes", "Harrell C, apparent", "Harrell C, optimism-corrected (200 bootstraps)",
              "Bootstrap calibration slope (shrinkage applied)", "Proportional hazards, global P", "Age non-linearity, P",
              "Adding race: likelihood-ratio P", "Temporal validation: development n / strokes (anchor <2023)",
              "Temporal validation: validation n / strokes (anchor ≥2023)", "Temporal validation: C at 3 years",
              "Temporal validation: calibration slope at 3 years", "Temporal validation: observed / expected 3-year risk",
              "Points score: Harrell C", "Points score range in cohort"),
  value = c(format(nrow(d), big.mark = ","), sum(d$y), f2(app_c, 3), f2(app_c - mean(boot$opt), 3), f2(shrink, 3), fp(ph),
            fp(nl$`Pr(>|Chi|)`[2]), fp(race_check$`Pr(>|Chi|)`[2]), sprintf("%s / %d", format(nrow(dv), big.mark = ","), sum(dv$y)),
            sprintf("%s / %d", format(nrow(vv), big.mark = ","), sum(vv$y)), f2(tv_c, 3), f2(tv_slope, 2), f2(oe, 2),
            f2(score_c, 3), sprintf("%d to %d", rng[1], rng[2])))
save_tab(val, "R2_validation", "Risk model validation (TRIPOD: development with internal and temporal validation)")
save_tab(rbind(cal_dev[, .(set = "Development, 5-year (shrunken model)", decile = dec, women, strokes, predicted_pct = f2(predicted), observed_pct = fci(observed, lo, hi))],
               cal_tv[, .(set = "Temporal validation, 3-year", decile = dec, women, strokes, predicted_pct = f2(predicted), observed_pct = fci(observed, lo, hi))],
               cal_score[, .(set = "Points score, 5-year", decile = dec, women, strokes, predicted_pct = f2(predicted), observed_pct = fci(observed, lo, hi))]),
         "R3_calibration", "Calibration: predicted vs observed (Kaplan-Meier) risk by decile")
save_tab(dca[, lapply(.SD, function(x) round(x, 4))], "R4_decision_curve", "Decision curve: net benefit of the model vs treat-all / treat-none (5-year)")
pts_tab <- rbind(age_bands[, .(item = paste("Age", lab), points)], bmi_bands[, .(item = paste("BMI", lab), points)],
                 hb_bands[, .(item = paste("Hemoglobin", lab, "g/dL"), points)], bin_pts[, .(item = label, points)])
save_tab(pts_tab, "R5_points", "Points score (1 point = stroke risk of 5 years of age)")
save_tab(score_risk[, .(score, women, strokes, predicted_5y_risk_pct = f2(risk5_pct))], "R6_score_to_risk", "Total points and predicted 5-year stroke risk")

spec <- list(model = "Cox, 5-year risk of stroke (excluding TIA), women 18-60 with benign uterine disease",
             horizon_years = H, baseline_survival_5y = S0, shrinkage = shrink, means = as.list(means), coef = as.list(coef(f)),
             predictors = as.list(LAB), points = list(age = age_bands, bmi = bmi_bands, hb = hb_bands, binary = bin_pts),
             score_coef = unname(coef(fs)), score_baseline_survival_5y = S0s, score_risk = score_risk[, .(score, risk5_pct)],
             validation = setNames(as.list(val$value), val$measure), cohort_5y_risk_pct = 100 * (1 - s_all),
             risk_distribution_pct = as.list(quantile(100 * d$risk5, c(.1, .25, .5, .75, .9, .99))))
write_json(spec, file.path(OUT, "risk_model_spec.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)
res_put("risk_model", list(validation = val, coefficients = coef_tab))
log_msg("Risk model", sprintf("C %.3f (optimism-corrected %.3f); temporal C %.3f; score C %.3f", app_c, app_c - mean(boot$opt), tv_c, score_c))
