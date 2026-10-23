# MED-HACK calculator model: the Cox model refitted with only the predictors that carry points (PI decision 2026-10-23),
# so the calculator asks only for the MED-HACK items. Same cohort, outcome, shrinkage and validation steps as R/15_risk_model.R.
source("R/00_setup.R")
co <- readRDS(file.path(DERIVED, "anchored_cohort.rds"))
set.seed(SEED)
d <- co[!is.na(grade) & !is.na(bmi) & age_anchor >= 18]      # same development sample as the full model
d[, `:=`(age = age_anchor, estrogen = as.integer(hormonal_type_cat %in% c("Combined OC", "Menopausal HT")), hb_def = pmax(13 - hb, 0), y = ev_stroke, t = py)]
PRED <- c("migraine_any", "estrogen", "dm", "htn", "haemoglobinopathy", "hb_def", "age", "cad", "ckd_pre")   # M E D H H A A C K
fml <- as.formula(paste("Surv(t, y) ~", paste(PRED, collapse = " + ")))
f <- coxph(fml, data = d, x = TRUE)
cidx <- function(lp, dd, h = Inf) concordance(Surv(pmin(dd$t, h), dd$y * (dd$t <= h)) ~ lp, reverse = TRUE)$concordance
app <- cidx(predict(f, type = "lp"), d)
boot <- rbindlist(lapply(1:200, function(b) { i <- sample(nrow(d), replace = TRUE); fb <- coxph(fml, data = d[i])
  lpo <- predict(fb, newdata = d, type = "lp")
  data.table(opt = cidx(predict(fb, newdata = d[i], type = "lp"), d[i]) - cidx(lpo, d), slope = coef(coxph(Surv(t, y) ~ lpo, data = d))) }))
shrink <- mean(boot$slope)
bh <- basehaz(f, centered = TRUE); S0 <- exp(-approx(bh$time, bh$hazard, xout = 5, rule = 2)$y)
risk <- 1 - S0^exp(shrink * drop(sweep(as.matrix(d[, ..PRED]), 2, f$means) %*% coef(f)))
d[, risk5 := risk]; d[, dec := cut(risk5, unique(quantile(risk5, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
cal <- d[, { s <- summary(survfit(Surv(t, y) ~ 1), times = 5, extend = TRUE); .(women = .N, strokes = sum(y), predicted_pct = 100 * mean(risk5), observed_pct = 100 * (1 - s$surv)) }, by = dec][order(dec)]
# temporal validation (develop <2023, validate >=2023, 3 years)
dv <- d[anchor < as.IDate("2023-01-01")]; vv <- d[anchor >= as.IDate("2023-01-01")]
fv <- coxph(fml, data = dv, x = TRUE); lpv <- drop(sweep(as.matrix(vv[, ..PRED]), 2, fv$means) %*% coef(fv))
tv_c <- cidx(lpv, vv, 3)
tv_slope <- coef(coxph(Surv(pmin(t, 3), y * (t <= 3)) ~ lpv, data = vv))
bhv <- basehaz(fv, centered = TRUE); S0v3 <- exp(-approx(bhv$time, bhv$hazard, xout = 3, rule = 2)$y)
pv <- 1 - S0v3^exp(lpv); obs3 <- 1 - summary(survfit(Surv(t, y) ~ 1, data = vv), times = 3, extend = TRUE)$surv
tv_oe <- obs3 / mean(pv)
full <- fromJSON(file.path(OUT, "risk_model_spec.json"))
val <- data.table(measure = c("Women / strokes", "Harrell C, apparent", "Harrell C, optimism-corrected", "Calibration slope (shrinkage)",
                              "Temporal validation C at 3 years", "Temporal validation calibration slope", "Temporal validation observed / expected 3-year risk", "Full 14-predictor model, optimism-corrected C (for comparison)"),
                  value = c(sprintf("%s / %d", format(nrow(d), big.mark = ","), sum(d$y)), f2(app, 3), f2(app - mean(boot$opt), 3), f2(shrink, 3), f2(tv_c, 3), f2(tv_slope, 2), f2(tv_oe, 2),
                            full$validation[["Harrell C, optimism-corrected (200 bootstraps)"]]))
save_tab(val, "R7_calculator_model_validation", "MED-HACK calculator model (point-carrying predictors only)")
save_tab(data.table(predictor = names(coef(f)), beta = round(coef(f), 6), hazard_ratio = fci(exp(coef(f)), exp(confint(f)[, 1]), exp(confint(f)[, 2])),
                    centre = round(f$means, 6)), "R8_calculator_model_coefficients")
save_tab(cal[, .(decile = dec, women, strokes, predicted_pct = f2(predicted_pct), observed_pct = f2(observed_pct))], "R9_calculator_model_calibration")
write_json(list(model = "MED-HACK calculator model: Cox, 5-year stroke risk, point-carrying predictors only", s0 = S0, shrink = shrink,
                coef = as.list(coef(f)), means = as.list(f$means), validation = setNames(as.list(val$value), val$measure),
                p99 = unname(quantile(100 * d$risk5, .99))),
           file.path(OUT, "calculator_model_spec.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10)
log_msg("Calculator model", sprintf("C %.3f (optimism-corrected %.3f), temporal %.3f", app, app - mean(boot$opt), tv_c))
