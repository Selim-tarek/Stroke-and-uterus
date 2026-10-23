# Part 5 / Part 6A: pre-specified candidate models on one cohort (baseline cohort), one outcome (acute-coded stroke) and one horizon (3 years),
# with full-process bootstrap validation (every data-dependent step repeated in each resample: fitting, spline knots, ridge penalty by CV).
source("post_audit/R/P00_setup.R")
md_open("P05_candidate_models.R", "P05 Candidate-model comparison (baseline cohort, acute-coded stroke, 3 years)")
b <- load_base(); require_cols(b, unique(c(unlist(CAND), "t", "y", "y_isch")), "baseline cohort")
MODELS <- list(
  C1_conventional_baseline = list(vars = CAND$C1_conventional_baseline, type = "cox", outcome = "y"),
  C2_current_nine_as_coded = list(vars = CAND$C2_current_nine_as_coded, type = "cox", outcome = "y"),
  C3_baseline_only         = list(vars = CAND$C3_baseline_only, type = "cox", outcome = "y"),
  C3s_baseline_only_age_spline = list(vars = CAND$C3_baseline_only, type = "cox", outcome = "y", spline = "age"),
  C4_reduced_baseline      = list(vars = CAND$C4_reduced_baseline, type = "cox", outcome = "y"),
  C5_baseline_only_ridge   = list(vars = CAND$C3_baseline_only, type = "ridge", outcome = "y"),
  C6_baseline_only_ischemic = list(vars = CAND$C3_baseline_only, type = "cox", outcome = "y_isch"))
# ---- design, fitting, Breslow baseline and prediction -------------------------------------------------------------------
mkX <- function(dd, m, knots = NULL) {
  X <- as.matrix(dd[, m$vars, with = FALSE]); storage.mode(X) <- "double"
  if (!is.null(m$spline)) { k <- if (is.null(knots)) quantile(dd$age, c(.1, .5, .9)) else knots
    X <- cbind(X[, colnames(X) != "age", drop = FALSE], Hmisc::rcspline.eval(dd$age, k, inclx = TRUE)); colnames(X)[(ncol(X) - 1):ncol(X)] <- c("age", "age'"); attr(X, "knots") <- k }
  X }
breslow <- function(lp, t, y) { o <- order(t); tt <- t[o]; yy <- y[o]; r <- rev(cumsum(rev(exp(lp[o]))))
  ut <- unique(tt[yy == 1]); h <- sapply(ut, function(s) sum(yy[tt == s]) / r[match(s, tt)]); list(time = ut, H = cumsum(h)) }
H0 <- function(bz, s) { i <- findInterval(s, bz$time); ifelse(i == 0, 0, bz$H[pmax(i, 1)]) }
fitm <- function(dd, m, nfolds = 10) {
  yv <- dd[[m$outcome]]; X <- mkX(dd, m)
  beta <- if (m$type == "cox") coef(coxph(Surv(dd$t, yv) ~ X)) else {
    cv <- cv.glmnet(X, Surv(dd$t, yv), family = "cox", alpha = 0, nfolds = nfolds); as.numeric(coef(cv, s = "lambda.min")) }
  names(beta) <- colnames(X); lp <- drop(X %*% beta)
  list(m = m, beta = beta, knots = attr(X, "knots"), bz = breslow(lp, dd$t, yv)) }
predlp <- function(fit, nd) drop(mkX(nd, fit$m, fit$knots) %*% fit$beta)
risk3 <- function(fit, nd, h = H) 1 - exp(-H0(fit$bz, h) * exp(predlp(fit, nd)))
evalm <- function(fit, nd) { x <- copy(nd)[, y := get(fit$m$outcome)]; lp <- predlp(fit, x); xt <- trunc_h(x, H)
  E <- sum(H0(fit$bz, xt$t) * exp(lp)); r <- 1 - exp(-H0(fit$bz, H) * exp(lp))
  c(c = cstat(lp, x, H)[["c"]], slope = unname(coef(coxph(Surv(xt$t, xt$y) ~ lp))), oe = sum(xt$y) / E, ipa = brier(r, x, H)[["ipa"]]) }

# ---- apparent performance and coefficients ---------------------------------------------------------------------------------
fits <- lapply(MODELS, fitm, dd = b)
app <- rbindlist(lapply(names(fits), function(n) { f <- fits[[n]]; x <- copy(b)[, y := get(f$m$outcome)]; lp <- predlp(f, x); cc <- cstat(lp, x, H); r <- risk3(f, x)
  fc <- flex_cal(r, x, H); e <- evalm(f, b)
  data.table(model = n, outcome = f$m$outcome, women = nrow(x), events_total = sum(x$y), events_3y = x[y == 1 & t <= H, .N], parameters = length(f$beta),
             events_per_parameter = round(sum(x$y) / length(f$beta), 1), apparent_c3 = sprintf("%.3f (%.3f–%.3f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]),
             apparent_oe3 = f2(e[["oe"]], 3), apparent_ici3_pp = f2(100 * fc$ici, 3), apparent_ipa3 = sprintf("%.2f%%", 100 * e[["ipa"]]), mean_risk3_pct = pct(mean(r), 2)) }))
coefs <- rbindlist(lapply(names(fits), function(n) { f <- fits[[n]]
  if (f$m$type == "cox") { x <- copy(b)[, y := get(f$m$outcome)]; X <- mkX(x, f$m); g <- coxph(Surv(x$t, x$y) ~ X); ci <- exp(confint(g))
    data.table(model = n, predictor = names(f$beta), hr = sprintf("%.2f (%.2f–%.2f)", exp(f$beta), ci[, 1], ci[, 2])) }
  else data.table(model = n, predictor = names(f$beta), hr = sprintf("%.2f (penalized; no CI)", exp(f$beta))) }))
aud_tab(dcast(coefs, predictor ~ model, value.var = "hr"), "P05a_coefficients", "Hazard ratios (95% CI, unshrunk; ridge: penalized point estimates, lambda by 10-fold CV). Age in C3s: restricted cubic spline, knots at the 10th/50th/90th percentiles (pre-specified Harrell default), coefficients not interpretable individually.")

# ---- full-process bootstrap -------------------------------------------------------------------------------------------------
B <- 200; cache <- file.path(DERIVED, "post_audit", "P05_bootstrap.rds")
if (file.exists(cache)) bt <- readRDS(cache) else {
  set.seed(SEED)
  bt <- rbindlist(lapply(1:B, function(i) { ix <- sample(nrow(b), replace = TRUE); db <- b[ix]
    rbindlist(lapply(names(MODELS), function(n) { fb <- fitm(db, MODELS[[n]], nfolds = 5); a <- evalm(fb, db); o <- evalm(fb, b)
      data.table(b = i, model = n, c_app = a[["c"]], c_test = o[["c"]], slope_test = o[["slope"]], oe_app = a[["oe"]], oe_test = o[["oe"]], ipa_app = a[["ipa"]], ipa_test = o[["ipa"]],
                 coef_name = names(fb$beta), coef = fb$beta) })) }))
  saveRDS(bt, cache) }
bs <- unique(bt[, .(b, model, c_app, c_test, slope_test, oe_app, oe_test, ipa_app, ipa_test)])
val <- bs[, .(optimism_c = mean(c_app - c_test), mc_se_c = sd(c_app - c_test) / sqrt(.N), calibration_slope = mean(slope_test),
              slope_interval = sprintf("%.2f–%.2f", quantile(slope_test, .025), quantile(slope_test, .975)), optimism_log_oe = mean(log(oe_app) - log(oe_test)),
              optimism_ipa = mean(ipa_app - ipa_test)), by = model]
val <- merge(app[, .(model, apparent_c3, apparent_ipa3)], val, by = "model")
val[, `:=`(corrected_c3 = f2(as.numeric(substr(apparent_c3, 1, 5)) - optimism_c, 3), corrected_ipa3 = sprintf("%.2f%%", as.numeric(sub("%", "", apparent_ipa3)) - 100 * optimism_ipa),
           corrected_citl_oe = f2(exp(-optimism_log_oe), 3))]
aud_tab(app, "P05b_apparent_performance", "Apparent performance at 3 years (in-sample; O/E by expected counts; ICI = mean |observed − predicted| from flexible calibration). C6 uses the ischemic outcome and is not comparable with the others.")
aud_tab(val[, .(model, apparent_c3, optimism_c = f2(optimism_c, 4), mc_se_c = f2(mc_se_c, 4), corrected_c3, calibration_slope = f2(calibration_slope, 3), slope_interval,
                corrected_citl_oe, apparent_ipa3, corrected_ipa3)],
        "P05c_full_process_validation", sprintf("Full-process bootstrap (%d resamples, patient level, seed %d). Calibration slope = mean slope of the bootstrap model's linear predictor in the original data (shrinkage needed); interval = 2.5th–97.5th percentile across resamples. corrected_citl_oe = exp(−optimism in log O/E).", B, SEED))
# paired difference in C vs the conventional comparator, and vs the current model (same resamples, test data)
pw <- dcast(bs, b ~ model, value.var = "c_test")
dC <- function(a, r) { z <- pw[[a]] - pw[[r]]; sprintf("%.4f (%.4f to %.4f)", mean(z), quantile(z, .025), quantile(z, .975)) }
aud_tab(data.table(comparison = c("C2 current nine − C1 conventional", "C3 baseline-only − C1 conventional", "C4 reduced − C1 conventional", "C5 ridge − C1 conventional",
                                  "C2 current nine − C3 baseline-only", "C3s spline − C3 linear age"),
                   delta_c3_test_95pct = c(dC("C2_current_nine_as_coded", "C1_conventional_baseline"), dC("C3_baseline_only", "C1_conventional_baseline"),
                                           dC("C4_reduced_baseline", "C1_conventional_baseline"), dC("C5_baseline_only_ridge", "C1_conventional_baseline"),
                                           dC("C2_current_nine_as_coded", "C3_baseline_only"), dC("C3s_baseline_only_age_spline", "C3_baseline_only"))),
        "P05d_paired_c_differences", "Difference in 3-year C of bootstrap-fitted models evaluated in the original data (optimism-adjusted comparison), 2.5th–97.5th percentiles")
# coefficient stability
st <- bt[!grepl("'", coef_name), .(boot_sd = sd(coef), pct_same_sign_as_full = pct(mean(sign(coef) == sign(median(coef))))), by = .(model, predictor = coef_name)]
aud_tab(dcast(st[model %in% c("C2_current_nine_as_coded", "C3_baseline_only", "C5_baseline_only_ridge")], predictor ~ model, value.var = "pct_same_sign_as_full"),
        "P05e_sign_stability", "Percentage of bootstrap resamples in which each coefficient has the same sign as its bootstrap median")
# individual prediction changes, current nine vs baseline-only
r2 <- risk3(fits$C2_current_nine_as_coded, b); r3 <- risk3(fits$C3_baseline_only, b)
bands <- function(r) cut(r, c(0, .005, .01, .02, .05, 1), c("<0.5%", "0.5–1%", "1–2%", "2–5%", "≥5%"), right = FALSE)
rc <- as.data.table(table(current_nine = bands(r2), baseline_only = bands(r3)))[N > 0]
aud_tab(dcast(rc, current_nine ~ baseline_only, value.var = "N", fill = 0), "P05f_reclassification", sprintf("3-year risk band, current nine (rows) vs baseline-only (columns). Same band for %s%% of women; ratio C2/C3 risk: median %.2f, 5th–95th percentile %.2f–%.2f.",
        pct(mean(bands(r2) == bands(r3))), median(r2 / r3), quantile(r2 / r3, .05), quantile(r2 / r3, .95)))
# DCA at 3 years: conventional, current nine, baseline-only (apparent; illustrative only — no clinical threshold is established)
th <- c(.005, .01, .015, .02, .03, .05)
aud_tab(rbindlist(lapply(th, function(p) data.table(threshold_pct = 100 * p, treat_all = round(nb_all(b, p, H), 5),
  C1 = round(nb_at(risk3(fits$C1_conventional_baseline, b), b, p, H), 5), C2 = round(nb_at(r2, b, p, H), 5), C3 = round(nb_at(r3, b, p, H), 5)))),
  "P05g_decision_curve_3y", "Apparent 3-year net benefit (illustrative; thresholds are not clinically justified decision points)")
res_put("P05", list(apparent = app, validation = val))
