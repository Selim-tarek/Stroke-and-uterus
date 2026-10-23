# Validation helpers for the audit (time-to-event; death and TIA treated as censoring, as in the published model).
fitc   <- function(dd, pred) coxph(fml_of(pred), data = dd, x = TRUE, model = TRUE)
lpof   <- function(f, nd, pred = names(coef(f))) drop(sweep(as.matrix(nd[, ..pred]), 2, f$means[pred]) %*% coef(f)[pred])
H0at   <- function(f, h) { bh <- basehaz(f, centered = TRUE); approx(bh$time, bh$hazard, xout = h, method = "constant", rule = 2, f = 0)$y }
riskh  <- function(f, nd, h, shrink = 1, pred = names(coef(f))) 1 - exp(-H0at(f, h) * exp(shrink * lpof(f, nd, pred)))
trunc_h <- function(dd, h) { x <- copy(dd); x[, `:=`(y = as.integer(y == 1 & t <= h), t = pmin(t, h))]; x }
cstat  <- function(score, dd, h = Inf) { cc <- harrell(score, dd, h); c(c = cc$concordance, se = sqrt(cc$var)) }
cal_slope <- function(lp, dd, h) { x <- trunc_h(dd, h); g <- coxph(Surv(t, y) ~ lp, data = cbind(x, lp = lp)); c(slope = unname(coef(g)), lo = confint(g)[1], hi = confint(g)[2]) }
# observed/expected by expected counts (sum of cumulative hazard to min(t, h)), Poisson CI on O
oe_counts <- function(f, dd, h, shrink = 1, pred = names(coef(f))) {
  x <- trunc_h(dd, h); bh <- basehaz(f, centered = TRUE)
  Hi <- approx(bh$time, bh$hazard, xout = x$t, method = "constant", rule = 2, f = 0)$y * exp(shrink * lpof(f, x, pred))
  O <- sum(x$y); E <- sum(Hi); ci <- c(qchisq(.025, 2 * O) / 2, qchisq(.975, 2 * O + 2) / 2)
  c(O = O, E = E, oe = O / E, lo = ci[1] / E, hi = ci[2] / E)
}
# observed/expected as KM risk at h / mean predicted risk (the published method)
oe_km <- function(pred, dd, h) { r <- km_risk(dd, h); m <- mean(pred); c(oe = r[["est"]] / m, lo = r[["lo"]] / m, hi = r[["hi"]] / m) }
# IPCW Brier score at h (censoring KM, deaths censored); index of prediction accuracy vs the KM null model
brier <- function(pred, dd, h) {
  G <- survfit(Surv(t, 1 - y) ~ 1, data = dd); Gf <- stepfun(G$time, c(1, G$surv))
  ev <- dd$y == 1 & dd$t <= h; al <- dd$t > h
  w <- numeric(nrow(dd)); w[ev] <- 1 / Gf(dd$t[ev] - 1e-8); w[al] <- 1 / Gf(h)
  bs <- mean(w * (as.numeric(ev) - pred)^2)
  p0 <- km_risk(dd, h)[["est"]]; bs0 <- mean(w * (as.numeric(ev) - p0)^2)
  c(brier = bs, brier_null = bs0, ipa = 1 - bs / bs0)
}
# flexible calibration (Austin, Harrell & van Klaveren 2020): Cox on RCS of cloglog(predicted risk)
flex_cal <- function(pred, dd, h, grid = NULL) {
  x <- trunc_h(dd, h); x[, cll := log(-log(1 - pred))]
  k <- quantile(x$cll, c(.1, .5, .9)); g <- coxph(Surv(t, y) ~ rms::rcs(cll, k), data = x)
  H <- H0at(g, h); obs <- 1 - exp(-H * exp(predict(g, type = "lp")))
  out <- list(ici = mean(abs(obs - pred)), e50 = median(abs(obs - pred)), e90 = unname(quantile(abs(obs - pred), .9)))
  if (!is.null(grid)) { nd <- data.frame(cll = log(-log(1 - grid))); out$curve <- 1 - exp(-H * exp(predict(g, newdata = nd, type = "lp"))) }
  out
}
# decision-curve net benefit (the published estimator: KM within the high-risk group)
nb_at <- function(pred, dd, pt, h) { hi <- pred >= pt; if (!any(hi)) return(0)
  s <- summary(survfit(Surv(t, y) ~ 1, data = dd[hi]), times = h, extend = TRUE)$surv
  (1 - s) * mean(hi) - s * mean(hi) * pt / (1 - pt) }
nb_all <- function(dd, pt, h) { s <- 1 - km_risk(dd, h)[["est"]]; (1 - s) - s * pt / (1 - pt) }
# MED-HACK development pipeline: 14-predictor model -> Sullivan points (B = beta_HTN/2 after heuristic shrinkage) -> keep predictors with non-zero points -> refit
heur_shrink <- function(f) { lr <- 2 * diff(f$loglik); max(0, (lr - length(coef(f))) / lr) }
pipeline_select <- function(dd) {
  f14 <- fitc(dd, PRED14); b <- coef(f14) * heur_shrink(f14); B <- b["htn"] / 2
  if (!is.finite(B) || B <= 0) return(PRED9)
  binp <- round(b[c("htn", "dm", "dyslipidemia", "smoke", "cad", "migraine_any", "ckd_pre", "haemoglobinopathy", "endo", "adeno", "estrogen")] / B)
  agep <- round(b["age"] * (c(40, 50, 57.5) - 30) / B); bmip <- round(b["bmi5"] * (c(27.5, 32.5, 39) / 5 - 22 / 5) / B)
  hbp <- round(b["hb_def"] * (c(2, 4, 6)) / B)
  c(names(binp)[binp != 0], if (any(agep != 0)) "age", if (any(bmip != 0)) "bmi5", if (any(hbp != 0)) "hb_def")
}
CONV <- c("age", "htn", "dm", "dyslipidemia", "smoke", "bmi5", "cad", "ckd_pre")   # conventional vascular risk factors available in the data
fmt_ci <- function(x, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), x[1], x[2], x[3])
