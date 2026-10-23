# Parts 1, 2, 5: recompute every disputed performance number from the patient-level files and reconcile populations, outcomes and estimators.
source("final_verification/R/V00_setup.R")
md_open("V01_reconcile.R", "V01 Numerical reconciliation")
co <- add_vars(load_cohort()); dev <- load_dev(); b <- load_base()
require_cols(b, unique(c(unlist(CAND), "y", "y_all", "y_acute", "t", "anchor", "index_date", "fibroids")), "baseline cohort")
stopifnot(all(complete.cases(b[, unique(unlist(CAND)), with = FALSE])))     # no missing predictor values in the baseline cohort

# ---- 1. event-count reconciliation: 345 (audit) -> 322 (post-audit) ---------------------------------------------------------------
s <- co[!is.na(hb) & age_anchor >= 18]
ev <- data.table(step = c("A. audit development sample: CBC, BMI recorded, age ≥ 18; outcome O1 (all coded strokes excl. TIA)",
                          "  minus procedure-anchored women (anchor before first coded diagnosis)",
                          "  = index-anchored, BMI recorded",
                          "  plus index-anchored women without BMI (BMI requirement dropped post-audit)",
                          "B. post-audit baseline cohort, outcome O1",
                          "  minus imaging-only events",
                          "  minus history/sequela-only events (no acute code)",
                          "C. post-audit baseline cohort, outcome O2 (acute-coded) — used in P05/P06/P08",
                          "  of which corroborated by dated positive imaging or NIHSS flowsheet within ±30 days (O6)"),
                 women = c(nrow(dev), dev[anchor < index_date, .N], dev[anchor == index_date, .N], s[anchor == index_date & is.na(bmi), .N], nrow(b), NA, NA, nrow(b), NA),
                 events = c(sum(dev$y), dev[anchor < index_date, sum(y)], dev[anchor == index_date, sum(y)], s[anchor == index_date & is.na(bmi), sum(ev_stroke)], sum(b$y_all),
                            b[y_all == 1 & imaging_only, .N], b[y_all == 1 & !imaging_only & !acute, .N], sum(b$y_acute), sum(b$y_corr)))
ev[, events := sc(events)]
aud_tab(ev, "V01a_event_reconciliation", "From 345 to 322 events (counts 1–10 shown as <11)")
stopifnot(sum(dev$y) - dev[anchor < index_date, sum(y)] + s[anchor == index_date & is.na(bmi), sum(ev_stroke)] == sum(b$y_all),
          sum(b$y_all) - b[y_all == 1 & (imaging_only | !acute), .N] == sum(b$y_acute))

# ---- 2. recompute temporal C for every candidate with survival::concordance directly (independent of cstat()) ---------------------
dv <- b[year(anchor) <= 2022]; vv <- b[year(anchor) >= 2023]
c_direct <- function(lp, dd, h) { yy <- dd$y * (dd$t <= h); tt <- pmin(dd$t, h); k <- concordance(Surv(tt, yy) ~ lp, reverse = TRUE); c(k$concordance, sqrt(k$var)) }
CAND$C1min_age_htn_cad_ckd <- c("age", "htn", "cad", "ckd_pre")    # added for verification: C4 minus Hb deficit (isolates Hb's contribution)
fits <- lapply(CAND, fitc2, dd = dv); lps <- lapply(fits, lpof, nd = vv)
tc <- rbindlist(lapply(names(CAND), function(n) data.table(model = n, predictors = paste(CAND[[n]], collapse = ", "), dev = "anchors 2018–2022", n_dev = nrow(dv), events_dev = sum(dv$y),
  validation = "anchors 2023+", n_val = nrow(vv), events_val_3y = vv[y == 1 & t <= H, .N], c3_temporal = ci_c(c_direct(lps[[n]], vv, H)),
  ipa3_temporal = sprintf("%.2f%%", 100 * brier(riskh(fits[[n]], vv, H), vv, H)[["ipa"]]), oe3_expected_counts = fmt_ci(oe_counts(fits[[n]], vv, H)[c("oe", "lo", "hi")], 2),
  oe3_km_over_mean = fmt_ci(oe_km(riskh(fits[[n]], vv, H), vv, H), 2))))
aud_tab(tc, "V01b_temporal_recomputed", "Temporal performance recomputed (baseline cohort, outcome O2, 3 years, Harrell C with follow-up truncated at 3 years, normal-approximation CI). O/E shown by both estimators used in the reports.")

# ---- 3. paired differences: temporal (fixed dev models, resample validation set) and development (proper sampling CI) ----------------
set.seed(SEED)
pb <- rbindlist(lapply(1:1000, function(i) { ix <- sample(nrow(vv), replace = TRUE); x <- vv[ix]; as.data.table(lapply(lps, function(l) c_direct(l[ix], x, H)[1])) }))
fd <- lapply(CAND, fitc2, dd = b); lpd <- lapply(fd, lpof, nd = b)
set.seed(SEED + 1)
pd <- rbindlist(lapply(1:1000, function(i) { ix <- sample(nrow(b), replace = TRUE); x <- b[ix]; as.data.table(lapply(lpd, function(l) c_direct(l[ix], x, H)[1])) }))
dci <- function(p, a, r) { z <- p[[a]] - p[[r]]; sprintf("%.3f (%.3f to %.3f)", mean(z), quantile(z, .025), quantile(z, .975)) }
P05 <- readRDS(file.path(DERIVED, "post_audit", "P05_bootstrap.rds")); bs <- unique(P05[, .(b, model, c_test)]); pw <- dcast(bs, b ~ model, value.var = "c_test")
p05 <- function(a, r) { z <- pw[[a]] - pw[[r]]; sprintf("%.3f [%.3f to %.3f]", mean(z), quantile(z, .025), quantile(z, .975)) }
cmp <- list(c("C2_current_nine_as_coded", "C3_baseline_only"), c("C2_current_nine_as_coded", "C4_reduced_baseline"), c("C3_baseline_only", "C1_conventional_baseline"),
            c("C4_reduced_baseline", "C1_conventional_baseline"), c("C2_current_nine_as_coded", "C1_conventional_baseline"), c("C4_reduced_baseline", "C1min_age_htn_cad_ckd"),
            c("C4_reduced_baseline", "C3_baseline_only"))
dt <- rbindlist(lapply(cmp, function(z) data.table(comparison = paste(z[1], "−", z[2]),
  dev_apparent_paired_95ci = dci(pd, z[1], z[2]),
  dev_P05d_fitting_interval = if (all(z %in% names(pw))) p05(z[1], z[2]) else "not in P05",
  temporal_paired_95ci = dci(pb, z[1], z[2]))))
aud_tab(dt, "V01c_paired_differences", "Paired ΔC at 3 years, same women and outcome within each column. dev_apparent: 1,000 paired resamples of the evaluation data with models fitted on all data (sampling CI; optimistic by ≈0.01). dev_P05d: the interval reported in P05d — spread of bootstrap-fitted models evaluated on the SAME original data; it reflects fitting variability only and is NOT a sampling confidence interval. temporal: 1,000 paired resamples of the 2023+ set, models fitted on 2018–2022.")

# ---- 4. follow-up and horizon numbers in both cohorts ---------------------------------------------------------------------------
fu <- rbindlist(lapply(list(list("audit development sample (O1)", dev), list("post-audit baseline cohort (O2)", b)), function(z) { x <- z[[2]]
  data.table(population = z[[1]], women = nrow(x), events = sum(x$y), median_followup = sprintf("%.2f (%.2f–%.2f)", median(x$t), quantile(x$t, .25), quantile(x$t, .75)),
             pct_observed_5y = pct(mean(x$t >= 5)), km_3y = { r <- km_risk(x, 3); sprintf("%s (%s–%s)", pct(r[1], 2), pct(r[2], 2), pct(r[3], 2)) }) }))
aud_tab(fu, "V01d_followup", "Follow-up and horizon statistics by population and outcome")
res_put("V01", list(events = ev, temporal = tc, diffs = dt))
