# Phase 6: revised time-updated analysis and landmark analysis.
source("R/00_setup.R"); source("R/model_fns.R")
N_CORES <- as.integer(Sys.getenv("N_CORES", min(4L, parallel::detectCores())))
e <- readRDS(file.path(DERIVED, "eligible.rds"))
hbl <- readRDS(file.path(DERIVED, "labs_long.rds"))$hb[src == "CBC", .(study_id, date, val)]
setorder(hbl, study_id, date)

# cohort: stroke-free at index (Hb at baseline NOT required), follow-up as in the primary analysis
w <- e[prevalent == 0 & timing_unknown == 0 & !is.na(fu_end0) & fu_end0 > index_date & !is.na(hormonal_type_cat)]
w[, `:=`(t_end = fu_end0,
         ev_stroke = as.integer(stroke_any == 1 & stroke_timing %in% 3 & stroke_type != 4 & !is.na(ev_date) & ev_date <= fu_end0),
         ev_any = as.integer(stroke_any == 1 & stroke_timing %in% 3 & !is.na(ev_date) & ev_date <= fu_end0))]
w[, py := as.numeric(t_end - index_date) / 365.25]
w[, grade := factor(fifelse(is.na(hb), "Untested", as.character(fifelse(hb >= 12, "None", fifelse(hb >= 10, "Mild", fifelse(hb >= 8, "Moderate", "Severe"))))),
                    levels = c(GRADES, "Untested"))]      # baseline status, used only as an imputation predictor
w <- prep_covs(w, "ev_stroke")

split_w <- function(W, L) {
  Le <- max(L, 1L)                         # a value is never effective on its own measurement day
  h <- merge(hbl, w[, .(study_id, index_date, t_end)], by = "study_id")
  h[, `:=`(s = date + Le, nxt = shift(date, -1) + Le), by = study_id]
  h[, e := pmin(fifelse(is.na(nxt), as.IDate("2200-01-01"), nxt), date + Le + W)]
  h[, `:=`(s = pmax(s, index_date), e = pmin(e, t_end))]
  cov <- h[e > s, .(study_id, s, e, val)]
  # gaps (no Hb in force) between index and end
  setorder(cov, study_id, s)
  bounds <- cov[, .(gs = c(w$index_date[match(study_id[1], w$study_id)], e), ge = c(s, w$t_end[match(study_id[1], w$study_id)])), by = study_id]
  gaps <- bounds[ge > gs, .(study_id, s = gs, e = ge, val = NA_real_)]
  none <- w[!study_id %in% cov$study_id, .(study_id, s = index_date, e = t_end, val = NA_real_)]
  iv <- rbind(cov, gaps, none)
  iv <- merge(iv, w[, .(study_id, t_end, ev_stroke, ev_any)], by = "study_id")
  iv[, `:=`(y_stroke = as.integer(ev_stroke == 1 & e == t_end), y_any = as.integer(ev_any == 1 & e == t_end),
            py = as.numeric(e - s) / 365.25)]
  cnt <- hbl[iv[, .(study_id, s, lo = s - 366L, hi = s - 1L)], on = .(study_id, date > lo, date <= hi), .N, by = .EACHI]
  iv[, tests_prior_y := cnt$N]
  iv[, tests_cat := cut(tests_prior_y, c(-1, 1, 3, 11, Inf), labels = c("0-1", "2-3", "4-11", "12+"))]
  iv[, tgrade := factor(fifelse(is.na(val), "No recent Hb", as.character(fifelse(val >= 12, "None", fifelse(val >= 10, "Mild", fifelse(val >= 8, "Moderate", "Severe"))))),
                        levels = c(GRADES, "No recent Hb"))]
  iv[, score := fifelse(is.na(val), 0, as.numeric(factor(as.character(tgrade), levels = GRADES)) - 1)]
  iv[, norecent := as.integer(is.na(val))]
  iv[py > 0, !c("t_end", "ev_stroke", "ev_any")]
}
# check: person-time and events are preserved by the split
chk <- split_w(365L, 0L)
stopifnot(abs(sum(chk$py) - sum(w$py)) < 1, sum(chk$y_stroke) == sum(w$ev_stroke))
log_msg("Time-updated", sprintf("cohort %s women, %.0f PY, %d strokes; split preserves person-time and events",
                                format(nrow(w), big.mark = ","), sum(w$py), sum(w$ev_stroke)))

imps <- impute_bmi(w, "ev_stroke")
covcols <- c("study_id", "age_index", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "ht", "dxgrp",
             "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y")
TG <- paste0("tgrade", c(GRADES[-1], "No recent Hb"))
pool_bv <- function(bl, Vl, terms) {               # Rubin's rules from stored coefficients / robust vcov
  B <- matrix(sapply(bl, function(b) b[terms]), nrow = length(terms)); m <- ncol(B)
  ubar <- Reduce(`+`, lapply(Vl, function(V) V[terms, terms, drop = FALSE])) / m
  bvar <- if (m > 1) cov(t(B)) else 0 * ubar
  q <- rowMeans(B); tv <- diag(ubar + (1 + 1 / m) * bvar); se <- sqrt(tv)
  r <- (1 + 1 / m) * diag(bvar) / diag(ubar); df <- (m - 1) * (1 + 1 / r)^2; df[!is.finite(df)] <- 1e6
  tq <- qt(0.975, df)
  data.table(term = terms, est = exp(q), lo = exp(q - tq * se), hi = exp(q + tq * se), p = 2 * pt(-abs(q / se), df))
}
fit_tu <- function(iv, y, extra = character(), label) {
  st <- parallel::mclapply(seq_along(imps), mc.cores = N_CORES, function(i) {
    d <- merge(iv, imps[[i]][, ..covcols], by = "study_id")
    f <- with_warnings(glm(reformulate(c("tgrade", PRIMARY_COVS, extra, "offset(log(py))"), y), poisson, data = d, model = FALSE, y = FALSE), paste(label, i))
    a <- list(b = coef(f), V = vcovCL(f, cluster = d$study_id, type = "HC0")); rm(f)
    ft <- glm(reformulate(c("score", "norecent", PRIMARY_COVS, extra, "offset(log(py))"), y), poisson, data = d, model = FALSE, y = FALSE)
    b <- list(b = coef(ft), V = vcovCL(ft, cluster = d$study_id, type = "HC0")); rm(ft, d); gc(verbose = FALSE)
    list(a = a, b = b)
  })
  a <- pool_bv(lapply(st, function(z) z$a$b), lapply(st, function(z) z$a$V), TG)
  b <- pool_bv(lapply(st, function(z) z$b$b), lapply(st, function(z) z$b$V), "score")
  ev_lv <- iv[, .(events = sum(get(y)), py = sum(py)), by = tgrade]
  out <- rbind(a[, .(level = sub("tgrade", "", term), est, lo, hi, p)], b[, .(level = "Per grade (trend, measured time)", est, lo, hi, p)])
  out[, `:=`(analysis = label, outcome = y, events = sum(iv[[y]]), person_years = sum(iv$py), women = uniqueN(iv$study_id))]
  list(tab = out, ev = ev_lv[, analysis := label])
}
specs <- list(
  list(W = 365L, L = 0L, extra = character(), label = "TU primary: 1-y carry-forward, no lag"),
  list(W = 365L, L = 30L, extra = character(), label = "TU: 30-day lag"),
  list(W = 365L, L = 90L, extra = character(), label = "TU: 90-day lag"),
  list(W = 365L, L = 0L, extra = "tests_cat", label = "TU: + Hb tests in previous year"),
  list(W = 365L, L = 30L, extra = "tests_cat", label = "TU: 30-day lag + Hb tests in previous year"),
  list(W = 1095L, L = 0L, extra = character(), label = "TU: 3-y carry-forward (original window)"))
res <- list(); evs <- list()
for (sp in specs) {
  iv <- split_w(sp$W, sp$L)
  for (y in c("y_stroke", "y_any")) {
    r <- fit_tu(iv, y, sp$extra, sp$label)
    res[[length(res) + 1]] <- r$tab; evs[[length(evs) + 1]] <- r$ev[, outcome := y]
  }
}
tu <- rbindlist(res)
tu[, outcome := fifelse(outcome == "y_stroke", "Stroke (excl. TIA)", "Stroke or TIA")]
save_tab(tu[, .(analysis, outcome, level, estimate = fci(est, lo, hi), p = fp(p), women, events, person_years = round(person_years))],
         "Table4_time_updated", "Time-updated anemia grade (most recent CBC Hb) and incident stroke; Poisson, SE clustered by woman, full adjustment (MI)")
ev_tab <- dcast(rbindlist(evs), analysis + outcome ~ tgrade, value.var = "events")
save_tab(ev_tab, "Table4a_time_updated_events", "Events by time-updated exposure level")
res_put("time_updated", tu)

# ---- landmark at 12 months (and 24) -------------------------------------------------------------------------
lm_fit <- function(months) {
  L <- as.integer(round(months * 30.4375))
  d <- w[!is.na(hb) & t_end > index_date + L & !(ev_stroke == 1 & ev_date <= index_date + L)]
  last <- merge(hbl, d[, .(study_id, index_date)], by = "study_id")[date > index_date & date <= index_date + L]
  setorder(last, study_id, -date); last <- last[, .SD[1], by = study_id][, .(study_id, last_hb = val)]
  low_any <- merge(hbl, d[, .(study_id, index_date)], by = "study_id")[date > index_date & date <= index_date + L, .(low = any(val < 12)), by = study_id]
  d <- merge(merge(d, last, by = "study_id", all.x = TRUE), low_any, by = "study_id", all.x = TRUE)
  d[, grp := factor(fifelse(hb >= 12, fifelse(low %in% TRUE, "New anemia in window", "Not anemic"),
                    fifelse(is.na(last_hb), "Anemic, no repeat Hb", fifelse(last_hb >= 12, "Anemic, resolved", "Anemic, persistent"))),
                    levels = c("Not anemic", "New anemia in window", "Anemic, resolved", "Anemic, persistent", "Anemic, no repeat Hb"))]
  d[, py := as.numeric(t_end - (index_date + L)) / 365.25]
  d <- d[!is.na(bmi)]
  f <- with_warnings(glm(reformulate(c("grp", PRIMARY_COVS, "offset(log(py))"), "ev_stroke"), poisson, data = d), paste("landmark", months))
  t <- robust_tab(f, vcovHC(f, type = "HC0"), paste0("grp", levels(d$grp)[-1]))
  # persistent vs resolved
  d2 <- copy(d); d2[, grp := relevel(grp, "Anemic, resolved")]
  f2 <- glm(reformulate(c("grp", PRIMARY_COVS, "offset(log(py))"), "ev_stroke"), poisson, data = d2)
  t2 <- robust_tab(f2, vcovHC(f2, type = "HC0"), "grpAnemic, persistent")
  g <- d[, .(women = .N, strokes = sum(ev_stroke), py = sum(py)), by = grp][order(grp)]
  g[, rate := 1000 * strokes / py]
  out <- merge(g, rbind(data.table(grp = "Not anemic", est = 1, lo = NA, hi = NA, p = NA), t[, .(grp = sub("grp", "", term), est, lo, hi, p)]), by = "grp", all.x = TRUE)
  out <- rbind(out, data.table(grp = "Persistent vs resolved", est = t2$est, lo = t2$lo, hi = t2$hi, p = t2$p), fill = TRUE)
  out[, landmark_months := months]
  out
}
lmk <- rbind(lm_fit(12), lm_fit(24))
save_tab(lmk[, .(landmark_months, group = grp, women, strokes, person_years = round(py), rate_per_1000py = f2(rate),
                 adjusted_RR = fifelse(grp == "Not anemic", "1.00 (reference)", fci(est, lo, hi)), p = fp(p))],
         "Table4b_landmark", "Landmark analysis: anemia status during the first 12/24 months after index and later stroke (complete case for BMI)")
res_put("landmark", lmk)
