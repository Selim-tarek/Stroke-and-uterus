# SAP addendum 1 (E1-E5): prediction, Hb trajectory, anemia type, effect modification, all-cause mortality.
# Cohort: anchored cohort from R/13_anchor_design.R (stroke-free at the earliest uterine evidence, CBC Hb within 3 y up to it).
source("R/00_setup.R"); source("R/model_fns.R")
co   <- readRDS(file.path(DERIVED, "anchored_cohort.rds"))
labs <- readRDS(file.path(DERIVED, "labs_long.rds"))
hb   <- labs$hb[src == "CBC", .(study_id, date, val)]
cm   <- readRDS(file.path(DERIVED, "comorbid_long.rds"))
sets <- readRDS(file.path(DERIVED, "ncc_sets.rds"))
cw <- prep_covs(co[!is.na(grade)], "ev_stroke")
imps <- impute_bmi(cw, "ev_stroke")
set.seed(SEED)
RUN <- strsplit(Sys.getenv("RUN_E", "E1,E2,E3,E4,E5"), ",")[[1]]

if ("E1" %in% RUN) {
# ============ E1. Does Hb improve stroke prediction? =============================================================
d1 <- copy(imps[[1]])
base_f <- Surv(py, ev_stroke) ~ ns(age_index, 3) + race4 + htn + dm + dyslipidemia + smoking3 + cad + migraine_any + bmi + ckd_pre
hb_f   <- update(base_f, . ~ . + ns(hb, 3))
cstats <- function(f, d) { lp <- predict(f, newdata = d, type = "lp")
  c(harrell = concordance(Surv(d$py, d$ev_stroke) ~ lp, reverse = TRUE)$concordance,
    uno5 = concordance(Surv(pmin(d$py, 5), d$ev_stroke * (d$py <= 5)) ~ lp, reverse = TRUE, timewt = "n/G2")$concordance) }
fA <- coxph(base_f, data = d1); fB <- coxph(hb_f, data = d1)
appA <- cstats(fA, d1); appB <- cstats(fB, d1)
B <- 200
boot <- rbindlist(lapply(seq_len(B), function(b) {
  i <- sample(nrow(d1), replace = TRUE); db <- d1[i]
  fa <- coxph(base_f, data = db); fb <- coxph(hb_f, data = db)
  data.table(optA = cstats(fa, db)["harrell"] - cstats(fa, d1)["harrell"], optB = cstats(fb, db)["harrell"] - cstats(fb, d1)["harrell"],
             dC = cstats(fb, db)["harrell"] - cstats(fa, db)["harrell"],
             optA_u = cstats(fa, db)["uno5"] - cstats(fa, d1)["uno5"], optB_u = cstats(fb, db)["uno5"] - cstats(fb, d1)["uno5"],
             dC_u = cstats(fb, db)["uno5"] - cstats(fa, db)["uno5"])
}))
lr <- anova(fA, fB)
e1 <- data.table(
  metric = c("Harrell C, base model", "Harrell C, base + Hb", "Harrell C difference (Hb added)",
             "Uno C (5 y), base model", "Uno C (5 y), base + Hb", "Uno C (5 y) difference",
             "Likelihood-ratio test for Hb (3 df)"),
  apparent = c(f2(appA[1], 3), f2(appB[1], 3), f2(appB[1] - appA[1], 3), f2(appA[2], 3), f2(appB[2], 3), f2(appB[2] - appA[2], 3),
               sprintf("chi2 %.2f", lr$Chisq[2])),
  optimism_corrected = c(f2(appA[1] - mean(boot$optA), 3), f2(appB[1] - mean(boot$optB), 3),
                         f2((appB[1] - mean(boot$optB)) - (appA[1] - mean(boot$optA)), 3),
                         f2(appA[2] - mean(boot$optA_u), 3), f2(appB[2] - mean(boot$optB_u), 3),
                         f2((appB[2] - mean(boot$optB_u)) - (appA[2] - mean(boot$optA_u)), 3), ""),
  bootstrap_95ci = c("", "", sprintf("%s to %s", f2(quantile(boot$dC, .025), 3), f2(quantile(boot$dC, .975), 3)), "", "",
                     sprintf("%s to %s", f2(quantile(boot$dC_u, .025), 3), f2(quantile(boot$dC_u, .975), 3)), sprintf("P = %s", fp(lr$`Pr(>|Chi|)`[2]))))
save_tab(e1, "E1_prediction", sprintf("Does Hb improve prediction of incident stroke? Cox models, n = %s, %d strokes; %d bootstrap resamples.",
                                      format(nrow(d1), big.mark = ","), sum(d1$ev_stroke), B))
# 5-year calibration by decile (base + Hb)
s5 <- function(f, d) { sf <- survfit(f, newdata = d); 1 - as.numeric(summary(sf, times = 5, extend = TRUE)$surv) }
d1[, pred5 := s5(fB, d1)]
d1[, dec := cut(pred5, unique(quantile(pred5, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
cal <- d1[, { km <- survfit(Surv(py, ev_stroke) ~ 1); s <- summary(km, times = 5, extend = TRUE)
  .(women = .N, strokes = sum(ev_stroke), predicted_5y_pct = 100 * mean(pred5), observed_5y_pct = 100 * (1 - s$surv),
    obs_lo = 100 * (1 - s$upper), obs_hi = 100 * (1 - s$lower)) }, by = dec][order(dec)]
save_tab(cal[, .(decile = dec, women, strokes, predicted_5y_pct = f2(predicted_5y_pct), observed_5y_pct = fci(observed_5y_pct, obs_lo, obs_hi))],
         "E1a_calibration", "Calibration of the base + Hb model: predicted vs observed (Kaplan-Meier) 5-year stroke risk by decile")
res_put("E1_prediction", list(table = e1, calibration = cal))

}
if ("E2" %in% RUN) {
# ============ E2. Hb trajectory before stroke (matched sets; stroke excl. TIA) ===================================
st <- co[, .(study_id, stroke_type, ev_any)]
case_sets <- merge(sets[case == 1], st, by = "study_id")[stroke_type != 4, set]
s2 <- sets[set %in% case_sets]
x <- merge(s2, hb, by = "study_id", allow.cartesian = TRUE)[date < ref & date >= ref - 1095]
x[, bin := floor(as.numeric(ref - date - 1) / 91.3125)]          # 0 = last 3 months ... 11 = 33-36 months before
pm <- x[, .(hb = mean(val)), by = .(set, study_id, case, bin)]
traj <- pm[, .(persons = .N, mean_hb = mean(hb), lo = mean(hb) - qt(.975, .N - 1) * sd(hb) / sqrt(.N),
               hi = mean(hb) + qt(.975, .N - 1) * sd(hb) / sqrt(.N), pct_anemic = 100 * mean(hb < 12)), by = .(case, bin)][order(case, bin)]
traj[, `:=`(group = fifelse(case == 1, "Stroke", "Matched controls"), months_before = -(bin + 1) * 3)]
tested <- s2[, .N, by = case][, .(case, n_total = N)]
traj <- merge(traj, tested, by = "case")
traj[, pct_tested := 100 * persons / n_total]
win_diff <- function(bins, label) {
  w <- x[bin %in% bins, .(hb = mean(val)), by = .(set, study_id, case)]
  f <- lm(hb ~ case + factor(set), data = w)
  t <- coeftest(f, vcovCL(f, cluster = w$set))["case", ]
  data.table(window = label, cases_with_hb = w[case == 1, .N], controls_with_hb = w[case == 0, .N],
             difference_g_dl = fci(t[1], t[1] - 1.96 * t[2], t[1] + 1.96 * t[2]), p = fp(t[4]))
}
e2 <- rbind(win_diff(4:11, "36 to 12 months before (long-standing)"), win_diff(0, "Last 3 months before"),
            win_diff(1:3, "12 to 3 months before"))
save_tab(e2, "E2_trajectory_differences", "Mean Hb difference, stroke cases minus matched controls (within matched sets; SE clustered by set)")
save_tab(traj[, .(group, months_before, persons, pct_tested = f2(pct_tested, 1), mean_hb = fci(mean_hb, lo, hi), pct_anemic = f2(pct_anemic, 1))],
         "E2a_trajectory_bins", "Hb in 3-month bins before the stroke (cases) or matched reference date (controls)")
p <- ggplot(traj, aes(months_before + 1.5, mean_hb, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(Stroke = "#eb6834", `Matched controls` = "#2a78d6"), name = NULL) +
  scale_fill_manual(values = c(Stroke = "#eb6834", `Matched controls` = "#2a78d6"), name = NULL) +
  geom_hline(yintercept = 12, colour = "#52514e", linetype = 2, linewidth = 0.3) +
  annotate("text", x = -35, y = 12.03, label = "Anemia threshold 12 g/dL", hjust = 0, vjust = 0, size = 3, colour = "#52514e") +
  scale_x_continuous(breaks = seq(-36, 0, 6)) +
  labs(x = "Months before stroke (cases) or matched date (controls)", y = "Mean hemoglobin, g/dL (95% CI)",
       title = "Hemoglobin in the 3 years before stroke",
       subtitle = sprintf("%d women with stroke after uterine diagnosis and %s matched controls (same age and time since diagnosis);\nmean of each woman's Hb values in each 3-month window, among women tested in that window",
                          length(case_sets), format(s2[case == 0, .N], big.mark = ","))) +
  theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8, colour = "#52514e"))
ggsave(file.path(OUT, "figures", "FigureE2_hb_trajectory.png"), p, width = 7.5, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "figures", "FigureE2_hb_trajectory.pdf"), p, width = 7.5, height = 4.8, bg = "white")
res_put("E2_trajectory", list(differences = e2, bins = traj))

}
if ("E3" %in% RUN) {
# ============ E3. Type of anemia at the anchor ======================================================================
an <- merge(hb, cw[, .(study_id, anchor)], by = "study_id")[date <= anchor & date >= anchor - 1095]
setorder(an, study_id, -date); an <- an[, .SD[1], by = study_id][, .(study_id, hb_date_a = date)]
near_idx <- function(tab, nm) {
  z <- merge(tab, merge(an, cw[, .(study_id, anchor)], by = "study_id"), by = "study_id")[abs(date - hb_date_a) <= 30 & date <= anchor]
  z[, gap := abs(as.integer(date - hb_date_a))]; setorder(z, study_id, gap)
  setnames(z[, .SD[1], by = study_id][, .(study_id, val)], "val", nm) }
ix <- Reduce(function(a, b) merge(a, b, by = "study_id", all = TRUE), list(an, near_idx(labs$mcv, "mcv_a"), near_idx(labs$rdw, "rdw_a")))
kk <- merge(cm[flag %in% c("anaemia_iron", "anaemia_blood_loss_acute")], cw[, .(study_id, anchor)], by = "study_id")[date <= anchor & date >= anchor - 365]
ix[, coded_ida := as.integer(study_id %in% kk$study_id)]
add_type <- function(z) {
  z <- merge(copy(z), ix[, .(study_id, mcv_a, rdw_a, coded_ida)], by = "study_id", all.x = TRUE)
  z[is.na(coded_ida), coded_ida := 0L]
  z[, cause := factor(fifelse(grade == "None", "No anemia",
                      fifelse(coded_ida == 1 | (!is.na(mcv_a) & mcv_a < 80), "Anemia with evident cause (iron deficiency / blood loss)",
                              "Anemia without evident cause")),
                      levels = c("No anemia", "Anemia with evident cause (iron deficiency / blood loss)", "Anemia without evident cause"))]
  z[, morph := factor(fifelse(grade == "None", "No anemia", fifelse(is.na(mcv_a), "Anemia, MCV missing",
                      fifelse(mcv_a < 80, "Microcytic anemia", fifelse(mcv_a <= 100, "Normocytic anemia", "Macrocytic anemia")))),
                      levels = c("No anemia", "Microcytic anemia", "Normocytic anemia", "Macrocytic anemia", "Anemia, MCV missing"))]
  z[, rdw_high := as.integer(rdw_a > 14.5)]
  z }
imps3 <- lapply(imps, add_type)
lvl_tab <- function(var, label, covs = PRIMARY_COVS) {
  levs <- levels(imps3[[1]][[var]])[-1]
  ev <- imps3[[1]][, .(women = .N, strokes = sum(ev_stroke)), by = var]
  ok <- levs[levs %in% ev[strokes >= 5, get(var)]]
  dd <- lapply(imps3, function(z) z[get(var) %in% c(levels(z[[var]])[1], ok)][, (var) := droplevels(get(var))])
  pf <- pois_fit(dd, "ev_stroke", var, covs, label)
  r <- pooled(pf, paste0(var, ok))
  out <- merge(ev[, .(level = as.character(get(var)), women, strokes)],
               rbind(data.table(level = levels(imps3[[1]][[var]])[1], est = 1, lo = NA, hi = NA, p = NA), r[, .(level = sub(var, "", term), est, lo, hi, p)]),
               by = "level", all.x = TRUE)
  out[, `:=`(analysis = label, estimate = fifelse(is.na(est), "not estimated (<5 strokes)", fifelse(is.na(lo), "1.00 (reference)", fci(est, lo, hi))))]
  out[order(match(level, levels(imps3[[1]][[var]])))]
}
e3 <- rbind(lvl_tab("cause", "Anemia by evident cause"), lvl_tab("morph", "Anemia by red-cell size (MCV)"), fill = TRUE)
# RDW, adjusted for Hb (continuous spline), among women with RDW
rdw_d <- lapply(imps3, function(z) z[!is.na(rdw_a)])
r_hi <- pooled(pois_fit(rdw_d, "ev_stroke", c("rdw_high", "ns(hb, 3)"), PRIMARY_COVS, "RDW high"), "rdw_high")
r_c  <- pooled(pois_fit(rdw_d, "ev_stroke", c("rdw_a", "ns(hb, 3)"), PRIMARY_COVS, "RDW per 1%"), "rdw_a")
e3 <- rbind(e3, data.table(analysis = "RDW-CV, adjusted for Hb", level = c("RDW-CV >14.5% vs ≤14.5%", "Per 1% higher RDW-CV"),
                           women = nrow(rdw_d[[1]]), strokes = sum(rdw_d[[1]]$ev_stroke), estimate = c(fci(r_hi$est, r_hi$lo, r_hi$hi), fci(r_c$est, r_c$lo, r_c$hi)),
                           p = c(r_hi$p, r_c$p)), fill = TRUE)
save_tab(e3[, .(analysis, level, women, strokes, adjusted_RR = estimate, p = fp(p))], "E3_anemia_type",
         "Type of anemia at the anchor and incident stroke (Poisson, full adjustment, MI)")
res_put("E3_anemia_type", e3)

}
if ("E4" %in% RUN) {
# ============ E4. Effect modification (any anemia) =====================================================================
imps4 <- lapply(imps, function(z) { z <- copy(z); z[, `:=`(anemia_any = as.integer(grade != "None"), fib_any = as.integer(fibroids == 1),
                                                            bleed = as.integer(uterine_bleeding == 1), age45 = as.integer(age_anchor >= 45),
                                                            black = as.integer(race4 == "Black"))]; z })
mods <- list(fib_any = c("Fibroids", "No fibroids (adenomyosis/endometriosis only)", "Any fibroids", NA),
             bleed = c("Heavy/abnormal bleeding", "No", "Yes", NA), age45 = c("Age at uterine diagnosis", "<45", "≥45", NA),
             black = c("Race", "Not Black", "Black", NA))
e4 <- rbindlist(lapply(names(mods), function(mv) {
  md <- mods[[mv]]; covs <- setdiff(PRIMARY_COVS, md[4])
  terms <- c("anemia_any", if (mv != "black") mv, paste0("anemia_any:", mv), covs, "offset(log(py))")   # Black is already a level of race4
  fits <- lapply(imps4, function(z) glm(reformulate(terms, "ev_stroke"), poisson, data = z))
  stopifnot(!anyNA(coef(fits[[1]])))
  pf <- pool_full(fits, vc_rob)
  b <- pf$b; V <- pf$V; it <- paste0("anemia_any:", mv)
  s0 <- c(b["anemia_any"], sqrt(V["anemia_any", "anemia_any"]))
  l1 <- b["anemia_any"] + b[it]; s1 <- sqrt(V["anemia_any", "anemia_any"] + V[it, it] + 2 * V["anemia_any", it])
  pint <- 2 * pnorm(-abs(b[it] / sqrt(V[it, it])))
  ev <- imps4[[1]][, .(women = .N, strokes = sum(ev_stroke), anemic_strokes = sum(ev_stroke * anemia_any)), by = mv]
  data.table(modifier = md[1], stratum = c(md[2], md[3]),
             women = ev[match(0:1, ev[[mv]]), women], strokes = ev[match(0:1, ev[[mv]]), strokes], strokes_in_anemic = ev[match(0:1, ev[[mv]]), anemic_strokes],
             any_anemia_RR = c(fci(exp(s0[1]), exp(s0[1] - 1.96 * s0[2]), exp(s0[1] + 1.96 * s0[2])), fci(exp(l1), exp(l1 - 1.96 * s1), exp(l1 + 1.96 * s1))),
             ratio_of_RRs = c(fci(exp(b[it]), exp(b[it] - 1.96 * sqrt(V[it, it])), exp(b[it] + 1.96 * sqrt(V[it, it]))), ""),
             p_interaction = c(pint, NA))
}))
e4[!is.na(p_interaction), q_BH := bh(p_interaction)]
save_tab(e4[, .(modifier, stratum, women, strokes, strokes_in_anemic, any_anemia_RR, ratio_of_RRs, p_interaction = fp(p_interaction), q_BH = fp(q_BH))],
         "E4_effect_modification", "Any anemia and incident stroke within subgroups (Poisson, full adjustment, MI); interaction by product term")
res_put("E4_effect_modification", e4)

}
if ("E5" %in% RUN) {
# ============ E5. All-cause mortality ============================================================================
imps5 <- lapply(imps, function(z) { z <- copy(z)
  z[, t_end_d := pmin(last, ADMIN_END, na.rm = TRUE)]
  z[, dead := as.integer(!is.na(death) & death <= t_end_d)]
  z[dead == 1, t_end_d := death]
  z[, py := as.numeric(t_end_d - anchor) / 365.25]
  z[py > 0] })
nd <- sum(imps5[[1]]$dead)
covs5 <- if (nd / 32 >= 10) PRIMARY_COVS else REDUCED_COVS
m5 <- grade_models(imps5, "dead", covs5, sprintf("All-cause death (%s)", if (identical(covs5, PRIMARY_COVS)) "full adjustment" else "reduced set, EPV<10 with full set"), cox = FALSE)
r5 <- imps5[[1]][, .(women = .N, deaths = sum(dead), py = sum(py)), by = grade][order(grade)]
r5[, `:=`(rate = 1000 * deaths / py, lo = 1000 * qchisq(.025, 2 * deaths) / 2 / py, hi = 1000 * qchisq(.975, 2 * (deaths + 1)) / 2 / py)]
save_tab(r5[, .(grade = GRADE_LABELS[as.character(grade)], women, deaths, person_years = round(py), death_rate_per_1000py = fci(rate, lo, hi))],
         "E5a_mortality_rates", "All-cause death after the uterine diagnosis by anemia at diagnosis (EHR-recorded deaths)")
save_tab(m5[, .(analysis, level, estimate = fci(est, lo, hi), p = fp(p), N, deaths = events, epv = round(epv, 1))], "E5_mortality_models",
         "Anemia at the uterine diagnosis and all-cause death (Poisson, robust SE, MI)")
res_put("E5_mortality", list(models = m5, rates = r5))
}
log_msg("Extensions", paste(RUN, collapse = ","), " complete")
