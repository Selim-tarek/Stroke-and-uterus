# Part 4: verify the drift claims and run one targeted test that can separate a calendar-PERIOD effect (coding / system / practice change that
# affects everyone after 2023) from a diagnosis-COHORT effect (women diagnosed from 2023 differ). Uses the pre-specified C4 predictors only.
source("final_verification/R/V00_setup.R")
md_open("V03_drift.R", "V03 Temporal calibration drift: verification and period-versus-cohort test")
b <- load_base(); b[, ay := year(anchor)]; CUT <- as.IDate("2023-01-01")
b[, exit_date := anchor + round(t * 365.25)]
# ---- 1. model developed with ADMINISTRATIVE CENSORING at 2022-12-31 (no calendar-2023+ outcome information used) ------------------
dv <- b[anchor < CUT]; dv[, `:=`(t0 = t, y0 = y)]
dv[exit_date >= CUT, `:=`(t = as.numeric(CUT - 1 - anchor) / 365.25, y = 0L)]
dv <- dv[t > 0]
f <- coxph(fml_of(CAND$C4_reduced_baseline), data = dv, x = TRUE, model = TRUE)
bh <- basehaz(f, centered = TRUE); Hf <- function(s) approx(bh$time, bh$hazard, xout = s, method = "constant", rule = 2, f = 0)$y
aud_log(sprintf("C4 developed on anchors < 2023 with follow-up censored at 2022-12-31: %d women, %d events, max follow-up %.2f y", nrow(dv), sum(dv$y), max(dv$t)))
# ---- 2. O/E in calendar 2023+ follow-up, by diagnosis cohort (left-truncated at the time since anchor on 2023-01-01) -----------------
win <- function(x) { x <- copy(x); x[, `:=`(t_in = pmax(0, as.numeric(CUT - anchor) / 365.25), t_out = t)]; x[t_out > t_in] }
oe_win <- function(x, tmax = Inf) { x <- copy(x); x[, t_cap := pmin(t_out, tmax)]; x <- x[t_cap > t_in]; x[, ev := as.integer(y == 1 & t <= t_cap)]
  # baseline hazard beyond the development follow-up range is held at its last value (no extrapolation of the hazard increase)
  E <- sum((Hf(x$t_cap) - Hf(x$t_in)) * exp(lpof(f, x))); O <- sum(x$ev); ci <- c(qchisq(.025, 2 * O) / 2, qchisq(.975, 2 * O + 2) / 2)
  data.table(women = nrow(x), observed = O, expected = round(E, 1), oe = sprintf("%.2f (%.2f–%.2f)", O / E, ci[1] / E, ci[2] / E)) }
old <- win(b[anchor < CUT]); new <- win(b[anchor >= CUT])
pc <- rbind(cbind(group = "anchored before 2023: calendar-2023+ follow-up within 0–3 y of the anchor", oe_win(old, tmax = 3)),
            cbind(group = "anchored 2023+: follow-up within 0–3 y of the anchor", oe_win(new, tmax = 3)))
aud_tab(pc, "V03a_period_vs_cohort", "Observed/expected strokes (acute-coded) in calendar-2023+ follow-up under C4 developed only on pre-2023 calendar time. Expected = Σ[H0(exit) − H0(entry)]·exp(lp), left-truncated. If O/E is raised in BOTH diagnosis cohorts, the drift is a calendar-period phenomenon rather than a property of women diagnosed recently.")
# ---- 3. Poisson model: calendar period vs diagnosis cohort, adjusted for time since anchor and attained age ---------------------------
sp <- as.data.table(survSplit(Surv(t, y) ~ ., data = b[, .(study_id, t, y, anchor, age, hb_def, htn, cad, ckd_pre, ay)], cut = c(30 / 365.25, 1, 3), episode = "tsa"))
sp[, `:=`(start_date = anchor + round(tstart * 365.25), stop_date = anchor + round(t * 365.25))]
# split each interval at 2023-01-01
pre <- copy(sp)[start_date < CUT][, `:=`(t_end = pmin(t, as.numeric(CUT - anchor) / 365.25))][, `:=`(y = as.integer(y == 1 & t <= t_end), t = t_end, period = "calendar < 2023")]
post <- copy(sp)[stop_date >= CUT][, `:=`(tstart = pmax(tstart, as.numeric(CUT - anchor) / 365.25), period = "calendar ≥ 2023")]
lx <- rbind(pre, post, fill = TRUE)[t > tstart]
lx[, `:=`(py = t - tstart, cohort = fifelse(ay < 2023, "anchored < 2023", "anchored ≥ 2023"), tsa = factor(tsa), att_age = age + tstart)]
g <- glm(y ~ period + cohort + tsa + att_age + hb_def + htn + cad + ckd_pre + offset(log(py)), family = poisson, data = lx)
rr <- function(k) sprintf("%.2f (%.2f–%.2f)", exp(coef(g)[k]), exp(confint.default(g)[k, 1]), exp(confint.default(g)[k, 2]))
aud_tab(data.table(term = c("calendar ≥ 2023 vs < 2023 (period)", "anchored ≥ 2023 vs < 2023 (cohort, within the same calendar period)"),
                   rate_ratio = c(rr("periodcalendar ≥ 2023"), rr("cohortanchored ≥ 2023"))),
        "V03b_poisson_period_cohort", "Poisson rate ratios (acute-coded stroke), mutually adjusted, plus time since anchor (0–30 d, 31 d–1 y, 1–3 y, >3 y), attained age and the C4 predictors. Period and cohort are only partly separable: all ≥2023-anchored person-time is in calendar ≥2023.")
rates <- lx[, .(events = sum(y), person_years = round(sum(py))), by = .(cohort, period)][order(cohort, period)]
rates[, rate_per_1000py := mapply(pois_ci, events, person_years)][, events := sc(events)]
aud_tab(rates, "V03c_rates_period_cohort", "Crude rates by diagnosis cohort and calendar period")
# ---- 4. outcome composition by calendar period of the event, and predictor-coding trends --------------------------------------------
e <- b[y_all == 1][, cal := fifelse(stroke_date < CUT, "event in calendar < 2023", "event in calendar ≥ 2023")]
aud_tab(e[, .(events_o1 = .N, pct_acute_code = pct(mean(acute & !imaging_only)), pct_corroborated = pct(mean(corroborated)), pct_ischemic = pct(mean(stroke_type == 1))), by = cal],
        "V03d_outcome_by_calendar", "Outcome composition by calendar period of the event (all O1 events, baseline cohort)")
tr <- rbindlist(lapply(c("htn", "dm", "migraine_any", "cad", "ckd_pre"), function(v) { gg <- glm(as.formula(paste(v, "~ ay + age + fibroids + endometriosis")), family = binomial, data = b)
  data.table(predictor = v, or_per_anchor_year = sprintf("%.3f (%.3f–%.3f)", exp(coef(gg)["ay"]), exp(confint.default(gg)["ay", 1]), exp(confint.default(gg)["ay", 2]))) }))
aud_tab(tr, "V03e_coding_trends", "Change in coded prevalence per anchor year, adjusted for age, fibroids and endometriosis (logistic)")
md_note("Not shown: O/E over ALL calendar-2023+ follow-up of pre-2023 anchors. Much of that follow-up lies more than 5 years after the anchor, beyond the development model's follow-up (censored at 2022-12-31, maximum 5.0 years). Expected counts there would rely on a baseline hazard held flat beyond its estimable range (an artefact, not evidence). The comparison above uses the same 0–3-year window since anchor for both cohorts.")
md_note("CKD is the only comorbidity whose date was verified from raw extracts (dx files 22–31, first code on or before the anchor). It shows no trend by anchor year, whereas hypertension, CAD and diabetes (master-file flags) fall 6–8% per year after adjustment. This is compatible with those flags partly reflecting codes recorded after the anchor, which accumulate more in older cohorts. It is also compatible with real case-mix change. These data cannot distinguish the two.")
res_put("V03", list(period_cohort = pc))
