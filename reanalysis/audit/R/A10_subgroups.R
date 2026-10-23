# Audit Q18: performance by clinically relevant subgroup (discrimination and calibration with uncertainty), and
# differential measurement/ascertainment that could produce unequal calibration.
source("audit/R/00_audit_setup.R"); source("audit/R/val_fns.R")
md_open("A10_subgroups.R", "A10 Subgroup performance and fairness")
d <- load_dev(); f9 <- fitc(d, PRED9); CUT <- as.IDate("2023-01-01")
d[, lp := lpof(f9, d)]
d[, `:=`(age_grp = cut(age, c(18, 35, 45, 61), c("18–34", "35–44", "45–60"), right = FALSE),
         hb_grp = cut(hb, c(0, 10, 12, 99), c("<10", "10–11.9", "≥12"), right = FALSE),
         race = as.character(race4), fibroid_any = fifelse(fibroids == 1, "any fibroids", "no fibroids"))]
SUBS <- list(age_grp = "age at anchor", dxgrp = "uterine diagnosis", fibroid_any = "fibroids", hb_grp = "hemoglobin, g/dL", htn = "hypertension",
             dm = "diabetes", estrogen = "estrogen-containing therapy", migraine_any = "migraine", race = "race (not a predictor)")
dv <- d[anchor < CUT]; vv <- d[anchor >= CUT]; fv <- fitc(dv, PRED9); vv[, lpv := lpof(fv, vv)]
one <- function(x, lpcol, f, h, shrink) { cc <- if (sum(x$y & x$t <= h) >= 5) cstat(x[[lpcol]], x, h) else c(NA, NA); o <- oe_counts(f, x, h, shrink)
  list(women = nrow(x), strokes = o[["O"]], c = if (is.na(cc[1])) "<5 events" else sprintf("%.2f (%.2f–%.2f)", cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2]),
       oe = fmt_ci(o[c("oe", "lo", "hi")], 2), flag = if (o[["O"]] < 20) "fewer than 20 events: imprecise" else "") }
tab <- rbindlist(lapply(names(SUBS), function(v) rbindlist(lapply(sort(unique(d[[v]])), function(l) {
  a <- one(d[get(v) == l], "lp", f9, 5, SPEC9$shrink); b <- one(vv[get(v) == l], "lpv", fv, 3, 1)
  data.table(subgroup = SUBS[[v]], level = as.character(l), dev_women = a$women, dev_strokes_5y = a$strokes, dev_c_5y = a$c, dev_oe_5y = a$oe,
             temporal_women = b$women, temporal_strokes_3y = b$strokes, temporal_c_3y = b$c, temporal_oe_3y = b$oe, note = trimws(paste(a$flag, if (nzchar(b$flag)) "(temporal)" else ""))) }))))
aud_tab(tab, "A10a_subgroup_performance", "Development: published nine-item model, apparent, 5-year O/E by expected counts (shrinkage applied). Temporal: model refitted on < 2023, evaluated on ≥ 2023 at 3 years. C within a subgroup defined by a predictor reflects only the other predictors.")
# differential measurement and ascertainment
co <- add_vars(load_cohort())
dm_ <- rbindlist(lapply(c("White", "Black", "Asian", "Other/unknown"), function(r) { x <- co[race4 == r]; e <- d[race4 == r & y == 1]
  data.table(race = r, anchored_women = nrow(x), pct_with_cbc_hb = pct(mean(!is.na(x$hb))), pct_smoking_documented = pct(mean(x$smoking3 != "Unknown")),
             pct_bmi_recorded = pct(mean(!is.na(x$bmi))), dev_strokes = nrow(e), pct_strokes_acute_code = pct(mean(grepl("acute", e$source_kinds))),
             pct_strokes_imaging_positive = pct(mean(e$stroke_confirmed_imaging %in% 1)), pct_hemoglobinopathy = pct(mean(x$haemoglobinopathy))) }))
aud_tab(dm_, "A10b_differential_measurement", "Measurement and outcome-ascertainment indicators by race (anchored cohort; strokes from the development sample)")
md_note("Hemoglobinopathy is coded in ", dm_[race == "Black", pct_hemoglobinopathy], "% of Black and ", dm_[race == "White", pct_hemoglobinopathy],
        "% of White women. Sickle cell trait alone is generally reported in roughly 7–9% of Black Americans (external figure, not from these data), so the flag almost certainly captures only diagnosed, clinically documented disease and under-ascertains trait; its meaning differs from 'carrier status'. ",
        "The race likelihood-ratio check in R/15 (P = 0.877) tests whether race adds discrimination, not whether the model is calibrated within each group.")
res_put("A10", list(subgroups = tab))
