# Audit Q3: cohort reconstruction, target population, time zero, censoring.
# Re-derives the anchor independently from the master file and reconciles every exclusion against
# output/tables/A3_anchored_flow.csv and the development sample of R/16_calculator_model.R.
source("audit/R/00_audit_setup.R")
md_open("A01_cohort_design.R", "A01 Cohort reconstruction and study design")
e  <- readRDS(file.path(DERIVED, "eligible.rds"))
mm <- readRDS(file.path(DERIVED, "master_derived.rds"))[eligible == 1, .(study_id, surgical_tx, surgical_tx_date, myomectomy_date, notes2 = notes)]
e  <- merge(e, mm, by = "study_id", suffixes = c("", ".m"))
co <- load_cohort(); d <- load_dev()

# ---- 1. independent anchor ------------------------------------------------------------------------------
PROC_RX <- "(myomectomy|uterine artery embolisation|hysterectomy|endometrial ablation) \\d{4}-\\d{2}-\\d{2}"
pr <- rbind(e[surgical_tx %in% 1:4 & !is.na(surgical_tx_date), .(study_id, date = surgical_tx_date)],
            e[!is.na(myomectomy_date), .(study_id, date = myomectomy_date)],
            e[!is.na(notes2), .(date = as.IDate(sub(".* ", "", regmatches(notes2, gregexpr(PROC_RX, notes2))[[1]]))), by = study_id])
pr <- pr[!is.na(date) & year(date) >= 1980, .(proc = min(date)), by = study_id]
e <- merge(e, pr, by = "study_id", all.x = TRUE)
e[, anchor2 := index_date]
e[!is.na(proc) & proc < index_date & (age_index - as.numeric(index_date - proc) / 365.25) >= 18, anchor2 := proc]
chk <- merge(co[, .(study_id, anchor)], e[, .(study_id, anchor2)], by = "study_id")
aud_log(sprintf("independent anchor matches saved anchor for %d of %d women", chk[anchor == anchor2, .N], nrow(chk)))

# ---- 2. flow with every exclusion broken down ------------------------------------------------------------
e[, prior := stroke_any == 1 & (is.na(stroke_date) | stroke_date <= anchor2)]
e[, why_prior := fifelse(!prior, NA_character_, fifelse(is.na(stroke_date), "undated stroke/TIA",
                  fifelse(stroke_type == 4, "TIA on or before anchor", "stroke on or before anchor")))]
s1 <- e[prior == FALSE]; s1h <- s1[!is.na(hormonal_type_cat)]
s1h[, t_end := pmin(last, death, ADMIN_END, na.rm = TRUE)]; s1h[stroke_any == 1, t_end := pmin(t_end, stroke_date)]
s1h[, nofu := is.na(t_end) | t_end <= anchor2]
s1h[, why_nofu := fifelse(!nofu, NA_character_, fifelse(is.na(t_end), "no last-contact date", fifelse(!is.na(death) & death <= anchor2, "death on/before anchor",
                    fifelse(t_end == anchor2, "last contact on the anchor date", "last contact before the anchor"))))]
s2 <- s1h[nofu == FALSE]
s3 <- co[!is.na(grade)]; s4 <- s3[!is.na(bmi)]; s5 <- s4[age_anchor >= 18]
flow <- rbindlist(list(
  data.table(step = "Eligible women (age 18–60 at first 2018+ diagnosis, no active malignancy)", n = nrow(e), reported = 39807),
  e[prior == TRUE, .(n = .N, reported = NA_real_), by = .(step = paste0("  excluded: ", why_prior))],
  data.table(step = "  excluded: stroke/TIA on/before anchor or undated (total)", n = e[prior == TRUE, .N], reported = 983),
  data.table(step = "  excluded: hormonal_type = 9 (unknown) — not shown in A3", n = s1[is.na(hormonal_type_cat), .N], reported = NA),
  s1h[nofu == TRUE, .(n = .N, reported = NA_real_), by = .(step = paste0("  excluded: ", why_nofu))],
  data.table(step = "  excluded: no follow-up after anchor (total)", n = s1h[nofu == TRUE, .N], reported = 491),
  data.table(step = "Anchored cohort, stroke-free at anchor", n = nrow(s2), reported = 38332),
  data.table(step = "  excluded: no CBC Hb in the 3 years up to the anchor", n = nrow(s2) - nrow(s3), reported = NA),
  data.table(step = "With CBC Hb", n = nrow(s3), reported = 24886),
  data.table(step = "  excluded: BMI missing", n = nrow(s3) - nrow(s4), reported = NA),
  data.table(step = "  excluded: age <18 at anchor", n = nrow(s4) - nrow(s5), reported = NA),
  data.table(step = "Development sample (R/16, and R/15 without the age filter)", n = nrow(s5), reported = 24006),
  data.table(step = "  strokes (excluding TIA) in development sample", n = sum(d$y), reported = 345)))
flow[, agrees := fifelse(is.na(reported), "", fifelse(n == reported, "yes", "NO"))]
aud_tab(flow, "A01a_cohort_flow", "Cohort flow re-derived independently; `reported` = value in output/tables/A3_anchored_flow.csv, R7 or the methods summary.")
md_note("Same anchor as the saved cohort for ", chk[anchor == anchor2, .N], " of ", nrow(chk), " women. ",
        "R/15 (14-predictor model) has no `age_anchor >= 18` filter; its sample is ", co[!is.na(grade) & !is.na(bmi), .N], " women, identical to R/16.")

# ---- 3. target population: diagnoses and overlap ----------------------------------------------------------
dx <- function(z, lab) z[, .(sample = lab, women = .N, fibroids = sum(fibroids == 1), adenomyosis = sum(adenomyosis == 1), endometriosis = sum(endometriosis == 1),
  fibroids_only = sum(dxgrp == "Fibroids only"), adenomyosis_only = sum(dxgrp == "Adenomyosis only"), endometriosis_only = sum(dxgrp == "Endometriosis only"),
  more_than_one = sum(dxgrp == ">1 condition"), strokes = sum(ev_stroke %in% 1))]
dxt <- rbind(dx(co, "Anchored cohort"), dx(d, "Development sample"), dx(d[y == 1], "Development sample, women with stroke"))
aud_tab(dxt, "A01b_diagnoses", "Uterine diagnoses (Codebook flags; any date in the 2018+ extract). Overlapping diagnoses are counted in each column.")
ov <- d[, .(women = .N, strokes = sum(y), person_years = round(sum(t)), rate_per_1000py = sprintf("%.2f", 1000 * sum(y) / sum(t))), by = .(fibroids, adenomyosis, endometriosis)][order(-women)]
aud_tab(ov, "A01c_diagnosis_overlap", "Development sample by combination of the three diagnoses")

# ---- 4. time zero ------------------------------------------------------------------------------------------
tz <- data.table(item = c("Development women with anchor = first coded diagnosis", "  of which first coded diagnosis in Jan–Mar 2018 (first quarter of the extract window)",
                          "Development women with anchor moved earlier by a procedure", "  strokes between the procedure anchor and the first coded diagnosis",
                          "  person-years between procedure anchor and first coded diagnosis", "Anchor before 2018"),
                 n = c(d[anchor == index_date, .N], d[anchor == index_date & index_date < as.IDate("2018-04-01"), .N], d[anchor < index_date, .N],
                       d[anchor < index_date & y == 1 & stroke_date <= index_date, .N], round(d[anchor < index_date, sum(as.numeric(pmin(index_date, t_end) - anchor)) / 365.25]),
                       d[anchor < as.IDate("2018-01-01"), .N]))
mon <- d[anchor == index_date, .N, by = .(m = format(index_date, "%Y-%m"))][order(m)]
tz <- rbind(tz, data.table(item = "  median women per calendar month of first coded diagnosis, 2018-04 to 2026-08", n = median(mon[m >= "2018-04" & m <= "2026-08", N])))
aud_tab(tz, "A01d_time_zero", "Time zero. The diagnosis extract window begins 2018-01-01; a pile-up of first diagnoses in early 2018 would indicate left truncation of prevalent disease. Compare the Jan–Mar 2018 count with 3 × the monthly median.")

# ---- 5. censoring ------------------------------------------------------------------------------------------
d[, end_reason := fifelse(ev_any == 1, fifelse(y == 1, "stroke (event)", "TIA (censored for stroke)"), fifelse(died == 1, "death",
                  fifelse(t_end == ADMIN_END, "administrative end 2026-09-25", "last recorded encounter")))]
cen <- d[, .(women = .N, median_years = round(median(t), 2)), by = end_reason][order(-women)]
aud_tab(cen, "A01e_end_of_follow_up", "Why follow-up ended (development sample)")
ck <- data.table(check = c("follow-up end later than 2026-09-25", "event date after the end of follow-up", "death date before the anchor",
                           "stroke date after last recorded encounter (event not counted)", "stroke date after death date",
                           "women whose last encounter is >1 year before 2026-09-25 (possible loss to follow-up)"),
                 n = c(d[t_end > ADMIN_END, .N], d[y == 1 & stroke_date > t_end, .N], d[!is.na(death) & death < anchor, .N],
                       co[!is.na(grade) & !is.na(bmi) & stroke_any == 1 & stroke_type != 4 & stroke_date > last, .N],
                       co[stroke_any == 1 & !is.na(death) & stroke_date > death, .N],
                       d[ev_any == 0 & died == 0 & t_end < ADMIN_END - 365, .N]))
aud_tab(ck, "A01f_censoring_checks", "Logical checks on follow-up dates")
res_put("A01", list(flow = flow, anchor_match = chk[anchor == anchor2, .N]))
