# Phase 2: cohort reconstruction, outcome audit, analysis dataset for the revised (incident) analysis.
source("R/00_setup.R")
m <- readRDS(file.path(DERIVED, "master.rds"))
labs <- readRDS(file.path(DERIVED, "labs_long.rds"))
cm <- readRDS(file.path(DERIVED, "comorbid_long.rds"))

# ---- covariates (Codebook coding; 9 = unknown -> NA except smoking, kept as a level) -------------
race4 <- function(r, b) fifelse(b %in% 1, "Black", fifelse(r %in% "White", "White",
                       fifelse(grepl("^asian", tolower(r)), "Asian", "Other/unknown")))
m[, `:=`(
  race4 = factor(race4(race, black_race), levels = c("White", "Black", "Asian", "Other/unknown")),
  smoking3 = factor(c("0" = "Never", "1" = "Ever", "2" = "Ever", "9" = "Unknown")[as.character(smoking)],
                    levels = c("Never", "Ever", "Unknown")),
  migraine_any = fifelse(migraine %in% c(1, 2, 9), 1L, fifelse(migraine == 0, 0L, NA_integer_)),
  hormonal_type_cat = factor(c("0" = "None", "1" = "Combined OC", "2" = "Progestin-only", "3" = "LNG-IUD",
                               "4" = "GnRH agonist/antagonist", "5" = "Menopausal HT", "6" = "Other/multiple")[as.character(hormonal_type)],
                             levels = c("None", "Combined OC", "Progestin-only", "LNG-IUD", "GnRH agonist/antagonist",
                                        "Menopausal HT", "Other/multiple")),
  dxgrp = factor(c("1" = "Fibroids only", "2" = "Adenomyosis only", "3" = "Endometriosis only", "4" = ">1 condition")[as.character(uterine_dx_group)],
                 levels = c("Fibroids only", "Adenomyosis only", "Endometriosis only", ">1 condition")))]
for (v in c("htn", "dm", "dyslipidemia", "cad", "afib", "chf", "vte_history", "thrombophilia", "uterine_bleeding"))
  m[get(v) == 9, (v) := NA]
m[, grade_cb := factor(c("0" = "None", "1" = "Mild", "2" = "Moderate", "3" = "Severe")[as.character(anemia)], levels = GRADES)]

# ---- follow-up dates --------------------------------------------------------
enc <- rbindlist(lapply(c(6, 8, 13), function(f) {
  x <- fread(dpath("enc", sprintf("MDE_Workflow_Results_for_Encounters_%d.csv", f)), colClasses = "character")
  data.table(cn = x$`Clinic Number`, f = f, death = as.IDate(substr(x$`Death Date`, 1, 10)),
             create = as.IDate(substr(x$`Create Date`, 1, 10)))}))
fu <- enc[, .(death = if (all(is.na(death))) as.IDate(NA) else min(death, na.rm = TRUE), last = create[f == 8][1]), by = cn]
fu[last > EXTRACT_CAP, last := EXTRACT_CAP]
m <- merge(m, fu[, .(mrn = as.numeric(cn), death, last)], by = "mrn", all.x = TRUE)

# ---- outcome audit (eligible) ------------------------------------------------
so <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "stroke_outcome"))
m <- merge(m, so[, .(mrn = as.numeric(MRN), source_kinds)], by = "mrn", all.x = TRUE)
m[, acute_code := as.integer(grepl("acute", source_kinds))]
e <- m[eligible == 1]
TYPE <- c("1" = "Ischemic stroke", "2" = "Intracerebral hemorrhage", "3" = "Subarachnoid hemorrhage",
          "4" = "TIA", "5" = "Cerebral venous thrombosis", "9" = "Stroke, unknown type")
TIMING <- c("1" = "Before index", "2" = "Same encounter", "3" = "After index", "9" = "Unknown (code 9)")
s <- e[stroke_any == 1]
s[, `:=`(type = TYPE[as.character(stroke_type)], timing = fifelse(is.na(stroke_timing), "Blank (undated)", TIMING[as.character(stroke_timing)]))]
out1 <- dcast(s[, .N, by = .(type, timing)], type ~ timing, value.var = "N", fill = 0)
out1[, Total := rowSums(.SD), .SDcols = -1]
out1 <- rbind(out1, cbind(data.table(type = "All"), out1[, lapply(.SD, sum), .SDcols = -1]))
save_tab(out1, "S4a_outcome_type_by_timing", "Stroke/TIA in eligible women (unique women): type by timing relative to index")
out2 <- s[, .N, by = .(timing, code_kinds = fifelse(is.na(source_kinds), "not in stroke_outcome sheet (imaging-derived)", source_kinds))][order(timing, -N)]
save_tab(out2, "S4b_outcome_code_kinds", "Kinds of stroke codes behind each woman's event (acute / history / sequela)")
# multiple event types recorded in notes
s[, multi := regmatches(notes, regexpr("multiple event types coded: [^;|]+", notes))[1], by = study_id]
s[, multi := fifelse(grepl("multiple event types", notes), sub(".*multiple event types coded: ([^;|]+).*", "\\1", notes), NA_character_)]
out3 <- s[!is.na(multi), .N, by = .(primary_type = type, all_types_coded = trimws(multi))][order(-N)]
save_tab(out3, "S4c_multiple_types", "Women with more than one stroke type coded (from notes); stroke_type holds the most specific")
res_put("outcome_audit", list(eligible_strokes = nrow(s), by_type = as.list(table(s$type)),
                              by_timing = as.list(table(s$timing)), multiple_types = s[!is.na(multi), .N],
                              dated = s[!is.na(stroke_date), .N], undated = s[is.na(stroke_date), .N],
                              after_index_acute_code = s[stroke_timing %in% 3, sum(acute_code)],
                              after_index_no_acute_code = s[stroke_timing %in% 3, sum(acute_code == 0)]))
chk_last <- e[!is.na(stroke_date) & stroke_date > last, .N]
log_msg("Follow-up", sprintf("last contact earlier than stroke date for %d eligible women; deaths recorded %d", chk_last, e[!is.na(death), .N]))

# ---- baseline laboratory exposure -------------------------------------------
hb <- labs$hb
base_hb <- function(src_ok, window) {
  x <- merge(hb[src %in% src_ok], e[, .(study_id, index_date)], by = "study_id")[date <= index_date & date >= index_date - window]
  x <- x[, .(val = median(val)), by = .(study_id, date)]       # CBC and POC on the same day -> median
  setorder(x, study_id, -date)
  x[, .SD[1], by = study_id][, .(study_id, hb = val, hb_date = date)]
}
b <- base_hb("CBC", 1095)
b1 <- base_hb("CBC", 365); setnames(b1, c("study_id", "hb_1y", "hb_1y_date"))
bp <- base_hb(c("CBC", "POC"), 1095); setnames(bp, c("study_id", "hb_poc", "hb_poc_date"))
ntest <- merge(hb[src == "CBC"], e[, .(study_id, index_date)], by = "study_id")[date <= index_date & date >= index_date - 730, .(hb_tests_2y = .N), by = study_id]
e <- Reduce(function(a, z) merge(a, z, by = "study_id", all.x = TRUE), list(e, b, b1, bp, ntest))
e[is.na(hb_tests_2y), hb_tests_2y := 0L]
grade <- function(x) factor(fifelse(x >= 12, "None", fifelse(x >= 10, "Mild", fifelse(x >= 8, "Moderate", "Severe"))), levels = GRADES)
e[, `:=`(grade = grade(hb), grade_1y = grade(hb_1y), grade_poc = grade(hb_poc))]
# red-cell indices on the CBC nearest the baseline Hb, within 30 d and not after index
near <- function(tab, nm) {
  x <- merge(tab, e[!is.na(hb_date), .(study_id, hb_date, index_date)], by = "study_id")[abs(date - hb_date) <= 30 & date <= index_date]
  x[, gap := abs(as.integer(date - hb_date))]; setorder(x, study_id, gap)
  setnames(x[, .SD[1], by = study_id][, .(study_id, val)], "val", nm)
}
for (z in list(list(labs$mcv, "mcv_lab"), list(labs$rdw, "rdw_lab"), list(labs$plt, "plt_lab")))
  e <- merge(e, near(z[[1]], z[[2]]), by = "study_id", all.x = TRUE)

# ---- baseline comorbidities (dated on or before index; recent = within 365 d before index) -------------
cmi <- merge(cm, e[, .(study_id, index_date)], by = "study_id")
flag_before <- function(flags, recent = FALSE, any_date = FALSE) {
  x <- cmi[flag %in% flags]
  if (!any_date) x <- x[date <= index_date]
  if (recent) x <- x[date >= index_date - 365]
  e$study_id %in% unique(x$study_id)
}
e[, `:=`(
  haemoglobinopathy = as.integer(flag_before(c("sickle_disease", "sickle_trait", "thal_minor", "thal_other", "other_haemoglobinopathy"), any_date = TRUE)),
  sickle_disease = as.integer(flag_before("sickle_disease", any_date = TRUE)),
  ckd_pre = as.integer(flag_before(c("ckd", "ckd_esrd"))),
  liver_pre = as.integer(flag_before(c("cirrhosis", "chronic_liver"))),
  alcohol_pre = as.integer(flag_before("alcohol_disorder")),
  ibd_pre = as.integer(flag_before("ibd")),
  hiv_pre = as.integer(flag_before("hiv")),
  malabs_pre = as.integer(flag_before("malabsorption")),
  gibleed_1y = as.integer(flag_before("gi_bleed", recent = TRUE)),
  preg_1y = as.integer(flag_before("pregnancy", recent = TRUE)),
  neoplastic_anaemia_pre = as.integer(flag_before("anaemia_neoplastic_chemo")),
  coded_ida_pre = as.integer(flag_before("anaemia_iron", recent = TRUE)),
  coded_bloodloss_pre = as.integer(flag_before("anaemia_blood_loss_acute", recent = TRUE)))]
e[, anaemia_cause_any := as.integer((haemoglobinopathy + ckd_pre + liver_pre + alcohol_pre + ibd_pre + hiv_pre + malabs_pre +
                                     gibleed_1y + preg_1y + neoplastic_anaemia_pre) > 0)]

# ---- incident cohort: flow ----------------------------------------------------------------
e[, prevalent := as.integer(stroke_any == 1 & stroke_timing %in% c(1, 2))]
e[, timing_unknown := as.integer(stroke_any == 1 & !(stroke_timing %in% c(1, 2, 3)))]
e[, fu_end0 := pmin(last, death, ADMIN_END, na.rm = TRUE)]
e[, ev_date := fifelse(stroke_any == 1 & stroke_timing %in% 3, stroke_date, as.IDate(NA))]
e[!is.na(ev_date) & ev_date <= ADMIN_END, fu_end0 := pmin(fu_end0, ev_date)]   # stroke ends follow-up for every outcome
flow <- data.table(step = "Records in stroke_analysis_master.csv", n = nrow(m))
for (r in list(c(1, "Excluded: gynecologic or other active malignancy (±6 months of index)"),
               c(6, "Excluded: age >60 at index"), c(3, "Excluded: age <18 at index"), c(5, "Excluded: age not derivable")))
  flow <- rbind(flow, data.table(step = r[2], n = m[exclusion_reason == as.numeric(r[1]), .N]))
flow <- rbind(flow, data.table(step = "Eligible women", n = nrow(e)))
c1 <- e[prevalent == 0 & timing_unknown == 0]
flow <- rbind(flow,
  data.table(step = "Excluded: stroke/TIA before index or at the index encounter", n = e[prevalent == 1, .N]),
  data.table(step = "Excluded: stroke/TIA with unknown timing (undated)", n = e[timing_unknown == 1, .N]),
  data.table(step = "Free of recorded stroke/TIA at index", n = nrow(c1)))
c2 <- c1[!is.na(hb)]
flow <- rbind(flow, data.table(step = "Excluded: no laboratory CBC Hb in the 3 years before or on the index date", n = nrow(c1) - nrow(c2)))
c3 <- c2[!is.na(fu_end0) & fu_end0 > index_date]
flow <- rbind(flow, data.table(step = "Excluded: no follow-up time after the index date", n = nrow(c2) - nrow(c3)),
                    data.table(step = "Incident-stroke cohort (revised primary analysis)", n = nrow(c3)))
save_tab(flow, "Fig1_flow_counts", "Cohort flow for the revised primary analysis")

ic <- c3
ic[, `:=`(t0 = index_date, t1 = fu_end0, py = as.numeric(fu_end0 - index_date) / 365.25)]
inc <- ic$stroke_any == 1 & ic$stroke_timing %in% 3 & !is.na(ic$ev_date) & ic$ev_date <= ic$t1
ic[, `:=`(
  ev_stroke = as.integer(inc & stroke_type != 4),
  ev_any = as.integer(inc),
  ev_isch = as.integer(inc & stroke_type == 1),
  ev_tia = as.integer(inc & stroke_type == 4),
  ev_ich = as.integer(inc & stroke_type == 2), ev_sah = as.integer(inc & stroke_type == 3),
  ev_cvt = as.integer(inc & stroke_type == 5), ev_unk = as.integer(inc & stroke_type == 9))]
ic[, died := as.integer(ev_any == 0 & !is.na(death) & death <= t1)]
ic[, ev_stroke_acute := as.integer(ev_stroke == 1 & acute_code == 1)]
lost <- e[stroke_any == 1 & stroke_timing %in% 3 & !is.na(hb) & study_id %in% c2$study_id & !(study_id %in% ic$study_id), .N]
log_msg("Cohort", sprintf("incident cohort %s women, %.0f person-years; strokes (excl. TIA) %d, stroke/TIA %d; deaths without stroke %d; after-index strokes lost because event on/after follow-up end check: %d",
                          format(nrow(ic), big.mark = ","), sum(ic$py), sum(ic$ev_stroke), sum(ic$ev_any), sum(ic$died), lost))
res_put("cohort", list(flow = setNames(as.list(flow$n), flow$step), person_years = sum(ic$py),
                       events = list(stroke = sum(ic$ev_stroke), stroke_or_tia = sum(ic$ev_any), ischemic = sum(ic$ev_isch),
                                     tia = sum(ic$ev_tia), ich = sum(ic$ev_ich), sah = sum(ic$ev_sah), cvt = sum(ic$ev_cvt),
                                     unknown = sum(ic$ev_unk), deaths_no_stroke = sum(ic$died)),
                       median_followup_years = median(ic$py)))
saveRDS(e, file.path(DERIVED, "eligible.rds"))
saveRDS(ic, file.path(DERIVED, "incident_cohort.rds"))
saveRDS(m, file.path(DERIVED, "master_derived.rds"))
