# Phase 1b: re-derive laboratory values from the raw dated extracts.
# Explicit inclusion rules per analyte; every dropped row is counted by reason.
source("R/00_setup.R")
m <- readRDS(file.path(DERIVED, "master.rds"))
link <- m[, .(cn = as.character(mrn), study_id)]

read_lab <- function(f) {
  x <- fread(dpath("labs", sprintf("MDE Workflow Results for Lab Tests (%d).csv", f)), colClasses = list(character = "Clinic Number"))
  setnames(x, c("Clinic Number", "Value (Numeric)", "Result Date", "LAB_UNIT_OF_MEASURE_TXT"), c("cn", "val", "rdate", "unit"))
  if (!"Lab Test Description" %in% names(x)) x[, `Lab Test Description` := NA_character_]
  if (!"Lab Test Code" %in% names(x)) x[, `Lab Test Code` := NA_character_]
  setnames(x, c("Lab Test Description", "Lab Test Code"), c("desc", "code"))
  x[, `:=`(unit = toupper(trimws(unit)), desc_n = toupper(gsub("\\.+", "", trimws(desc))),
           date = as.IDate(substr(rdate, 1, 10)), code = as.character(code))]
  x[, .(cn, desc, desc_n, code, unit, val, date)]
}

steps <- list()
clean <- function(x, analyte, keep_desc, unit_ok, lo, hi) {
  keep_desc[is.na(keep_desc)] <- FALSE; unit_ok[is.na(unit_ok)] <- FALSE
  s <- data.table(analyte = analyte, step = "raw rows", rows = nrow(x))
  x <- x[keep_desc]; s <- rbind(s, data.table(analyte, "test description/type matches inclusion rule", nrow(x)), use.names = FALSE)
  x <- x[unit_ok[keep_desc]]; s <- rbind(s, data.table(analyte, "unit accepted", nrow(x)), use.names = FALSE)
  x <- x[!is.na(val) & val != 9999999]; s <- rbind(s, data.table(analyte, "numeric value, not sentinel 9999999", nrow(x)), use.names = FALSE)
  x <- x[val >= lo & val <= hi]; s <- rbind(s, data.table(analyte, sprintf("plausible range %s–%s", lo, hi), nrow(x)), use.names = FALSE)
  x <- x[!is.na(date) & year(date) >= 1990 & date <= EXTRACT_CAP]; s <- rbind(s, data.table(analyte, "valid date (1990 to 2026-09-27)", nrow(x)), use.names = FALSE)
  x <- merge(x, link, by = "cn"); s <- rbind(s, data.table(analyte, "linked to a cohort woman", nrow(x)), use.names = FALSE)
  steps[[analyte]] <<- s
  x
}
daily <- function(x, by_extra = NULL) x[, .(val = median(val), n_same_day = .N), by = c("study_id", "date", by_extra)]

# ---- Hemoglobin (Lab Tests 2) ----------------------------------------------
h <- read_lab(2)
d <- h$desc_n
not_blood_total <- grepl("QL|, U$|^U |URINE|UPHB|A1C|RETIC|PLASMA|, P$|^HEMOGLOBIN (A|A2|F)|HB A$|HGB A$|HEMOGLOBIN A,|VARIANT|O2|OXY|METHEM|CARBOXY|COMMENT|SOURCE|^CBC|DIFFERENTIAL|HBC TOTAL", d)
gas <- grepl("ARTERIAL|VENOUS|GAS|ABG|VBG|PABG|^T ?T?HB|^T HGB|TOTAL HGB|HGB TOTAL|^THB", d)
poc <- grepl("POC", d)
cbc_name <- grepl("^(HEMOGLOBIN|HGB|HB)(, B|\\(B\\)| BLD-MCNC| ?MHS)?$|^EXTI (NON-GHS )?HEMOGLOBIN(, B)? - |^EXTM HEMOGLOBIN$", d)
h[, hb_class := fifelse(not_blood_total, "excluded: not whole-blood total Hb",
                fifelse(gas, "excluded: blood gas / co-oximetry",
                fifelse(poc, "POC", fifelse(cbc_name, "CBC", "excluded: unrecognised description"))))]
desc_tab <- h[, .(rows = .N, median_value = median(val[val < 1e6], na.rm = TRUE)), by = .(description = desc, class = hb_class)][order(class, -rows)]
save_tab(desc_tab, "S3c_hb_descriptions", "Every Hb test description in Lab Tests 2 and its classification")
g_dl <- h$unit %in% c("G/DL", "GDL", "GM/DL", "GMDL", "G(HB)/DL", "")
hb_cbc <- clean(h, "Hb (CBC, primary)", h$hb_class == "CBC", g_dl, 3, 20)
hb_poc <- clean(h, "Hb (POC, sensitivity S4)", h$hb_class == "POC", g_dl, 3, 20)
hb <- rbind(daily(hb_cbc)[, src := "CBC"], daily(hb_poc)[, src := "POC"])
rm(h, hb_cbc, hb_poc); gc()

# ---- MCV (Lab Tests 3) ------------------------------------------------------
x <- read_lab(3)
mcv <- daily(clean(x, "MCV", grepl("MEAN CORPUSCULAR VOLUME|^MCV$", x$desc_n), x$unit %in% c("FL", ""), 50, 130))
# ---- MCH / MCHC (Lab Tests 4) -------------------------------------------------
x <- read_lab(4)
mchc <- daily(clean(x, "MCHC", grepl("MCHC", x$desc_n), x$unit %in% c("G/DL", "GM/DL", "GDL", ""), 20, 45))
mch <- daily(clean(x, "MCH", grepl("^MCH$|MCH \\(MEAN CORPUSCULAR HEMOGLOBIN\\)", x$desc_n) & !grepl("MCHC", x$desc_n),
                   x$unit %in% c("PG", ""), 10, 60))
# ---- RDW (Lab Tests 7; 6 checked as duplicate) --------------------------------
x <- read_lab(7)
x6 <- read_lab(6)
k <- c("cn", "code", "date", "val")
dup6 <- nrow(fintersect(unique(x6[, ..k]), unique(x[, ..k])))
log_msg("Labs", sprintf("Lab Tests 6: %s of %s distinct rows also in Lab Tests 7 -> file 6 not used",
                        format(dup6, big.mark = ","), format(nrow(unique(x6[, ..k])), big.mark = ",")))
res_put("lab6_in_lab7", list(rows_in_7 = dup6, rows_6 = nrow(unique(x6[, ..k])))); rm(x6)
code_med <- x[, .(cm = median(val, na.rm = TRUE)), by = code]
x[code_med, cm := i.cm, on = "code"]
is_cv <- (x$unit == "%") | (x$unit == "" & x$cm <= 30)
rdw <- daily(clean(x, "RDW-CV", is_cv, rep(TRUE, nrow(x)), 8, 40))
# ---- Platelets (Lab Tests 8) ----------------------------------------------------
x <- read_lab(8)
plt_desc <- grepl("PLATELET", x$desc_n) & !grepl("ESTIMATE|\\bEST\\b|LARGE|COMMENT|UNIT INFO|MEAN|MPV|VOLUME|CLUMP|AGG|FUNCTION|ANTIBOD|SMEAR|MORPH|IMMATURE|IPF", x$desc_n)
plt_unit <- !grepl("%|FL", x$unit)   # all count units in the extract are x10^3/uL-equivalent (checked: medians 200-270)
plt <- daily(clean(x, "Platelets", plt_desc, plt_unit, 5, 1500))
rm(x); gc()

lab_steps <- rbindlist(steps)
save_tab(lab_steps, "S3d_lab_cleaning_steps", "Laboratory cleaning: rows remaining after each step")
saveRDS(list(hb = hb, mcv = mcv, mch = mch, mchc = mchc, rdw = rdw, plt = plt), file.path(DERIVED, "labs_long.rds"))

# ---- Validation against the Codebook hgb (rebuilt from a different export) -----
# Codebook rule: closest to index within +-1095 d, ties to pre-index, POC included, blood gas excluded.
e <- m[eligible == 1, .(study_id, index_date, hgb, hgb_date)]
cand <- merge(hb, e[, .(study_id, index_date)], by = "study_id")[abs(date - index_date) <= 1095]
cand[, `:=`(dist = abs(as.integer(date - index_date)), pre = as.integer(date <= index_date))]
setorder(cand, study_id, dist, -pre)
cb_like <- cand[, .(hb_re = val[1], date_re = date[1]), by = study_id]
v <- merge(e, cb_like, by = "study_id", all.x = TRUE)
val <- data.table(
  check = c("Eligible women with Codebook hgb", "Eligible women with re-derived Hb (Codebook rule)",
            "Both present", "Same date", "Same date and value within 0.05 g/dL",
            "Codebook hgb present, no re-derived value", "Re-derived value present, Codebook hgb blank",
            "Anemia (<12) classification agrees (both present)"),
  n = c(sum(!is.na(v$hgb)), sum(!is.na(v$hb_re)), v[!is.na(hgb) & !is.na(hb_re), .N],
        v[!is.na(hgb) & !is.na(hb_re) & hgb_date == date_re, .N],
        v[!is.na(hgb) & !is.na(hb_re) & hgb_date == date_re & abs(hgb - hb_re) <= 0.05, .N],
        v[!is.na(hgb) & is.na(hb_re), .N], v[is.na(hgb) & !is.na(hb_re), .N],
        v[!is.na(hgb) & !is.na(hb_re), sum((hgb < 12) == (hb_re < 12))]))
save_tab(val, "S3e_hb_validation", "Codebook hgb vs Hb re-derived from Lab Tests 2 with the Codebook rule")
res_put("hb_validation", setNames(as.list(val$n), val$check))
log_msg("Labs", paste(sprintf("%s: %s", val$check, format(val$n, big.mark = ",")), collapse = "; "))
