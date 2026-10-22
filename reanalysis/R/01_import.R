# Phase 1a: import the analysis master and Codebook, validate structure,
# build the variable dictionary. Patient-level output -> derived/ (git-ignored).
source("R/00_setup.R")
cat("# Run log\n\n", file = LOG_FILE)
cat("# Deviations from the original analysis and from the SAP\n\n| Item | Reason |\n|---|---|\n", file = DEV_FILE)
cat("# Model warnings (captured, not suppressed)\n\n| Model | Warnings |\n|---|---|\n", file = WARN_FILE)
if (file.exists(RES_FILE)) file.remove(RES_FILE)

m <- fread(dpath("stroke_analysis_master.csv"), na.strings = c("", "NA"))
log_msg("Import", sprintf("stroke_analysis_master.csv: %s rows x %d columns", format(nrow(m), big.mark = ","), ncol(m)))

# ---- structural checks ----------------------------------------------------
chk <- data.table(
  check = c("rows", "unique study_id", "unique MRN (link key)", "missing MRN",
            "eligible==1 equals exclusion_reason==0"),
  value = c(nrow(m), uniqueN(m$study_id), uniqueN(m$mrn), sum(is.na(m$mrn)),
            as.integer(all((m$eligible == 1) == (m$exclusion_reason == 0)))))
# eligibility re-derived from the Codebook rule (age 18-60 at index, no active malignancy)
m[, elig_rederived := as.integer(!is.na(age_index) & age_index >= 18 & age_index <= 60 & malignancy_active == 0)]
chk <- rbind(chk, data.table(check = c("eligible == re-derived rule (age 18-60, no active malignancy)",
                                       "re-derived eligible count"),
                             value = c(sum(m$eligible == m$elig_rederived), sum(m$elig_rederived))))
save_tab(chk, "S3a_structure_checks")

# ---- Codebook and variable dictionary ------------------------------------
cb <- as.data.table(read_excel(dpath("stroke.xlsx"), sheet = "Codebook", .name_repair = "minimal"))
cb <- cb[!is.na(Variable)]
setnames(cb, c("num", "domain", "variable", "label", "type", "values", "required", "notes"))
role <- c(study_id = "identifier", eligible = "selection", exclusion_reason = "selection", index_date = "time origin",
          age_index = "confounder", race = "confounder", bmi = "confounder", smoking = "confounder",
          htn = "confounder", dm = "confounder", dyslipidemia = "confounder", cad = "confounder",
          migraine = "confounder", hormonal_type = "confounder", uterine_dx_group = "confounder",
          uterine_bleeding = "exposure determinant (sensitivity)", thrombophilia = "undated (sensitivity)",
          afib = "undated (sensitivity)", chf = "undated (sensitivity)", vte_history = "undated (sensitivity)",
          hgb = "exposure (original)", hgb_date = "exposure timing (original)", hgb_days = "exposure timing (original)",
          anemia = "exposure (original)", mcv = "secondary exposure (Codebook; undated)", platelets = "secondary exposure (Codebook; undated)",
          ferritin = "exploratory exposure", ferritin_date = "exploratory exposure timing", iron_deficiency = "exploratory exposure",
          stroke_any = "outcome", stroke_type = "outcome", stroke_date = "outcome timing", stroke_timing = "outcome timing",
          malignancy_active = "selection", lab_date = "NOT USED (back-filled with d-dimer date)", ldl = "NOT USED (implausible distribution)")
e <- m[eligible == 1]
dict <- rbindlist(lapply(setdiff(names(m), c("mrn", "dob", "elig_rederived")), function(v) {
  x <- e[[v]]
  cbr <- cb[variable == v]
  vals <- if (is.numeric(x) && uniqueN(x, na.rm = TRUE) <= 12)
    paste(sprintf("%s: %s", names(table(x)), format(as.integer(table(x)), big.mark = ",")), collapse = "; ") else
    if (is.numeric(x)) sprintf("median %s (IQR %s–%s)", f2(median(x, na.rm = TRUE), 1),
                               f2(quantile(x, .25, na.rm = TRUE), 1), f2(quantile(x, .75, na.rm = TRUE), 1)) else
    sprintf("%d distinct values", uniqueN(x))
  data.table(variable = v,
             codebook_label = if (nrow(cbr)) cbr$label[1] else "(not in Codebook)",
             codebook_type = if (nrow(cbr)) cbr$type[1] else "",
             codebook_values_units = if (nrow(cbr)) cbr$values[1] else "",
             observed_type = class(x)[1],
             observed_eligible = vals,
             missing_blank_eligible = sum(is.na(x)),
             missing_pct = round(100 * mean(is.na(x)), 2),
             code9_eligible = if (is.numeric(x)) sum(x == 9, na.rm = TRUE) else NA_integer_,
             role = if (!is.na(role[v])) role[v] else "descriptive / not used")
}))
save_tab(dict, "variable_dictionary", "Variable dictionary, eligible women (N = 39,807)")
log_msg("Dictionary", sprintf("%d variables documented (MRN and DOB omitted)", nrow(dict)))

# ---- date fields ----------------------------------------------------------
for (v in c("index_date", "hgb_date", "ferritin_date", "stroke_date", "lab_date", "surgical_tx_date", "myomectomy_date"))
  m[, (v) := as.IDate(get(v))]
dchk <- rbindlist(lapply(c("index_date", "hgb_date", "ferritin_date", "stroke_date", "lab_date"), function(v)
  data.table(field = v, n_eligible = sum(!is.na(e[[v]])), min = as.character(min(m[eligible == 1][[v]], na.rm = TRUE)),
             max = as.character(max(m[eligible == 1][[v]], na.rm = TRUE)),
             sentinel_1000 = sum(year(m[eligible == 1][[v]]) < 1990, na.rm = TRUE))))
save_tab(dchk, "S3b_date_fields")

m[, dob := NULL]
saveRDS(m, file.path(DERIVED, "master.rds"))
res_put("import", list(rows = nrow(m), eligible = sum(m$eligible == 1),
                       eligible_rederived_agree = sum(m$eligible == m$elig_rederived)))
