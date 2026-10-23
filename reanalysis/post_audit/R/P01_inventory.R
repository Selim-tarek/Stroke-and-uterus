# Stage 1: input inventory, integrity check of the published model, and decision timeline from git.
source("post_audit/R/P00_setup.R")
md_open("P01_inventory.R", "P01 Input inventory, integrity check and decision timeline")

# ---- 1. files actually present --------------------------------------------------------------------------------
fl <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
inv <- data.table(file = sub(paste0("^", DATA_DIR, "/"), "", fl), size_mb = round(file.size(fl) / 1e6, 1))
aud_tab(inv, "P01a_data_files", "Every file under data/ (read-only source)")
lab_desc <- rbindlist(lapply(c(2, 3, 4, 7, 8), function(f) { x <- fread(dpath("labs", sprintf("MDE Workflow Results for Lab Tests (%d).csv", f)), select = 3, nrows = 2e5)
  data.table(file = paste0("Lab Tests (", f, ")"), top_descriptions = paste(head(names(sort(table(x[[1]]), decreasing = TRUE)), 4), collapse = "; ")) }))
aud_tab(lab_desc, "P01b_lab_contents", "Third column of each lab extract (first 200,000 rows). Files 3 and 7 use a different column layout (MCV and RDW values / dates appear in that position; see R/02_labs.R). All six extracts are CBC components only: no HbA1c, glucose, lipids, creatinine or BMI.")
sh <- excel_sheets(dpath("stroke.xlsx"))
dated <- c(smoking_source = "statement_date", radiology_stroke = "earliest_positive_date", nihss_flowsheet = "first_recorded", stroke_outcome = "earliest_event_date",
           myomectomy_procedures = "procedure_date", ddimer_source = "result_date", hgb_ferritin_update = "result_date")
aud_tab(data.table(sheet = names(dated), date_field = dated, present = names(dated) %in% sh),
        "P01c_dated_workbook_sheets", sprintf("Per-patient dated sheets in stroke.xlsx (of %d sheets). Term lists (migraine_terms, cad_terms, hormonal_meds, ...) hold counts only, no patient dates.", length(sh)))
mx <- fread(dpath("meds", "MDE_Workflow_Results_for_Medications_Administered_9.csv"), select = 15)[[1]]
md_note("Medication extract present: file 9 only; generic names containing 'estradiol' or 'estrogen' in ", sum(grepl("estradiol|estrogen", mx, ignore.case = TRUE)),
        " of ", format(length(mx), big.mark = ","), " rows (iron-containing combined OCs). The extract used to derive `hormonal_type` (Codebook: Administered Medications) is not in the project.")

# ---- 2. integrity: the published nine-item model reproduces from the current derived data ------------------------
d <- load_dev(); require_cols(d, c(PRED9, "t", "y", "anchor", "index_date"), "development sample")
f <- coxph(fml_of(PRED9), data = d)
chk <- data.table(item = c("women", "strokes", "max |beta − spec|", "S0(5) − spec"),
                  value = c(nrow(d), sum(d$y), signif(max(abs(coef(f)[PRED9] - unlist(SPEC9$coef[PRED9]))), 3), signif(s0_at(f, 5) - SPEC9$s0, 3)))
aud_tab(chk, "P01d_integrity", "Published nine-item model refitted from derived/anchored_cohort.rds (reported: 24,006 women, 345 strokes)")
if (nrow(d) != 24006 || sum(d$y) != 345 || max(abs(coef(f)[PRED9] - unlist(SPEC9$coef[PRED9]))) > 1e-6) stop("Derived data no longer reproduce the published model; stop and investigate.")

# ---- 3. decision timeline (git) ---------------------------------------------------------------------------------
g <- system2("git", shQuote(c("log", "--reverse", "--date=format:%Y-%m-%d %H:%M", "--format=%h|%ad|%s", "--", "SAP.md", "R/14_extensions.R", "R/15_risk_model.R",
                      "R/16_calculator_model.R", "calculator/bud_stroke_calculator.html", "audit/AUDIT_REPORT.md")), stdout = TRUE)
gt <- as.data.table(tstrsplit(g, "|", fixed = TRUE)); setnames(gt, c("commit", "time", "message"))
aud_tab(gt, "P01e_git_timeline", "Commits touching the SAP, risk-model scripts, calculator and audit, in order")
res_put("P01", list(integrity = chk))
