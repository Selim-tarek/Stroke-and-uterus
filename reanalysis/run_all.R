# Run the whole reanalysis from the source files. From the reanalysis/ directory:
#   Rscript run_all.R
# Source data are read from ../data (override with STROKE_DATA_DIR). Patient-level intermediates go to derived/ (git-ignored).
steps <- c("R/01_import.R", "R/02_labs.R", "R/03_comorbid.R", "R/04_cohort.R", "R/05_missing.R", "R/06_reproduce.R",
           "R/07_primary.R", "R/08_secondary.R", "R/09_timeupdated.R", "R/10_sensitivity.R", "R/11_tables_figures.R", "R/12_manuscript.R", "R/13_anchor_design.R")
if (!file.exists("output/original_reported.csv"))
  stop("Run first: python3 tools/export_original_results.py ../outputs/_state.pkl output/original_reported.csv")
for (s in steps) {
  message("\n==== ", s, " ====")
  t <- system.time(status <- system2("Rscript", s))
  if (status != 0) stop("Step failed: ", s)
  message(sprintf("%s finished in %.0f s", s, t[["elapsed"]]))
}
writeLines(capture.output(sessionInfo()), "output/logs/sessionInfo.txt")
