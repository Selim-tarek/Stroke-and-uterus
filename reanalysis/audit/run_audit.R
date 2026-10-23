# Run the full statistical audit of the MED-HACK calculator. From reanalysis/:  Rscript audit/run_audit.R
# Reads derived/ and output/ produced by run_all.R; writes only to audit/output/ and derived/audit/ (nothing in the original analysis is changed).
# Requires node + playwright (global) for A05; A08 (~10 min) and A11 (~5 min) are the slow steps. Delete derived/audit/A08_bootstrap.rds to recompute the A08 bootstrap.
steps <- sprintf("audit/R/%s.R", c("A01_cohort_design", "A02_predictors_missingness", "A03_outcome", "A04_followup_competing", "A05_reproduce_calculator",
                                   "A06_diagnostics_age", "A07_points", "A08_internal_validation", "A09_temporal_validation", "A10_subgroups", "A11_multiple_imputation"))
for (s in steps) {
  message("\n==== ", s, " ====")
  t <- system.time(status <- system2("Rscript", s))
  if (status != 0) stop("Step failed: ", s)
  message(sprintf("%s finished in %.0f s", s, t[["elapsed"]]))
}
writeLines(capture.output(sessionInfo()), "audit/output/sessionInfo.txt")
