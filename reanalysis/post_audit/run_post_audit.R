# Post-audit investigation of the MED-HACK model. From reanalysis/:  Rscript post_audit/run_post_audit.R
# Requires derived/ from run_all.R. Writes only to post_audit/output/ and derived/post_audit/ (patient level, git-ignored).
# Order matters: P02 builds derived/post_audit/baseline_cohort.rds used by P03–P06. P05 (~25 min) caches its bootstrap in
# derived/post_audit/P05_bootstrap.rds; delete that file to recompute it.
steps <- sprintf("post_audit/R/%s.R", c("P01_inventory", "P02_baseline_reconstruction", "P03_outcome_definitions", "P04_horizon_censoring",
                                        "P05_candidate_models", "P06_temporal_drift", "P07_calculator_support", "P08_fibroid_subgroups"))
for (s in steps) {
  message("\n==== ", s, " ====")
  t <- system.time(status <- system2("Rscript", s))
  if (status != 0) stop("Step failed: ", s)
  message(sprintf("%s finished in %.0f s", s, t[["elapsed"]]))
}
