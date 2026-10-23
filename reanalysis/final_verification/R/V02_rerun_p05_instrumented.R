# Part 3: re-run the P05 full-process bootstrap from scratch (same seed) to (a) test that the cached bootstrap reproduces and
# (b) record, for each resample, the number of events among hbp_pre-positive women — the original cache did not store resample indices.
# The original script is NOT modified: its text is read, four documented substitutions are applied, and the copy is run.
src <- readLines("post_audit/R/P05_candidate_models.R")
subs <- list(
  c('cache <- file.path(DERIVED, "post_audit", "P05_bootstrap.rds")', 'cache <- file.path(DERIVED, "final_verification", "P05_bootstrap_rerun.rds")'),
  c('md_open("P05_candidate_models.R", "P05 Candidate-model comparison (baseline cohort, acute-coded stroke, 3 years)")',
    'AUD <- "final_verification/output/p05_rerun"; OUT <- AUD; dir.create(file.path(AUD, "tables"), FALSE, TRUE); LOG_FILE <- file.path(AUD, "run_log.md"); RES_FILE <- file.path(AUD, "results.json"); md_open("P05_candidate_models.R", "P05 re-run (verification copy)")'),
  c('data.table(b = i, model = n, c_app', 'data.table(b = i, model = n, hbp_pos_events = sum(db[[MODELS[[n]]$outcome]][db$hbp_pre == 1]), hbp_pos_women = sum(db$hbp_pre), c_app'))
for (s in subs) { k <- grep(s[1], src, fixed = TRUE); if (length(k) != 1) stop("substitution target not found exactly once: ", s[1]); src[k] <- sub(s[1], s[2], src[k], fixed = TRUE) }
f <- file.path(tempdir(), "P05_verification_copy.R"); writeLines(src, f)
dir.create("derived/final_verification", FALSE, TRUE)
if (file.exists("derived/final_verification/P05_bootstrap_rerun.rds")) file.remove("derived/final_verification/P05_bootstrap_rerun.rds")
source(f, echo = FALSE)
