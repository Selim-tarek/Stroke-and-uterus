# Phase 5b: secondary incident outcomes (same cohort and time origin as the primary analysis).
source("R/00_setup.R"); source("R/model_fns.R")
ic <- readRDS(file.path(DERIVED, "incident_cohort.rds"))[!is.na(hormonal_type_cat)]
OUTS <- c(ev_isch = "Ischemic stroke", ev_any = "Stroke or TIA (Codebook stroke_any)", ev_tia = "TIA",
          ev_ich = "Intracerebral hemorrhage", ev_sah = "Subarachnoid hemorrhage", ev_cvt = "Cerebral venous thrombosis",
          ev_unk = "Stroke, unknown type")
rows <- list(); cnt <- list()
for (y in names(OUTS)) {
  d <- prep_covs(ic, y)
  cnt[[y]] <- d[, .(outcome = OUTS[y], events = sum(get(y))), by = grade]
  n_ev <- sum(d[[y]])
  if (n_ev < 10) { rows[[y]] <- data.table(outcome = OUTS[y], level = "not modeled", estimate = sprintf("%d events (<10): counts only", n_ev), events = n_ev); next }
  covs <- if (n_ev / 32 >= 10) PRIMARY_COVS else if (n_ev / 16 >= 10) REDUCED_COVS else c("age_index", "htn")
  spec <- if (n_ev / 32 >= 10) "full adjustment" else if (n_ev / 16 >= 10) "reduced pre-declared set (EPV<10 with full set)" else "age + hypertension (few events)"
  dl <- if (spec == "full adjustment" || spec %like% "reduced") impute_bmi(d, y, m = 20) else list(d)
  # categorical estimates only if every anemia level has >= 5 events
  min_ev <- min(d[, sum(get(y)), by = grade]$V1)
  r <- grade_models(dl, y, covs, paste0("Secondary: ", OUTS[y]), cox = FALSE)
  if (min_ev < 5) r <- r[!(level %in% GRADES[-1])]
  rows[[y]] <- r[, .(outcome = OUTS[y], level, estimate = fci(est, lo, hi), p_raw = p, adjustment = spec, events = n_ev, epv = round(epv, 1),
                     note = fifelse(min_ev < 5, "grade categories not estimated (<5 events in a level)", ""))]
}
sec <- rbindlist(rows, fill = TRUE)
sec[level == "Per grade (trend)", q_BH := bh(p_raw)]
sec[, `:=`(p = fp(p_raw), q_BH = fp(q_BH))]
sec[epv < 10, note := trimws(paste(note, "EPV<10: unstable, exploratory"))]
save_tab(sec[, !"p_raw"], "Table3_secondary_outcomes", "Secondary incident outcomes. Benjamini-Hochberg q-values across the per-grade trend tests.")
cnts <- dcast(rbindlist(cnt), outcome ~ grade, value.var = "events")
save_tab(cnts, "Table3a_secondary_event_counts", "Incident events by baseline anemia grade, per outcome")
res_put("secondary", sec)
