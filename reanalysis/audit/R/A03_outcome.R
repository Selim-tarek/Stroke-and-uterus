# Audit Q5: outcome algorithm, event provenance, subtype feasibility, sensitivity to a stricter outcome,
# and a chart-review sampling frame (study_id only, written to derived/audit, git-ignored).
source("audit/R/00_audit_setup.R")
md_open("A03_outcome.R", "A03 Stroke outcome")
d <- load_dev()
d[, gap := as.numeric(stroke_date - anchor)]
TYPE <- c("1" = "ischemic", "2" = "intracerebral hemorrhage", "3" = "subarachnoid hemorrhage", "4" = "TIA", "5" = "cerebral venous thrombosis", "9" = "unspecified")
d[, `:=`(type = TYPE[as.character(stroke_type)], imaging_only = grepl("CASE FOUND ON IMAGING ONLY", notes),
         has_acute = grepl("acute", source_kinds), multi_type = grepl("multiple event types coded", notes))]
ev <- d[y == 1]

alg <- data.table(step = c("Source", "Event types", "TIA", "Imaging-only cases", "Event date", "Timing", "One record per woman", "Not done"),
  rule = c("Problem List + Diagnosis extracts, any date in the record, then positive brain-imaging reports (Codebook: stroke_any). Code lists are in the workbook sheet stroke_outcome; the raw stroke diagnosis extracts are not in data/.",
           "stroke_type from the most specific code: 1 ischemic, 2 ICH, 3 SAH, 4 TIA, 5 CVT, 9 unspecified (188 women with >1 type, flagged in notes).",
           "Excluded from the outcome (R/13_anchor_design.R:55), but a TIA ENDS follow-up for stroke (R/13:52, censored).",
           "Included if infarct on imaging without a code, including chronic/incidental infarcts and hedged wording (Codebook stroke_any note). Dated by the scan.",
           "Earliest coded diagnosis or problem date (Codebook stroke_date); not chart-verified. Acute, history and sequela codes all count.",
           "Incident if stroke_date > anchor (R/13:50); stroke on or before anchor, or undated, excluded.",
           "Single stroke_any/stroke_type/stroke_date per woman, so recurrent events cannot be counted and the first event is the earliest code of any kind.",
           "No chart adjudication; no source-document review was possible in this audit (no notes or charts in the project)."))
aud_tab(alg, "A03a_algorithm", "Outcome algorithm, from the Codebook notes and the scripts")

prov <- rbind(ev[, .(category = "all strokes", n = .N)],
  ev[, .N, by = .(category = paste("type:", type))][order(-N)][, .(category, n = N)],
  data.table(category = c("has an acute-event code", "history/sequela codes only (no acute code, not imaging-only)", "imaging-only case (no code)",
                          "  of which imaging text says chronic/remote/old", "brain imaging reviewed and positive (stroke_confirmed_imaging = 1)", "brain imaging reviewed and negative (= 0)", "no brain imaging (blank)",
                          "more than one event type coded", "event within 30 days after anchor", "event 31–90 days after anchor", "event 91–365 days after anchor", "event > 1 year after anchor"),
             n = c(ev[has_acute == TRUE, .N], ev[has_acute == FALSE & imaging_only == FALSE, .N], ev[imaging_only == TRUE, .N],
                   ev[imaging_only == TRUE & grepl("chronic|remote|old ", notes, ignore.case = TRUE), .N], ev[stroke_confirmed_imaging %in% 1, .N],
                   ev[stroke_confirmed_imaging %in% 0, .N], ev[is.na(stroke_confirmed_imaging), .N], ev[multi_type == TRUE, .N],
                   ev[gap <= 30, .N], ev[gap > 30 & gap <= 90, .N], ev[gap > 90 & gap <= 365, .N], ev[gap > 365, .N])))
prov[, pct := pct(n / nrow(ev))]
aud_tab(prov, "A03b_event_provenance", "The 345 development-sample strokes by subtype, code provenance, imaging and timing")

# subtypes: events, person-time and feasibility (events per candidate parameter for the nine-item model)
sub <- d[, .(events = sum(y)), by = type][!is.na(type) & type != "TIA"]
sub <- rbind(sub, data.table(type = "TIA (censoring event)", events = d[ev_tia == 1, .N]))
sub[, `:=`(rate_per_1000py = sapply(events, function(k) pois_ci(k, sum(d$t))), events_per_parameter_9 = round(events / 9, 1))]
aud_tab(sub[order(-events)], "A03c_subtypes", sprintf("Events by subtype over %s person-years (development sample). Fewer than ~10 events per parameter makes a subtype-specific model unstable.", format(round(sum(d$t)), big.mark = ",")))

# sensitivity: stricter outcome definitions (other events censored at their date)
fit_def <- function(ok, lab) {
  dd <- copy(d); dd[, y := as.integer(y == 1 & ok)]
  f <- coxph(fml_of(PRED9), data = dd)
  data.table(definition = lab, strokes = sum(dd$y), predictor = names(coef(f)), hr = sprintf("%.2f (%.2f–%.2f)", exp(coef(f)), exp(confint(f)[, 1]), exp(confint(f)[, 2])),
             c = f2(concordance(f)$concordance, 3))
}
defs <- rbind(fit_def(rep(TRUE, nrow(d)), "A. as published"),
              fit_def(d$has_acute, "B. acute-event code required"),
              fit_def(d$has_acute & d$gap > 30, "C. acute code and > 30 days after anchor"),
              fit_def(d$has_acute & d$gap > 30 & d$stroke_type == 1, "D. as C, ischemic only"))
w <- dcast(defs, predictor ~ definition, value.var = "hr")
w <- rbind(w, dcast(unique(defs[, .(definition, predictor = "Harrell C (apparent)", c)]), predictor ~ definition, value.var = "c"),
           dcast(unique(defs[, .(definition, predictor = "strokes", s = as.character(strokes))]), predictor ~ definition, value.var = "s"))
aud_tab(w, "A03d_strict_outcome", "Nine-item Cox model (unshrunk HR, 95% CI) under stricter outcome definitions")

# chart-review sampling frame: all non-acute / imaging-only / ≤30-day events plus a random sample of the rest; PPV precision
set.seed(SEED)
must <- ev[has_acute == FALSE | imaging_only == TRUE | gap <= 30, study_id]
rest <- setdiff(ev$study_id, must)
samp <- data.table(study_id = c(must, sample(rest, 100 - length(must))))
samp[, stratum := fifelse(study_id %in% must, "certain review (non-acute, imaging-only or ≤30 d)", "random sample of remaining events")]
fwrite(samp, file.path(DERIVED, "audit", "chart_review_sample_events.csv"))
ppv_hw <- function(p, n) 1.96 * sqrt(p * (1 - p) / n)
aud_tab(data.table(stratum = c("certain review", "random sample of remaining events", "total"), n = c(length(must), 100 - length(must), 100),
                   events_in_stratum = c(length(must), length(rest), nrow(ev))),
        "A03e_chart_review_frame", sprintf("Proposed blinded chart-review sample (study_id list in derived/audit/chart_review_sample_events.csv, not committed). Estimate PPV with weights = stratum size / sampled. For a PPV near 0.85, n = 100 gives ±%.2f; sensitivity needs a sample of non-events and is not feasible at a 1.4%% event rate without thousands of reviews.", ppv_hw(.85, 100)))
res_put("A03", list(provenance = prov, strict = defs))
