# Phase 8: figures and the editable table workbook / Word file, generated from the saved outputs only.
source("R/00_setup.R")
suppressPackageStartupMessages({ library(officer); library(flextable) })
TAB <- function(n) fread(file.path(OUT, "tables", paste0(n, ".csv")), encoding = "UTF-8")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
theme_pub <- theme_minimal(base_size = 10, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK), plot.title = element_text(face = "bold", colour = INK),
        strip.text = element_text(face = "bold", hjust = 0, colour = INK), plot.background = element_rect(fill = "white", colour = NA))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT, "figures", paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

# ---- Figure 1: cohort flow ---------------------------------------------------------------------------------
fl <- TAB("Fig1_flow_counts")
main <- fl[!step %like% "^Excluded"]; exc <- fl[step %like% "^Excluded"]
nb <- function(x) format(x, big.mark = ",")
boxes <- data.table(
  y = c(6, 4, 2, 0),
  label = c(sprintf("Records in the analysis file\nN = %s", nb(fl[1, n])),
            sprintf("Eligible women (18\u201360 y, no active malignancy)\nN = %s", nb(main[step == "Eligible women", n])),
            sprintf("Free of recorded stroke/TIA at index\nN = %s", nb(main[step %like% "^Free", n])),
            sprintf("Incident-stroke cohort\nN = %s\n%s incident strokes (excl. TIA)", nb(main[step %like% "^Incident", n]),
                    nb(fromJSON(RES_FILE)$cohort$events$stroke))))
side <- data.table(y = c(5, 3, 1),
  label = c(paste(sprintf("%s: %s", sub("Excluded: ", "", exc$step[1:4]), nb(exc$n[1:4])), collapse = "\n"),
            paste(sprintf("%s: %s", sub("Excluded: ", "", exc$step[5:6]), nb(exc$n[5:6])), collapse = "\n"),
            paste(sprintf("%s: %s", sub("Excluded: ", "", exc$step[7:8]), nb(exc$n[7:8])), collapse = "\n")))
p1 <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 5.55, yend = 0.65), colour = INK2, linewidth = 0.4, arrow = arrow(length = unit(2, "mm"))) +
  geom_segment(data = side, aes(x = 0, xend = 0.9, y = y, yend = y), colour = INK2, linewidth = 0.4, arrow = arrow(length = unit(2, "mm"))) +
  geom_label(data = boxes, aes(x = 0, y = y, label = label), size = 3, label.r = unit(1, "mm"), label.padding = unit(2.5, "mm"), fill = "white", colour = INK) +
  geom_label(data = side, aes(x = 0.95, y = y, label = label), hjust = 0, size = 2.6, label.r = unit(1, "mm"), fill = "#f6f6f4", colour = INK2) +
  xlim(-1.4, 4.6) + ylim(-0.7, 6.7) + theme_void()
save_fig(p1, "Figure1_flow", 8.5, 6.5)

# ---- Figure 2: continuous Hb (spline) with Hb distribution ----------------------------------------------------------
sp <- readRDS(file.path(DERIVED, "spline_primary.rds"))
cur <- sp$curve
hist_dt <- data.table(hb = sp$hb[sp$hb >= sp$range[1] & sp$hb <= sp$range[2]])
top <- ggplot(cur, aes(hb, rr)) +
  geom_hline(yintercept = 1, colour = INK2, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = BLUE, alpha = 0.18) +
  geom_line(colour = BLUE, linewidth = 0.8) +
  geom_vline(xintercept = c(8, 10, 12), colour = GRID, linetype = 2, linewidth = 0.3) +
  annotate("point", x = 13, y = 1, colour = INK, size = 2) +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  labs(x = NULL, y = "Adjusted rate ratio vs Hb 13 g/dL (log scale)",
       title = "Baseline hemoglobin and incident stroke",
       subtitle = sprintf("Restricted cubic spline (knots %s g/dL); shaded band 95%% CI;\noverall P = %s, non-linearity P = %s",
                          paste(f2(sp$knots, 1), collapse = ", "), fp(fromJSON(RES_FILE)$primary$spline_p_overall), fp(fromJSON(RES_FILE)$primary$spline_p_nonlinear))) +
  scale_x_continuous(limits = sp$range, expand = c(0, 0)) + theme_pub + theme(plot.margin = margin(5, 15, 0, 5), plot.subtitle = element_text(colour = INK2, size = 8))
bot <- ggplot(hist_dt, aes(hb)) + geom_histogram(binwidth = 0.25, boundary = 0, fill = "#9ec5f4", colour = "white", linewidth = 0.2) +
  scale_x_continuous(limits = sp$range, expand = c(0, 0)) + scale_y_continuous(labels = function(x) format(x, big.mark = ",")) + labs(x = "Baseline hemoglobin, g/dL (1st\u201399th percentile shown)", y = "Women") + theme_pub
g2 <- gridExtra_arrange <- NULL
if (requireNamespace("patchwork", quietly = TRUE)) {
  p2 <- patchwork::wrap_plots(top, bot, ncol = 1, heights = c(3, 1))
  save_fig(p2, "Figure2_spline", 7, 6)
} else {
  png(file.path(OUT, "figures", "Figure2_spline.png"), width = 7, height = 6, units = "in", res = 300)
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(4, 1)))
  print(top, vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1)); print(bot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1))
  dev.off()
  pdf(file.path(OUT, "figures", "Figure2_spline.pdf"), width = 7, height = 6)
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(4, 1)))
  print(top, vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1)); print(bot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1))
  dev.off()
}

# ---- Figure 3: forest plot, separated by design and effect measure ---------------------------------------------------
parse_ci <- function(s) { m <- regmatches(s, regexec("([0-9.]+) \\(([0-9.]+)\u2013([0-9.]+)\\)", s)); t(sapply(m, function(z) if (length(z) == 4) as.numeric(z[2:4]) else rep(NA, 3))) }
t2 <- TAB("Table2_primary_models")[analysis %like% "^Primary" & measure %like% "Rate" & level %in% c("Mild", "Moderate", "Severe", "Per grade (trend)")]
t3 <- TAB("Table3_secondary_outcomes")[outcome %in% c("Ischemic stroke", "TIA") & level %in% c("Mild", "Moderate", "Severe", "Per grade (trend)")]
t4 <- TAB("Table4_time_updated")[outcome == "Stroke (excl. TIA)" & analysis %in% c("TU primary: 1-y carry-forward, no lag", "TU: 30-day lag") & level %in% c("Mild", "Moderate", "Severe", "Per grade (trend, measured time)")]
xs <- TAB("S2a_reproduction_comparison")[Analysis == "Cross-sectional logistic, adjusted" & Outcome %like% "Any stroke" & Level %in% c("Mild", "Moderate", "Severe")]
fp_dt <- rbind(
  data.table(panel = "A. Incident stroke after index (primary)\nRate ratio, baseline Hb", series = "Stroke (excl. TIA)", level = t2$level, txt = t2$estimate),
  data.table(panel = "B. Secondary incident outcomes\nRate ratio, baseline Hb", series = t3$outcome, level = t3$level, txt = t3$estimate),
  data.table(panel = "C. Time-updated Hb (most recent value)\nRate ratio", series = fifelse(t4$analysis %like% "30-day", "30-day lag", "No lag"),
             level = sub(", measured time", "", t4$level), txt = t4$estimate),
  data.table(panel = "D. Original cross-sectional analysis\nOdds ratio, any stroke ever (includes prevalent)", series = "Original model (reproduced)",
             level = xs$Level, txt = xs$`Reproduced in R`))
fp_dt[, c("est", "lo", "hi") := as.data.table(parse_ci(txt))]
fp_dt[, multi := uniqueN(series) > 1, by = panel]
fp_dt[, lab := fifelse(multi, paste0(level, " \u2014 ", series), level)]
fp_dt[, lab := factor(lab, levels = rev(unique(lab)))]
fp_dt[, series := factor(series, levels = unique(series))]
cols <- setNames(c(BLUE, BLUE, ORANGE, BLUE, ORANGE, INK2), c("Stroke (excl. TIA)", "Ischemic stroke", "TIA", "No lag", "30-day lag", "Original model (reproduced)"))
shp <- setNames(c(16, 16, 17, 16, 17, 15), names(cols))
p3 <- ggplot(fp_dt[!is.na(est)], aes(est, lab, colour = series, shape = series)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6) +
  geom_point(size = 2.4) +
  geom_text(aes(x = 13, label = txt), hjust = 1, size = 2.5, colour = INK2) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4, 6), limits = c(0.45, 13.2)) +
  scale_colour_manual(values = cols, guide = "none") + scale_shape_manual(values = shp, guide = "none") +
  labs(x = "Estimate vs no anemia (Hb \u226512 g/dL), log scale, 95% CI", y = NULL) +
  theme_pub
save_fig(p3, "Figure3_forest", 8.5, 10)

# ---- Supplementary figure: cumulative incidence by baseline grade --------------------------------------------------
aj <- readRDS(file.path(DERIVED, "aj_primary.rds"))
k <- which(aj$states == "stroke")
cdt <- rbindlist(lapply(seq_along(aj$strata), function(i) {
  idx <- if (i == 1) 1:aj$strata[1] else (sum(aj$strata[1:(i - 1)]) + 1):sum(aj$strata[1:i])
  data.table(grade = sub("grade=", "", names(aj$strata)[i]), time = c(0, aj$time[idx]), cif = c(0, aj$pstate[idx, k]))}))
cdt[, grade := factor(GRADE_LABELS[grade], levels = GRADE_LABELS)]
pS <- ggplot(cdt[time <= 6], aes(time, 100 * cif, linetype = grade)) + geom_step(colour = BLUE, linewidth = 0.6) +
  labs(x = "Years after index", y = "Cumulative incidence of stroke, %", linetype = NULL,
       title = "Cumulative incidence of stroke by baseline anemia grade", subtitle = "Aalen\u2013Johansen; death as competing event; crude") +
  theme_pub + theme(legend.position = "bottom", plot.subtitle = element_text(colour = INK2, size = 8))
save_fig(pS, "FigureS1_cumulative_incidence", 7, 4.5)

# ---- editable tables: Excel workbook + Word -------------------------------------------------------------------------
order <- c("Table1_baseline_by_grade", "Table2_primary_models", "Table2a_rates_by_grade", "Table2b_continuous_hb", "Table3_secondary_outcomes",
           "Table3a_secondary_event_counts", "Table4_time_updated", "Table4a_time_updated_events", "Table4b_landmark", "Table5_sensitivity",
           "S1a_tested_vs_untested", "S1b_rates_tested_vs_untested", "S1c_missingness_incident_cohort", "S1d_missing_patterns",
           "S2a_reproduction_comparison", "S2b_ph_test_primary", "S3a_structure_checks", "S3b_date_fields", "S3c_hb_descriptions",
           "S3d_lab_cleaning_steps", "S3e_hb_validation", "S3f_comorbidity_flags", "S4a_outcome_type_by_timing", "S4b_outcome_code_kinds",
           "S4c_multiple_types", "Fig1_flow_counts", "variable_dictionary")
order <- order[file.exists(file.path(OUT, "tables", paste0(order, ".csv")))]
wb <- createWorkbook(); modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
hs <- createStyle(textDecoration = "bold", fgFill = "#2a78d6", fontColour = "white", wrapText = TRUE)
for (n in order) { s <- substr(n, 1, 31); addWorksheet(wb, s); writeData(wb, s, TAB(n), headerStyle = hs); setColWidths(wb, s, 1:60, "auto"); freezePane(wb, s, firstRow = TRUE) }
saveWorkbook(wb, file.path(OUT, "reanalysis_tables.xlsx"), overwrite = TRUE)
doc <- read_docx()
for (n in order[!order %in% c("S3c_hb_descriptions", "variable_dictionary", "S3f_comorbidity_flags")]) {
  d <- TAB(n)
  doc <- body_add_par(doc, n, style = "heading 2")
  ft <- flextable(d) |> fontsize(size = 7, part = "all") |> font(fontname = "Arial", part = "all") |> bold(part = "header") |> autofit()
  doc <- body_add_flextable(doc, ft)
}
print(doc, target = file.path(OUT, "reanalysis_tables.docx"))
log_msg("Outputs", sprintf("%d tables -> reanalysis_tables.xlsx / .docx; figures 1-3 + S1", length(order)))
