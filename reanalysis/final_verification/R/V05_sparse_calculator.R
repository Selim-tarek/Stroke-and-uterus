# Parts 3 and 7: denominators for the sparse predictors, and facts about the live calculator read directly from its source.
source("final_verification/R/V00_setup.R")
md_open("V05_sparse_calculator.R", "V05 Sparse predictors and calculator facts")
dev <- load_dev(); b <- load_base()
sp <- data.table(predictor = c("estrogen (as coded), published development sample", "hemoglobinopathy, any date, published development sample", "hemoglobinopathy before anchor, baseline cohort",
                               "CAD, baseline cohort", "CKD, baseline cohort"),
                 women = c(sum(dev$estrogen), sum(dev$haemoglobinopathy), sum(b$hbp_pre), sum(b$cad), sum(b$ckd_pre)),
                 events = c(dev[estrogen == 1, sum(y)], dev[haemoglobinopathy == 1, sum(y)], b[hbp_pre == 1, sum(y)], b[cad == 1, sum(y)], b[ckd_pre == 1, sum(y)]),
                 outcome = c("O1", "O1", "O2", "O2", "O2"))
sp[, events := sc(events)]
aud_tab(sp, "V05a_sparse_denominators", "Women and events among predictor-positive women (counts 1–10 shown as <11)")
h <- readLines("calculator/standalone/index.html", warn = FALSE); live_same <- identical(tools::md5sum("calculator/standalone/index.html")[[1]], tools::md5sum("calculator/bud_stroke_calculator.html")[[1]])
hit <- function(rx) { k <- grep(rx, h, ignore.case = TRUE); if (length(k)) paste(k, collapse = ", ") else "absent" }
cf <- data.table(item = c("'Research use only' text (line numbers)", "'not externally validated' text", "mentions predictor timing / leakage", "mentions underprediction in recent patients",
                          "horizon shown", "inputs required before a result is shown", "'unknown' option for yes/no items", "upper cap or suppression of extreme risks",
                          "warning above the development 99th percentile", "points table rows shown"),
                 finding = c(hit("Research use only"), hit("not (yet )?externally validated"), hit("leak|ascertain|recorded after"), hit("underpredict|underestimat"),
                             if (length(grep("5-year", h))) "5-year" else "?", "all 9 (collect(): complete = FALSE until every item answered; age 18–60, Hb 3–20 g/dL)",
                             if (length(grep("Not known|unknown", h, ignore.case = TRUE))) "present" else "absent (Yes/No buttons only)",
                             "none (risk computed for any valid input)", hit("p99"), "0–13 (ROWS array); points total can reach 19"))
aud_tab(cf, "V05b_calculator_source", sprintf("Facts read from calculator/standalone/index.html (byte-identical to the deployed page on 2026-10-23, audit A05; identical to calculator/bud_stroke_calculator.html: %s)", live_same))
res_put("V05", list(sparse = sp, calculator = cf))
