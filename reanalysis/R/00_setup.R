# Shared setup: paths, packages, logging, helpers.
# Run every script from the reanalysis/ directory (run_all.R does this).

invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))   # labels contain en-dashes and \u2265
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(sandwich)
  library(lmtest)
  library(survival)
  library(splines)
  library(mice)
  library(ggplot2)
  library(jsonlite)
  library(openxlsx)
  library(brglm2)
})

DATA_DIR <- Sys.getenv("STROKE_DATA_DIR", "../data")   # relative to reanalysis/
OUT      <- "output"
DERIVED  <- "derived"                                   # patient-level, git-ignored
for (d in c(DERIVED, file.path(OUT, c("tables", "figures", "logs")))) dir.create(d, FALSE, TRUE)

SEED      <- 20261022L
ADMIN_END <- as.IDate("2026-09-25")   # last date covered by the diagnosis extracts
EXTRACT_CAP <- as.IDate("2026-09-27") # latest real encounter date; later dates are scheduled visits
GRADES <- c("None", "Mild", "Moderate", "Severe")
GRADE_LABELS <- c(None = "No anemia (Hb ≥12)", Mild = "Mild (10.0–11.9)",
                  Moderate = "Moderate (8.0–9.9)", Severe = "Severe (<8.0)")

dpath <- function(...) file.path(DATA_DIR, ...)

# ---- logging --------------------------------------------------------------
LOG_FILE <- file.path(OUT, "logs", "run_log.md")
log_msg <- function(section, ...) {
  txt <- paste0("- **", section, "** — ", paste0(...))
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
  message(section, ": ", paste0(...))
}
DEV_FILE <- file.path(OUT, "logs", "deviations.md")
deviation <- function(what, why) {
  cat(sprintf("| %s | %s |\n", what, why), file = DEV_FILE, append = TRUE)
}
WARN_FILE <- file.path(OUT, "logs", "model_warnings.md")
# capture warnings from a model call without suppressing them
with_warnings <- function(expr, label) {
  w <- character()
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  if (length(w)) cat(sprintf("| %s | %s |\n", label, paste(unique(w), collapse = "; ")),
                     file = WARN_FILE, append = TRUE)
  attr(val, "warnings") <- unique(w)
  val
}

# ---- results store (machine-readable) -------------------------------------
RES_FILE <- file.path(OUT, "results.json")
res_put <- function(key, value) {
  r <- if (file.exists(RES_FILE)) fromJSON(RES_FILE, simplifyVector = FALSE) else list()
  r[[key]] <- value
  write_json(r, RES_FILE, auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

save_tab <- function(dt, name, caption = "") {
  fwrite(dt, file.path(OUT, "tables", paste0(name, ".csv")), bom = TRUE)
  attr(dt, "caption") <- caption
  invisible(dt)
}

# ---- formatting -----------------------------------------------------------
f2 <- function(x, d = 2) ifelse(is.na(x), "", formatC(round(x, d), format = "f", digits = d))
fci <- function(est, lo, hi, d = 2) ifelse(is.na(est), "", sprintf("%s (%s–%s)", f2(est, d), f2(lo, d), f2(hi, d)))
fp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
n_pct <- function(n, d) sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / d)

# exp(coef) with robust SE for selected terms
robust_tab <- function(fit, vc, terms = NULL) {
  ct <- coeftest(fit, vcov. = vc)
  ct <- ct[, , drop = FALSE]
  dt <- data.table(term = rownames(ct), b = ct[, 1], se = ct[, 2], p = ct[, 4])
  if (!is.null(terms)) dt <- dt[term %in% terms]
  dt[, `:=`(est = exp(b), lo = exp(b - 1.959964 * se), hi = exp(b + 1.959964 * se))]
  dt[]
}

# Rubin's rules for a list of fits with robust vcov
pool_robust <- function(fits, vcfun, terms) {
  B <- sapply(fits, function(f) coef(f)[terms])
  V <- lapply(fits, function(f) vcfun(f)[terms, terms, drop = FALSE])
  B <- matrix(B, nrow = length(terms))
  m <- ncol(B)
  qbar <- rowMeans(B)
  ubar <- Reduce(`+`, V) / m
  bvar <- if (m > 1) cov(t(B)) else matrix(0, length(terms), length(terms))
  tvar <- ubar + (1 + 1 / m) * bvar
  se <- sqrt(diag(tvar))
  r <- (1 + 1 / m) * diag(bvar) / diag(ubar)
  df <- (m - 1) * (1 + 1 / r)^2
  df[!is.finite(df)] <- 1e6
  p <- 2 * pt(-abs(qbar / se), df)
  tq <- qt(0.975, df)
  data.table(term = terms, b = qbar, se = se, p = p,
             est = exp(qbar), lo = exp(qbar - tq * se), hi = exp(qbar + tq * se))
}

evalue <- function(rr, lo = NA, hi = NA) {          # VanderWeele & Ding, ratio scale
  f <- function(x) { if (is.na(x)) return(NA); if (x < 1) x <- 1 / x; x + sqrt(x * (x - 1)) }
  ci <- if (is.na(lo)) NA else if (lo > 1) lo else if (hi < 1) hi else 1
  c(point = f(rr), ci = if (identical(ci, 1)) 1 else f(ci))
}

bh <- function(p) p.adjust(p, method = "BH")

# ---- descriptive table: n (%) of known values, median (IQR); missing shown separately ----------
VAR_LABELS <- c(age_index = "Age at index, years", race4 = "Race", bmi = "BMI, kg/m²", smoking3 = "Smoking",
  htn = "Hypertension", dm = "Diabetes", dyslipidemia = "Dyslipidemia", cad = "Coronary artery disease (before index)",
  migraine_any = "Migraine (any)", hormonal_type_cat = "Hormonal therapy type", dxgrp = "Uterine diagnosis group",
  uterine_bleeding = "Heavy/abnormal uterine bleeding", thrombophilia = "Thrombophilia (undated)",
  afib = "Atrial fibrillation (undated)", chf = "Heart failure (undated)", vte_history = "VTE (undated)",
  haemoglobinopathy = "Hemoglobinopathy (any date)", ckd_pre = "CKD (before index)", liver_pre = "Chronic liver disease (before index)",
  alcohol_pre = "Alcohol use disorder (before index)", ibd_pre = "IBD (before index)", hiv_pre = "HIV (before index)",
  malabs_pre = "Malabsorption (before index)", gibleed_1y = "GI bleeding (year before index)", preg_1y = "Pregnancy (year before index)",
  anaemia_cause_any = "Any anemia-causing condition", hb = "Baseline Hb, g/dL", mcv_lab = "MCV, fL", rdw_lab = "RDW-CV, %",
  plt_lab = "Platelets, ×10³/µL", hb_tests_2y = "Hb tests in 2 y before index", py = "Follow-up, years",
  ev_stroke = "Incident stroke (excl. TIA)", ev_any = "Incident stroke or TIA", died = "Death without stroke")
describe <- function(d, group, vars) {
  lv <- levels(factor(d[[group]]))
  rows <- list()
  for (v in vars) {
    x <- d[[v]]; lab <- if (!is.na(VAR_LABELS[v])) VAR_LABELS[v] else v
    pfun <- function() tryCatch(if (is.numeric(x) && uniqueN(x, na.rm = TRUE) > 2)
      kruskal.test(x ~ factor(d[[group]]))$p.value else chisq.test(table(x, d[[group]]))$p.value, error = function(e) NA)
    if (is.numeric(x) && uniqueN(x, na.rm = TRUE) > 2) {
      cells <- c(All = sprintf("%s (%s–%s)", f2(median(x, na.rm = TRUE), 1), f2(quantile(x, .25, na.rm = TRUE), 1), f2(quantile(x, .75, na.rm = TRUE), 1)),
                 sapply(lv, function(g) { y <- x[d[[group]] == g]; sprintf("%s (%s–%s)", f2(median(y, na.rm = TRUE), 1), f2(quantile(y, .25, na.rm = TRUE), 1), f2(quantile(y, .75, na.rm = TRUE), 1)) }))
      rows[[length(rows) + 1]] <- data.table(Characteristic = lab, Level = "median (IQR)", t(cells), P = fp(pfun()))
    } else {
      levs <- if (is.factor(x)) levels(x) else if (uniqueN(x, na.rm = TRUE) <= 2) 1 else sort(unique(na.omit(x)))
      p <- fp(pfun())
      for (L in levs) {
        cells <- c(All = n_pct(sum(x == L, na.rm = TRUE), sum(!is.na(x))),
                   sapply(lv, function(g) { y <- x[d[[group]] == g]; n_pct(sum(y == L, na.rm = TRUE), sum(!is.na(y))) }))
        rows[[length(rows) + 1]] <- data.table(Characteristic = lab, Level = if (identical(levs, 1)) "yes" else as.character(L), t(cells), P = p)
        p <- ""
      }
    }
    if (any(is.na(x))) {
      cells <- c(All = n_pct(sum(is.na(x)), length(x)), sapply(lv, function(g) { y <- x[d[[group]] == g]; n_pct(sum(is.na(y)), length(y)) }))
      rows[[length(rows) + 1]] <- data.table(Characteristic = lab, Level = "missing", t(cells), P = "")
    }
  }
  out <- rbindlist(rows, use.names = FALSE)
  setnames(out, c("Characteristic", "Level", "All", lv, "P"))
  out
}
