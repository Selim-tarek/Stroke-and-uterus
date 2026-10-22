# Shared setup: paths, packages, logging, helpers.
# Run every script from the reanalysis/ directory (run_all.R does this).

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
