# Fallbacks for janitor::make_clean_names() and binom::binom.confint(), for
# environments where those two CRAN packages are unavailable. Sourced by
# stroke_prevalence.R only when the real packages are missing; identical
# behaviour for the single function used from each.

if (!requireNamespace("janitor", quietly = TRUE)) {
  make_clean_names <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x[x == ""] <- "x"
    x <- ifelse(grepl("^[0-9]", x), paste0("x", x), x)
    make.unique(x, sep = "_")
  }
}

if (!requireNamespace("binom", quietly = TRUE)) {
  # Wilson score interval, matching binom.confint(method = "wilson").
  binom.confint <- function(x, n, conf.level = 0.95, method = "wilson") {
    if (!identical(method, "wilson"))
      stop("shim implements the Wilson interval only")
    x <- as.numeric(x); n <- as.numeric(n)
    z  <- qnorm(1 - (1 - conf.level) / 2)
    p  <- x / n
    z2 <- z^2
    centre <- (x + z2 / 2) / (n + z2)
    half   <- (z / (n + z2)) * sqrt(x * (n - x) / n + z2 / 4)
    data.frame(method = method, x = x, n = n, mean = p,
               lower = pmax(0, centre - half),
               upper = pmin(1, centre + half))
  }
}
