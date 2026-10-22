# Model functions shared by the revised primary, secondary and sensitivity analyses.
PRIMARY_COVS <- c("ns(age_index, 3)", "race4", "bmi", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "ht",
                  "dxgrp", "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre",
                  "gibleed_1y", "preg_1y")
REDUCED_COVS <- c("ns(age_index, 3)", "race4", "bmi", "htn", "dm", "dyslipidemia", "smoking3", "ckd_pre")
GTERMS <- paste0("grade", GRADES[-1])

prep_covs <- function(d, y) {
  d <- copy(d)
  d[, ibd_hiv_pre := as.integer(ibd_pre == 1 | hiv_pre == 1)]
  ev <- d[!is.na(hormonal_type_cat), .(ev = sum(get(y))), by = hormonal_type_cat]
  small <- setdiff(c(as.character(ev[ev < 5, hormonal_type_cat]),
                     setdiff(levels(d$hormonal_type_cat), as.character(ev$hormonal_type_cat))), c("None", "Other/multiple"))
  h <- as.character(d$hormonal_type_cat); h[h %in% small] <- "Other/multiple"
  d[, ht := factor(h, levels = setdiff(levels(d$hormonal_type_cat), small))]
  attr(d, "ht_merged") <- small
  d
}

# multiple imputation of BMI (pmm) with outcome and Nelson-Aalen hazard in the imputation model
impute_bmi <- function(d, y, m = as.integer(Sys.getenv("N_IMP", "20"))) {
  d <- copy(d)
  sf <- survfit(Surv(d$py, d[[y]]) ~ 1)                     # Nelson-Aalen cumulative hazard at each woman's time
  d[, na_haz := stepfun(sf$time, c(0, sf$cumhaz))(py)]
  vars <- c("bmi", "grade", "age_index", "race4", "smoking3", "htn", "dm", "dyslipidemia", "cad", "migraine_any", "ht", "dxgrp",
            "haemoglobinopathy", "ckd_pre", "liver_pre", "alcohol_pre", "ibd_hiv_pre", "malabs_pre", "gibleed_1y", "preg_1y",
            "na_haz", y)
  sub <- as.data.frame(d[, ..vars])
  meth <- setNames(rep("", ncol(sub)), names(sub)); meth["bmi"] <- "pmm"
  imp <- mice(sub, m = m, method = meth, seed = SEED, printFlag = FALSE)
  if (!is.null(imp$loggedEvents)) {
    le <- as.data.table(imp$loggedEvents)[, .N, by = .(meth, out)]
    cat(sprintf("| mice (%s) | logged events: %s |\n", y, paste(sprintf("%s %s (x%d)", le$meth, le$out, le$N), collapse = "; ")),
        file = WARN_FILE, append = TRUE)
  }
  lapply(seq_len(m), function(i) { z <- copy(d); z[, bmi := complete(imp, i)$bmi]; z })
}

# Poisson rate model (robust SE), pooled over imputations when `dl` is a list
pois_fit <- function(dl, y, expo, covs, label) {
  if (is.data.frame(dl)) dl <- list(dl)
  fml <- reformulate(c(expo, covs, "offset(log(py))"), y)
  fits <- lapply(seq_along(dl), function(i) with_warnings(glm(fml, poisson, data = dl[[i]]), paste0(label, " [imp ", i, "]")))
  list(fits = fits, n = nrow(dl[[1]]), events = sum(dl[[1]][[y]]), df = length(coef(fits[[1]])) - 1,
       converged = all(sapply(fits, `[[`, "converged")))
}
cox_fit <- function(dl, y, expo, covs, label) {
  if (is.data.frame(dl)) dl <- list(dl)
  fml <- as.formula(paste0("Surv(py, ", y, ") ~ ", paste(c(expo, covs), collapse = " + ")))
  fits <- lapply(seq_along(dl), function(i) with_warnings(coxph(fml, data = dl[[i]], ties = "efron"), paste0(label, " [imp ", i, "]")))
  list(fits = fits, n = nrow(dl[[1]]), events = sum(dl[[1]][[y]]))
}
vc_rob <- function(f) vcovHC(f, type = "HC0")
pooled <- function(obj, terms, vcfun = vc_rob) pool_robust(obj$fits, vcfun, terms)

# full pooled coefficient vector and total variance (for linear combinations, e.g. spline contrasts)
pool_full <- function(fits, vcfun) {
  B <- sapply(fits, coef); m <- ncol(B)
  ubar <- Reduce(`+`, lapply(fits, vcfun)) / m
  bvar <- if (m > 1) cov(t(B)) else 0 * ubar
  list(b = rowMeans(B), V = ubar + (1 + 1 / m) * bvar)
}

# standard set of estimates for one outcome/specification: categories + trend (+ any anemia)
grade_models <- function(dl, y, covs, label, cox = TRUE) {
  if (is.data.frame(dl)) dl <- list(dl)
  pf <- pois_fit(dl, y, "grade", covs, label)
  cat_t <- pooled(pf, GTERMS)
  dl2 <- lapply(dl, function(z) { z <- copy(z); z[, score := as.numeric(grade) - 1]; z[, anemia_any := as.integer(grade != "None")]; z })
  tr <- pooled(pois_fit(dl2, y, "score", covs, paste(label, "trend")), "score")
  an <- pooled(pois_fit(dl2, y, "anemia_any", covs, paste(label, "any anemia")), "anemia_any")
  out <- rbind(cat_t[, .(level = sub("grade", "", term), est, lo, hi, p)],
               tr[, .(level = "Per grade (trend)", est, lo, hi, p)],
               an[, .(level = "Any anemia (Hb <12)", est, lo, hi, p)])
  out[, measure := "Rate ratio (Poisson, robust SE)"]
  if (cox) {
    cf <- cox_fit(dl, y, "grade", covs, paste(label, "Cox"))
    ct <- pooled(cf, GTERMS, vcov)
    ctr <- pooled(cox_fit(dl2, y, "score", covs, paste(label, "Cox trend")), "score", vcov)
    out <- rbind(out, rbind(ct[, .(level = sub("grade", "", term), est, lo, hi, p)], ctr[, .(level = "Per grade (trend)", est, lo, hi, p)])[, measure := "Hazard ratio (Cox)"])
  }
  out[, `:=`(analysis = label, outcome = y, N = pf$n, events = pf$events, df = pf$df, epv = pf$events / pf$df,
             converged = pf$converged, imputations = length(dl))]
  out[]
}
