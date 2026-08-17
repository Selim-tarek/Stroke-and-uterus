# =============================================================================
# Prevalence of stroke among patients with non-cancerous uterine pathology
#
# Cross-sectional analysis of the "Data Entry" sheet.
#   - overall + subgroup prevalence with Wilson 95% CIs
#   - unadjusted and adjusted odds ratios for every risk factor
#
# Coding follows the workbook's Codebook tab. Where the data departs from the
# codebook (fibroid_count, race) the departure is handled explicitly below and
# reported in the audit sheets.
#
# install.packages(c("readxl","janitor","lubridate","dplyr","tidyr","stringr",
#                    "binom","broom","purrr","ggplot2","scales","forcats",
#                    "openxlsx","rlang","brglm2","car"))
# =============================================================================

library(readxl);   library(lubridate)
suppressWarnings(suppressMessages({
  if (requireNamespace("janitor", quietly = TRUE)) library(janitor)
  if (requireNamespace("binom",  quietly = TRUE)) library(binom)
}))
source("compat_shims.R")
library(dplyr);    library(tidyr);   library(stringr)
library(broom);    library(purrr)
library(ggplot2);  library(scales);  library(forcats)
library(openxlsx); library(rlang)

# ------------------------- USER SETTINGS -------------------------------------
input_file <- Sys.getenv("STROKE_XLSX", "stroke.xlsx")   # <<-- path to your workbook
sheet_name <- "Data Entry"

outdir <- "results"
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

min_epv <- 10   # events per variable below which a model is flagged as unstable

# =============================================================================
# CODEBOOK VALUE LABELS
# =============================================================================
labs_yesno   <- c("0" = "No", "1" = "Yes", "9" = "Unknown")

code_labels <- list(
  uterine_dx_group = c("1" = "Fibroids only", "2" = "Adenomyosis only",
                       "3" = "Endometriosis only", "4" = "More than one",
                       "9" = "Unknown"),
  hormonal_type    = c("0" = "None", "1" = "Combined OC", "2" = "Progestin-only",
                       "3" = "LNG-IUD", "4" = "GnRH agonist/antagonist",
                       "5" = "Menopausal HT", "6" = "Other", "9" = "Unknown"),
  surgical_tx      = c("0" = "None", "1" = "Myomectomy", "2" = "Hysterectomy",
                       "3" = "Uterine artery embolisation",
                       "4" = "Endometrial ablation", "5" = "Other", "9" = "Unknown"),
  smoking          = c("0" = "Never", "1" = "Current", "2" = "Former", "9" = "Unknown"),
  migraine         = c("0" = "No", "1" = "Without aura", "2" = "With aura", "9" = "Unknown"),
  antithrombotic   = c("0" = "None", "1" = "Antiplatelet", "2" = "Anticoagulant",
                       "3" = "Both", "9" = "Unknown"),
  stroke_type      = c("1" = "Ischaemic stroke", "2" = "Intracerebral haemorrhage",
                       "3" = "Subarachnoid haemorrhage", "4" = "TIA",
                       "5" = "Cerebral venous thrombosis", "9" = "Unknown"),
  stroke_timing    = c("1" = "Before index", "2" = "Same encounter",
                       "3" = "After index", "9" = "Unknown"),
  stroke_etiology  = c("1" = "Cardioembolic", "2" = "Large-artery atherosclerosis",
                       "3" = "Small-vessel occlusion", "4" = "Intracranial atherosclerosis",
                       "5" = "Cryptogenic/ESUS", "6" = "Hypercoagulable state",
                       "7" = "Paradoxical embolus", "8" = "Other determined",
                       "9" = "Unknown"),
  vascular_territory = c("1" = "Anterior circulation", "2" = "Posterior circulation",
                         "3" = "Both", "9" = "Unknown"),
  num_infarct      = c("1" = "One", "2" = "Two", "3" = "Three or more", "9" = "Unknown"),
  exclusion_reason = c("0" = "Not excluded", "1" = "Active malignancy",
                       "2" = "Pathology not confirmed", "3" = "Age <18",
                       "4" = "No uterus documented", "5" = "Insufficient records",
                       "6" = "Other (age >60)")
)

# True 0/1 flags: 9 becomes NA (in this dataset 9 is essentially absent for these).
binary_vars <- c("fibroids", "adenomyosis", "endometriosis", "uterine_bleeding",
                 "hormonal_tx", "htn", "dm", "dyslipidemia", "cad", "afib", "chf",
                 "vte_history", "thrombophilia", "malignancy_active",
                 "malignancy_ever", "ddimer_elevated", "pfo", "cardiac_thrombus",
                 "stroke_confirmed_imaging")

# Multi-level coded fields. "9 = Unknown" is kept as an explicit factor level
# rather than NA: smoking is 73% unknown in this dataset, and converting that to
# NA would delete three quarters of the cohort from every adjusted model.
factor_vars <- c("uterine_dx_group", "hormonal_type", "surgical_tx", "smoking",
                 "migraine", "antithrombotic")

# Genuinely continuous. Blanks are missing; 9 / 99 are REAL values here and must
# never be treated as missing codes (an HDL or platelet count of 99 is ordinary).
cont_vars <- c("age_index", "bmi", "fibroid_max_cm", "hgb", "mcv", "ferritin",
               "platelets", "ldl", "hdl", "trig", "hba1c", "inr", "ddimer")

# Plausibility ranges taken from the Codebook's "permitted values" column.
plausible <- list(
  age_index = c(18, 100), bmi = c(10, 80), fibroid_max_cm = c(0.1, 40),
  hgb = c(3, 20), mcv = c(50, 130), ferritin = c(1, 2000),
  platelets = c(5, 1500), ldl = c(10, 400), hdl = c(5, 150),
  trig = c(10, 2000), hba1c = c(3, 20), inr = c(0.5, 10), ddimer = c(10, 20000)
)

# Bug 4a: eligibility is age 18-60 at index, and the observed eligible range is
# exactly 18-60. The old bands ran to 69, so the 433 patients aged exactly 60
# were binned as "60-69" -- a label implying ineligible patients had leaked into
# the cohort. They had not; only the bin definition was wrong. Bands now stop at
# the eligibility ceiling.
age_breaks <- c(-Inf, 29, 39, 49, 60)
age_labels <- c("18-29", "30-39", "40-49", "50-60")

# ---- Adjustment sets --------------------------------------------------------
# Chosen from the observed missingness in THIS dataset. The clinical set uses
# only variables that are complete or near-complete, so the adjusted models keep
# essentially the whole cohort. The lab set is exploratory: ferritin (58%),
# hba1c (57%), inr (69%) and d-dimer (85%) are missing often enough that
# requiring them all at once leaves a small, selected subsample.
adj_sets <- list(
  core = c("age_index", "bmi"),

  clinical = c("age_index", "bmi", "htn", "dm", "dyslipidemia", "cad", "afib",
               "chf", "smoking", "migraine", "vte_history", "thrombophilia",
               "antithrombotic", "uterine_bleeding", "hormonal_tx",
               "surgical_tx", "uterine_dx_group"),

  # Bug 3: hormonal_tx, surgical_tx and antithrombotic are ascertained at or
  # after index_date, but 44.9% of strokes occurred BEFORE index_date. For those
  # patients these are post-outcome variables, and conditioning on them can
  # induce bias rather than remove it. This variant drops all three.
  clinical_no_treatment = c("age_index", "bmi", "htn", "dm", "dyslipidemia",
                            "cad", "afib", "chf", "smoking", "migraine",
                            "vte_history", "thrombophilia", "uterine_bleeding",
                            "uterine_dx_group"),

  full_with_labs = c("age_index", "bmi", "htn", "dm", "dyslipidemia", "cad",
                     "afib", "chf", "smoking", "migraine", "vte_history",
                     "thrombophilia", "antithrombotic", "uterine_bleeding",
                     "hormonal_tx", "surgical_tx", "uterine_dx_group",
                     "hgb", "mcv", "platelets", "ldl")
)

primary_adjustment <- "clinical"

# ---- Collinear exposure family (Bug 1) --------------------------------------
# uterine_dx_group is a deterministic recoding of the three binary pathology
# flags: knowing fibroids/adenomyosis/endometriosis fixes the group exactly.
# Including the group as a covariate while one of the flags is the exposure puts
# the same information in the model twice. The symptom is a CI that widens while
# the point estimate barely moves -- variance inflation, not confounding
# control. These variables must never share a model.
uterine_flags  <- c("fibroids", "adenomyosis", "endometriosis")
uterine_family <- c(uterine_flags, "uterine_dx_group", "fibroid_count_cat",
                    "fibroid_max_cm")

# Builds the covariate vector for one exposure. The two branches are mutually
# exclusive by construction, and the result never contains the exposure itself.
build_covs <- function(exposure, set_name) {
  covs <- setdiff(adj_sets[[set_name]], adj_exclude)

  if (exposure %in% uterine_family) {
    # Exposure is part of the uterine family -> remove the whole family.
    covs <- setdiff(covs, c("uterine_dx_group", uterine_flags))
  } else {
    # Otherwise keep the group as the single summary of uterine pathology,
    # and never add the individual flags alongside it.
    covs <- setdiff(covs, uterine_flags)
  }

  covs <- setdiff(covs, exposure)

  # Fails loudly rather than silently fitting a model containing its own
  # exposure as a covariate.
  if (exposure %in% covs)
    stop("Exposure '", exposure, "' appears in its own covariate set for '",
         set_name, "'.")
  if (exposure %in% uterine_family &&
      any(c("uterine_dx_group", uterine_flags) %in% covs))
    stop("Collinear uterine terms retained for exposure '", exposure, "'.")

  intersect(covs, names(model_df))
}

# race and malignancy_active never enter an adjusted model:
# race is free-text with dozens of unstandardised levels, and active malignancy
# is an exclusion criterion (so it is constant in the analysis population).
adj_exclude <- c("race", "race_group", "malignancy_active")

# =============================================================================
# HELPERS
# =============================================================================
wb <- createWorkbook(); .used <- character()

make_sheet_name <- function(base) {
  s <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", base), 1, 31)
  o <- s; i <- 1
  while (s %in% .used) { s <- substr(paste0(o, "_", i), 1, 31); i <- i + 1 }
  .used <<- c(.used, s); s
}

save_result <- function(df, name) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  sh <- make_sheet_name(name)
  addWorksheet(wb, sh); writeData(wb, sh, df); freezePane(wb, sh, firstRow = TRUE)
  invisible(TRUE)
}

plot_files <- character()
add_plot <- function(p, file, h = 4) {
  pf <- file.path(outdir, "plots", file)
  ggsave(pf, p, width = 8, height = h, dpi = 150, limitsize = FALSE)
  plot_files <<- c(plot_files, pf); invisible(pf)
}

blanks <- c("", "na", "n/a", "nd", "not done", "unknown", "missing", ".", "-", "null")

clean_chr <- function(x) {
  s <- str_squish(as.character(x))
  s[str_to_lower(s) %in% blanks] <- NA_character_
  s
}

# Continuous parse. Blanks -> NA; 9 and 99 are kept as real values.
to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  s <- clean_chr(x)
  suppressWarnings(as.numeric(str_replace_all(s, "[^0-9.eE+-]", "")))
}

# 0/1 flags: 9 (and 99) mean "unknown" -> NA.
to_bin01 <- function(x) {
  s <- clean_chr(x)
  suppressWarnings(n <- as.numeric(s))
  out <- rep(NA_real_, length(s))
  out[!is.na(n) & n == 1] <- 1
  out[!is.na(n) & n == 0] <- 0
  txt <- str_to_lower(s)
  out[is.na(out) & txt %in% c("yes", "y", "true", "t", "present")] <- 1
  out[is.na(out) & txt %in% c("no", "n", "false", "f", "absent", "none")] <- 0
  out
}

# Coded factor, labelled from the codebook, "9 = Unknown" retained as a level.
to_labelled <- function(x, map, unknown_as_level = TRUE) {
  s <- clean_chr(x)
  lab <- unname(map[s])
  lab[is.na(lab) & !is.na(s)] <- s[is.na(lab) & !is.na(s)]   # unmapped -> keep raw
  if (!unknown_as_level) lab[lab %in% "Unknown"] <- NA_character_
  lvls <- unique(c(unname(map), sort(setdiff(lab, c(unname(map), NA)))))
  factor(lab, levels = lvls[lvls %in% lab])
}

# Excel date cells arrive as serial numbers when the sheet is read as text.
to_date_flex <- function(x) {
  s <- clean_chr(x)
  out <- as.Date(rep(NA_real_, length(s)), origin = "1970-01-01")
  ser <- !is.na(s) & str_detect(s, "^[0-9]{5}(\\.[0-9]+)?$")
  if (any(ser)) out[ser] <- as.Date(as.numeric(s[ser]), origin = "1899-12-30")
  rest <- !is.na(s) & !ser
  if (any(rest))
    out[rest] <- as.Date(parse_date_time(s[rest], orders = c("ymd", "mdy", "dmy"),
                                         truncated = 3))
  out
}

# binom.confint() errors when n = 0, so both wrappers guard for it.
prev_wilson <- function(events, n) {
  events <- as.integer(events); n <- as.integer(n)
  if (is.na(n) || is.na(events) || n == 0)
    return(tibble(n = n, events = events, prev = NA_real_,
                  lower = NA_real_, upper = NA_real_))
  ci <- binom.confint(events, n, method = "wilson")
  tibble(n = n, events = events, prev = ci$mean, lower = ci$lower, upper = ci$upper)
}

prev_wilson_vec <- function(events, n) {
  out <- tibble(prev = rep(NA_real_, length(n)), lower = NA_real_, upper = NA_real_)
  ok <- !is.na(n) & !is.na(events) & n > 0
  if (any(ok)) {
    ci <- binom.confint(events[ok], n[ok], method = "wilson")
    out$prev[ok] <- ci$mean; out$lower[ok] <- ci$lower; out$upper[ok] <- ci$upper
  }
  out
}

# glm -> quasibinomial -> Firth. brglmFit is a fitting METHOD for glm(), not a
# formula-level function, so it has to be passed through method =.
safe_glm_fit <- function(formula, data, try_firth = TRUE) {
  f1 <- try(glm(formula, data = data, family = binomial("logit"),
                control = glm.control(maxit = 50)), silent = TRUE)
  if (!inherits(f1, "try-error") && isTRUE(f1$converged))
    return(list(fit = f1, method = "logistic"))
  if (try_firth && requireNamespace("brglm2", quietly = TRUE)) {
    f2 <- try(glm(formula, data = data, family = binomial("logit"),
                  method = brglm2::brglmFit), silent = TRUE)
    if (!inherits(f2, "try-error")) return(list(fit = f2, method = "firth"))
  }
  warning("Model failed: ", deparse1(formula))
  NULL
}

tidy_safe <- function(fit) {
  td <- try(broom::tidy(fit, conf.int = TRUE), silent = TRUE)
  if (!inherits(td, "try-error")) return(as_tibble(td))
  cm <- try(coef(summary(fit)), silent = TRUE)
  if (inherits(cm, "try-error") || is.null(cm))
    return(tibble(term = character(), estimate = double(), std.error = double(),
                  p.value = double(), conf.low = double(), conf.high = double()))
  d <- as.data.frame(cm); est <- d[[1]]; se <- d[[2]]
  tibble(term = rownames(d), estimate = est, std.error = se,
         p.value = if (ncol(d) >= 4) d[[4]] else NA_real_,
         conf.low = est - 1.96 * se, conf.high = est + 1.96 * se)
}

# ---- Bug 2: sparse factor levels -------------------------------------------
# A level carrying almost no patients or almost no events causes complete or
# quasi-complete separation: the coefficient runs off to +/-Inf and the Wald CI
# explodes (the observed case was OR 0.00 with an upper bound of 2.5e11). Such
# levels are set to NA before fitting so the affected rows drop out, and every
# drop is logged.
min_level_n      <- 10
min_level_events <- 5
rare_level_log <- list()

drop_rare_levels <- function(d, vars, context) {
  for (v in intersect(vars, names(d))) {
    if (!is.factor(d[[v]])) next
    tb <- d %>% filter(!is.na(.data[[v]])) %>%
      group_by(.lvl = droplevels(.data[[v]])) %>%
      summarise(n = n(), events = sum(stroke_flag == 1, na.rm = TRUE),
                .groups = "drop")
    bad <- tb %>% filter(n < min_level_n | events < min_level_events)
    if (!nrow(bad)) next

    rare_level_log[[length(rare_level_log) + 1]] <<- bad %>%
      transmute(context = context, variable = v, level = as.character(.lvl),
                n, events,
                reason = ifelse(n < min_level_n,
                                paste0("n < ", min_level_n),
                                paste0("events < ", min_level_events)))

    d[[v]][as.character(d[[v]]) %in% as.character(bad$.lvl)] <- NA
    d[[v]] <- droplevels(d[[v]])
  }
  d
}

# ---- Bug 1: collinearity diagnostics ---------------------------------------
# Reported for every fitted model so a repeat of the uterine_dx_group problem
# shows up as a number rather than as a quietly inflated confidence interval.
# GVIF^(1/(2*df)) is the scale-invariant form for multi-level factors; squaring
# it puts it back on the usual "VIF > 5 is suspicious" scale. The condition
# number is a whole-design backstop that needs no extra package.
model_collinearity <- function(fit) {
  out <- list(max_vif = NA_real_, max_vif_term = NA_character_,
              condition_number = NA_real_)
  cn <- try({
    mm <- model.matrix(fit)
    keep <- apply(mm, 2, function(z) sd(z) > 0)
    kappa(scale(mm[, keep, drop = FALSE]), exact = FALSE)
  }, silent = TRUE)
  if (!inherits(cn, "try-error")) out$condition_number <- round(as.numeric(cn), 1)

  if (requireNamespace("car", quietly = TRUE)) {
    v <- try(car::vif(fit), silent = TRUE)
    if (!inherits(v, "try-error") && length(v)) {
      vals <- if (is.matrix(v)) v[, ncol(v)]^2 else v   # GVIF^(1/(2df)) -> VIF scale
      if (length(vals)) {
        out$max_vif <- round(max(vals, na.rm = TRUE), 2)
        out$max_vif_term <- names(vals)[which.max(vals)]
      }
    }
  }
  out
}

# Drop covariates that are constant in a given complete-case subset.
usable_covs <- function(data, covs) {
  covs[map_lgl(covs, function(cv) {
    x <- data[[cv]]
    if (is.null(x)) return(FALSE)
    if (is.factor(x)) nlevels(droplevels(x)) >= 2 else length(unique(x[!is.na(x)])) >= 2
  })]
}

# =============================================================================
# READ
# =============================================================================
message("Reading ", input_file, " [", sheet_name, "] ...")

# The sheet's declared dimensions are ~803,000 x 61 because Excel counts
# formatted-but-empty trailing rows, while only ~53,000 rows hold patients.
# Reading the full declared range means ~49 million cells, and scanning every
# row to find the header is slower still -- so the header is located from a
# small probe and the data is then read in a bounded window that grows only if
# it actually fills up.
probe <- read_excel(input_file, sheet = sheet_name, col_names = FALSE,
                    col_types = "text", n_max = 30, .name_repair = "minimal")

hdr_row <- which(apply(probe, 1, function(r)
  any(grepl("study_id|eligible|stroke_any", tolower(as.character(r))))))[1]
if (is.na(hdr_row))
  stop("Could not locate the header row in the first 30 rows of sheet '",
       sheet_name, "'.")

hdr_names <- make_clean_names(as.character(unlist(probe[hdr_row, ])))
message("Header found on row ", hdr_row, " (", length(hdr_names), " columns).")

read_window <- function(n) {
  suppressWarnings(
    read_excel(input_file, sheet = sheet_name, skip = hdr_row, col_names = FALSE,
               col_types = "text", n_max = n, .name_repair = "minimal")
  )
}

n_try <- 100000
repeat {
  raw <- read_window(n_try)
  names(raw) <- hdr_names[seq_len(ncol(raw))]
  raw <- raw %>% filter(!is.na(study_id), str_squish(study_id) != "")
  # If the window came back short of its limit, every patient row was inside it.
  if (nrow(raw) < n_try) break
  n_try <- n_try * 2
  message("Window filled; re-reading with n_max = ", n_try, " ...")
}

message("Patient rows: ", nrow(raw))

# =============================================================================
# RECODE
# =============================================================================
df <- raw

for (v in intersect(cont_vars, names(df)))   df[[v]] <- to_num(df[[v]])
for (v in intersect(binary_vars, names(df))) df[[v]] <- to_bin01(df[[v]])
for (v in intersect(names(code_labels), names(df)))
  df[[v]] <- to_labelled(df[[v]], code_labels[[v]])

for (v in intersect(c("dob", "index_date", "surgical_tx_date", "lab_date",
                      "stroke_date"), names(df)))
  df[[v]] <- to_date_flex(df[[v]])

df$eligible <- to_bin01(df$eligible)
df$stroke_flag      <- to_bin01(df$stroke_any)
df$stroke_flag_sens <- to_bin01(df$stroke_any_incl_imaging)

# ---- fibroid_count ----------------------------------------------------------
# The Codebook defines this as coded (0=None, 1=Single, 2=2-4, 3=>=5, 9=Unknown)
# but the data holds raw counts (0-8) mixed with free text ("2+", "Multiple").
# Treating it as a plain number would read "3" as three fibroids in some rows and
# as the ">=5" category in others, so it is rebuilt as an explicit ordinal and
# the ambiguity is recorded in the audit.
if ("fibroid_count" %in% names(df)) {
  fc_raw <- clean_chr(df$fibroid_count)
  fc_num <- suppressWarnings(as.numeric(fc_raw))
  df$fibroid_count_cat <- factor(
    case_when(
      is.na(fc_raw)                              ~ NA_character_,
      fc_raw == "9"                              ~ NA_character_,   # coded Unknown
      !is.na(fc_num) & fc_num == 0               ~ "None",
      !is.na(fc_num) & fc_num == 1               ~ "Single",
      !is.na(fc_num) & fc_num >= 2 & fc_num <= 4 ~ "2-4",
      !is.na(fc_num) & fc_num >= 5               ~ ">=5",
      str_detect(str_to_lower(fc_raw), "multiple|\\+|many") ~ "Multiple (unspecified)",
      TRUE                                       ~ NA_character_),
    # Bug 2: "Unknown" held 3 patients and produced OR 0.00 with an upper CI of
    # ~2.5e11 -- complete separation, not a finding. An explicitly unknown count
    # carries no information about fibroid burden, so it is missing data rather
    # than a category, and is coded NA.
    levels = c("None", "Single", "2-4", ">=5", "Multiple (unspecified)"))

  save_result(
    tibble(raw_value = fc_raw) %>% count(raw_value, name = "n") %>% arrange(desc(n)),
    "audit_fibroid_count_raw")
}

# ---- race -------------------------------------------------------------------
# Free text in this dataset (42,837 "White", plus "African American",
# "Asian Indian", "Choose Not to Disclose", ...), not the codebook's 1-9 codes.
# Collapsed for reporting only; race does not enter any model.
if ("race" %in% names(df)) {
  r <- str_to_lower(clean_chr(df$race))
  df$race_group <- factor(case_when(
    is.na(r)                                              ~ NA_character_,
    str_detect(r, "white|caucasian")                      ~ "White",
    str_detect(r, "black|african")                        ~ "Black or African American",
    str_detect(r, "asian|filipino|chinese|indian|japanese|korean|vietnamese") ~ "Asian",
    str_detect(r, "american indian|alaska")               ~ "American Indian/Alaska Native",
    str_detect(r, "hawaii|pacific")                       ~ "Native Hawaiian/Pacific Islander",
    str_detect(r, "more than one|multi|two or more")      ~ "More than one",
    str_detect(r, "choose not|declin|refus|unknown|unable") ~ "Unknown / not disclosed",
    TRUE                                                  ~ "Other"),
    levels = c("White", "Black or African American", "Asian",
               "American Indian/Alaska Native", "Native Hawaiian/Pacific Islander",
               "More than one", "Other", "Unknown / not disclosed"))
  save_result(count(tibble(race = clean_chr(df$race)), race, name = "n") %>%
                arrange(desc(n)), "audit_race_raw")
}

if ("age_index" %in% names(df))
  df$age_group <- cut(df$age_index, breaks = age_breaks, labels = age_labels, right = TRUE)

# ---- implausible continuous values -> NA, with an audit trail ---------------
# Bug 4b: this audit counts every row in the file, but 12_missingness counts only
# the analysis population. Both numbers were right and they disagreed because the
# denominators differ -- e.g. all 360 out-of-range ages (345 under 18, 15 over
# 100) belong to INELIGIBLE rows, so age_index is genuinely 0% missing among the
# eligible. Both denominators are now reported side by side.
is_elig <- !is.na(df$eligible) & df$eligible == 1 &
  (is.na(df$malignancy_active) | df$malignancy_active == 0)

range_audit <- map_dfr(intersect(names(plausible), names(df)), function(v) {
  rg <- plausible[[v]]
  bad <- !is.na(df[[v]]) & (df[[v]] < rg[1] | df[[v]] > rg[2])
  df[[v]][bad] <<- NA_real_
  tibble(variable = v, allowed_min = rg[1], allowed_max = rg[2],
         n_set_to_na_all_rows = sum(bad),
         n_set_to_na_in_analysis_population = sum(bad & is_elig))
})
save_result(range_audit, "audit_out_of_range")

# =============================================================================
# ANALYSIS POPULATION
# =============================================================================
n_rows <- nrow(df)
n_elig <- sum(df$eligible == 1, na.rm = TRUE)

df_elig <- df %>% filter(eligible == 1,
                         is.na(malignancy_active) | malignancy_active == 0)

save_result(
  tibble(step = c("patient rows on Data Entry",
                  "eligible = 1",
                  "and no active malignancy (analysis population)",
                  "with a recorded primary outcome",
                  "stroke events (stroke_any = 1)"),
         n = c(n_rows, n_elig, nrow(df_elig),
               sum(!is.na(df_elig$stroke_flag)),
               sum(df_elig$stroke_flag == 1, na.rm = TRUE))),
  "01_flow")

if ("exclusion_reason" %in% names(df))
  save_result(df %>% filter(eligible == 0) %>%
                count(exclusion_reason, name = "n") %>% arrange(desc(n)),
              "02_exclusions")

# =============================================================================
# OVERALL PREVALENCE
# =============================================================================
n_obs <- sum(!is.na(df_elig$stroke_flag))
n_ev  <- sum(df_elig$stroke_flag == 1, na.rm = TRUE)

overall <- bind_rows(
  prev_wilson(n_ev, n_obs) %>%
    mutate(definition = "Primary: stroke_any, patients with a recorded outcome"),
  prev_wilson(sum(df_elig$stroke_flag_sens == 1, na.rm = TRUE),
              sum(!is.na(df_elig$stroke_flag_sens))) %>%
    mutate(definition = "Sensitivity: including imaging-only infarcts"),
  prev_wilson(n_ev, nrow(df_elig)) %>%
    mutate(definition = "Primary, missing outcome counted as no stroke")
) %>%
  mutate(prevalence_pct = 100 * prev, lower_pct = 100 * lower, upper_pct = 100 * upper,
         reported = sprintf("%.2f%% (95%% CI %.2f-%.2f)",
                            prevalence_pct, lower_pct, upper_pct)) %>%
  select(definition, n, events, prevalence_pct, lower_pct, upper_pct, reported)

save_result(overall, "03_overall_prevalence")
message("\n--- Overall prevalence ---"); print(as.data.frame(overall[, c("definition", "n", "events", "reported")]))

# Characteristics of the strokes themselves.
for (v in intersect(c("stroke_type", "stroke_timing", "stroke_etiology",
                      "vascular_territory", "num_infarct",
                      "stroke_confirmed_imaging"), names(df_elig))) {
  tb <- df_elig %>% filter(stroke_flag == 1) %>%
    mutate(level = as.character(.data[[v]])) %>%
    count(level, name = "n") %>%
    mutate(variable = v, pct_of_strokes = round(100 * n / sum(n), 1)) %>%
    select(variable, level, n, pct_of_strokes) %>% arrange(desc(n))
  if (nrow(tb)) save_result(tb, paste0("04_strokes_by_", v))
}

# =============================================================================
# SUBGROUP PREVALENCE
# =============================================================================
subgroup_vars <- intersect(
  c("age_group", "race_group", "uterine_dx_group", "fibroids", "adenomyosis",
    "endometriosis", "fibroid_count_cat", "uterine_bleeding", "hormonal_tx",
    "hormonal_type", "surgical_tx", "htn", "dm", "dyslipidemia", "cad", "afib",
    "chf", "smoking", "migraine", "vte_history", "thrombophilia",
    "antithrombotic", "malignancy_ever", "ddimer_elevated"),
  names(df_elig))

label_level <- function(v, x) {
  if (v %in% binary_vars) return(unname(labs_yesno[as.character(x)]))
  as.character(x)
}

subgroups <- list()
message("\nSubgroup prevalence over ", length(subgroup_vars), " variables ...")
for (v in subgroup_vars) {
  message("  ", v)
  tb <- df_elig %>%
    mutate(.l = label_level(v, .data[[v]])) %>%
    group_by(.l) %>%
    summarise(n = n(),
              observed = sum(!is.na(stroke_flag)),
              events   = sum(stroke_flag == 1, na.rm = TRUE), .groups = "drop") %>%
    arrange(is.na(.l))

  if (!nrow(tb)) next

  # Unadjusted test across observed levels.
  cmp <- df_elig %>% filter(!is.na(.data[[v]]), !is.na(stroke_flag))
  p_un <- NA_real_
  if (nrow(cmp)) {
    tab <- table(droplevels(factor(cmp[[v]])), cmp$stroke_flag)
    if (all(dim(tab) >= 2))
      p_un <- tryCatch({
        ex <- suppressWarnings(chisq.test(tab)$expected)
        if (any(ex < 5)) fisher.test(tab, simulate.p.value = TRUE, B = 2000)$p.value
        else chisq.test(tab)$p.value
      }, error = function(e) NA_real_)
  }

  tb2 <- bind_cols(tb, prev_wilson_vec(tb$events, tb$observed)) %>%
    transmute(variable = v,
              level = ifelse(is.na(.l), "(missing)", .l),
              n, observed, events,
              prevalence_pct = round(100 * prev, 2),
              lower_pct = round(100 * lower, 2),
              upper_pct = round(100 * upper, 2),
              reported = sprintf("%.2f%% (%.2f-%.2f)", 100 * prev, 100 * lower, 100 * upper),
              p_unadjusted = p_un)

  save_result(tb2, paste0("prev_", v))
  subgroups[[v]] <- tb2

  add_plot(
    ggplot(tb2, aes(fct_inorder(level), prevalence_pct / 100)) +
      geom_col(fill = "steelblue") +
      geom_errorbar(aes(ymin = lower_pct / 100, ymax = upper_pct / 100), width = .25) +
      scale_y_continuous(labels = percent_format(accuracy = .1)) +
      labs(title = paste("Stroke prevalence by", v), x = NULL, y = "Prevalence (95% CI)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)),
    paste0("prev_", v, ".png"))
}
if (length(subgroups)) save_result(bind_rows(subgroups), "05_prevalence_all_subgroups")

# Continuous variables by stroke status.
save_result(
  map_dfr(intersect(cont_vars, names(df_elig)), function(v) {
    d <- df_elig %>% filter(!is.na(stroke_flag))
    g <- split(d[[v]], d$stroke_flag)
    tibble(variable = v,
           n_nonmissing = sum(!is.na(df_elig[[v]])),
           pct_missing  = round(100 * mean(is.na(df_elig[[v]])), 1),
           median_no_stroke = median(g[["0"]], na.rm = TRUE),
           median_stroke    = median(g[["1"]], na.rm = TRUE),
           mean_no_stroke   = round(mean(g[["0"]], na.rm = TRUE), 2),
           mean_stroke      = round(mean(g[["1"]], na.rm = TRUE), 2),
           sd_no_stroke     = round(sd(g[["0"]], na.rm = TRUE), 2),
           sd_stroke        = round(sd(g[["1"]], na.rm = TRUE), 2),
           p_wilcoxon = tryCatch(wilcox.test(d[[v]] ~ d$stroke_flag)$p.value,
                                 error = function(e) NA_real_))
  }), "06_continuous_by_stroke")

# =============================================================================
# LOGISTIC REGRESSION -- UNADJUSTED AND ADJUSTED ODDS RATIOS
# =============================================================================
exposures <- setdiff(
  unique(c(intersect(c("uterine_dx_group", "fibroids", "adenomyosis",
                       "endometriosis", "fibroid_count_cat", "uterine_bleeding",
                       "hormonal_tx", "hormonal_type", "surgical_tx", "htn", "dm",
                       "dyslipidemia", "cad", "afib", "chf", "smoking", "migraine",
                       "vte_history", "thrombophilia", "antithrombotic",
                       "malignancy_ever", "ddimer_elevated"), names(df_elig)),
            intersect(cont_vars, names(df_elig)))),
  adj_exclude)

# Binary exposures are relabelled Yes/No so the OR rows read properly.
model_df <- df_elig
for (v in intersect(binary_vars, names(model_df)))
  model_df[[v]] <- factor(unname(labs_yesno[as.character(model_df[[v]])]),
                          levels = c("No", "Yes"))

# Cached analysis frame, so a re-run that only touches the modelling code does
# not have to re-read and re-clean the workbook. Set STROKE_PREP_ONLY=1 to stop
# here after writing it.
saveRDS(list(model_df = model_df, df_elig = df_elig, df = df,
             exposures = exposures),
        file.path(outdir, "prepared_data.rds"))
if (nzchar(Sys.getenv("STROKE_PREP_ONLY"))) {
  message("STROKE_PREP_ONLY set - stopping after data preparation.")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# TABLE 1 -- cohort characteristics, overall and by stroke status.
# Continuous variables as median (IQR) with a Wilcoxon test; categorical as
# n (%) within each column, with a chi-square / Fisher test on the first row.
# ---------------------------------------------------------------------------
t1_vars <- intersect(
  c("age_index", "age_group", "bmi", "race_group", "uterine_dx_group",
    "fibroids", "adenomyosis", "endometriosis", "fibroid_count_cat",
    "uterine_bleeding", "hormonal_tx", "hormonal_type", "surgical_tx",
    "htn", "dm", "dyslipidemia", "cad", "afib", "chf", "smoking", "migraine",
    "vte_history", "thrombophilia", "antithrombotic", "malignancy_ever",
    "hgb", "mcv", "ferritin", "platelets", "ldl", "hdl", "trig", "hba1c",
    "inr", "ddimer_elevated"),
  names(model_df))

t1_base <- model_df %>% filter(!is.na(stroke_flag))
n_no <- sum(t1_base$stroke_flag == 0); n_yes <- sum(t1_base$stroke_flag == 1)

fmt_p1 <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

table1 <- map_dfr(t1_vars, function(v) {
  x <- t1_base[[v]]

  if (is.numeric(x) && !all(na.omit(x) %in% c(0, 1))) {
    med <- function(z) {
      z <- z[!is.na(z)]
      if (!length(z)) return("")
      q <- quantile(z, c(.25, .5, .75))
      sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3])
    }
    p <- tryCatch(wilcox.test(x ~ t1_base$stroke_flag)$p.value,
                  error = function(e) NA_real_)
    tibble(Variable = v, Level = "median (IQR)",
           Overall = med(x),
           `No stroke` = med(x[t1_base$stroke_flag == 0]),
           Stroke = med(x[t1_base$stroke_flag == 1]),
           `Missing (n)` = sum(is.na(x)),
           P = fmt_p1(p))
  } else {
    f <- droplevels(factor(x))
    p <- NA_real_
    tb <- table(f, t1_base$stroke_flag)
    if (all(dim(tb) >= 2))
      p <- tryCatch({
        ex <- suppressWarnings(chisq.test(tb)$expected)
        if (any(ex < 5)) fisher.test(tb, simulate.p.value = TRUE, B = 2000)$p.value
        else chisq.test(tb)$p.value
      }, error = function(e) NA_real_)

    pct <- function(k, d) if (d > 0) sprintf("%d (%.1f%%)", k, 100 * k / d) else ""
    lv <- levels(f)
    tibble(
      Variable = v,
      Level = lv,
      Overall   = map_chr(lv, ~ pct(sum(f == .x, na.rm = TRUE), sum(!is.na(f)))),
      `No stroke` = map_chr(lv, ~ pct(sum(f == .x & t1_base$stroke_flag == 0, na.rm = TRUE), n_no)),
      Stroke    = map_chr(lv, ~ pct(sum(f == .x & t1_base$stroke_flag == 1, na.rm = TRUE), n_yes)),
      `Missing (n)` = c(sum(is.na(f)), rep(NA_integer_, length(lv) - 1)),
      P = c(fmt_p1(p), rep("", length(lv) - 1)))
  }
})

save_result(
  bind_rows(
    tibble(Variable = "N", Level = "", Overall = as.character(nrow(t1_base)),
           `No stroke` = as.character(n_no), Stroke = as.character(n_yes),
           `Missing (n)` = NA_integer_, P = ""),
    table1),
  "TABLE_1_characteristics")

extract_terms <- function(td, exposure, sd_val = NA_real_) {
  sel <- td %>% filter(term != "(Intercept)", startsWith(term, exposure))
  if (!nrow(sel)) return(NULL)
  sel %>% transmute(
    exposure = exposure,
    level = ifelse(term == exposure,
                   if (is.na(sd_val)) "per 1 unit" else sprintf("per 1 SD (%.3g)", sd_val),
                   sub(paste0("^", exposure), "", term)),
    OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high),
    p.value)
}

# ---- unadjusted -------------------------------------------------------------
unadj <- map_dfr(exposures, function(ex) {
  d <- model_df %>% select(stroke_flag, all_of(ex)) %>%
    drop_rare_levels(ex, context = paste0("unadjusted: ", ex)) %>%
    drop_na() %>% mutate(across(where(is.factor), droplevels))
  if (!nrow(d)) return(NULL)
  if (is.factor(d[[ex]]) && nlevels(d[[ex]]) < 2) return(NULL)
  ev <- sum(d$stroke_flag == 1)
  if (ev < 5) return(NULL)
  r <- safe_glm_fit(as.formula(paste("stroke_flag ~", ex)), d)
  if (is.null(r)) return(NULL)
  out <- extract_terms(tidy_safe(r$fit), ex)
  if (is.null(out)) return(NULL)
  out %>% mutate(n_complete = nrow(d), events = ev, method = r$method)
})

if (nrow(unadj))
  save_result(unadj %>% mutate(reported = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high)),
              "07_unadjusted_OR")

# ---- adjusted, one exposure at a time, per adjustment set -------------------
# `population` selects the analysis sample: "all" (primary) or "incident"
# (Bug 3 sensitivity -- strokes recorded before index_date removed).
fit_adjusted <- function(exposure, set_name, data = model_df, population = "all") {
  covs <- build_covs(exposure, set_name)      # Bug 1: mutually exclusive by construction

  d <- data %>%
    select(all_of(intersect(c("stroke_flag", exposure, covs), names(data)))) %>%
    drop_rare_levels(c(exposure, covs),
                     context = paste(population, set_name, exposure, sep = " / ")) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  n_cc <- nrow(d); ev <- if (n_cc) sum(d$stroke_flag == 1) else 0
  ok_exp <- n_cc > 0 && (if (is.factor(d[[exposure]])) nlevels(d[[exposure]]) >= 2
                         else length(unique(d[[exposure]])) >= 2)
  keep <- if (n_cc) usable_covs(d, intersect(covs, names(d))) else character()
  npar <- length(keep) + 1
  epv  <- if (npar) ev / npar else NA_real_

  meta <- tibble(exposure = exposure, adjustment = set_name, population = population,
                 n_complete = n_cc, events = ev, n_covariates = length(keep),
                 covariates = paste(keep, collapse = ", "),
                 events_per_variable = round(epv, 1),
                 stable = !is.na(epv) & epv >= min_epv,
                 status = "ok", method = NA_character_,
                 max_vif = NA_real_, max_vif_term = NA_character_,
                 condition_number = NA_real_)

  if (!ok_exp) { meta$status <- "exposure constant after complete-case"; return(list(meta = meta)) }
  if (ev < 5)  { meta$status <- "fewer than 5 events"; return(list(meta = meta)) }

  r <- safe_glm_fit(as.formula(paste("stroke_flag ~",
                                     paste(c(exposure, keep), collapse = " + "))), d)
  if (is.null(r)) { meta$status <- "did not converge"; return(list(meta = meta)) }
  meta$method <- r$method

  cd <- model_collinearity(r$fit)
  meta$max_vif <- cd$max_vif
  meta$max_vif_term <- cd$max_vif_term
  meta$condition_number <- cd$condition_number
  if (!is.na(cd$max_vif) && cd$max_vif > 5)
    meta$status <- paste0("ok (VIF ", cd$max_vif, " on ", cd$max_vif_term, ")")

  or <- extract_terms(tidy_safe(r$fit), exposure)
  if (!is.null(or))
    or <- or %>% mutate(adjustment = set_name, population = population,
                        n_complete = n_cc, events = ev,
                        method = r$method, stable = meta$stable)

  # Continuous exposures also get a per-SD estimate.
  if (!is.factor(d[[exposure]])) {
    s <- sd(d[[exposure]], na.rm = TRUE)
    if (!is.na(s) && s > 0) {
      d2 <- d; d2[[exposure]] <- (d2[[exposure]] - mean(d2[[exposure]])) / s
      r2 <- safe_glm_fit(as.formula(paste("stroke_flag ~",
                                          paste(c(exposure, keep), collapse = " + "))), d2)
      if (!is.null(r2)) {
        o2 <- extract_terms(tidy_safe(r2$fit), exposure, sd_val = s)
        if (!is.null(o2))
          or <- bind_rows(or, o2 %>% mutate(adjustment = set_name, population = population,
                                            n_complete = n_cc, events = ev,
                                            method = r2$method, stable = meta$stable))
      }
    }
  }
  list(meta = meta, or = or)
}

# ---- Bug 3: incident-stroke sensitivity population --------------------------
# 44.9% of strokes (691/1538) are recorded as occurring BEFORE index_date, so
# treatment variables measured at/after index cannot have preceded them.
# Patients whose stroke predates index are removed; everyone without a stroke
# stays in the denominator. Events with blank or unknown timing (189 + 13) are
# retained -- they are not positively known to precede index -- and that choice
# is recorded in the flow sheet.
incident_df <- model_df %>%
  filter(is.na(stroke_flag) | stroke_flag == 0 |
           is.na(stroke_timing) | as.character(stroke_timing) != "Before index")

save_result(
  tibble(population = c("primary (all eligible)", "incident-stroke sensitivity"),
         n = c(nrow(model_df), nrow(incident_df)),
         strokes = c(sum(model_df$stroke_flag == 1, na.rm = TRUE),
                     sum(incident_df$stroke_flag == 1, na.rm = TRUE)),
         note = c("", "strokes timed 'Before index' removed; unknown/blank timing retained")),
  "13_sensitivity_population")

# ---- run the grid -----------------------------------------------------------
# Primary population gets every adjustment set; the sensitivity population is
# run on the two clinical sets only, which is where the temporality issue bites.
grid <- bind_rows(
  expand_grid(exposure = exposures, set = names(adj_sets), population = "all"),
  expand_grid(exposure = exposures,
              set = c("clinical", "clinical_no_treatment"), population = "incident")
)

meta_list <- list(); adj_list <- list()
message("\nFitting ", nrow(grid), " adjusted models (this is the slow part) ...")

for (i in seq_len(nrow(grid))) {
  ex <- grid$exposure[i]; sn <- grid$set[i]; pp <- grid$population[i]
  message(sprintf("  [%3d/%d] %-22s %-22s %s", i, nrow(grid), ex, sn, pp))
  o <- fit_adjusted(ex, sn, data = if (pp == "incident") incident_df else model_df,
                    population = pp)
  key <- paste(ex, sn, pp)
  meta_list[[key]] <- o$meta
  if (!is.null(o$or)) adj_list[[key]] <- o$or
}

save_result(bind_rows(meta_list), "09_model_diagnostics")
if (length(rare_level_log))
  save_result(bind_rows(rare_level_log) %>% distinct(), "14_dropped_sparse_levels")

if (length(adj_list)) {
  adj <- bind_rows(adj_list) %>%
    mutate(reported = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high))
  save_result(adj, "08_adjusted_OR_all_sets")
  for (sn in names(adj_sets))
    save_result(filter(adj, adjustment == sn, population == "all"),
                paste0("08_adjOR_", sn))
  save_result(filter(adj, population == "incident"), "08_adjOR_incident_sens")

  # ---------------------------------------------------------------------------
  # TABLE 2 -- the single publication table.
  # One row per exposure level: denominator, events, prevalence, unadjusted OR
  # and adjusted OR side by side, with explicit reference rows. Reference levels
  # have no coefficient in the model, so they are added here rather than being
  # silently absent (which is what makes a raw coefficient dump hard to read).
  # ---------------------------------------------------------------------------
  if (nrow(unadj)) {

    descriptive <- map_dfr(exposures, function(ex) {
      d <- model_df %>% filter(!is.na(.data[[ex]]), !is.na(stroke_flag))
      if (!nrow(d)) return(NULL)

      if (is.factor(d[[ex]])) {
        d$.l <- droplevels(d[[ex]])
        tb <- d %>% group_by(.l) %>%
          summarise(n = n(), events = sum(stroke_flag == 1), .groups = "drop")
        ci <- prev_wilson_vec(tb$events, tb$n)
        tibble(exposure = ex, level = as.character(tb$.l),
               n = tb$n, events = tb$events,
               prevalence = sprintf("%.2f%% (%.2f-%.2f)",
                                    100 * ci$prev, 100 * ci$lower, 100 * ci$upper),
               is_ref = as.character(tb$.l) == levels(d$.l)[1])
      } else {
        q <- quantile(d[[ex]], c(.25, .5, .75), na.rm = TRUE)
        tibble(exposure = ex, level = "per 1 unit",
               n = nrow(d), events = sum(d$stroke_flag == 1),
               prevalence = sprintf("median %.1f (IQR %.1f-%.1f)", q[2], q[1], q[3]),
               is_ref = FALSE)
      }
    })

    fmt_or <- function(e, l, u) ifelse(is.na(e), "", sprintf("%.2f (%.2f-%.2f)", e, l, u))
    fmt_p  <- function(p) ifelse(is.na(p), "",
                                 ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

    table2 <- descriptive %>%
      left_join(unadj %>% select(exposure, level, uOR = OR, uL = OR_low,
                                 uH = OR_high, up = p.value),
                by = c("exposure", "level")) %>%
      left_join(adj %>% filter(adjustment == primary_adjustment,
                               population == "all") %>%
                  select(exposure, level, aOR = OR, aL = OR_low, aH = OR_high,
                         ap = p.value, n_model = n_complete, stable),
                by = c("exposure", "level")) %>%
      left_join(adj %>% filter(adjustment == "clinical_no_treatment",
                               population == "incident") %>%
                  select(exposure, level, sOR = OR, sL = OR_low, sH = OR_high),
                by = c("exposure", "level")) %>%
      transmute(
        Variable = exposure,
        Level = level,
        N = n,
        Strokes = events,
        `Prevalence (95% CI)` = prevalence,
        `Unadjusted OR (95% CI)` = ifelse(is_ref, "1.00 (reference)", fmt_or(uOR, uL, uH)),
        `P (unadjusted)` = ifelse(is_ref, "", fmt_p(up)),
        `Adjusted OR (95% CI)` = ifelse(is_ref, "1.00 (reference)", fmt_or(aOR, aL, aH)),
        `P (adjusted)` = ifelse(is_ref, "", fmt_p(ap)),
        # Bug 3 sensitivity: incident strokes only, treatment covariates removed.
        `Sensitivity OR (95% CI)` = ifelse(is_ref, "1.00 (reference)", fmt_or(sOR, sL, sH)),
        `N in adjusted model` = n_model,
        Note = case_when(is_ref ~ "",
                         is.na(stable) ~ "not estimated",
                         !stable ~ paste0("unstable (<", min_epv, " events/variable)"),
                         TRUE ~ ""))

    save_result(table2, "TABLE_2_main_results")
    message("\nTable 2 written: ", nrow(table2), " rows.")
  }

  fp <- adj %>% filter(adjustment == primary_adjustment, population == "all", stable,
                       is.finite(OR), is.finite(OR_low), is.finite(OR_high),
                       OR_high < 50, OR_low > 0.01)
  if (nrow(fp))
    add_plot(
      ggplot(fp %>% mutate(lab = paste0(exposure, ": ", level)),
             aes(OR, fct_rev(fct_inorder(lab)))) +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
        geom_point() +
        geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = .25) +
        scale_x_log10() +
        labs(title = paste0("Adjusted odds of stroke (", primary_adjustment, " model)"),
             x = "Odds ratio (log scale)", y = NULL) +
        theme_minimal(),
      "forest_adjusted.png", h = max(4, .3 * nrow(fp)))
}

# ---- single multivariable model per adjustment set --------------------------
for (sn in names(adj_sets)) {
  covs <- intersect(setdiff(adj_sets[[sn]], adj_exclude), names(model_df))
  d <- model_df %>% select(all_of(c("stroke_flag", covs))) %>% drop_na() %>%
    mutate(across(where(is.factor), droplevels))
  ev <- if (nrow(d)) sum(d$stroke_flag == 1) else 0
  keep <- if (nrow(d)) usable_covs(d, covs) else character()

  save_result(tibble(adjustment = sn, n_complete = nrow(d), events = ev,
                     n_covariates = length(keep),
                     events_per_variable = ifelse(length(keep), round(ev / length(keep), 1), NA)),
              paste0("11_mv_summary_", sn))

  if (!length(keep) || ev < 5) next
  r <- safe_glm_fit(as.formula(paste("stroke_flag ~", paste(keep, collapse = " + "))), d)
  if (is.null(r)) next
  cdm <- model_collinearity(r$fit)
  save_result(tidy_safe(r$fit) %>% filter(term != "(Intercept)") %>%
                transmute(term, OR = exp(estimate), OR_low = exp(conf.low),
                          OR_high = exp(conf.high), p.value,
                          reported = sprintf("%.2f (%.2f-%.2f)",
                                             exp(estimate), exp(conf.low), exp(conf.high)),
                          n_complete = nrow(d), events = ev, method = r$method,
                          max_vif = cdm$max_vif, condition_number = cdm$condition_number),
              paste0("11_multivariable_", sn))
}

# =============================================================================
# BEFORE / AFTER -- effect of the Bug 1 fix on the three pathology exposures
# Refits the OLD specification (uterine_dx_group retained as a covariate while
# a pathology flag is the exposure) next to the corrected one, so the variance
# inflation is visible rather than asserted.
# =============================================================================
message("\nBefore/after comparison for the pathology exposures ...")

fit_spec <- function(exposure, set_name, covs) {
  covs <- setdiff(intersect(covs, names(model_df)), exposure)
  d <- model_df %>%
    select(all_of(c("stroke_flag", exposure, covs))) %>%
    drop_na() %>% mutate(across(where(is.factor), droplevels))
  if (!nrow(d)) return(NULL)
  keep <- usable_covs(d, covs)
  r <- safe_glm_fit(as.formula(paste("stroke_flag ~",
                                     paste(c(exposure, keep), collapse = " + "))), d)
  if (is.null(r)) return(NULL)
  o <- extract_terms(tidy_safe(r$fit), exposure)
  if (is.null(o) || !nrow(o)) return(NULL)
  cd <- model_collinearity(r$fit)
  o %>% slice(n()) %>%
    transmute(exposure, level, OR, OR_low, OR_high,
              ci_width = OR_high - OR_low,
              n = nrow(d), events = sum(d$stroke_flag == 1),
              max_vif = cd$max_vif, condition_number = cd$condition_number)
}

ba_total <- length(uterine_flags) * length(adj_sets); ba_i <- 0
message("  (", ba_total, " paired refits, two models each)")

before_after <- map_dfr(uterine_flags, function(ex) {
  map_dfr(names(adj_sets), function(sn) {
    ba_i <<- ba_i + 1
    message(sprintf("  [%2d/%d] %-16s %s", ba_i, ba_total, ex, sn))
    old_covs <- setdiff(adj_sets[[sn]], adj_exclude)          # uterine_dx_group kept
    new_covs <- build_covs(ex, sn)                            # collinear family removed
    b <- fit_spec(ex, sn, old_covs); a <- fit_spec(ex, sn, new_covs)
    if (is.null(b) && is.null(a)) return(NULL)
    tibble(
      exposure = ex, adjustment = sn,
      n = if (!is.null(a)) a$n else NA_integer_,
      events = if (!is.null(a)) a$events else NA_integer_,
      before_OR = if (!is.null(b)) sprintf("%.2f (%.2f-%.2f)", b$OR, b$OR_low, b$OR_high) else "",
      after_OR  = if (!is.null(a)) sprintf("%.2f (%.2f-%.2f)", a$OR, a$OR_low, a$OR_high) else "",
      before_ci_width = if (!is.null(b)) round(b$ci_width, 3) else NA_real_,
      after_ci_width  = if (!is.null(a)) round(a$ci_width, 3) else NA_real_,
      ci_narrowed_by_pct = if (!is.null(a) && !is.null(b) && b$ci_width > 0)
        round(100 * (b$ci_width - a$ci_width) / b$ci_width, 1) else NA_real_,
      before_max_vif = if (!is.null(b)) b$max_vif else NA_real_,
      after_max_vif  = if (!is.null(a)) a$max_vif else NA_real_)
  })
})

if (nrow(before_after)) {
  save_result(before_after, "15_before_after_collinearity_fix")
  message("\n--- Bug 1 before/after (uterine_dx_group removed from covariates) ---")
  print(as.data.frame(before_after %>%
    select(exposure, adjustment, n, events, before_OR, after_OR, ci_narrowed_by_pct)),
    row.names = FALSE)
}

# =============================================================================
# DEDICATED MODEL -- multivariable logistic regression for ISCHAEMIC stroke
# among women aged 18-60 with non-cancerous uterine pathology.
#
# Separate from the exposure-at-a-time grid above: this is a single model in
# which every covariate is entered simultaneously, and an OR is reported for
# each of them rather than for one exposure of interest.
#
# Outcome definition (primary): stroke_any = 1 AND stroke_type = "Ischaemic
# stroke". Patients whose only stroke was haemorrhagic, SAH, TIA or CVT are
# a competing outcome, not controls, so they are REMOVED from the primary
# analysis; a sensitivity model keeps them as non-cases. Strokes with an
# unknown/blank type are also removed from the primary model (they cannot be
# classified) and counted in the flow table.
#
# Uterine pathology enters as the three binary flags. uterine_dx_group is a
# deterministic recoding of those flags, so it is excluded here for the same
# collinearity reason enforced by build_covs().
# =============================================================================
message("\nFitting dedicated ischaemic-stroke multivariable model ...")

isch_type_label <- "Ischaemic stroke"

isch_covs_full <- intersect(
  c("age_index", "bmi", "fibroids", "adenomyosis", "endometriosis",
    "uterine_bleeding", "htn", "dm", "dyslipidemia", "cad", "afib", "chf",
    "smoking", "migraine", "vte_history", "thrombophilia", "antithrombotic",
    "hormonal_tx", "surgical_tx"),
  names(model_df))

# Temporality-safe variant: drops the three variables ascertained at or after
# index_date (same rationale as clinical_no_treatment above).
isch_covs_notx <- setdiff(isch_covs_full,
                          c("antithrombotic", "hormonal_tx", "surgical_tx"))

# ---- build the outcome ------------------------------------------------------
mk_isch <- function(d, competing = c("exclude", "control")) {
  competing <- match.arg(competing)
  stype <- if ("stroke_type" %in% names(d)) as.character(d$stroke_type) else NA_character_
  is_case  <- !is.na(d$stroke_flag) & d$stroke_flag == 1 &
    !is.na(stype) & stype == isch_type_label
  is_other <- !is.na(d$stroke_flag) & d$stroke_flag == 1 & !is_case
  d$stroke_flag <- ifelse(is_case, 1L,
                          ifelse(is_other,
                                 if (competing == "exclude") NA_integer_ else 0L,
                                 ifelse(!is.na(d$stroke_flag) & d$stroke_flag == 0,
                                        0L, NA_integer_)))
  d[!is.na(d$stroke_flag), , drop = FALSE]
}

# Age band is an eligibility criterion, but it is re-applied explicitly so the
# model population matches the stated 18-60 definition even if a stray row slips
# past the eligible flag.
age_ok <- !is.na(model_df$age_index) & model_df$age_index >= 18 &
  model_df$age_index <= 60
isch_base <- model_df[age_ok, , drop = FALSE]

# Concordance (C-statistic / AUC) from the rank of fitted probabilities. Ties
# contribute a half, which is the Mann-Whitney definition.
c_stat <- function(y, p) {
  y <- as.integer(y); ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p)
  round((sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 3)
}

isch_flow <- tibble(
  step = c("Eligible, no active malignancy",
           "Aged 18-60 at index",
           "Any stroke",
           "  ischaemic (cases)",
           "  other stroke type (competing, removed from primary)",
           "  unknown/blank stroke type (removed from primary)",
           "Primary model population (cases + stroke-free controls)"),
  n = c(
    nrow(model_df),
    nrow(isch_base),
    sum(isch_base$stroke_flag == 1, na.rm = TRUE),
    sum(isch_base$stroke_flag == 1 &
          !is.na(isch_base$stroke_type) &
          as.character(isch_base$stroke_type) == isch_type_label, na.rm = TRUE),
    sum(isch_base$stroke_flag == 1 &
          !is.na(isch_base$stroke_type) &
          as.character(isch_base$stroke_type) != isch_type_label, na.rm = TRUE),
    sum(isch_base$stroke_flag == 1 & is.na(isch_base$stroke_type), na.rm = TRUE),
    nrow(mk_isch(isch_base, "exclude"))))
save_result(isch_flow, "16_ischemic_flow")

# ---- fit one specification --------------------------------------------------
fit_isch <- function(data, covs, label) {
  d <- data %>%
    select(all_of(intersect(c("stroke_flag", covs), names(data)))) %>%
    drop_rare_levels(covs, context = paste0("ischaemic model: ", label)) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  keep <- usable_covs(d, intersect(covs, names(d)))
  if (!length(keep) || !nrow(d)) return(NULL)

  ev <- sum(d$stroke_flag == 1)
  npar <- sum(map_int(keep, function(v)
    if (is.factor(d[[v]])) nlevels(d[[v]]) - 1L else 1L))

  r <- safe_glm_fit(as.formula(paste("stroke_flag ~",
                                     paste(keep, collapse = " + "))), d)
  if (is.null(r)) return(NULL)

  td <- tidy_safe(r$fit) %>% filter(term != "(Intercept)")
  cd <- model_collinearity(r$fit)

  # Map each coefficient back to its variable and level.
  var_of <- function(tm) {
    hit <- keep[map_lgl(keep, ~ startsWith(tm, .x))]
    if (!length(hit)) return(NA_character_)
    hit[which.max(nchar(hit))]
  }
  td$Variable <- map_chr(td$term, var_of)
  td$Level <- map2_chr(td$term, td$Variable, function(tm, v) {
    if (is.na(v)) return(tm)
    if (identical(tm, v)) "per 1 unit" else substring(tm, nchar(v) + 1)
  })

  # Reference level for each factor, so the table reads like a paper table.
  refs <- map_dfr(keep, function(v) {
    if (!is.factor(d[[v]])) return(NULL)
    tibble(Variable = v, Level = levels(d[[v]])[1], is_ref = TRUE)
  })

  res <- td %>%
    transmute(Variable, Level,
              OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high),
              p.value, is_ref = FALSE) %>%
    bind_rows(refs %>% mutate(OR = 1, OR_low = NA_real_, OR_high = NA_real_,
                              p.value = NA_real_)) %>%
    arrange(match(Variable, keep), desc(is_ref))

  # n and events behind each level, for the reader who wants the denominator.
  lvl_n <- function(v, l) {
    if (!is.factor(d[[v]])) return(c(nrow(d), sum(d$stroke_flag == 1)))
    idx <- as.character(d[[v]]) == l
    c(sum(idx), sum(idx & d$stroke_flag == 1))
  }
  cnt <- map2(res$Variable, res$Level, lvl_n)
  res$n <- map_int(cnt, ~ as.integer(.x[1]))
  res$events <- map_int(cnt, ~ as.integer(.x[2]))

  list(
    table = res %>% transmute(
      Model = label, Variable, Level, n, events,
      `Adjusted OR (95% CI)` = ifelse(is_ref, "1.00 (reference)",
                                      sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high)),
      OR, OR_low, OR_high,
      P = fmt_p1(p.value)),
    meta = tibble(
      Model = label, n = nrow(d), events = ev,
      prevalence = sprintf("%.2f%%", 100 * ev / nrow(d)),
      n_covariates = length(keep), n_parameters = npar,
      events_per_variable = round(ev / npar, 1),
      stable = ev / npar >= min_epv,
      method = r$method,
      c_statistic = c_stat(d$stroke_flag, fitted(r$fit)),
      AIC = round(AIC(r$fit), 1),
      max_vif = cd$max_vif, max_vif_term = cd$max_vif_term,
      condition_number = cd$condition_number,
      covariates = paste(keep, collapse = ", ")),
    fit = r$fit, data = d, keep = keep)
}

# ---- unadjusted counterparts, same population -------------------------------
isch_unadj <- function(data, covs, label) {
  map_dfr(covs, function(v) {
    d <- data %>% select(all_of(c("stroke_flag", v))) %>%
      drop_rare_levels(v, context = paste0("ischaemic unadjusted: ", v)) %>%
      drop_na() %>% mutate(across(where(is.factor), droplevels))
    if (!nrow(d) || !length(usable_covs(d, v))) return(NULL)
    if (sum(d$stroke_flag == 1) < 5) return(NULL)
    r <- safe_glm_fit(as.formula(paste("stroke_flag ~", v)), d)
    if (is.null(r)) return(NULL)
    tidy_safe(r$fit) %>% filter(term != "(Intercept)") %>%
      transmute(Variable = v,
                Level = ifelse(term == v, "per 1 unit", substring(term, nchar(v) + 1)),
                `Unadjusted OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)",
                                                   exp(estimate), exp(conf.low), exp(conf.high)),
                `P (unadjusted)` = fmt_p1(p.value))
  })
}

isch_primary_df <- mk_isch(isch_base, "exclude")

m_main <- fit_isch(isch_primary_df, isch_covs_full, "Primary (all covariates)")
m_notx <- fit_isch(isch_primary_df, isch_covs_notx,
                   "Sensitivity: treatment covariates dropped")
m_ctrl <- fit_isch(mk_isch(isch_base, "control"), isch_covs_full,
                   "Sensitivity: other stroke types as non-cases")

isch_incident <- mk_isch(
  isch_base %>% filter(is.na(stroke_flag) | stroke_flag == 0 |
                         is.na(stroke_timing) |
                         as.character(stroke_timing) != "Before index"),
  "exclude")
m_inc <- fit_isch(isch_incident, isch_covs_notx,
                  "Sensitivity: incident strokes only, no treatment covariates")

if (!is.null(m_main)) {
  ua <- isch_unadj(isch_primary_df, m_main$keep, "primary")
  tab3 <- m_main$table %>%
    left_join(ua, by = c("Variable", "Level")) %>%
    mutate(`Unadjusted OR (95% CI)` = ifelse(`Adjusted OR (95% CI)` == "1.00 (reference)",
                                             "1.00 (reference)",
                                             coalesce(`Unadjusted OR (95% CI)`, "")),
           `P (unadjusted)` = coalesce(`P (unadjusted)`, "")) %>%
    select(Variable, Level, n, events,
           `Unadjusted OR (95% CI)`, `P (unadjusted)`,
           `Adjusted OR (95% CI)`, `P (adjusted)` = P)

  hdr <- tibble(Variable = "MODEL", Level = m_main$meta$Model,
                n = m_main$meta$n, events = m_main$meta$events,
                `Unadjusted OR (95% CI)` = "", `P (unadjusted)` = "",
                `Adjusted OR (95% CI)` = sprintf("C-statistic %.3f; EPV %.1f",
                                                 m_main$meta$c_statistic,
                                                 m_main$meta$events_per_variable),
                `P (adjusted)` = "")

  save_result(bind_rows(hdr, tab3), "TABLE_3_ischemic_multivariable")

  message("\n--- Ischaemic stroke, multivariable model (women 18-60) ---")
  message("n = ", m_main$meta$n, ", ischaemic strokes = ", m_main$meta$events,
          " (", m_main$meta$prevalence, "), C = ", m_main$meta$c_statistic)
  print(as.data.frame(tab3 %>% select(Variable, Level, n, events,
                                      `Adjusted OR (95% CI)`, `P (adjusted)`)),
        row.names = FALSE)
}

isch_all <- compact(list(m_main, m_notx, m_ctrl, m_inc))
if (length(isch_all)) {
  save_result(map_dfr(isch_all, ~ .x$table), "16_ischemic_all_models")
  save_result(map_dfr(isch_all, ~ .x$meta),  "16_ischemic_diagnostics")
}

# =============================================================================
# AUDIT + SAVE
# =============================================================================
df_out <- df_elig %>% select(-any_of(c("mrn", "dob", "notes")))

save_result(tibble(variable = names(df_out),
                   n_nonmissing = map_int(names(df_out), ~ sum(!is.na(df_out[[.x]]))),
                   n_missing    = map_int(names(df_out), ~ sum(is.na(df_out[[.x]]))),
                   pct_missing  = round(100 * map_dbl(names(df_out),
                                                      ~ mean(is.na(df_out[[.x]]))), 1)),
            "12_missingness")

if (length(plot_files)) {
  sp <- make_sheet_name("Plots"); addWorksheet(wb, sp); r <- 1
  for (pf in plot_files) {
    insertImage(wb, sp, pf, startRow = r, startCol = 1, width = 20, height = 10, units = "cm")
    r <- r + 22
  }
}

out_xlsx <- file.path(outdir, "stroke_prevalence_results.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("\nDone.")
message("Workbook: ", normalizePath(out_xlsx))
message("Plots:    ", normalizePath(file.path(outdir, "plots")))
message("\nMain results table = sheet '10_MAIN_unadj_vs_adj'.")
message("Check '09_model_diagnostics' -- the 'stable' column flags models with <",
        min_epv, " events per variable.")
