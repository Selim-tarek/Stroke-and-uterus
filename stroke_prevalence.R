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
#                    "openxlsx","rlang","brglm2"))
# =============================================================================

library(readxl);   library(janitor); library(lubridate)
library(dplyr);    library(tidyr);   library(stringr)
library(binom);    library(broom);   library(purrr)
library(ggplot2);  library(scales);  library(forcats)
library(openxlsx); library(rlang)

# ------------------------- USER SETTINGS -------------------------------------
input_file <- "stroke.xlsx"     # <<-- path to your workbook
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

age_breaks <- c(-Inf, 39, 49, 59, 69, Inf)
age_labels <- c("<40", "40-49", "50-59", "60-69", "70+")

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

  full_with_labs = c("age_index", "bmi", "htn", "dm", "dyslipidemia", "cad",
                     "afib", "chf", "smoking", "migraine", "vte_history",
                     "thrombophilia", "antithrombotic", "uterine_bleeding",
                     "hormonal_tx", "surgical_tx", "uterine_dx_group",
                     "hgb", "mcv", "platelets", "ldl")
)

primary_adjustment <- "clinical"

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
      fc_raw == "9"                              ~ "Unknown",
      !is.na(fc_num) & fc_num == 0               ~ "None",
      !is.na(fc_num) & fc_num == 1               ~ "Single",
      !is.na(fc_num) & fc_num >= 2 & fc_num <= 4 ~ "2-4",
      !is.na(fc_num) & fc_num >= 5               ~ ">=5",
      str_detect(str_to_lower(fc_raw), "multiple|\\+|many") ~ "Multiple (unspecified)",
      TRUE                                       ~ "Unknown"),
    levels = c("None", "Single", "2-4", ">=5", "Multiple (unspecified)", "Unknown"))

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
range_audit <- map_dfr(intersect(names(plausible), names(df)), function(v) {
  rg <- plausible[[v]]
  bad <- !is.na(df[[v]]) & (df[[v]] < rg[1] | df[[v]] > rg[2])
  df[[v]][bad] <<- NA_real_
  tibble(variable = v, allowed_min = rg[1], allowed_max = rg[2],
         n_set_to_na = sum(bad))
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
  d <- model_df %>% select(stroke_flag, all_of(ex)) %>% drop_na() %>%
    mutate(across(where(is.factor), droplevels))
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
fit_adjusted <- function(exposure, set_name) {
  covs <- setdiff(adj_sets[[set_name]], c(exposure, adj_exclude))
  covs <- intersect(covs, names(model_df))
  d <- model_df %>%
    select(all_of(intersect(c("stroke_flag", exposure, covs), names(model_df)))) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  n_cc <- nrow(d); ev <- if (n_cc) sum(d$stroke_flag == 1) else 0
  ok_exp <- n_cc > 0 && (if (is.factor(d[[exposure]])) nlevels(d[[exposure]]) >= 2
                         else length(unique(d[[exposure]])) >= 2)
  keep <- if (n_cc) usable_covs(d, intersect(covs, names(d))) else character()
  npar <- length(keep) + 1
  epv  <- if (npar) ev / npar else NA_real_

  meta <- tibble(exposure = exposure, adjustment = set_name, n_complete = n_cc,
                 events = ev, n_covariates = length(keep),
                 events_per_variable = round(epv, 1),
                 stable = !is.na(epv) & epv >= min_epv,
                 status = "ok", method = NA_character_)

  if (!ok_exp) { meta$status <- "exposure constant after complete-case"; return(list(meta = meta)) }
  if (ev < 5)  { meta$status <- "fewer than 5 events"; return(list(meta = meta)) }

  r <- safe_glm_fit(as.formula(paste("stroke_flag ~",
                                     paste(c(exposure, keep), collapse = " + "))), d)
  if (is.null(r)) { meta$status <- "did not converge"; return(list(meta = meta)) }
  meta$method <- r$method

  or <- extract_terms(tidy_safe(r$fit), exposure)
  if (!is.null(or))
    or <- or %>% mutate(adjustment = set_name, n_complete = n_cc, events = ev,
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
          or <- bind_rows(or, o2 %>% mutate(adjustment = set_name, n_complete = n_cc,
                                            events = ev, method = r2$method,
                                            stable = meta$stable))
      }
    }
  }
  list(meta = meta, or = or)
}

meta_list <- list(); adj_list <- list()
n_fits <- length(exposures) * length(adj_sets); i_fit <- 0
message("\nFitting ", n_fits, " adjusted models (this is the slow part) ...")

for (ex in exposures) for (sn in names(adj_sets)) {
  i_fit <- i_fit + 1
  message(sprintf("  [%2d/%d] %-22s %s", i_fit, n_fits, ex, sn))
  o <- fit_adjusted(ex, sn)
  meta_list[[paste(ex, sn)]] <- o$meta
  if (!is.null(o$or)) adj_list[[paste(ex, sn)]] <- o$or
}

save_result(bind_rows(meta_list), "09_model_diagnostics")

if (length(adj_list)) {
  adj <- bind_rows(adj_list) %>%
    mutate(reported = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high))
  save_result(adj, "08_adjusted_OR_all_sets")
  for (sn in names(adj_sets))
    save_result(filter(adj, adjustment == sn), paste0("08_adjOR_", sn))

  # Side-by-side unadjusted vs primary adjusted -- the main results table.
  if (nrow(unadj)) {
    main <- unadj %>%
      transmute(exposure, level,
                unadjusted_OR = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high),
                p_unadjusted = p.value) %>%
      full_join(
        adj %>% filter(adjustment == primary_adjustment) %>%
          transmute(exposure, level,
                    adjusted_OR = reported, p_adjusted = p.value,
                    n_model = n_complete, events_model = events, stable),
        by = c("exposure", "level"))
    save_result(main, "10_MAIN_unadj_vs_adj")
  }

  fp <- adj %>% filter(adjustment == primary_adjustment, stable,
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
  save_result(tidy_safe(r$fit) %>% filter(term != "(Intercept)") %>%
                transmute(term, OR = exp(estimate), OR_low = exp(conf.low),
                          OR_high = exp(conf.high), p.value,
                          reported = sprintf("%.2f (%.2f-%.2f)",
                                             exp(estimate), exp(conf.low), exp(conf.high)),
                          n_complete = nrow(d), events = ev, method = r$method),
              paste0("11_multivariable_", sn))
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
