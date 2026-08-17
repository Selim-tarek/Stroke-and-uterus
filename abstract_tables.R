# =============================================================================
# ISC 2027 abstract -- Table 2 gaps and Table 3 (multivariable model).
#
# Reads the cached analysis frame written by stroke_prevalence.R, so run that
# first (STROKE_PREP_ONLY=1 is enough).
#
# Table 3 follows the abstract's predictor list, which differs from the
# pipeline's clinical adjustment set: it adds race and haemoglobin, collapses
# smoking to ever/never, restricts hormone therapy to oestrogen-containing
# preparations, and omits diabetes.
#   Rscript abstract_tables.R
# =============================================================================
suppressMessages({
  library(dplyr); library(purrr); library(stringr); library(openxlsx)
})

cache <- readRDS(file.path("results", "prepared_data.rds"))
d <- cache$model_df

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
pct <- function(k, d) {
  p <- 100 * k / d
  map2_chr(k, p, ~ sprintf("%s (%s%%)", format(.x, big.mark = ","),
                           formatC(.y, format = "f",
                                   digits = if (.y < 1) 1 else 0)))
}

# ---- outcome: ischaemic stroke vs stroke-free -------------------------------
# Other stroke types are a competing outcome, not controls, so they leave the
# analysis rather than joining the denominator.
stype <- as.character(d$stroke_type)
d$isch <- case_when(
  !is.na(d$stroke_flag) & d$stroke_flag == 1 & !is.na(stype) & stype == "Ischaemic stroke" ~ 1L,
  !is.na(d$stroke_flag) & d$stroke_flag == 1                                               ~ NA_integer_,
  !is.na(d$stroke_flag) & d$stroke_flag == 0                                               ~ 0L,
  TRUE                                                                                     ~ NA_integer_)

# =============================================================================
# TABLE 2 -- the rows left blank: uterine pathology among ischaemic strokes
# =============================================================================
n_isch <- sum(d$isch == 1, na.rm = TRUE)

t2 <- d %>%
  filter(!is.na(uterine_dx_group)) %>%
  group_by(Characteristic = as.character(uterine_dx_group)) %>%
  summarise(group_n = n(), isch = sum(isch == 1, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    `n (%)`               = pct(isch, n_isch),
    `Prevalence in group` = sprintf("%.1f%%", 100 * isch / group_n)) %>%
  arrange(match(Characteristic, c("Fibroids only", "Endometriosis only",
                                  "Adenomyosis only", "More than one")))

cat("\n=== Table 2: uterine pathology among the", n_isch, "ischaemic strokes ===\n")
print(as.data.frame(t2), row.names = FALSE)

# =============================================================================
# TABLE 3 -- multivariable logistic regression, abstract's predictor list
# =============================================================================
m <- d %>% filter(!is.na(isch))

m$age10 <- m$age_index / 10

m$race_cat <- factor(
  case_when(
    is.na(m$race_group)                                            ~ NA_character_,
    as.character(m$race_group) == "White"                          ~ "White",
    as.character(m$race_group) == "Black or African American"      ~ "Black/African American",
    as.character(m$race_group) == "Asian"                          ~ "Asian",
    as.character(m$race_group) == "Native Hawaiian/Pacific Islander" ~ "Native Hawaiian or Pacific Islander",
    TRUE                                                           ~ "Other / unknown"),
  levels = c("White", "Black/African American", "Asian",
             "Native Hawaiian or Pacific Islander", "Other / unknown"))

# "History of smoking (current or former)". Unknown is 71% of the cohort and is
# kept as its own level rather than dropped or folded into never.
m$smoking_hx <- factor(
  case_when(
    is.na(m$smoking)                                       ~ NA_character_,
    as.character(m$smoking) == "Never"                     ~ "Never",
    as.character(m$smoking) %in% c("Current", "Former")    ~ "Ever (current or former)",
    TRUE                                                   ~ "Unknown"),
  levels = c("Never", "Ever (current or former)", "Unknown"))

# Oestrogen-containing preparations only: combined OC and menopausal HT.
# Progestin-only, LNG-IUD, GnRH agonists and "Other" are not oestrogen-based.
# The codebook notes hormonal_type undercaptures combined OC, so this is a
# floor on exposure, not a complete ascertainment.
m$estrogen_hrt <- factor(
  case_when(
    is.na(m$hormonal_type)                                                 ~ NA_character_,
    as.character(m$hormonal_type) %in% c("Combined OC", "Menopausal HT")   ~ "Yes",
    TRUE                                                                   ~ "No"),
  levels = c("No", "Yes"))

preds <- c("age10", "race_cat", "htn", "dyslipidemia", "cad", "afib", "chf",
           "smoking_hx", "migraine", "hgb", "estrogen_hrt", "uterine_dx_group")

fit_df <- m %>% select(isch, all_of(preds))

# Sparse levels (n < 10 or fewer than 5 ischaemic strokes) cannot be estimated:
# they produce separation, not a finding. They are set to NA and reported in the
# table as "not estimable" so the reader sees them rather than wondering where
# they went.
sparse <- list()
for (v in preds) {
  if (!is.factor(fit_df[[v]])) next
  tb <- fit_df %>% filter(!is.na(.data[[v]])) %>%
    group_by(.lvl = droplevels(.data[[v]])) %>%
    summarise(n = n(), ev = sum(isch == 1), .groups = "drop")
  bad <- tb %>% filter(n < 10 | ev < 5)
  if (!nrow(bad)) next
  sparse[[length(sparse) + 1]] <- bad %>%
    transmute(variable = v, level = as.character(.lvl), n, events = ev)
  fit_df[[v]][as.character(fit_df[[v]]) %in% as.character(bad$.lvl)] <- NA
  fit_df[[v]] <- droplevels(fit_df[[v]])
}
sparse <- if (length(sparse)) bind_rows(sparse) else tibble()

cc <- fit_df %>% tidyr::drop_na() %>% mutate(across(where(is.factor), droplevels))
fit <- glm(isch ~ ., data = cc, family = binomial("logit"))

td <- broom::tidy(fit, conf.int = TRUE) %>% filter(term != "(Intercept)")
var_of <- function(tm) { h <- preds[map_lgl(preds, ~ startsWith(tm, .x))]
                         if (length(h)) h[which.max(nchar(h))] else NA_character_ }
td$Variable <- map_chr(td$term, var_of)
td$Level <- map2_chr(td$term, td$Variable,
                     ~ if (identical(.x, .y)) "" else substring(.x, nchar(.y) + 1))

lvl_n <- function(v, l) {
  if (!is.factor(cc[[v]])) return(c(nrow(cc), sum(cc$isch == 1)))
  i <- as.character(cc[[v]]) == l
  c(sum(i), sum(i & cc$isch == 1))
}
counts <- map2(td$Variable, td$Level, lvl_n)

table3 <- tibble(
  Variable = td$Variable, Level = td$Level,
  n = map_int(counts, ~ as.integer(.x[1])),
  events = map_int(counts, ~ as.integer(.x[2])),
  aOR = exp(td$estimate), lo = exp(td$conf.low), hi = exp(td$conf.high),
  p = td$p.value) %>%
  mutate(`Adjusted odds ratio` = sprintf("%.2f", aOR),
         `95% CI` = sprintf("%.2f-%.2f", lo, hi),
         `p-value` = fmt_p(p))

cstat <- local({
  p <- fitted(fit); y <- cc$isch
  n1 <- sum(y == 1); n0 <- sum(y == 0); r <- rank(p)
  round((sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 3)
})

cat("\n=== Table 3: n =", nrow(cc), " ischaemic strokes =", sum(cc$isch == 1),
    " C =", cstat, "===\n")
print(as.data.frame(table3 %>% select(Variable, Level, n, events,
                                      `Adjusted odds ratio`, `95% CI`, `p-value`)),
      row.names = FALSE)
if (nrow(sparse)) { cat("\nNot estimable (sparse):\n"); print(as.data.frame(sparse), row.names = FALSE) }

# ---- display labels, and a CSV the document generator reads ----------------
lab <- c(age10 = "Age, per 10 years", race_cat = "Race",
         htn = "Hypertension", dyslipidemia = "Dyslipidemia",
         cad = "Coronary artery disease", afib = "AFib",
         chf = "Congestive heart failure",
         smoking_hx = "History of smoking", migraine = "Migraine",
         hgb = "Hemoglobin, per 1 g/dL",
         estrogen_hrt = "Use of estrogen-containing hormone therapy",
         uterine_dx_group = "Uterine pathology")

# Rows the model could not estimate are carried into the table as explicit
# "not estimable" entries rather than disappearing.
sparse_rows <- if (nrow(sparse)) sparse %>%
  transmute(Variable = lab[variable], Level = level, n = n, events = events,
            `Adjusted odds ratio` = "-", `95% CI` = "-",
            `p-value` = "not estimable") else tibble()

# Reference levels, printed so the reader can see what each OR is against.
# Binary yes/no predictors are shown as a single row instead, so their "No"
# reference is left out here.
ref_rows <- map_dfr(preds, function(v) {
  if (!is.factor(cc[[v]])) return(NULL)
  lv <- levels(cc[[v]])
  if (length(lv) == 2 && all(lv %in% c("No", "Yes"))) return(NULL)
  i <- as.character(cc[[v]]) == lv[1]
  tibble(Variable = lab[[v]], Level = lv[1], n = sum(i), events = sum(i & cc$isch == 1),
         `Adjusted odds ratio` = "1.00", `95% CI` = "reference", `p-value` = "")
})

docx_t3 <- table3 %>%
  transmute(Variable = lab[Variable], Level, n, events,
            `Adjusted odds ratio`, `95% CI`, `p-value`) %>%
  bind_rows(sparse_rows, ref_rows) %>%
  mutate(.is_ref = `95% CI` == "reference") %>%
  arrange(match(Variable, lab), desc(.is_ref)) %>%
  select(-.is_ref)

write.csv(docx_t3, file.path("results", "docx_table3.csv"), row.names = FALSE)
write.csv(t2 %>% transmute(Characteristic, `n (%)`, `Prevalence in group`),
          file.path("results", "docx_table2.csv"), row.names = FALSE)
write.csv(tibble(n = nrow(cc), events = sum(cc$isch == 1), c_statistic = cstat),
          file.path("results", "docx_model.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Table2_gaps");  writeData(wb, "Table2_gaps", t2)
addWorksheet(wb, "Table3");       writeData(wb, "Table3", table3)
addWorksheet(wb, "Table3_model"); writeData(wb, "Table3_model", tibble(
  n = nrow(cc), ischaemic_strokes = sum(cc$isch == 1), c_statistic = cstat,
  AIC = round(AIC(fit), 1), predictors = paste(preds, collapse = ", ")))
if (nrow(sparse)) { addWorksheet(wb, "Table3_not_estimable"); writeData(wb, "Table3_not_estimable", sparse) }
saveWorkbook(wb, file.path("results", "abstract_tables.xlsx"), overwrite = TRUE)
cat("\nWrote results/abstract_tables.xlsx\n")
