# stroke_prev_pipeline_no_mice_exclude_race_malignancy_nomice.R
# No mice. Excludes 'race' and 'malignancy_active' from all adjusted logistic/Poisson models.
# Outputs to "nomice/".
#
# Install missing packages if needed:
# install.packages(c("readxl","janitor","lubridate","dplyr","tidyr","stringr","binom","broom","purrr","sandwich","lmtest","ggplot2","scales","forcats","openxlsx","readr","brglm2"))
#
library(readxl); library(janitor); library(lubridate)
library(dplyr); library(tidyr); library(stringr)
library(binom)
library(broom); library(purrr); library(sandwich); library(lmtest)
library(ggplot2); library(scales); library(forcats)
library(openxlsx); library(readr)

# ------------------------- USER SETTINGS -------------------------
input_file <- "//mfad.mfroot.org/ResearchFL/MPL-Collaboration/Dr.Lin_projects/Non-Cx uterine pathology and stroke/stroke.xlsx"  # <<-- set your actual .xlsx path
sheet_name <- "Data Entry"

outdir <- "nomice"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# Variables lists
adj_covs_default <- c(
  "age_index","race","bmi","fibroids","adenomyosis","endometriosis","uterine_dx_group",
  "fibroid_count","fibroid_max_cm","uterine_bleeding","hormonal_tx","hormonal_type","surgical_tx",
  "htn","dm","dyslipidemia","cad","afib","chf","smoking","migraine","vte_history","thrombophilia",
  "antithrombotic","malignancy_active","malignancy_ever",
  "hgb","mcv","ferritin","platelets","ldl","hdl","trig","hba1c","inr","ddimer_elevated","ddimer"
)

# Categorical & continuous lists (used for prevalence/modeling)
cat_vars <- c("htn","dm","dyslipidemia","cad","afib","chf","smoking","migraine",
              "vte_history","thrombophilia","antithrombotic","malignancy_active","malignancy_ever",
              "ddimer_elevated",
              "fibroids","adenomyosis","endometriosis","uterine_dx_group","race")

cont_vars <- c("hgb","mcv","ferritin","platelets","ldl","hdl","trig","hba1c","inr","ddimer",
               "age_index","bmi","fibroid_max_cm","fibroid_count")

clinical_cutpoints <- list(
  hgb = c(10, 12, 16), mcv = c(80, 100), ferritin = c(30, 300),
  platelets = c(100, 150, 400), ldl = c(70, 100, 130, 160),
  hdl = c(40, 60), trig = c(150, 200, 500), hba1c = c(5.7, 6.5, 8.0), inr = c(0.8, 1.2, 1.5)
)

# Exclude these from adjustment models
adj_exclude <- c("race", "malignancy_active")

age_breaks <- c(-Inf, 39, 49, 59, 69, Inf)
age_labels <- c("<40","40-49","50-59","60-69","70+")

# ------------------------- WORKBOOK HELPERS -------------------------
wb <- createWorkbook(); .used_sheet_names <- character()
make_sheet_name <- function(base){
  s <- gsub("[\\[\\]\\*\\?/\\\\]", "_", base); s <- substr(s,1,31)
  orig <- s; i <- 1
  while (s %in% .used_sheet_names){ s <- substr(paste0(orig,"_",i),1,31); i<-i+1 }
  .used_sheet_names <<- c(.used_sheet_names, s); s
}
save_result <- function(df, base_name){
  sheet <- make_sheet_name(base_name); addWorksheet(wb, sheet)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  writeData(wb, sheet = sheet, x = df); invisible(TRUE)
}
plot_files <- character()

# ------------------------- HELPERS -------------------------
num_missing_codes <- c(9,99); char_missing_codes <- c("9","99","NA","")
recode_missing <- function(x){
  if (is.numeric(x)){ x[x %in% num_missing_codes] <- NA; return(x) }
  x <- as.character(x); x[x %in% char_missing_codes] <- NA; x
}
to_bin01_codebook <- function(x){
  x0 <- recode_missing(x); suppressWarnings(xn <- as.numeric(x0))
  if (!all(is.na(xn))) return(ifelse(is.na(xn), NA_real_, ifelse(xn==1,1, ifelse(xn==0,0, NA_real_))))
  xl <- tolower(as.character(x0))
  ifelse(xl %in% c("1","yes","y","true","t"),1, ifelse(xl %in% c("0","no","n","false","f"),0, NA_real_))
}
recode_smoking <- function(x){
  x0 <- recode_missing(x); suppressWarnings(xn<-as.numeric(x0))
  tibble(smoking_cat = case_when(is.na(xn) ~ NA_character_, xn==0 ~ "Never", xn==1 ~ "Current", xn==2 ~ "Former", TRUE ~ as.character(x0)),
         smoking_current = ifelse(xn==1,1, ifelse(xn %in% c(0,2),0, NA_real_)),
         smoking_former = ifelse(xn==2,1, ifelse(xn %in% c(0,1),0, NA_real_)))
}

# Wilson CI, safe when n == 0 or inputs are missing (binom.confint errors on n = 0).
prev_wilson <- function(events, n){
  events <- as.integer(events); n <- as.integer(n)
  if (is.na(n) || n == 0 || is.na(events)) {
    return(tibble(n = n, events = events, prev = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  ci <- binom.confint(events, n, method = "wilson")
  tibble(n = n, events = events, prev = ci$mean, lower = ci$lower, upper = ci$upper)
}

# Vectorised Wilson for grouped tables; rows with n = 0 come back as NA rather than erroring.
prev_wilson_vec <- function(events, n){
  out <- tibble(prev = rep(NA_real_, length(n)), lower = NA_real_, upper = NA_real_)
  ok <- !is.na(n) & n > 0 & !is.na(events)
  if (any(ok)) {
    ci <- binom.confint(events[ok], n[ok], method = "wilson")
    out$prev[ok] <- ci$mean; out$lower[ok] <- ci$lower; out$upper[ok] <- ci$upper
  }
  out
}

# Safe fitting helpers
safe_glm_fit <- function(formula, data, family = binomial(link = "logit"), try_firth = TRUE){
  fit_try <- try(glm(formula, data=data, family=family, control=glm.control(maxit=50, epsilon=1e-8)), silent=TRUE)
  if (!inherits(fit_try,"try-error") && isTRUE(fit_try$converged)) return(list(fit=fit_try, method="glm"))
  fit_q <- try(glm(formula, data=data, family=quasibinomial(link="logit"), control=glm.control(maxit=50, epsilon=1e-8)), silent=TRUE)
  if (!inherits(fit_q,"try-error") && isTRUE(fit_q$converged)) return(list(fit=fit_q, method="quasibinomial"))
  if (try_firth && requireNamespace("brglm2", quietly=TRUE)){
    # brglmFit is a fitting *method* for glm(), not a formula-level entry point.
    fit_f <- try(glm(formula, data = data, family = binomial(link = "logit"),
                     method = brglm2::brglmFit), silent = TRUE)
    if (!inherits(fit_f,"try-error")) return(list(fit=fit_f, method="firth"))
  }
  warning("safe_glm_fit failed for formula: ", deparse(formula))
  return(NULL)
}
tidy_safe <- function(fit){ td <- try(broom::tidy(fit, conf.int = TRUE), silent=TRUE); if(!inherits(td,"try-error")) return(as_tibble(td)); cm<-try(coef(summary(fit)), silent=TRUE); if(!inherits(cm,"try-error") && !is.null(cm)){ cm_df<-as.data.frame(cm); terms<-rownames(cm_df); est<-cm_df[,1]; se<-if(ncol(cm_df)>=2) cm_df[,2] else sqrt(diag(vcov(fit))); pval<-if(ncol(cm_df)>=4) cm_df[,4] else rep(NA_real_, length(est)); conf_low<-est-1.96*se; conf_high<-est+1.96*se; return(tibble(term=terms, estimate=as.numeric(est), std.error=as.numeric(se), statistic=NA_real_, p.value=as.numeric(pval), conf.low=conf_low, conf.high=conf_high)) }; tibble(term=character(), estimate=double(), std.error=double(), statistic=double(), p.value=double(), conf.low=double(), conf.high=double()) }

# ------------------------- READ + RECODE -------------------------
raw <- read_excel(input_file, sheet = sheet_name, col_names = FALSE, na = c("", "NA")) %>% as_tibble()
header_row <- which(apply(raw,1,function(r) any(grepl("study_id|mrn|dob|eligible|stroke_any", tolower(as.character(r))))))[1]
if (is.na(header_row)) stop("Cannot find header row")
header <- as.character(unlist(raw[header_row,])); colnames(raw) <- janitor::make_clean_names(header)
raw <- raw %>% slice((header_row+1):n())

possible_dates <- c("dob","date_of_birth","index_date","surgical_tx_date","lab_date","abstraction_date","stroke_date")
for (dc in intersect(possible_dates, names(raw))) raw[[dc]] <- lubridate::parse_date_time(as.character(raw[[dc]]), orders = c("ymd","mdy","dmy"), truncated = 3) %>% as.Date()

df <- raw %>%
  mutate(across(intersect(c("age_index","bmi","fibroid_max_cm","hgb","mcv","ferritin","platelets","ldl","hdl","trig","hba1c","inr","ddimer","fibroid_count"), names(.)), ~ suppressWarnings(as.numeric(recode_missing(.))))) %>%
  mutate(
    eligible = suppressWarnings(as.numeric(recode_missing(eligible))),
    fibroids = if ("fibroids" %in% names(.)) to_bin01_codebook(fibroids) else NA_real_,
    uterine_dx_group = if ("uterine_dx_group" %in% names(.)) as.character(uterine_dx_group) else NA_character_,
    smoking = if ("smoking" %in% names(.)) raw$smoking else NA_character_
  ) %>%
  bind_cols(if ("smoking" %in% names(raw)) recode_smoking(raw$smoking) else tibble(smoking_cat = NA_character_, smoking_current = NA_real_, smoking_former = NA_real_)) %>%
  mutate(
    stroke_flag = if ("stroke_any" %in% names(raw)) to_bin01_codebook(raw$stroke_any) else NA_real_,
    stroke_flag_sens = if ("stroke_any_incl_imaging" %in% names(raw)) to_bin01_codebook(raw$stroke_any_incl_imaging) else NA_real_,
    race = if ("race" %in% names(.)) as.character(race) else NA_character_,
    malignancy_active = if ("malignancy_active" %in% names(.)) to_bin01_codebook(malignancy_active) else NA_real_,
    age_group = if ("age_index" %in% names(.)) cut(as.numeric(age_index), breaks = age_breaks, labels = age_labels, right=TRUE) else NA
  )

cat_to_factor <- intersect(c("uterine_dx_group","hormonal_type","surgical_tx","antithrombotic","race","smoking_cat","smoking"), names(df))
for (cname in cat_to_factor) df[[cname]] <- as.factor(df[[cname]])

df_eligible <- df %>% filter(eligible == 1 & (is.na(malignancy_active) | malignancy_active == 0))

# ------------------------- Prevalence & summaries (race/malignancy still included here if present) -------------------------
n_total <- nrow(df_eligible)
n_obs_stroke <- sum(!is.na(df_eligible$stroke_flag))
n_events <- sum(df_eligible$stroke_flag == 1, na.rm = TRUE)
save_result(tibble(total_rows = n_total, observed_outcome_rows = n_obs_stroke, events = n_events), "flow_overview")

# prevalence overall
save_result(prev_wilson(n_events, n_obs_stroke) %>% mutate(prev_pct = prev*100, lower_pct = lower*100, upper_pct = upper*100), "overall_prevalence_primary_complete_case")
save_result(prev_wilson(sum(df_eligible$stroke_flag_sens == 1, na.rm=TRUE), sum(!is.na(df_eligible$stroke_flag_sens))) %>% mutate(prev_pct = prev*100), "overall_prevalence_incl_imaging")
save_result(prev_wilson(n_events, nrow(df_eligible)) %>% mutate(prev_pct = prev*100), "overall_prevalence_assume_missing_no")

# prevalence by categorical vars
for (v in intersect(cat_vars, names(df_eligible))){
  var_sym <- rlang::sym(v)
  tbl <- df_eligible %>% group_by(!!var_sym) %>% summarise(n = n(), observed = sum(!is.na(stroke_flag)), events = sum(stroke_flag==1, na.rm=TRUE), .groups="drop")
  if (nrow(tbl)>0){
    tbl2 <- bind_cols(tbl, prev_wilson_vec(tbl$events, tbl$observed)) %>% mutate(prev_pct = prev*100)
    save_result(tbl2, paste0("prevalence_by_", v))
    p <- ggplot(tbl2, aes(x = fct_inorder(as.character(!!var_sym)), y = prev)) + geom_col(fill="steelblue") + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.25) + scale_y_continuous(labels = percent_format(accuracy=0.1)) + labs(title=paste("Prevalence by", v), x=v, y="Prevalence (95% CI)") + theme_minimal()
    pf <- file.path(outdir, "plots", paste0("prevalence_by_", v, ".png")); ggsave(pf, p, width=7, height=4); plot_files <- c(plot_files, pf)
  }
}

# continuous summaries
cont_summary <- list()
for (v in intersect(cont_vars, names(df_eligible))){
  df_eligible[[v]] <- suppressWarnings(as.numeric(recode_missing(df_eligible[[v]])))
  n_nonmiss <- sum(!is.na(df_eligible[[v]]))
  cont_summary[[v]] <- tibble(variable=v, n_nonmiss = n_nonmiss, n_missing = sum(is.na(df_eligible[[v]])), median = ifelse(n_nonmiss>0, median(df_eligible[[v]], na.rm=TRUE), NA_real_), mean = ifelse(n_nonmiss>0, mean(df_eligible[[v]], na.rm=TRUE), NA_real_), sd = ifelse(n_nonmiss>0, sd(df_eligible[[v]], na.rm=TRUE), NA_real_))
  save_result(cont_summary[[v]], paste0("summary_", v))
}

if (length(cont_summary)>0) save_result(bind_rows(cont_summary), "continuous_summary_all")

# ------------------------- ADJUSTED MODELS (exclude race & malignancy_active) -------------------------
# Prepare adj_covs after excluding race & malignancy_active
adj_covs_filtered <- setdiff(adj_covs_default, adj_exclude)
adj_covs <- intersect(adj_covs_filtered, names(df_eligible))

# record model complete-case counts
model_cc_counts <- list()

# Continuous per-variable adjusted models
cont_results_unit <- list(); cont_results_sd <- list()
# reduced core covariates (exclude race & malignancy_active)
core_covs <- intersect(c("age_index","bmi","htn","dm"), names(df_eligible))

for (v in intersect(cont_vars, names(df_eligible))){
  # core complete-case counts
  vars_core <- intersect(c("stroke_flag", v, core_covs), names(df_eligible))
  ncc_core <- df_eligible %>% select(all_of(vars_core)) %>% drop_na() %>% nrow()
  ev_core  <- if(ncc_core>0) df_eligible %>% select(all_of(vars_core)) %>% drop_na() %>% summarise(sum = sum(stroke_flag==1, na.rm=TRUE)) %>% pull(sum) else 0
  model_cc_counts[[paste0("cont_core_",v)]] <- tibble(variable=v, model="core_covs", n_complete=ncc_core, events=ev_core)

  # full covariates present (without race & malignancy_active)
  model_covs_present <- intersect(setdiff(adj_covs, v), names(df_eligible))
  model_vars_here <- intersect(c("stroke_flag", v, model_covs_present), names(df_eligible))
  df_model <- df_eligible %>% select(all_of(model_vars_here)) %>% mutate(across(all_of(model_covs_present), ~ if(is.character(.) ) as.factor(.) else .)) %>% drop_na()
  ncc_full <- nrow(df_model); ev_full <- if(ncc_full>0) sum(df_model$stroke_flag==1, na.rm=TRUE) else 0
  model_cc_counts[[paste0("cont_full_",v)]] <- tibble(variable=v, model="full_covs", n_complete=ncc_full, events=ev_full)

  if (ncc_full>0){
    f_unit <- as.formula(paste("stroke_flag ~", paste(c(v, model_covs_present), collapse=" + ")))
    res <- safe_glm_fit(f_unit, data=df_model, try_firth=TRUE)
    if (!is.null(res)){
      td <- tidy_safe(res$fit)
      row_v <- which(td$term == v | grepl(paste0("^",v), td$term))
      if (length(row_v)>=1 && !is.na(td$estimate[row_v[1]])){
        est <- td$estimate[row_v[1]]; se <- td$std.error[row_v[1]]; pval <- td$p.value[row_v[1]]
        cont_results_unit[[v]] <- tibble(variable=v, type="per_unit", OR=exp(est), OR_low=exp(est-1.96*se), OR_high=exp(est+1.96*se), p.value=pval, n_complete=ncc_full)
      } else cont_results_unit[[v]] <- tibble(variable=v, type="per_unit", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
    } else cont_results_unit[[v]] <- tibble(variable=v, type="per_unit", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)

    sd_v <- sd(df_model[[v]], na.rm=TRUE)
    if (!is.na(sd_v) && sd_v>0){
      df_model[[paste0(v,"_z")]] <- (df_model[[v]] - mean(df_model[[v]], na.rm=TRUE))/sd_v
      f_sd <- as.formula(paste("stroke_flag ~", paste(c(paste0(v,"_z"), model_covs_present), collapse=" + ")))
      res2 <- safe_glm_fit(f_sd, data=df_model, try_firth=TRUE)
      if (!is.null(res2)){
        td2 <- tidy_safe(res2$fit)
        row_z <- which(td2$term == paste0(v,"_z") | grepl(paste0("scale\\(",v,"\\)"), td2$term))
        if (length(row_z)>=1 && !is.na(td2$estimate[row_z[1]])){
          est2 <- td2$estimate[row_z[1]]; se2 <- td2$std.error[row_z[1]]; pval2 <- td2$p.value[row_z[1]]
          cont_results_sd[[v]] <- tibble(variable=v, type="per_sd", OR=exp(est2), OR_low=exp(est2-1.96*se2), OR_high=exp(est2+1.96*se2), p.value=pval2, n_complete=ncc_full)
        } else cont_results_sd[[v]] <- tibble(variable=v, type="per_sd", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
      } else cont_results_sd[[v]] <- tibble(variable=v, type="per_sd", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
    } else cont_results_sd[[v]] <- tibble(variable=v, type="per_sd", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
  } else {
    cont_results_unit[[v]] <- tibble(variable=v, type="per_unit", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
    cont_results_sd[[v]] <- tibble(variable=v, type="per_sd", OR=NA_real_, OR_low=NA_real_, OR_high=NA_real_, p.value=NA_real_, n_complete=ncc_full)
  }
}

if (length(cont_results_unit)>0) save_result(bind_rows(cont_results_unit), "adjusted_ORs_continuous_per_unit")
if (length(cont_results_sd)>0) save_result(bind_rows(cont_results_sd), "adjusted_ORs_continuous_per_SD")

# Categorical adjusted models (exclude race & malignancy_active from covariates)
cat_model_results <- list()
for (v_raw in cat_vars){
  var_use <- if (v_raw=="smoking" && "smoking_cat" %in% names(df_eligible)) "smoking_cat" else v_raw
  if (! var_use %in% names(df_eligible)){ message("Skipping categorical not found: ", var_use); next }
  df_eligible[[var_use]] <- as.factor(df_eligible[[var_use]])
  # covariates for categorical model: adj_covs but exclude race & malignancy_active and the exposure var itself
  model_covs_present <- intersect(setdiff(adj_covs, c(var_use)), names(df_eligible))
  model_vars_here <- intersect(c("stroke_flag", var_use, model_covs_present), names(df_eligible))
  df_model <- df_eligible %>% select(all_of(model_vars_here)) %>% mutate(across(all_of(model_covs_present), ~ if(is.character(.)) as.factor(.) else .)) %>% drop_na()
  model_cc_counts[[paste0("cat_full_",var_use)]] <- tibble(variable=var_use, model="full_covs_excl_race_malignancy", n_complete=nrow(df_model), events = if(nrow(df_model)>0) sum(df_model$stroke_flag==1, na.rm=TRUE) else 0)
  if (nrow(df_model) < 10) { message("Too few cases for categorical model:", var_use); next }
  # A factor that collapsed to a single level in the complete-case subset cannot be fit.
  if (nlevels(droplevels(df_model[[var_use]])) < 2) { message("Exposure constant after complete-case: ", var_use); next }
  df_model <- df_model %>% mutate(across(where(is.factor), droplevels))
  f_cat <- as.formula(paste("stroke_flag ~", paste(c(var_use, model_covs_present), collapse=" + ")))
  res <- safe_glm_fit(f_cat, data=df_model, try_firth=TRUE)
  if (is.null(res)){ message("Categorical model failed:", var_use); next }
  td <- tidy_safe(res$fit)
  sel <- td %>% filter(term != "(Intercept)" & grepl(paste0("^", var_use), term))
  if (nrow(sel)==0) sel <- td %>% filter(term != "(Intercept)" & grepl(var_use, term))
  if (nrow(sel)==0) { message("No factor terms for ", var_use); next }
  sel2 <- sel %>% mutate(level = sub(paste0("^",var_use),"",term), level = ifelse(level=="", term, level), OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high), variable = var_use) %>%
    select(variable, term, level, estimate, std.error, p.value, OR, OR_low, OR_high)
  cat_model_results[[var_use]] <- sel2
}
if (length(cat_model_results)>0) save_result(bind_rows(cat_model_results), "adjusted_ORs_categorical_models")

# Full multivariable model excluding race & malignancy_active
full_model_covs <- intersect(setdiff(adj_covs_default, adj_exclude), names(df_eligible))
model_vars_full <- intersect(c("stroke_flag", full_model_covs), names(df_eligible))
df_model_full <- df_eligible %>% select(all_of(model_vars_full)) %>% mutate(across(all_of(full_model_covs), ~ if(is.character(.) || is.logical(.)) as.factor(.) else .)) %>% drop_na()
save_result(tibble(full_model_n_complete = nrow(df_model_full), full_model_events = if(nrow(df_model_full)>0) sum(df_model_full$stroke_flag==1, na.rm=TRUE) else 0), "full_model_complete_case_summary")

if (nrow(df_model_full) > 0){
  f_full <- as.formula(paste("stroke_flag ~", paste(setdiff(names(df_model_full), "stroke_flag"), collapse=" + ")))
  res_full <- safe_glm_fit(f_full, data=df_model_full, try_firth=TRUE)
  if (!is.null(res_full)){
    td_full <- tidy_safe(res_full$fit) %>% mutate(OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high))
    save_result(td_full %>% select(term, estimate, std.error, p.value, OR, OR_low, OR_high), "adjusted_ORs_logistic_fullmodel_excl_race_malignancy")
  } else {
    save_result(tibble(message="Full model failed with glm, quasibinomial and Firth"), "adjusted_ORs_logistic_fullmodel_excl_race_malignancy_error")
  }
} else {
  save_result(tibble(message="No complete cases for full model after excluding race/malignancy_active"), "adjusted_ORs_logistic_fullmodel_excl_race_malignancy_error")
}

# Reduced core model (age_index, bmi, htn, dm; excludes race & malignancy_active)
core_covs <- intersect(c("age_index","bmi","htn","dm"), names(df_eligible))
if (length(core_covs)>0){
  model_vars_core <- intersect(c("stroke_flag", core_covs), names(df_eligible))
  df_model_core <- df_eligible %>% select(all_of(model_vars_core)) %>% mutate(across(all_of(core_covs), ~ if(is.character(.) ) as.factor(.) else .)) %>% drop_na()
  save_result(tibble(core_model_n_complete = nrow(df_model_core), core_model_events = if(nrow(df_model_core)>0) sum(df_model_core$stroke_flag==1, na.rm=TRUE) else 0), "core_model_complete_case_summary")
  if (nrow(df_model_core)>0){
    f_core <- as.formula(paste("stroke_flag ~", paste(setdiff(names(df_model_core), "stroke_flag"), collapse=" + ")))
    res_core <- safe_glm_fit(f_core, data=df_model_core, try_firth=TRUE)
    if (!is.null(res_core)){
      td_core <- tidy_safe(res_core$fit) %>% mutate(OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high))
      save_result(td_core %>% select(term, estimate, std.error, p.value, OR, OR_low, OR_high), "adjusted_ORs_logistic_reduced_coremodel_excl_race_malignancy")
    } else {
      save_result(tibble(message="Reduced core model failed"), "adjusted_ORs_logistic_reduced_coremodel_excl_race_malignancy_error")
    }
  } else save_result(tibble(message="No complete cases for reduced core model"), "adjusted_ORs_logistic_reduced_coremodel_excl_race_malignancy_error")
} else save_result(tibble(message="No core covariates present"), "adjusted_ORs_logistic_reduced_coremodel_excl_race_malignancy_error")

# Save model cc counts and audit
if (length(model_cc_counts)>0) save_result(bind_rows(model_cc_counts), "model_complete_case_counts")
to_drop <- intersect(c("mrn","dob"), names(df_eligible)); df_anon <- df_eligible %>% select(-any_of(to_drop)); save_result(df_anon, "analysis_dataset_deidentified")
audit_vars <- intersect(unique(c(cat_vars, cont_vars, adj_covs_default)), names(df_anon))
audit_tbl <- tibble(variable=audit_vars, n_nonmissing = map_int(audit_vars, ~ sum(!is.na(df_anon[[.x]]))), n_missing = map_int(audit_vars, ~ sum(is.na(df_anon[[.x]]))))
save_result(audit_tbl, "audit_nonmissing_counts")

# embed plots if any
if (length(plot_files)>0){
  sheet_plots <- make_sheet_name("Plots"); addWorksheet(wb, sheet_plots); start_row <- 1
  for (pf in plot_files){ insertImage(wb, sheet = sheet_plots, file = pf, startRow = start_row, startCol = 1, width = 18, height = 8, units="cm"); start_row <- start_row + 25 }
}

out_xlsx <- file.path(outdir, "combined_results_no_mice_excl_race_malignancy.xlsx")
if (file.exists(out_xlsx)) file.remove(out_xlsx)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Done. Workbook: ", normalizePath(out_xlsx))
message("Plots saved to: ", normalizePath(file.path(outdir, "plots")))
