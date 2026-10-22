# Phase 1c: dated comorbidity flags from the Diagnosis extracts (22-31).
# Code lists are the project's prefix rules (ICD-10/ICD-9), ported verbatim from analysis/comorbid_codes.py
# so that the comparison isolates implementation errors; HIC rows by description, file 31 (description only) for pregnancy.
source("R/00_setup.R")
m <- readRDS(file.path(DERIVED, "master.rds"))
link <- m[, .(cn = as.character(mrn), study_id)]

RULES <- list(
  sickle_disease = c("D570", "D571", "D572", "D574", "D578", "2826"),
  sickle_trait = c("D573", "2825"),
  thal_minor = c("D563", "28246", "D564"),
  thal_other = c("D560", "D561", "D562", "D565", "D568", "D569", "2824"),
  other_haemoglobinopathy = c("D580", "D581", "D582", "D588", "D589", "D55", "2827", "2820", "2821", "2822", "2823", "2829"),
  ckd_esrd = c("N186", "N185", "I120", "I1311", "I132", "Z992", "Z940", "Z4901", "Z4902", "Z4931", "Z4932", "5856", "5855",
               "40301", "40311", "40391", "40402", "40403", "40412", "40413", "40492", "40493", "V4511", "V420", "V560", "V568"),
  ckd = c("N18", "I12", "I13", "E0822", "E0922", "E1022", "E1122", "E1322", "D631", "28521", "585", "403", "404"),
  cirrhosis = c("K703", "K704", "K717", "K72", "K74", "K766", "K767", "I85", "5712", "5715", "5716", "5722", "5723", "5724",
                "4560", "4561", "4562", "5728"),
  chronic_liver = c("K70", "K73", "K75", "K760", "B18", "5710", "5711", "5713", "5714", "5718", "5719", "0702", "0703", "07044", "07054"),
  alcohol_disorder = c("F101", "F102", "F1092", "F1093", "F1094", "F1095", "F1096", "F1098", "F1099", "K292", "K852", "K860",
                       "G621", "I426", "G312", "303", "3050", "291", "5353", "4255", "V113"),
  gi_bleed = c("K920", "K921", "K922", "K625", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266", "K270", "K272",
               "K274", "K276", "K280", "K282", "K284", "K286", "K31811", "K5521", "K5701", "K5711", "K5713", "K5721", "K5731",
               "K5733", "K5741", "K5751", "K5753", "K5781", "K5791", "K5793", "5780", "5781", "5789", "5693", "5310", "5312",
               "5314", "5316", "5320", "5322", "5324", "5326", "5330", "5332", "5334", "5336", "5340", "5342", "5344", "5346"),
  ibd = c("K50", "K51", "555", "556"),
  malabsorption = c("K900", "K904", "K909", "5790", "Z9884", "V4586"),
  hiv = c("B20", "Z21", "042", "V08", "O987"),
  menopause = c("Z780", "N951", "E2831", "E2839", "E894", "6272", "6274", "2563", "V4981"),
  pregnancy = c("O", "Z331", "Z34", "Z3A", "Z37", "Z39", "V22", "V23", "V24", "V27", as.character(630:679)))
ANAEMIA_RULES <- list(
  anaemia_pregnancy = c("O990", "O9081", "6482"), anaemia_iron = c("D50", "280"), anaemia_b12_folate = c("D51", "D52", "281"),
  anaemia_nutritional = "D53", anaemia_blood_loss_acute = c("D62", "2851"),
  anaemia_neoplastic_chemo = c("D630", "D6481", "28522", "2853"),
  anaemia_chronic_disease = c("D631", "D638", "28521", "28529"),
  anaemia_haemolytic_aplastic = c("D59", "D60", "D61", "283", "284"),
  anaemia_unspecified = c("D649", "D6489", "D648", "2859", "2858"))
EXCLUDE <- c("7948", "5738", "Z0283", "Z811", "F1090", "F1091", "Z7141", "Z114", "R75", "Z206", "Z2981", "79571", "V016",
             "D75A", "V780", "V781")
HIC_RULES <- c(sickle_disease = "sickle.?cell (disease|anemia|anaemia)|hb.?ss", sickle_trait = "sickle.?cell trait",
  thal_minor = "thalass?emia (minor|trait)", thal_other = "thalass?emia",
  ckd_esrd = "dialysis|end.?stage renal|kidney transplant|renal transplant", ckd = "chronic (renal|kidney)",
  cirrhosis = "cirrhosis|portal hypertension|varices",
  chronic_liver = "fatty,? liver|liver, fatty|steatohepatitis|chronic hepatitis|alcoholic liver",
  alcohol_disorder = "alcoholi(c|sm)|alcohol (abuse|dependence|addiction)|abuse, alcohol",
  gi_bleed = "melena|hematemesis|gastrointestinal (hemorrhage|bleed)|rectal (bleed|hemorrhage)",
  ibd = "crohn|ulcerative colitis", malabsorption = "celiac|sprue|malabsorption|bariatric",
  pregnancy = "deliver|pregnan|postpartum|puerper|cesarean|labor|gestation|abortion|miscarriage",
  hiv = "hiv test positive|human immunodeficiency virus (disease|infection|positive)|\\baids\\b",
  menopause = "^menopause",
  anaemia_iron = "anemia, (iron deficiency|blood loss|microcytic)|iron deficiency, secondary to bleeding",
  anaemia_unspecified = "^anemia")
PREG_INCLUDE <- "pregnan|gestation|deliver|postpartum|post-partum|puerper|prenatal|antepartum|trimester|labou?r|cesarean|caesarean|abortion|miscarriage|ectopic|mother|multigravida|primigravida|obstetric|eclampsia|placent|fetal|fetus|lactation|breastfeed|\\bpreg\\b|chldbrth|childbirth|gravidarum|peripartum|puerp"
PREG_EXCLUDE <- "^personal history|not pregnant|result unknown|result negative|unconfirmed|family planning|avoid pregnancy|^pregnancy test$|^pregnancy examination or test$|contracepti(ve|on) management|^history of"
CURRENT_PREG <- "^(pregnancy|supervision|high.?risk|encounter for supervision)"

files <- Sys.glob(dpath("dx", "MDE_Workflow_Results_for_Diagnosis_*.csv"))
dx <- rbindlist(lapply(files, function(f) {
  x <- fread(f, colClasses = "character"); setnames(x, c("cn", "system", "code", "desc", "date")); x[, file := basename(f)] }))
f31 <- dpath("dx", "_cache", "MDE_Workflow_Results_for_Diagnosis_31.csv")   # cached CSV of the delivered .xlsb
if (file.exists(f31)) {
  x <- fread(f31, colClasses = "character"); setnames(x, c("cn", "desc", "date"))
  dx <- rbind(dx, x[, .(cn, system = "DESC_ONLY", code = "", desc, date, file = "Diagnosis_31 (description only)")])
} else log_msg("Comorbid", "Diagnosis_31 cache not found; pregnancy from description-only extract not available")
dx[, date := as.IDate(substr(date, 1, 10))]
dx <- unique(dx, by = c("cn", "system", "code", "desc", "date"))
n_raw <- nrow(dx)

clean_code <- function(c) gsub("[^A-Z0-9]", "", toupper(c))
first_match <- function(codes, rules) {
  out <- rep(NA_character_, length(codes))
  for (fl in names(rules)) {
    pat <- paste0("^(", paste(rules[[fl]], collapse = "|"), ")")
    hit <- is.na(out) & grepl(pat, codes)
    out[hit] <- fl
  }
  out
}
# ICD rows: one row per individual code
icd <- dx[!system %in% c("HIC", "DESC_ONLY")]
icd <- icd[, .(code1 = clean_code(trimws(unlist(strsplit(code, ","))))), by = .(cn, system, code, desc, date, file)]
uc <- unique(icd$code1)
excl <- grepl(paste0("^(", paste(EXCLUDE, collapse = "|"), ")"), uc)
cls <- data.table(code1 = uc, flag = ifelse(excl, NA, first_match(uc, RULES)), aflag = first_match(uc, ANAEMIA_RULES))
icd <- merge(icd, cls, by = "code1")
icd[flag == "menopause" & grepl("perimenopaus|hot flash|hot flush", desc, ignore.case = TRUE), flag := NA]
long <- rbind(icd[!is.na(flag), .(cn, date, flag)], icd[!is.na(aflag), .(cn, date, flag = aflag)])
# HIC rows
hic <- dx[system == "HIC"]
hic[, flag := NA_character_]
for (fl in names(HIC_RULES)) hic[is.na(flag) & grepl(HIC_RULES[[fl]], tolower(desc)), flag := fl]
long <- rbind(long, hic[!is.na(flag), .(cn, date, flag)])
# description-only pregnancy rows
dsc <- dx[system == "DESC_ONLY"]
dsc[, keep := grepl(PREG_INCLUDE, desc, ignore.case = TRUE) & !grepl(PREG_EXCLUDE, desc, ignore.case = TRUE) &
          !(grepl("personal history|\\bprior\\b", desc, ignore.case = TRUE) & !grepl(CURRENT_PREG, desc, ignore.case = TRUE))]
long <- rbind(long, dsc[keep == TRUE, .(cn, date, flag = "pregnancy")])
long <- merge(unique(long), link, by = "cn")[, .(study_id, flag, date)]
long <- long[!is.na(date) & year(date) >= 1900]
saveRDS(long, file.path(DERIVED, "comorbid_long.rds"))

summ <- long[, .(women_any_date = uniqueN(study_id), rows = .N, first = min(date), last = max(date)), by = flag][order(flag)]
save_tab(summ, "S3f_comorbidity_flags", "Dated comorbidity flags from the diagnosis extracts (all records)")
log_msg("Comorbid", sprintf("%s distinct diagnosis rows; %d flags; ICD codes unmatched to any rule are not counted",
                            format(n_raw, big.mark = ","), nrow(summ)))
