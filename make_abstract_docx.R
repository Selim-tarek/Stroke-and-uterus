# =============================================================================
# Writes the ISC 2027 abstract tables into a formatted Word document.
# Numbers come from results/docx_table*.csv (written by abstract_tables.R), so
# nothing here is transcribed by hand.
#   Rscript abstract_tables.R && Rscript make_abstract_docx.R
# Emits abstract_tables.json for the docx generator.
# =============================================================================
suppressMessages(library(jsonlite))

t2 <- read.csv("results/docx_table2.csv", check.names = FALSE)
t3 <- read.csv("results/docx_table3.csv", check.names = FALSE)
mm <- read.csv("results/docx_model.csv", check.names = FALSE)

write_json(list(table2 = t2, table3 = t3, model = mm),
           "results/abstract_tables.json", auto_unbox = TRUE, pretty = TRUE)
cat("Wrote results/abstract_tables.json\n")
