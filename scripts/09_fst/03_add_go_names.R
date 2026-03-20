#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(GO.db)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage:\n",
    "Rscript 03_add_go_names.R <input_go_table.tsv> <output_go_table.tsv>\n"
  )
}

input_file <- args[1]
output_file <- args[2]

cat("Reading GO table...\n")
go <- read_tsv(input_file, show_col_types = FALSE)

if (!("ID" %in% colnames(go))) {
  stop("Input file must contain an 'ID' column.")
}

if (!("Description" %in% colnames(go))) {
  stop("Input file must contain a 'Description' column.")
}

cat("Fetching GO term names from GO.db...\n")

get_go_name <- function(go_id) {
  term <- tryCatch(GO.db::Term(go_id), error = function(e) NA_character_)
  if (length(term) == 0 || is.na(term)) {
    return(NA_character_)
  }
  term
}

go$go_name <- vapply(go$ID, get_go_name, FUN.VALUE = character(1))
go$Description <- ifelse(
  !is.na(go$go_name) & go$go_name != "",
  go$go_name,
  go$Description
)
go$go_name <- NULL

cat("Writing output...\n")
write_tsv(go, output_file)

cat("Done\n")
cat("Output:\n")
cat(output_file, "\n")
