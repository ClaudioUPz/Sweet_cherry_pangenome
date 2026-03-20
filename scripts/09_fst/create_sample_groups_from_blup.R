#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript create_sample_groups_from_blup.R <BLUP.xlsx> <output.tsv>\n"
  )
}

blup_path <- args[1]
out_tsv   <- args[2]

message("Reading BLUP spreadsheet from: ", blup_path)

# Read the full spreadsheet without predefined column names
blup_raw <- read_excel(blup_path, col_names = FALSE)

message("Raw spreadsheet dimensions: ", nrow(blup_raw), " rows x ", ncol(blup_raw), " columns")

# Detect the header row containing 'Accession code'
header_row_candidates <- which(
  apply(blup_raw, 1, function(r) any(r == "Accession code"))
)

if (length(header_row_candidates) == 0) {
  stop("No row containing 'Accession code' was found in the spreadsheet.")
}

header_row <- header_row_candidates[1]
message("Detected header row: ", header_row)

# Extract header values
header_vec <- as.character(blup_raw[header_row, ])

# Replace missing header names with placeholder names
header_vec[is.na(header_vec)] <- paste0("COL", which(is.na(header_vec)))

# Extract data below the detected header row
if (header_row == nrow(blup_raw)) {
  stop("The detected header row is the last row of the spreadsheet; no data found below it.")
}

blup <- blup_raw[(header_row + 1):nrow(blup_raw), ]
colnames(blup) <- header_vec

message("Detected columns:")
print(colnames(blup)[1:min(20, ncol(blup))])

# Check required columns
if (!all(c("Accession code", "Status") %in% colnames(blup))) {
  stop("Required columns 'Accession code' and/or 'Status' were not found after redefining the header.")
}

# Extract accession IDs and assign normalized group labels
sample_groups <- blup %>%
  transmute(
    sample_id  = `Accession code`,
    status_raw = Status
  ) %>%
  mutate(
    group = case_when(
      status_raw == "Landrace"        ~ "landrace",
      status_raw == "Early selection" ~ "early_selection",
      status_raw == "Modern breeding" ~ "modern_breeding",
      TRUE                            ~ NA_character_
    )
  ) %>%
  select(sample_id, group)

message("Group summary:")
print(table(sample_groups$group, useNA = "ifany"))

na_count <- sum(is.na(sample_groups$group))
if (na_count > 0) {
  warning(na_count, " samples were assigned group = NA. Check the Status column if this was not expected.")
}

write.table(
  sample_groups,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Sample group table written to: ", out_tsv)
