#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript calculate_fst_snps_landrace_vs_modern.R ",
    "<snps.genotypes.gt.tsv> <snps.sample_order.txt> <sample_groups_checked.tsv> <output_fst.tsv>\n"
  )
}

geno_file <- args[1]
sample_order_file <- args[2]
groups_file <- args[3]
out_file <- args[4]

cat("Genotype table:", geno_file, "\n")
cat("Sample order file:", sample_order_file, "\n")
cat("Sample group table:", groups_file, "\n")
cat("Output FST file:", out_file, "\n\n")

#------------------------------------------------------------
# 1. Read sample order
#------------------------------------------------------------
samples_order <- read_lines(sample_order_file)
samples_order <- samples_order[samples_order != ""]

cat("Number of samples in VCF order file:", length(samples_order), "\n")

#------------------------------------------------------------
# 2. Read and normalize group labels
#------------------------------------------------------------
groups <- read_tsv(groups_file, show_col_types = FALSE)
colnames(groups)[1:2] <- c("SampleID", "Status")

groups <- groups %>%
  mutate(Status = case_when(
    Status %in% c("Landrace", "landrace") ~ "landrace",
    Status %in% c("Modern breeding", "modern_breeding", "modern") ~ "modern",
    Status %in% c("Early selection", "early_selection", "early") ~ "early",
    TRUE ~ NA_character_
  ))

cat("Group summary after normalization:\n")
print(table(groups$Status, useNA = "ifany"))
cat("\n")

groups_in_vcf <- groups %>%
  filter(SampleID %in% samples_order)

cat("Samples with known group and present in VCF:\n")
print(table(groups_in_vcf$Status, useNA = "ifany"))
cat("\n")

groups_lm <- groups_in_vcf %>%
  filter(Status %in% c("landrace", "modern"))

n_land <- sum(groups_lm$Status == "landrace")
n_mod  <- sum(groups_lm$Status == "modern")

cat("Landrace samples:", n_land, "\n")
cat("Modern breeding samples:", n_mod, "\n\n")

if (n_land < 2 || n_mod < 2) {
  stop("Too few samples in one or both populations for FST calculation.")
}

#------------------------------------------------------------
# 3. Read genotype header and clean sample names
#------------------------------------------------------------
clean_colname <- function(x) {
  x <- sub("^#", "", x)
  x <- sub("^\\[[0-9]+\\]", "", x)
  x <- sub("^\\s+", "", x)
  x <- sub("\\s+$", "", x)
  x <- sub(":GT$", "", x)
  x
}

header_line <- readLines(geno_file, n = 1)
raw_cols <- strsplit(header_line, "\t")[[1]]
clean_cols <- vapply(raw_cols, clean_colname, FUN.VALUE = character(1))

cat("First cleaned column names:\n")
print(head(clean_cols, 10))
cat("\n")

#------------------------------------------------------------
# 4. Read genotype table
#------------------------------------------------------------
geno <- read_tsv(
  geno_file,
  skip = 1,
  col_names = clean_cols,
  show_col_types = FALSE
)

cat("Example genotype data columns:\n")
print(head(geno[, 1:min(4, ncol(geno))]))
cat("\n")

samples_in_geno <- colnames(geno)[-(1:2)]
cat("Number of samples in genotype table:", length(samples_in_geno), "\n")

#------------------------------------------------------------
# 5. Align sample order and group labels
#------------------------------------------------------------
common_samples <- intersect(samples_order, samples_in_geno)
cat("Samples shared between genotype table and sample order file:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No shared samples between genotype table and sample order file.")
}

samples_use <- samples_order[samples_order %in% samples_in_geno]
geno <- geno[, c("CHROM", "POS", samples_use)]

land_samples <- groups_lm %>%
  filter(Status == "landrace") %>%
  pull(SampleID)

modern_samples <- groups_lm %>%
  filter(Status == "modern") %>%
  pull(SampleID)

land_samples <- intersect(land_samples, samples_use)
modern_samples <- intersect(modern_samples, samples_use)

cat("Landrace samples used:", length(land_samples), "\n")
cat("Modern samples used:", length(modern_samples), "\n\n")

if (length(land_samples) < 2 || length(modern_samples) < 2) {
  stop("Too few samples in one or both populations after intersecting with genotype data.")
}

#------------------------------------------------------------
# 6. Recode genotypes as ALT allele counts
#------------------------------------------------------------
geno_mat <- as.matrix(geno[, samples_use])

geno_alt <- matrix(NA_real_, nrow = nrow(geno_mat), ncol = ncol(geno_mat))
colnames(geno_alt) <- samples_use

geno_alt[geno_mat == "0/0"] <- 0
geno_alt[geno_mat == "0/1" | geno_mat == "1/0"] <- 1
geno_alt[geno_mat == "1/1"] <- 2

idx_land <- match(land_samples, samples_use)
idx_mod  <- match(modern_samples, samples_use)

alt_land <- geno_alt[, idx_land, drop = FALSE]
alt_mod  <- geno_alt[, idx_mod, drop = FALSE]

#------------------------------------------------------------
# 7. Compute per-SNP Hudson FST
#------------------------------------------------------------
n1 <- 2 * rowSums(!is.na(alt_land))
n2 <- 2 * rowSums(!is.na(alt_mod))

alt1 <- rowSums(alt_land, na.rm = TRUE)
alt2 <- rowSums(alt_mod, na.rm = TRUE)

p1 <- alt1 / n1
p2 <- alt2 / n2

min_alleles <- 4
valid <- n1 >= min_alleles & n2 >= min_alleles

pi1 <- 2 * p1 * (1 - p1)
pi2 <- 2 * p2 * (1 - p2)
pi_within <- (pi1 + pi2) / 2

pi_between <- p1 * (1 - p2) + p2 * (1 - p1)

fst_hudson <- (pi_between - pi_within) / pi_between
fst_hudson[!valid | pi_between == 0 | is.na(pi_between)] <- NA_real_

#------------------------------------------------------------
# 8. Save results
#------------------------------------------------------------
res <- tibble(
  CHROM = geno$CHROM,
  POS   = geno$POS,
  N_alleles_landrace = n1,
  N_alleles_modern   = n2,
  p_landrace = p1,
  p_modern   = p2,
  Fst_Hudson = fst_hudson
)

cat("Summary of Hudson FST values:\n")
print(summary(res$Fst_Hudson))
cat("\n")

write_tsv(res, out_file)
cat("Results written to:", out_file, "\n")
