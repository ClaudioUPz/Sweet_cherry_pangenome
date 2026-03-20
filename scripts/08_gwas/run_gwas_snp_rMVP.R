#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(bigmemory)
  library(rMVP)
  library(readr)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ragg)
})

#============================================================
# 0) Parameters
#============================================================

pheno_xlsx <- "Description_accessions_short_read_and_BLUPs.xlsx"
trait_name <- "Fru_Firm_G"
marker_label <- "SNP"

geno_file <- "ref_and_bam_bai/mvp_snps/snp_matrix_num.txt"
map_file_numeric <- "ref_and_bam_bai/mvp_snps/snp_map.numeric.txt"
map_file_text <- "ref_and_bam_bai/mvp_snps/snp_map.txt"

id_map_out <- paste0("id_mapping_", marker_label, "_min20.csv")

#============================================================
# 1) Phenotype data
#============================================================

phe_blup <- readxl::read_excel(pheno_xlsx, sheet = 1, skip = 1)

phe_blup <- phe_blup %>%
  dplyr::filter(Species == "P. avium")

meta_cols  <- c("Accession code", "Accession name", "Country of origin",
                "Status", "Species", "Parent 1", "Parent 2")
trait_cols <- setdiff(names(phe_blup), meta_cols)

phe_blup <- phe_blup %>%
  dplyr::filter(if_any(all_of(trait_cols), ~ !is.na(.)))

if (!("Accession code" %in% names(phe_blup))) {
  stop("Column 'Accession code' not found in phenotype file.")
}
if (!(trait_name %in% names(phe_blup))) {
  stop(paste0("Trait '", trait_name, "' not found in phenotype file."))
}

phe <- data.frame(
  Taxa  = as.character(phe_blup[["Accession code"]]),
  Trait = suppressWarnings(as.numeric(phe_blup[[trait_name]])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

phe <- phe[is.finite(phe$Trait), , drop = FALSE]
phe_input <- phe

#============================================================
# 2) Genotype and map files
#============================================================

geno <- read.table(
  geno_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

map_num <- read.table(
  map_file_numeric,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

map_txt <- read.table(
  map_file_text,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#============================================================
# 2.1) Harmonize sample names between genotype and phenotype
#============================================================

corrections <- c(
  "AZ_223_5_1" = "AZ_223/5_1",
  "KZ_073B1"   = "KZ_073B_1",
  "Rouge_du_Roussillon" = "Rouge du Roussillon",
  "Rustiques_des_Pyrenees" = "Rustiques_des_Pyrénées",
  "UZ_01B1"    = "UZ_01B_1",
  "Goldwich"   = "Goldrich"
)

taxa_g_raw <- colnames(geno)[-1]
taxa_g_rec <- dplyr::recode(taxa_g_raw, !!!corrections)
taxa_g_simple <- sub("__.*$", "", taxa_g_rec)
colnames(geno)[-1] <- taxa_g_simple

cat("Example genotype sample names:\n")
print(head(data.frame(
  raw = taxa_g_raw,
  recoded = taxa_g_rec,
  simplified = taxa_g_simple
), 10))

cat("Example phenotype sample names:\n")
print(head(phe_input$Taxa, 10))

taxa_g <- colnames(geno)[-1]
taxa_p <- phe_input$Taxa
common <- intersect(taxa_g, taxa_p)

cat("Number of shared individuals:", length(common), "\n")

if (length(common) == 0) {
  cat("\n[DEBUG] No shared sample IDs found.\n")
  cat("Genotype examples:\n")
  print(head(taxa_g, 10))
  cat("Phenotype examples:\n")
  print(head(taxa_p, 10))
  stop("No shared individuals between genotype and phenotype tables.")
}

geno <- geno[, c("ID", common)]
phe  <- phe_input[match(common, phe_input$Taxa), , drop = FALSE]

stopifnot(ncol(geno) - 1 == nrow(phe))
stopifnot(all(phe$Taxa == common))

#============================================================
# 3) Map file for rMVP
#============================================================

if (!all(c("SNP", "Chrom", "Position") %in% names(map_num))) {
  stop("Numeric map file must contain columns: SNP, Chrom, Position")
}
if (!all(c("SNP", "Chrom", "Position") %in% names(map_txt))) {
  stop("Text map file must contain columns: SNP, Chrom, Position")
}

map <- map_num
map$Chrom_txt <- map_txt$Chrom[match(map_num$SNP, map_txt$SNP)]

cat("Example chromosome names from text map:\n")
print(head(map$Chrom_txt, 10))

map$ChromNum <- suppressWarnings(
  as.integer(sub("^chr([0-9]+).*", "\\1", as.character(map$Chrom_txt)))
)

cat("Chromosome numeric summary:\n")
print(summary(map$ChromNum))
cat("Number of missing chromosome values:", sum(is.na(map$ChromNum)), "\n")

if (!any(!is.na(map$ChromNum))) {
  stop("Chromosome numbers could not be parsed from text map.")
}

map$Position <- suppressWarnings(as.integer(map$Position))

map_mvp <- data.frame(
  SNP      = map$SNP,
  Chrom    = map$ChromNum,
  Position = map$Position
)

map_mvp <- map_mvp[match(geno$ID, map_mvp$SNP), , drop = FALSE]

keep_non_na <- !is.na(map_mvp$SNP)
if (any(!keep_non_na)) {
  map_mvp <- map_mvp[keep_non_na, , drop = FALSE]
  geno    <- geno[keep_non_na, , drop = FALSE]
}

dup <- duplicated(map_mvp$SNP)
if (any(dup)) {
  map_mvp <- map_mvp[!dup, , drop = FALSE]
  geno    <- geno[!dup, , drop = FALSE]
}

stopifnot(nrow(map_mvp) == nrow(geno))
stopifnot(all(is.finite(map_mvp$Chrom)))
stopifnot(all(!is.na(map_mvp$Position)))
stopifnot(all(geno$ID == map_mvp$SNP))

#============================================================
# 4) Genotype matrix and filters
#============================================================

G <- as.matrix(geno[, -1])
rownames(G) <- geno$ID
storage.mode(G) <- "numeric"

stopifnot(all(phe$Taxa == colnames(G)))
stopifnot(ncol(G) == nrow(phe))

callrate <- rowMeans(!is.na(G))
pA  <- rowMeans(G / 2, na.rm = TRUE)
maf <- pmin(pA, 1 - pA)

keep <- (callrate >= 0.8) & (maf >= 0.01)
G <- G[keep, , drop = FALSE]
map_mvp <- map_mvp[keep, , drop = FALSE]

dup_mask <- duplicated(map_mvp$SNP) | duplicated(map_mvp$SNP, fromLast = TRUE)
cat("Duplicated marker IDs detected:", sum(dup_mask), "\n")

if (any(dup_mask)) {
  idx_grp <- ave(seq_len(nrow(map_mvp)), map_mvp$SNP, FUN = seq_along)
  new_ids <- ifelse(dup_mask, paste0(map_mvp$SNP, ":v", idx_grp), map_mvp$SNP)

  id_map <- data.frame(
    SNP_new  = new_ids,
    SNP_old  = map_mvp$SNP,
    Chrom    = map_mvp$Chrom,
    Position = map_mvp$Position,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(id_map, id_map_out, row.names = FALSE)

  map_mvp$SNP <- new_ids
  rownames(G) <- new_ids
} else {
  id_map <- data.frame(
    SNP_new  = map_mvp$SNP,
    SNP_old  = map_mvp$SNP,
    Chrom    = map_mvp$Chrom,
    Position = map_mvp$Position,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(id_map, id_map_out, row.names = FALSE)
}

stopifnot(!anyDuplicated(map_mvp$SNP))
stopifnot(all(map_mvp$SNP == rownames(G)))

for (i in seq_len(nrow(G))) {
  mu <- mean(G[i, ], na.rm = TRUE)
  if (is.finite(mu)) G[i, is.na(G[i, ])] <- mu
}

stopifnot(all(phe$Taxa == colnames(G)))
stopifnot(all(map_mvp$SNP == rownames(G)))

#============================================================
# 5) big.matrix, PCA, and kinship matrix
#============================================================

reset_bm <- function(bpath, stem) {
  try({
    if (exists("geno_bm")) rm(geno_bm, envir = .GlobalEnv)
    if (exists("geno_bm_ind")) rm(geno_bm_ind, envir = .GlobalEnv)
    gc()
    fdesc <- file.path(bpath, paste0(stem, ".desc"))
    fbin  <- file.path(bpath, paste0(stem, ".bin"))
    if (file.exists(fdesc)) file.remove(fdesc)
    if (file.exists(fbin)) file.remove(fbin)
    Sys.sleep(0.2)
  }, silent = TRUE)
}

bpath <- "mvp_cache"
dir.create(bpath, showWarnings = FALSE, recursive = TRUE)

reset_bm(bpath, "geno3")
reset_bm(bpath, "geno_ind")

if (file.exists(file.path(bpath, "geno3.bin"))) file.remove(file.path(bpath, "geno3.bin"))
if (file.exists(file.path(bpath, "geno3.desc"))) file.remove(file.path(bpath, "geno3.desc"))
if (file.exists(file.path(bpath, "geno_ind.bin"))) file.remove(file.path(bpath, "geno_ind.bin"))
if (file.exists(file.path(bpath, "geno_ind.desc"))) file.remove(file.path(bpath, "geno_ind.desc"))

geno_bm <- bigmemory::as.big.matrix(
  G,
  type = "double",
  backingfile = "geno3.bin",
  descriptorfile = "geno3.desc",
  backingpath = bpath
)

G_t <- t(G)
geno_bm_ind <- bigmemory::as.big.matrix(
  G_t,
  type = "double",
  backingfile = "geno_ind.bin",
  descriptorfile = "geno_ind.desc",
  backingpath = bpath
)

K <- MVP.K.VanRaden(geno_bm_ind)

PC_res <- MVP.PCA(M = geno_bm_ind, pcs.keep = 3)
PC <- if (is.list(PC_res)) PC_res$eigenVectors else as.matrix(PC_res)

samples <- rownames(G_t)
rownames(PC) <- samples
colnames(PC) <- paste0("PC", seq_len(ncol(PC)))

PC <- PC[match(phe$Taxa, rownames(PC)), , drop = FALSE]
stopifnot(all(rownames(PC) == phe$Taxa))

cat("Individuals in phenotype/genotype:", nrow(phe), ncol(geno_bm), "\n")
cat("Markers before filtering:", nrow(map), "\n")
cat("Markers after filtering:", nrow(map_mvp), "\n")
cat("Trait variance:", var(phe$Trait), "\n")
cat("Trait range:", range(phe$Trait), "\n")

#============================================================
# 6) GWAS with rMVP
#============================================================

glm_out <- paste0("mvp_out_glm_", trait_name, "_", marker_label, "_min20")
mlm_out <- paste0("mvp_out_mlm_", trait_name, "_", marker_label, "_min20")

dir.create(glm_out, showWarnings = FALSE, recursive = TRUE)
dir.create(mlm_out, showWarnings = FALSE, recursive = TRUE)

res_glm <- MVP(
  phe         = phe,
  geno        = geno_bm,
  map         = map_mvp,
  method      = "GLM",
  CV.GLM      = PC,
  file.output = FALSE,
  outpath     = paste0(glm_out, "/"),
  ncpus       = 1
)

res_mlm <- MVP(
  phe         = phe,
  geno        = geno_bm,
  map         = map_mvp,
  K           = K,
  CV.MLM      = PC,
  method      = "MLM",
  file.output = FALSE,
  outpath     = paste0(mlm_out, "/"),
  ncpus       = 1
)

#============================================================
# 7) Manhattan plots
#============================================================

bonf_thr <- 0.05 / nrow(map_mvp)

out_png <- "plots_gwas_png"
dir.create(out_png, showWarnings = FALSE, recursive = TRUE)

ragg::agg_png(
  filename = file.path(out_png, paste0("Manhattan_GLM_", trait_name, "_", marker_label, ".png")),
  width = 8,
  height = 4,
  units = "in",
  res = 400
)

MVP.Report(
  MVP           = res_glm,
  plot.type     = "m",
  LOG10         = TRUE,
  threshold     = bonf_thr,
  threshold.lwd = 2.5,
  threshold.col = "red",
  cex           = 0.5,
  cex.axis      = 1.3,
  cex.lab       = 1.4,
  lwd.axis      = 2,
  xlab          = "Chromosome",
  ylab          = expression(-log[10](p)),
  box           = TRUE,
  file.output   = FALSE,
  memo          = paste0("GLM_", trait_name, "_", marker_label)
)

dev.off()

ragg::agg_png(
  filename = file.path(out_png, paste0("Manhattan_MLM_", trait_name, "_", marker_label, ".png")),
  width = 8,
  height = 4,
  units = "in",
  res = 400
)

MVP.Report(
  MVP           = res_mlm,
  plot.type     = "m",
  LOG10         = TRUE,
  threshold     = bonf_thr,
  threshold.lwd = 2.5,
  threshold.col = "red",
  cex           = 0.5,
  cex.axis      = 1.3,
  cex.lab       = 1.4,
  lwd.axis      = 2,
  xlab          = "Chromosome",
  ylab          = expression(-log[10](p)),
  box           = TRUE,
  file.output   = FALSE,
  memo          = paste0("MLM_", trait_name, "_", marker_label)
)

dev.off()

#============================================================
# 8) Export significant hits and genomic regions
#============================================================

sig_file  <- file.path(mlm_out, "Trait.MLM_signals.csv")
imap_file <- id_map_out

if (file.exists(sig_file)) {
  sig  <- read.csv(sig_file, check.names = FALSE)
  imap <- read.csv(imap_file, check.names = FALSE)

  imap <- imap %>% select(SNP_new, SNP_old)

  ann <- sig %>%
    rename(SNP_new = SNP, p = `Trait.MLM`) %>%
    left_join(imap, by = "SNP_new") %>%
    mutate(
      CHROM_vcf = paste0("chr", Chrom, "__REFERENCE"),
      POS_vcf   = Position
    )

  out_ann <- paste0("hits_", trait_name, "_", marker_label, "_min20_annot.csv")
  write.csv(ann, out_ann, row.names = FALSE)

  regions <- ann %>%
    transmute(CHROM_vcf, START = POS_vcf, END = POS_vcf)

  out_reg <- paste0("hits_", trait_name, "_", marker_label, "_min20.regions.tsv")
  write.table(
    regions,
    out_reg,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
} else {
  cat("[WARNING] Signal file not found. Skipping annotated hit export.\n")
}

#============================================================
# 9) Diagnostic plots and genomic inflation
#============================================================

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

p_hist <- ggplot(phe, aes(x = Trait)) +
  geom_histogram(bins = 40, fill = "grey80", color = "grey40") +
  geom_density(linewidth = 1) +
  labs(
    x = "BLUP phenotype",
    y = "Density",
    title = paste0("BLUP distribution: ", trait_name, " (", marker_label, ")")
  ) +
  theme_bw()

ggsave(
  paste0("plots/BLUP_distribution_", trait_name, "_", marker_label, ".png"),
  p_hist, width = 6, height = 4, dpi = 300
)

p_qq <- ggplot(phe, aes(sample = Trait)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = paste0("BLUP QQ plot: ", trait_name, " (", marker_label, ")")) +
  theme_bw()

ggsave(
  paste0("plots/BLUP_QQ_", trait_name, "_", marker_label, ".png"),
  p_qq, width = 5, height = 5, dpi = 300
)

mlm_all_file <- file.path(mlm_out, "Trait.MLM.csv")

if (file.exists(mlm_all_file)) {
  mlm_all <- read.csv(mlm_all_file)

  median_p <- median(mlm_all$Trait.MLM, na.rm = TRUE)
  cat("Median p-value:", median_p, "\n")

  chisq <- qchisq(1 - mlm_all$Trait.MLM, 1)
  lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  cat("Genomic inflation factor (lambda):", lambda_gc, "\n")

  obs <- -log10(sort(mlm_all$Trait.MLM))
  exp <- -log10(ppoints(length(obs)))

  png(
    paste0("plots/QQ_MLM_", trait_name, "_", marker_label, ".png"),
    width = 1400,
    height = 1400,
    res = 200
  )
  plot(
    exp, obs,
    pch = 19, cex = 0.4,
    xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)",
    main = paste0(
      "QQ plot MLM (lambda=", round(lambda_gc, 3),
      ") - ", trait_name, " (", marker_label, ")"
    )
  )
  abline(0, 1, lwd = 2)
  dev.off()

  phe_pcs <- cbind(phe, PC)
  cor_res <- cor(phe_pcs$Trait, phe_pcs[, grep("^PC", names(phe_pcs))])
  print(cor_res)
} else {
  cat("[WARNING] MLM results file not found. Skipping lambda and QQ plot.\n")
}
