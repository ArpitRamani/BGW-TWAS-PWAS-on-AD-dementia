#!/usr/bin/env Rscript

# --------------------------------------------------------------------------
# BGW association test (cis + trans)
#
# Reads per-gene BGW-xWAS weights, genotype dosages, and GWAS summary
# statistics, then computes a burden-style Z statistic and TWAS/PWAS
# p-value for each gene in a gene list.
#
# Usage:
#   Rscript association_test.R <gene_list_file> <chr> <num_cores>
#
# Data paths are read from config.R. Resolution order:
#   1. environment variable CONFIG_R (set this in your SLURM / shell config)
#   2. config.R located one directory above this script
#
# config.R must define: WEIGHT_DIR, GENO_DIR, SS_DIR, RESULT_DIR
# (see config.example.R)
# --------------------------------------------------------------------------

args <- commandArgs(TRUE)

gene_list_file <- args[[1]]
chr            <- args[[2]]
num_cores      <- as.numeric(args[[3]])

suppressMessages({
  library(data.table)
  library(tidyverse)
  library(foreach)
  library(doParallel)
})

# ---- locate and source config ----
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}

config_path <- Sys.getenv("CONFIG_R", unset = "")
if (!nzchar(config_path)) {
  config_path <- file.path(get_script_dir(), "..", "config.R")
}
if (!file.exists(config_path)) {
  stop(sprintf("config.R not found. Set CONFIG_R or create config.R (looked for: %s)", config_path))
}
source(config_path)

for (v in c("WEIGHT_DIR", "GENO_DIR", "SS_DIR", "RESULT_DIR")) {
  if (!exists(v)) stop(sprintf("config.R is missing required variable: %s", v))
}
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

gene_list <- scan(gene_list_file, what = character())

make_id <- function(chrom, pos, ref, alt) {
  chrom <- sub("^[Cc][Hh][Rr]", "", as.character(chrom))
  paste(chrom, pos, ref, alt, sep = ":")
}

assoc_test_func <- function(gene_name) {

  weight_file  <- file.path(WEIGHT_DIR, gene_name, paste0(gene_name, "_BGW_xQTL_weights.txt"))
  geno_file    <- file.path(GENO_DIR,   paste0(gene_name, ".geno.txt"))
  GWAS_SS_file <- file.path(SS_DIR,     paste0(gene_name, ".SS.txt"))

  na_row <- data.frame(GeneName         = gene_name,
                       ZStat            = NA, PVAL_twas        = NA,
                       N_with_weight    = NA, Ncis_with_weight = NA,
                       N_tested         = NA, Ncis_tested      = NA,
                       sum_cpp_cis      = NA, sum_cpp_trans    = NA)

  if (!file.exists(weight_file) || !file.exists(geno_file) || !file.exists(GWAS_SS_file)) {
    return(na_row)
  }

  weight <- tryCatch({
    first_line <- readLines(weight_file, n = 1)
    has_header <- grepl("^[A-Za-z#]", trimws(first_line))
    read.table(weight_file, header = has_header)
  }, error = function(e) NULL)

  if (is.null(weight) || nrow(weight) == 0 || ncol(weight) < 13) return(na_row)

  colnames(weight)[1:13] <- c('CHROM','POS','ID','REF','ALT','MAF',
                               'Trans','CPP','Beta','mBeta','ChisqTest','Pval_svt','Rank')
  weight$Trans <- as.integer(as.character(weight$Trans))
  weight$ID    <- make_id(weight$CHROM, weight$POS, weight$REF, weight$ALT)
  if (any(duplicated(weight$ID))) weight <- weight[!duplicated(weight$ID), ]

  nsnp_weight <- nrow(weight)
  ncis_weight <- sum(weight$Trans == 0, na.rm = TRUE)

  SS <- tryCatch({
    ss <- read.table(GWAS_SS_file, header = TRUE, comment.char = "#", fill = TRUE)
    if (ncol(ss) < 7) stop()
    colnames(ss)[1:8] <- c("CHROM","POS","rsID","REF","ALT","Beta_gwas","Beta_sd_gwas","PVAL_gwas")
    ss$ID     <- make_id(ss$CHROM, ss$POS, ss$REF, ss$ALT)
    ss$Zscore <- ss$Beta_gwas / ss$Beta_sd_gwas
    ss$Zscore[ss$Zscore == -Inf] <- -50
    ss$Zscore[ss$Zscore ==  Inf] <-  50
    ss
  }, error = function(e) NULL)

  if (is.null(SS) || nrow(SS) == 0) {
    na_row$N_with_weight    <- nsnp_weight
    na_row$Ncis_with_weight <- ncis_weight
    return(na_row)
  }

  dosage <- tryCatch(as.data.frame(fread(geno_file, header = TRUE)), error = function(e) NULL)

  if (is.null(dosage) || nrow(dosage) == 0) {
    na_row$N_with_weight    <- nsnp_weight
    na_row$Ncis_with_weight <- ncis_weight
    return(na_row)
  }

  colnames(dosage)[1] <- "CHROM"
  dosage$ID <- sub("^[Cc][Hh][Rr]", "", dosage$ID)

  fixed_meta  <- c("CHROM", "POS", "ID", "REF", "ALT")
  sample_cols <- setdiff(colnames(dosage), fixed_meta)

  if (length(sample_cols) == 0) {
    na_row$N_with_weight    <- nsnp_weight
    na_row$Ncis_with_weight <- ncis_weight
    return(na_row)
  }

  if (any(duplicated(dosage$ID))) dosage <- dosage[!duplicated(dosage$ID), ]

  dosage <- semi_join(dosage, weight, by = "ID")
  weight <- semi_join(weight, dosage, by = "ID")

  weight1 <- semi_join(weight, SS, by = "ID")
  SS1     <- semi_join(SS,     weight, by = "ID")

  SS_flip        <- SS
  SS_flip$ID     <- make_id(SS_flip$CHROM, SS_flip$POS, SS_flip$ALT, SS_flip$REF)
  SS_flip$Zscore <- -1 * SS_flip$Zscore
  weight2 <- semi_join(weight,  SS_flip, by = "ID")
  SS2     <- semi_join(SS_flip, weight,  by = "ID")

  weight <- rbind(weight1, weight2)
  SS     <- rbind(SS1, SS2)
  SS$ID  <- make_id(SS$CHROM, SS$POS, SS$REF, SS$ALT)

  weight <- weight[order(weight$ID), ]
  SS     <- SS[order(SS$ID), ]
  dosage <- semi_join(dosage, weight, by = "ID")
  dosage <- dosage[order(dosage$ID), ]

  nsnp_used     <- nrow(weight)
  ncis_used     <- sum(weight$Trans == 0, na.rm = TRUE)
  sum_cpp_cis   <- sum(weight$CPP[weight$Trans == 0],  na.rm = TRUE)
  sum_cpp_trans <- sum(weight$CPP[weight$Trans != 0],  na.rm = TRUE)

  stat <- NA; pval <- NA

  if (nrow(dosage) > 1 &&
      nrow(weight) == nrow(dosage) &&
      nrow(SS)     == nrow(weight) &&
      nrow(SS)     == nrow(dosage)) {

    dosage2 <- t(as.matrix(dosage[, sample_cols]))
    covmat  <- cov(dosage2, use = "complete.obs")
    zscore  <- SS$Zscore
    sigma_l <- sqrt(diag(covmat))
    w       <- weight$CPP * weight$Beta

    if (length(zscore) == length(sigma_l) && length(sigma_l) == length(w)) {
      numerator   <- as.numeric(crossprod(zscore * sigma_l, w))
      denominator <- sqrt(as.numeric(crossprod(w, covmat %*% w)))

      if (is.finite(denominator) && denominator != 0) {
        stat <- numerator / denominator
        pval <- 2 * pnorm(-abs(stat))
      }
    }
  }

  return(data.frame(GeneName         = gene_name,
                    ZStat            = stat,
                    PVAL_twas        = pval,
                    N_with_weight    = nsnp_weight,
                    Ncis_with_weight = ncis_weight,
                    N_tested         = nsnp_used,
                    Ncis_tested      = ncis_used,
                    sum_cpp_cis      = sum_cpp_cis,
                    sum_cpp_trans    = sum_cpp_trans))
}

registerDoParallel(cores = num_cores)
assoc_out_dt <- foreach(i = seq_along(gene_list), .combine = rbind) %dopar% {
  assoc_test_func(gene_list[i])
}

out_file <- file.path(RESULT_DIR, paste0("CHR_", chr, "_assoc.tsv"))
write_tsv(assoc_out_dt, file = out_file)
