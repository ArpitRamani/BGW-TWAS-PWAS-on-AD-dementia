library(data.table)
library(dplyr)

input5_base  <- "/home/jyang51/YangLabData/aramani/GIFT2/input1"
input_base   <- "/home/jyang51/YangLabData/aramani/GIFT2/input1"
geno_base    <- "/home/jyang51/YangLabData/aramani/BTARS/output"
weights_base <- "/home/jyang51/YangLabData/jyang/BGW_TPWAS_AD/wkdir"
region_file  <- file.path(input5_base, "region_from_matched.txt")
gemma_dir    <- "/projects/YangLabData/qliu/GIFT/Final_input_revise/AD"
prune_dir        <- "/projects/YangLabData/qliu/GIFT/ROSMAP_prune/prune/0.5"
PRUNE_REGION_IDS <- c(3, 9, 11, 13)

dir.create(input_base, showWarnings = FALSE, recursive = TRUE)

regions <- fread(region_file)
regions[, CHROM   := gsub("^chr", "", gsub("^CHR", "", as.character(CHROM)))]
regions[, txStart := suppressWarnings(as.integer(txStart))]
regions[, txEnd   := suppressWarnings(as.integer(txEnd))]
region_ids <- sort(unique(regions$region))


load_pruned_ids <- function(chr, prune_dir) {
  prune_file <- file.path(prune_dir, paste0("CHR", chr, ".vcf.gz"))
  if (!file.exists(prune_file)) {
    warning(paste0("no pruned vcf for chr", chr))
    return(NULL)
  }
  vcf <- fread(cmd = paste0("zcat ", prune_file, " | grep -v '^##'"),
               header = TRUE, sep = "\t", check.names = FALSE)
  id_col <- grep("^ID$", colnames(vcf), value = TRUE)
  ids <- if (length(id_col) == 0) vcf[[3]] else vcf[[id_col]]
  unique(ids)
}


build_geno <- function(region_dt, GWAS_stacked, geno_base) {

  gwas_dedup <- GWAS_stacked[!duplicated(GWAS_stacked$ps), ]
  gwas_dedup$geno_id      <- paste(gwas_dedup$chr, gwas_dedup$ps, gwas_dedup$allelel0, gwas_dedup$allele1, sep = ":")
  gwas_dedup$geno_id_flip <- paste(gwas_dedup$chr, gwas_dedup$ps, gwas_dedup$allele1, gwas_dedup$allelel0, sep = ":")
  keep_ids <- unique(c(gwas_dedup$geno_id, gwas_dedup$geno_id_flip))

  combined_geno <- data.frame()
  all_snp_info  <- data.frame()

  for (gene_name in region_dt$GeneName) {
    geno_file <- file.path(geno_base, paste0(gene_name, ".geno.txt"))
    if (!file.exists(geno_file)) {
      warning(paste0("no geno: ", gene_name))
      next
    }
    gene_geno <- read.delim(geno_file, header = TRUE, check.names = FALSE,
                            stringsAsFactors = FALSE, comment.char = "")
    snp_info  <- gene_geno[, 1:5]
    colnames(snp_info) <- c("CHROM", "POS", "ID", "REF", "ALT")
    geno_data <- gene_geno[, 6:ncol(gene_geno)]

    if (nrow(combined_geno) == 0) {
      combined_geno <- geno_data
      all_snp_info  <- snp_info
    } else {
      combined_geno <- rbind(combined_geno, geno_data)
      all_snp_info  <- rbind(all_snp_info,  snp_info)
    }
  }

  if (nrow(combined_geno) == 0) stop("no geno files found")

  all_snp_info$key <- paste(all_snp_info$CHROM, all_snp_info$POS, all_snp_info$REF, all_snp_info$ALT, sep = ":")
  dup_idx       <- !duplicated(all_snp_info$key)
  all_snp_info  <- all_snp_info[dup_idx, ]
  combined_geno <- combined_geno[dup_idx, ]

  keep_snps      <- all_snp_info$ID %in% keep_ids
  final_snp_info <- all_snp_info[keep_snps, c("CHROM", "POS", "ID", "REF", "ALT")]
  final_geno     <- combined_geno[keep_snps, ]

  if (nrow(final_snp_info) == 0) stop("no snps after geno filter")

  final_snp_info$CHROM <- as.numeric(final_snp_info$CHROM)
  final_snp_info$POS   <- as.numeric(final_snp_info$POS)
  sort_idx       <- order(final_snp_info$CHROM, final_snp_info$POS)
  final_snp_info <- final_snp_info[sort_idx, ]
  final_geno     <- final_geno[sort_idx, ]

  cbind(final_snp_info, final_geno)
}


drop_snps <- function(mask, GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, region_dir) {
  GWAS_stacked <- GWAS_stacked[mask, ]
  SNP_stacked  <- SNP_stacked[mask]

  gene_boundaries <- cumsum(c(0, pindex$V1))
  new_pindex <- integer(nrow(pindex))
  for (gi in seq_len(nrow(pindex))) {
    gene_rows      <- seq(gene_boundaries[gi] + 1, gene_boundaries[gi + 1])
    new_pindex[gi] <- sum(mask[gene_rows])
  }
  pindex <- data.table(V1 = new_pindex)

  write.table(GWAS_stacked, file.path(region_dir, "GWAS.txt"),    sep = "\t", quote = FALSE, col.names = TRUE,  row.names = FALSE)
  write.table(SNP_stacked,  file.path(region_dir, "snplist.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  fwrite(pindex,            file.path(region_dir, "pindex.txt"),  sep = "\t", quote = FALSE, col.names = FALSE)

  eqtl_files <- list.files(eqtl_out_dir, pattern = "\\.txt$", full.names = TRUE)
  for (ef in eqtl_files) {
    eq <- read.delim(ef, stringsAsFactors = FALSE)
    eq <- eq[eq$ps %in% GWAS_stacked$ps, ]
    write.table(eq, ef, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }

  list(GWAS_stacked = GWAS_stacked, SNP_stacked = SNP_stacked, pindex = pindex)
}


for (rid in region_ids) {

  region_dt <- regions[regions$region == rid]
  region_dt <- region_dt[order(GeneName)]

  input5_region_dir <- file.path(input5_base, paste0("region", rid))
  input9_region_dir <- file.path(input_base,  paste0("region", rid))
  dir.create(input9_region_dir, showWarnings = FALSE, recursive = TRUE)

  chr_vals <- unique(region_dt$CHROM)
  chr_vals <- chr_vals[!is.na(chr_vals) & chr_vals != ""]
  if (length(chr_vals) != 1) stop(paste0("region ", rid, ": multi-chr"))
  chr <- chr_vals[1]

  # --- FIX: collect all chromosomes needed from weights (cis + trans) ---
  all_weight_chrs <- c()
  for (gene in region_dt$GeneName) {
    wf <- file.path(weights_base, gene, paste0(gene, "_BGW_xQTL_weights.txt"))
    if (!file.exists(wf)) next
    w <- fread(wf)
    if ("#CHR" %in% colnames(w)) setnames(w, "#CHR", "CHR")
    w <- w[w$mBeta != 0, ]
    all_weight_chrs <- unique(c(all_weight_chrs, as.character(w$CHR)))
  }
  message(sprintf("[region %s] loading GEMMA for chromosomes: %s",
                  rid, paste(sort(as.integer(all_weight_chrs)), collapse = ", ")))

  # Load GEMMA for all relevant chromosomes and stack into one table
  gemma <- rbindlist(lapply(all_weight_chrs, function(c) {
    f <- file.path(gemma_dir, paste0("CHR", c, "_GEMMA.txt"))
    if (!file.exists(f)) {
      warning(paste0("no GEMMA file for chr", c))
      return(NULL)
    }
    g <- fread(f)
    g[, chr := gsub("^chr|^CHR", "", as.character(chr))]
    g[, ps  := suppressWarnings(as.integer(ps))]
    g
  }))
  # --- END FIX ---

  eqtl_out_dir <- file.path(input9_region_dir, "final_input", "eQTL")
  dir.create(eqtl_out_dir, showWarnings = FALSE, recursive = TRUE)

  GWAS_stacked <- data.frame()
  SNP_stacked  <- c()
  pindex       <- data.table(V1 = integer())

  for (j in seq_len(nrow(region_dt))) {

    gene <- region_dt$GeneName[j]

    weights_file <- file.path(weights_base, gene, paste0(gene, "_BGW_xQTL_weights.txt"))
    if (!file.exists(weights_file)) stop(paste0("no weights: ", gene))

    weights <- fread(weights_file)
    if ("#CHR" %in% colnames(weights)) setnames(weights, "#CHR", "CHR")

    weights <- weights[weights$mBeta != 0, ]
    if (nrow(weights) == 0) stop(paste0("empty weights: ", gene))
    weights[, POS := as.integer(POS)]

    gemma_sub <- as.data.frame(gemma[ps %in% weights$POS])
    weights   <- as.data.frame(weights)

    gwas_rows    <- list()
    weights_rows <- list()
    is_flipped   <- c()

    for (k in seq_len(nrow(gemma_sub))) {
      pos <- gemma_sub$ps[k]
      ref <- gemma_sub$allelel0[k]
      alt <- gemma_sub$allele1[k]

      w  <- weights[weights$POS == pos & weights$REF == ref & weights$ALT == alt, ]
      fl <- FALSE
      if (nrow(w) == 0) {
        w  <- weights[weights$POS == pos & weights$REF == alt & weights$ALT == ref, ]
        fl <- TRUE
      }
      if (nrow(w) == 0) next

      gwas_rows    <- c(gwas_rows,    list(gemma_sub[k, ]))
      weights_rows <- c(weights_rows, list(w[1, ]))
      is_flipped   <- c(is_flipped,   fl)
    }

    if (length(gwas_rows) == 0) stop(paste0("no gwas/weights overlap: ", gene))

    gwas_matched    <- do.call(rbind, gwas_rows)
    weights_matched <- do.call(rbind, weights_rows)
    rownames(gwas_matched)    <- NULL
    rownames(weights_matched) <- NULL

    eqtl_fmt        <- gwas_matched
    eqtl_fmt$n_obs  <- 931
    eqtl_fmt$af     <- weights_matched$MAF
    eqtl_fmt$beta   <- ifelse(is_flipped, -weights_matched$mBeta, weights_matched$mBeta)
    eqtl_fmt$se     <- abs(weights_matched$mBeta / sqrt(weights_matched$ChisqTest))
    eqtl_fmt$p_wald <- weights_matched$Pval_svt

    for (n in seq_len(nrow(eqtl_fmt))) {
      if (eqtl_fmt$allelel0[n] != weights_matched$REF[n]) {
        tmp                  <- eqtl_fmt$allelel0[n]
        eqtl_fmt$allelel0[n] <- eqtl_fmt$allele1[n]
        eqtl_fmt$allele1[n]  <- tmp
      }
    }

    write.table(eqtl_fmt, file.path(eqtl_out_dir, paste0("eQTL", gene, ".txt")),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    GWAS_stacked <- rbind(GWAS_stacked, gwas_matched)
    SNP_stacked  <- c(SNP_stacked, gwas_matched$rs)
    pindex       <- rbind(pindex, data.table(V1 = nrow(gwas_matched)))
  }

  if (sum(pindex$V1) != length(SNP_stacked)) stop(paste0("pindex/snplist mismatch: region ", rid))

  n_snps_before_prune <- nrow(GWAS_stacked)
  if (rid %in% PRUNE_REGION_IDS) {
    message(sprintf("[region %s] %d snps, pruning", rid, n_snps_before_prune))

    pruned_ids <- load_pruned_ids(chr, prune_dir)

    if (!is.null(pruned_ids)) {
      gwas_key      <- paste(GWAS_stacked$chr, GWAS_stacked$ps, GWAS_stacked$allelel0, GWAS_stacked$allele1, sep = ":")
      gwas_key_flip <- paste(GWAS_stacked$chr, GWAS_stacked$ps, GWAS_stacked$allele1, GWAS_stacked$allelel0, sep = ":")
      prune_mask    <- gwas_key %in% pruned_ids | gwas_key_flip %in% pruned_ids
      n_kept        <- sum(prune_mask)
      message(sprintf("[region %s] pruned: %d -> %d", rid, n_snps_before_prune, n_kept))

      if (n_kept == 0) {
        warning(sprintf("[region %s] all snps pruned, skipping", rid))
      } else {
        res          <- drop_snps(prune_mask, GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, input9_region_dir)
        GWAS_stacked <- res$GWAS_stacked
        SNP_stacked  <- res$SNP_stacked
        pindex       <- res$pindex
      }
    } else {
      warning(sprintf("[region %s] no pruned vcf, skipping", rid))
    }
  }

  write.table(GWAS_stacked, file.path(input9_region_dir, "GWAS.txt"),    sep = "\t", quote = FALSE, col.names = TRUE,  row.names = FALSE)
  write.table(SNP_stacked,  file.path(input9_region_dir, "snplist.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  fwrite(pindex,            file.path(input9_region_dir, "pindex.txt"),  sep = "\t", quote = FALSE, col.names = FALSE)

  geno_out  <- build_geno(region_dt, GWAS_stacked, geno_base)
  geno_file <- file.path(input9_region_dir, paste0("region", rid, ".geno.txt"))
  write.table(geno_out, geno_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

  expr_files <- list.files(input9_region_dir, pattern = "_expr\\.txt$", full.names = TRUE)
  EXPR <- NULL
  for (f in expr_files) {
    expr <- tryCatch(read.delim(f, header = FALSE), error = function(e) NULL)
    if (is.null(expr) || nrow(expr) == 0) next
    expr <- t(expr)
    colnames(expr) <- expr[1, ]
    expr <- as.data.frame(t(expr[-1, ]))
    expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(x)))
    expr[is.na(expr)] <- 0
    EXPR <- if (is.null(EXPR)) expr else bind_rows(EXPR, expr)
  }
  if (is.null(EXPR) || nrow(EXPR) == 0) {
    EXPR_corr <- matrix(0, nrow = nrow(region_dt), ncol = nrow(region_dt))
  } else {
    EXPR_corr <- cor(t(EXPR), use = "pairwise.complete.obs")
  }

  geno_pos     <- as.integer(geno_out[, "POS"])
  stacked_pos  <- as.integer(GWAS_stacked[, "ps"])
  geno_idx     <- match(stacked_pos, geno_pos)
  matched_mask <- !is.na(geno_idx)

  if (any(!matched_mask)) {
    res          <- drop_snps(!is.na(geno_idx), GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, input9_region_dir)
    GWAS_stacked <- res$GWAS_stacked
    SNP_stacked  <- res$SNP_stacked
    pindex       <- res$pindex
    geno_idx     <- geno_idx[matched_mask]
  }

  geno_mat  <- round(as.matrix(geno_out[geno_idx, 6:ncol(geno_out)]))
  row_vars  <- apply(geno_mat, 1, var, na.rm = TRUE)
  mono_mask <- row_vars == 0

  if (any(mono_mask)) {
    res          <- drop_snps(!mono_mask, GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, input9_region_dir)
    GWAS_stacked <- res$GWAS_stacked
    SNP_stacked  <- res$SNP_stacked
    pindex       <- res$pindex
    geno_mat     <- geno_mat[!mono_mask, ]
  }

  ld_corr <- cor(t(geno_mat), use = "pairwise.complete.obs")
  ld_corr[is.na(ld_corr)] <- 0

  final_input_dir <- file.path(input9_region_dir, "final_input")
  dir.create(final_input_dir, showWarnings = FALSE, recursive = TRUE)

  write.table(EXPR_corr, file.path(final_input_dir, "R_matrix.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(ld_corr,   file.path(final_input_dir, "GWASLD.txt"),   sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(ld_corr,   file.path(final_input_dir, "LD_eQTL.txt"),  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

  file.copy(file.path(input9_region_dir, "GWAS.txt"),    file.path(final_input_dir, "GWAS.txt"),    overwrite = TRUE)
  file.copy(file.path(input9_region_dir, "snplist.txt"), file.path(final_input_dir, "snplist.txt"), overwrite = TRUE)
  file.copy(file.path(input9_region_dir, "pindex.txt"),  file.path(final_input_dir, "pindex.txt"),  overwrite = TRUE)
}


out_run_dir <- file.path(input_base, "GIFT_running")
out_log_dir <- file.path(out_run_dir, "logs")
dir.create(out_run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_log_dir, showWarnings = FALSE, recursive = TRUE)

Whole <- "#!/bin/bash"

for (i in region_ids) {

  region_name <- paste0("region", i)

  R_script <- c(
    paste0('region <- "', region_name, '"'),
    "library(GIFT)",
    paste0('input_base <- "', input_base, '"'),
    'final_input_dir <- paste0(input_base, "/", region, "/final_input")',
    'eQTLfilelocation <- paste0(final_input_dir, "/eQTL/")',
    'GWASfile   <- paste0(final_input_dir, "/GWAS.txt")',
    'eQTLLDfile <- paste0(final_input_dir, "/LD_eQTL.txt")',
    'GWASLDfile <- paste0(final_input_dir, "/GWASLD.txt")',
    'snplist <- read.table(paste0(final_input_dir, "/snplist.txt"))$V1',
    'pindex  <- read.table(paste0(final_input_dir, "/pindex.txt"))$V1',
    'convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)',
    'gene      <- convert$gene',
    'Zscore1   <- convert$Zscore1',
    'Zscore2   <- convert$Zscore2',
    'LDmatrix1 <- convert$LDmatrix1',
    'LDmatrix2 <- convert$LDmatrix2',
    'n1 <- 931',
    'GWASresult <- read.table(GWASfile, header=TRUE)',
    'n2 <- mean(GWASresult[,5])',
    'R  <- as.matrix(read.table(paste0(final_input_dir, "/R_matrix.txt")))',
    'result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, n1, n2, gene, pindex, R=R, maxiter=100, tol=1e-3, pleio=0, ncores=4, in_sample_LD=FALSE, filter=TRUE, split=5)',
    'dir.create(paste0(final_input_dir, "/output"), showWarnings=FALSE)',
    'result$p[result$p == 0] <- .Machine$double.xmin',
    'write.table(result, paste0(final_input_dir, "/output/GIFT_", region, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")'
  )

  rfile  <- file.path(out_run_dir, paste0("GIFT_", region_name, ".R"))
  shfile <- file.path(out_run_dir, paste0("GIFT_", region_name, ".sh"))

  write.table(R_script, rfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

  bash_script <- c(
    "#!/bin/bash",
    paste0("#SBATCH --job-name=GIFT_", region_name),
    "#SBATCH --ntasks=1",
    "#SBATCH --nodes=1",
    "#SBATCH --cpus-per-task=4",
    "#SBATCH --mem=32G",
    "#SBATCH --partition=yanglab,day-long-cpu,week-long-cpu,month-long-cpu",
    paste0("#SBATCH --output=", out_log_dir, "/GIFT_", region_name, ".out.txt"),
    paste0("#SBATCH --error=",  out_log_dir, "/GIFT_", region_name, ".err.txt"),
    'source "/home/qliu259/miniconda3/etc/profile.d/conda.sh"',
    'conda activate "/projects/YangLabData/qliu/Conda/r_env/"',
    paste0('Rscript "', rfile, '"')
  )

  write.table(bash_script, shfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  Whole <- c(Whole, paste0('sbatch "', shfile, '"'))
}

submit_all <- file.path(out_run_dir, "submit_all_regions.sh")
write.table(Whole, submit_all, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cat("submit jobs: bash", submit_all, "\n")
system(paste("bash", submit_all))
