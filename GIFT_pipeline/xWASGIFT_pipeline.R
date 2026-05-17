#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
} else {
  stop("config.R not found. Copy config.example.R to config.R and edit paths.")
}

MERGED_GC_FILE       <- file.path(XWAS_DIR, "merged_gc_adjusted.tsv")
SIG_FILE             <- file.path(XWAS_DIR, "merged_gc_adjusted_significant.tsv")
ALL_MATCHED_FILE     <- file.path(XWAS_DIR, "ALL_GENES_matched_combined_significant.txt")
MISSING_GENES_FILE   <- file.path(XWAS_DIR, "missing_genes_in_matched.txt")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript xwas_gift_pipeline.R <STEP>")
step <- as.integer(args[1])
if (is.na(step) || step < 1 || step > 5) stop("STEP must be 1..5")

dir.create(XWAS_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(INPUT_BASE,  showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(MATCHED_DIR, showWarnings = FALSE, recursive = TRUE)


weights_path_for <- function(gene) {
  p1 <- file.path(WEIGHTS_BASE, gene, paste0(gene, "_BGW_xQTL_weights.txt"))
  if (file.exists(p1)) return(p1)
  if (nzchar(WEIGHTS_FALLBACK)) {
    p2 <- file.path(WEIGHTS_FALLBACK, gene, paste0(gene, "_BGW_xQTL_weights.txt"))
    if (file.exists(p2)) return(p2)
  }
  NA_character_
}
ss_path_for <- function(gene) {
  for (d in SS_DIRS) {
    p <- file.path(d, paste0(gene, SS_SUFFIX))
    if (file.exists(p)) return(p)
  }
  NA_character_
}
geno_path_for <- function(gene) {
  for (d in GENO_BASES) {
    p <- file.path(d, paste0(gene, GENO_SUFFIX))
    if (file.exists(p)) return(p)
  }
  NA_character_
}


if (step == 1) {

  files <- list.files(ASSOC_DIR,
                      pattern = paste0("^", ASSOC_PREFIX, "[0-9]+", ASSOC_SUFFIX, "$"),
                      full.names = TRUE)
  if (length(files) == 0) stop("no assoc files in ", ASSOC_DIR)

  all <- rbindlist(lapply(files, function(f) {
    chr <- as.integer(gsub(paste0(ASSOC_PREFIX, "|", ASSOC_SUFFIX), "", basename(f)))
    d <- fread(f)
    d[, CHROM := chr]
    d
  }), fill = TRUE)

  if (!ASSOC_GENE_COL %in% colnames(all)) stop("missing gene col: ", ASSOC_GENE_COL)
  if (!ASSOC_PVAL_COL %in% colnames(all)) stop("missing pval col: ", ASSOC_PVAL_COL)
  if (!ASSOC_STAT_COL %in% colnames(all)) stop("missing stat col: ", ASSOC_STAT_COL)

  setnames(all, ASSOC_GENE_COL, "GeneName")
  setnames(all, ASSOC_PVAL_COL, "Pvalue")
  setnames(all, ASSOC_STAT_COL, "RawStat")

  all[, Pvalue := as.numeric(Pvalue)]
  all[, RawStat := as.numeric(RawStat)]
  all <- all[!is.na(Pvalue) & Pvalue > 0 & Pvalue <= 1]
  all <- all[!is.na(RawStat)]

  if (ASSOC_STAT_TYPE == "z") {
    all[, ChiSquare := RawStat^2]
  } else {
    all[, ChiSquare := RawStat]
  }

  median_chi <- median(all$ChiSquare, na.rm = TRUE)
  lambda     <- median_chi / qchisq(0.5, df = 1)

  all[, ChiSquare_adj := ChiSquare / lambda]
  all[, Pvalue_adj    := pchisq(ChiSquare_adj, df = 1, lower.tail = FALSE)]
  all[, Pvalue        := Pvalue_adj]

  fwrite(all, MERGED_GC_FILE, sep = "\t")

  sig <- all[Pvalue < P_CUTOFF]
  fwrite(sig, SIG_FILE, sep = "\t")
}


if (step == 2) {

  if (!file.exists(SIG_FILE)) stop("missing ", SIG_FILE)
  sig <- fread(SIG_FILE)
  if (nrow(sig) == 0) stop("no significant genes")

  genes <- unique(sig$GeneName)

  missing <- list()

  for (gene in genes) {

    weights_file <- weights_path_for(gene)
    ss_file      <- ss_path_for(gene)
    if (is.na(ss_file)) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "no_ss")
      next
    }
    out_file <- file.path(MATCHED_DIR, paste0(gene, "_matched.txt"))

    if (file.exists(out_file)) next

    if (is.na(weights_file)) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "no_weights")
      next
    }

    weights <- tryCatch(fread(weights_file), error = function(e) NULL)
    ss      <- tryCatch(fread(ss_file),      error = function(e) NULL)
    if (is.null(weights) || is.null(ss) || nrow(weights) == 0 || nrow(ss) == 0) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "empty_input")
      next
    }

    if ("#CHR" %in% colnames(weights)) setnames(weights, "#CHR", "CHR")
    weights <- weights[mBeta != 0]
    if (nrow(weights) == 0) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "no_nonzero_weights")
      next
    }
    weights[, key      := paste(CHR, POS, REF, ALT, sep = ":")]
    weights[, key_flip := paste(CHR, POS, ALT, REF, sep = ":")]

    setnames(ss, colnames(ss), gsub("^#", "", colnames(ss)))
    ss_cols <- colnames(ss)
    chr_col <- ss_cols[grep("^(#?CHR|chr|CHROM)$", ss_cols, ignore.case = TRUE)][1]
    pos_col <- ss_cols[grep("^(POS|ps|BP|position)$", ss_cols, ignore.case = TRUE)][1]
    ref_col <- ss_cols[grep("^(REF|allele0|allelel0|a0)$", ss_cols, ignore.case = TRUE)][1]
    alt_col <- ss_cols[grep("^(ALT|allele1|a1)$", ss_cols, ignore.case = TRUE)][1]

    if (any(is.na(c(chr_col, pos_col, ref_col, alt_col)))) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "ss_col_detect_failed")
      next
    }

    setnames(ss, c(chr_col, pos_col, ref_col, alt_col), c("CHR_ss", "POS_ss", "REF_ss", "ALT_ss"))
    ss[, CHR_ss := gsub("^chr|^CHR", "", as.character(CHR_ss))]
    ss[, POS_ss := suppressWarnings(as.integer(POS_ss))]
    ss[, key_ss := paste(CHR_ss, POS_ss, REF_ss, ALT_ss, sep = ":")]

    matched_strict <- merge(weights, ss, by.x = "key",      by.y = "key_ss", suffixes = c("", ".ss"))
    matched_flip   <- merge(weights, ss, by.x = "key_flip", by.y = "key_ss", suffixes = c("", ".ss"))
    if (nrow(matched_flip) > 0) matched_flip[, flipped := TRUE]
    if (nrow(matched_strict) > 0) matched_strict[, flipped := FALSE]
    matched <- rbindlist(list(matched_strict, matched_flip), use.names = TRUE, fill = TRUE)

    if (nrow(matched) == 0) {
      missing[[length(missing) + 1]] <- list(gene = gene, reason = "no_overlap")
      next
    }

    fwrite(matched, out_file, sep = "\t")
  }

  if (length(missing) > 0) {
    miss_dt <- rbindlist(missing)
    fwrite(miss_dt, MISSING_GENES_FILE, sep = "\t")
  }
}


if (step == 3) {

  files <- list.files(MATCHED_DIR, pattern = "_matched\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("no matched files")

  combined <- rbindlist(lapply(files, function(f) {
    gene <- sub("_matched\\.txt$", "", basename(f))
    d <- fread(f)
    d[, GeneName := gene]
    d
  }), fill = TRUE)

  fwrite(combined, ALL_MATCHED_FILE, sep = "\t")
}


if (step == 4) {

  if (!file.exists(ALL_MATCHED_FILE)) stop("missing ", ALL_MATCHED_FILE)
  if (!file.exists(EXPR_FILE)) stop("missing expression file: ", EXPR_FILE)

  combined <- fread(ALL_MATCHED_FILE)
  expr <- fread(EXPR_FILE, select = c("CHROM", "GeneStart", "GeneEnd", "GeneName"))
  expr[, CHROM := gsub("^chr|^CHR", "", as.character(CHROM))]
  expr[, GeneStart := suppressWarnings(as.integer(GeneStart))]
  expr[, GeneEnd   := suppressWarnings(as.integer(GeneEnd))]

  unique_genes <- unique(combined$GeneName)
  gene_coords  <- expr[GeneName %in% unique_genes]

  setorder(gene_coords, CHROM, GeneStart)
  gene_coords[, txStart := pmax(GeneStart - 1e6, 0)]
  gene_coords[, txEnd   := GeneEnd + 1e6]

  gene_coords[, region := 0L]
  rid <- 0L
  cur_chr <- ""
  cur_end <- -Inf

  for (i in seq_len(nrow(gene_coords))) {
    chr <- gene_coords$CHROM[i]
    s   <- gene_coords$txStart[i]
    e   <- gene_coords$txEnd[i]
    if (chr != cur_chr || s > cur_end) {
      rid <- rid + 1L
      cur_chr <- chr
      cur_end <- e
    } else {
      cur_end <- max(cur_end, e)
    }
    gene_coords[i, region := rid]
  }

  region_counts <- gene_coords[, .N, by = region]
  keep_regions  <- region_counts[N > 1, region]
  gene_coords   <- gene_coords[region %in% keep_regions]

  gene_coords[, region := match(region, sort(unique(region)))]

  fwrite(gene_coords, REGION_FILE, sep = "\t")
}


if (step == 5) {

  if (!file.exists(REGION_FILE)) stop("missing ", REGION_FILE)
  if (!file.exists(EXPR_FILE)) stop("missing expression file: ", EXPR_FILE)

  regions <- fread(REGION_FILE)
  regions[, CHROM := gsub("^chr|^CHR", "", as.character(CHROM))]
  regions <- as.data.frame(regions)
  region_ids <- sort(unique(regions$region))

  options(scipen = 999)

  expr <- read.delim(EXPR_FILE, check.names = FALSE)
  meta_cols <- c("CHROM", "GeneStart", "GeneEnd", "GeneName")
  if ("TargetID" %in% colnames(expr)) meta_cols <- c(meta_cols[1:3], "TargetID", "GeneName")
  n_meta <- length(meta_cols)

  for (i in region_ids) {
    region <- regions[regions$region == i, ]
    region <- region[order(region$GeneName), ]

    for (j in seq_len(nrow(region))) {
      gene_name <- region$GeneName[j]
      gene_expr <- expr[which(expr[, "GeneName"] == gene_name), ]
      if (nrow(gene_expr) == 0) next
      gene_expr <- t(gene_expr[, (n_meta + 1):ncol(gene_expr)])
      gene_expr <- cbind(rownames(gene_expr), gene_expr)
      gene_expr[, 1] <- gsub("\\.", "-", gene_expr[, 1])
      gene_expr <- na.omit(gene_expr)

      out_file <- file.path(INPUT_BASE, paste0("region", i),
                            paste0(gene_name, "_expr.txt"))
      dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
      write.table(gene_expr, out_file, sep = "\t", quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
  }

  eQTL <- c()
  for (i in region_ids) {
    region <- regions[regions$region == i, ]
    region <- region[order(region$GeneName), ]
    for (j in seq_len(nrow(region))) {
      gene <- region$GeneName[j]
      region_dir <- file.path(INPUT_BASE, paste0("region", i))
      gene_dir   <- file.path(region_dir, gene)
      expr_path  <- file.path(region_dir, paste0(gene, "_expr.txt"))
      vcf_path   <- file.path(region_dir, "temp", paste0(gene, ".vcf.gz"))
      eQTL <- c(eQTL,
        paste0("cd ", region_dir),
        paste0("mkdir -p ", gene_dir),
        paste0("cd ", gene_dir),
        paste0(ESTEP_MCMC, " -vcf ", vcf_path,
               " -p ", expr_path,
               " -maf 0.01 -o eQTL_SumData -LDwindow 1 -GTfield GT -saveSS -zipSS")
      )
    }
  }
  eqtl_bash_file <- file.path(INPUT_BASE, "eQTL_bash.txt")
  write.table(eQTL, eqtl_bash_file, sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)

  build_geno <- function(region_dt, GWAS_stacked) {

    gwas_dedup <- GWAS_stacked[!duplicated(GWAS_stacked$ps), ]
    gwas_dedup$geno_id      <- paste(gwas_dedup$chr, gwas_dedup$ps, gwas_dedup$allelel0, gwas_dedup$allele1, sep = ":")
    gwas_dedup$geno_id_flip <- paste(gwas_dedup$chr, gwas_dedup$ps, gwas_dedup$allele1, gwas_dedup$allelel0, sep = ":")
    keep_ids <- unique(c(gwas_dedup$geno_id, gwas_dedup$geno_id_flip))

    combined_geno <- data.frame()
    all_snp_info  <- data.frame()

    for (gene_name in region_dt$GeneName) {
      geno_file <- geno_path_for(gene_name)
      if (is.na(geno_file)) {
        warning(paste0("no geno: ", gene_name)); next
      }
      gene_geno <- read.delim(geno_file, header = TRUE, check.names = FALSE,
                              stringsAsFactors = FALSE, comment.char = "")
      snp_info  <- gene_geno[, 1:5]
      colnames(snp_info) <- c("CHROM", "POS", "ID", "REF", "ALT")
      geno_data <- gene_geno[, 6:ncol(gene_geno)]
      if (nrow(combined_geno) == 0) {
        combined_geno <- geno_data; all_snp_info <- snp_info
      } else {
        combined_geno <- rbind(combined_geno, geno_data)
        all_snp_info  <- rbind(all_snp_info, snp_info)
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
      gene_rows <- seq(gene_boundaries[gi] + 1, gene_boundaries[gi + 1])
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

  load_pruned_ids <- function(chr, prune_dir) {
    prune_file <- file.path(prune_dir, paste0("CHR", chr, ".vcf.gz"))
    if (!file.exists(prune_file)) {
      warning(paste0("no pruned vcf for chr", chr)); return(NULL)
    }
    vcf <- fread(cmd = paste0("zcat ", prune_file, " | grep -v '^##'"),
                 header = TRUE, sep = "\t", check.names = FALSE)
    id_col <- grep("^ID$", colnames(vcf), value = TRUE)
    ids <- if (length(id_col) == 0) vcf[[3]] else vcf[[id_col]]
    unique(ids)
  }


  for (rid in region_ids) {

    region_dt <- regions[regions$region == rid, ]
    region_dt <- region_dt[order(region_dt$GeneName), ]

    region_dir <- file.path(INPUT_BASE, paste0("region", rid))
    dir.create(region_dir, showWarnings = FALSE, recursive = TRUE)

    chr_vals <- unique(region_dt$CHROM)
    chr_vals <- chr_vals[!is.na(chr_vals) & chr_vals != ""]
    if (length(chr_vals) != 1) stop(paste0("region ", rid, ": multi-chr"))
    chr <- chr_vals[1]

    gemma_file <- file.path(GEMMA_DIR, paste0("CHR", chr, "_GEMMA.txt"))
    if (!file.exists(gemma_file)) stop(paste0("missing GEMMA: ", gemma_file))
    gemma <- fread(gemma_file)
    gemma[, chr := gsub("^chr|^CHR", "", as.character(chr))]
    gemma[, ps  := suppressWarnings(as.integer(ps))]

    eqtl_out_dir <- file.path(region_dir, "final_input", "eQTL")
    dir.create(eqtl_out_dir, showWarnings = FALSE, recursive = TRUE)

    GWAS_stacked <- data.frame()
    SNP_stacked  <- c()
    pindex       <- data.table(V1 = integer())

    for (j in seq_len(nrow(region_dt))) {

      gene <- region_dt$GeneName[j]
      weights_file <- weights_path_for(gene)
      if (is.na(weights_file)) stop(paste0("no weights: ", gene))

      weights <- fread(weights_file)
      if ("#CHR" %in% colnames(weights)) setnames(weights, "#CHR", "CHR")
      weights <- weights[weights$mBeta != 0, ]
      if (CIS_TRANS_FILTER == "cis")   weights <- weights[weights$Trans == 0, ]
      if (CIS_TRANS_FILTER == "trans") weights <- weights[weights$Trans == 1, ]
      if (nrow(weights) == 0) stop(paste0("empty weights post-filter: ", gene))
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
      eqtl_fmt$n_obs  <- N_SAMPLES
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

    if (sum(pindex$V1) != length(SNP_stacked))
      stop(paste0("pindex/snplist mismatch: region ", rid))

    if (DO_PRUNE && rid %in% PRUNE_REGION_IDS) {
      pruned_ids <- load_pruned_ids(chr, PRUNE_DIR)
      if (!is.null(pruned_ids)) {
        gwas_key      <- paste(GWAS_stacked$chr, GWAS_stacked$ps, GWAS_stacked$allelel0, GWAS_stacked$allele1, sep = ":")
        gwas_key_flip <- paste(GWAS_stacked$chr, GWAS_stacked$ps, GWAS_stacked$allele1, GWAS_stacked$allelel0, sep = ":")
        prune_mask    <- gwas_key %in% pruned_ids | gwas_key_flip %in% pruned_ids
        n_kept        <- sum(prune_mask)
        if (n_kept == 0) {
          warning(sprintf("[region %d] all snps pruned, skipping", rid))
        } else {
          res          <- drop_snps(prune_mask, GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, region_dir)
          GWAS_stacked <- res$GWAS_stacked
          SNP_stacked  <- res$SNP_stacked
          pindex       <- res$pindex
        }
      }
      zero_gene_idx <- which(pindex$V1 == 0)
      if (length(zero_gene_idx) > 0) {
        kept_gene_idx <- which(pindex$V1 > 0)
        dropped_genes <- region_dt$GeneName[zero_gene_idx]

        for (g in dropped_genes) {
          ef <- file.path(eqtl_out_dir, paste0("eQTL", g, ".txt"))
          if (file.exists(ef)) file.remove(ef)
          xf <- file.path(region_dir, paste0(g, "_expr.txt"))
          if (file.exists(xf)) file.remove(xf)
        }

        pindex    <- data.table(V1 = pindex$V1[kept_gene_idx])
        region_dt <- region_dt[kept_gene_idx, ]
      }
    }

    write.table(GWAS_stacked, file.path(region_dir, "GWAS.txt"),    sep = "\t", quote = FALSE, col.names = TRUE,  row.names = FALSE)
    write.table(SNP_stacked,  file.path(region_dir, "snplist.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    fwrite(pindex,            file.path(region_dir, "pindex.txt"),  sep = "\t", quote = FALSE, col.names = FALSE)

    geno_out  <- build_geno(region_dt, GWAS_stacked)
    geno_file <- file.path(region_dir, paste0("region", rid, ".geno.txt"))
    write.table(geno_out, geno_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    expr_files <- list.files(region_dir, pattern = "_expr\\.txt$", full.names = TRUE)
    EXPR <- NULL
    for (f in expr_files) {
      expr_dt <- tryCatch(read.delim(f, header = FALSE), error = function(e) NULL)
      if (is.null(expr_dt) || nrow(expr_dt) == 0) next
      expr_dt <- t(expr_dt)
      colnames(expr_dt) <- expr_dt[1, ]
      expr_dt <- as.data.frame(t(expr_dt[-1, ]))
      expr_dt[] <- lapply(expr_dt, function(x) suppressWarnings(as.numeric(x)))
      expr_dt[is.na(expr_dt)] <- 0
      EXPR <- if (is.null(EXPR)) expr_dt else bind_rows(EXPR, expr_dt)
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
      res          <- drop_snps(!is.na(geno_idx), GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, region_dir)
      GWAS_stacked <- res$GWAS_stacked
      SNP_stacked  <- res$SNP_stacked
      pindex       <- res$pindex
      geno_idx     <- geno_idx[matched_mask]
    }

    geno_mat  <- round(as.matrix(geno_out[geno_idx, 6:ncol(geno_out)]))
    row_vars  <- apply(geno_mat, 1, var, na.rm = TRUE)
    mono_mask <- row_vars == 0
    if (any(mono_mask)) {
      res          <- drop_snps(!mono_mask, GWAS_stacked, SNP_stacked, pindex, eqtl_out_dir, region_dir)
      GWAS_stacked <- res$GWAS_stacked
      SNP_stacked  <- res$SNP_stacked
      pindex       <- res$pindex
      geno_mat     <- geno_mat[!mono_mask, ]
    }

    ld_corr <- cor(t(geno_mat), use = "pairwise.complete.obs")
    ld_corr[is.na(ld_corr)] <- 0

    final_input_dir <- file.path(region_dir, "final_input")
    dir.create(final_input_dir, showWarnings = FALSE, recursive = TRUE)
    write.table(EXPR_corr, file.path(final_input_dir, "R_matrix.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(ld_corr,   file.path(final_input_dir, "GWASLD.txt"),   sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(ld_corr,   file.path(final_input_dir, "LD_eQTL.txt"),  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    file.copy(file.path(region_dir, "GWAS.txt"),    file.path(final_input_dir, "GWAS.txt"),    overwrite = TRUE)
    file.copy(file.path(region_dir, "snplist.txt"), file.path(final_input_dir, "snplist.txt"), overwrite = TRUE)
    file.copy(file.path(region_dir, "pindex.txt"),  file.path(final_input_dir, "pindex.txt"),  overwrite = TRUE)
  }

  out_run_dir <- file.path(INPUT_BASE, "GIFT_running")
  out_log_dir <- file.path(out_run_dir, "logs")
  dir.create(out_run_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_log_dir, showWarnings = FALSE, recursive = TRUE)

  Whole <- "#!/bin/bash"

  for (i in region_ids) {
    region_name <- paste0("region", i)
    R_script <- c(
      paste0('region <- "', region_name, '"'),
      "library(GIFT)",
      paste0('input_base <- "', INPUT_BASE, '"'),
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
      paste0('n1 <- ', N_SAMPLES),
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
      paste0("#SBATCH --job-name=xwasGIFT_", region_name),
      "#SBATCH --ntasks=1",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=4",
      "#SBATCH --mem=32G",
      paste0("#SBATCH --partition=", SLURM_PARTITIONS),
      paste0("#SBATCH --output=", out_log_dir, "/GIFT_", region_name, ".out.txt"),
      paste0("#SBATCH --error=",  out_log_dir, "/GIFT_", region_name, ".err.txt"),
      CONDA_INIT,
      CONDA_ENV,
      paste0('Rscript "', rfile, '"')
    )
    write.table(bash_script, shfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    Whole <- c(Whole, paste0('sbatch "', shfile, '"'))
  }

  submit_all <- file.path(out_run_dir, "submit_all_regions.sh")
  write.table(Whole, submit_all, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
