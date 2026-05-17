#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

matched_file <- "/home/jyang51/YangLabData/aramani/TWAS2/Results/ALL_GENES_matched_combined_significant.txt"
gene_coords_file <- "/home/jyang51/YangLabData/jyang/BGW_TPWAS_AD/Data/RNAseq_BulkBrain/ROSMAP_expr_TIGAR_format_WGS_IDs_b38_2023_chr1-22.tsv"
output_sorted <- "/home/jyang51/YangLabData/aramani/GIFT/ALL_GENES_matched_sorted.txt"
region_file <- "/home/jyang51/YangLabData/aramani/GIFT2/input1/region_from_matched.txt"
base_region_dir <- "/home/jyang51/YangLabData/aramani/GIFT2/input1"

matched_data <- fread(matched_file)
matched_sorted <- matched_data %>% arrange(CHROM, POS)
fwrite(matched_sorted, output_sorted, sep = "\t", quote = FALSE)
cat("Rows in matched_sorted:", nrow(matched_sorted), "\n")

cat("Reading gene coordinates from expression file...\n")
gene_coords <- fread(gene_coords_file, select = c("CHROM", "GeneStart", "GeneEnd", "GeneName"))

genes_with_snps <- unique(matched_sorted$GeneName)
cat("Total genes with SNPs:", length(genes_with_snps), "\n")

gene_boundaries <- gene_coords %>%
  filter(GeneName %in% genes_with_snps) %>%
  rename(txStart = GeneStart, txEnd = GeneEnd) %>%
  mutate(start = txStart - 1e6, end = txEnd + 1e6)

cat("Genes with both coordinates and SNP data:", nrow(gene_boundaries), "\n")

setDT(gene_boundaries)
gene_boundaries <- gene_boundaries[order(CHROM, start)]

merge_intervals <- function(dt) {
  if (nrow(dt) == 0) return(data.table())
  if (nrow(dt) == 1) return(dt[, .(start, end)])

  merged <- list()
  current_start <- dt$start[1]
  current_end <- dt$end[1]

  for (i in 2:nrow(dt)) {
    if (dt$start[i] <= current_end) {
      current_end <- max(current_end, dt$end[i])
    } else {
      merged[[length(merged) + 1]] <- list(start = current_start, end = current_end)
      current_start <- dt$start[i]
      current_end <- dt$end[i]
    }
  }
  merged[[length(merged) + 1]] <- list(start = current_start, end = current_end)
  rbindlist(merged)
}

merged_regions <- gene_boundaries[, merge_intervals(.SD), by = CHROM]
merged_regions[, region := .I]
cat("Total merged regions:", nrow(merged_regions), "\n")

setkey(merged_regions, CHROM, start, end)
setkey(gene_boundaries, CHROM, start, end)

gene_region_map <- foverlaps(
  gene_boundaries, merged_regions,
  nomatch = 0
)[, .(GeneName, CHROM, txStart, txEnd, region)]

cat("Rows after foverlaps:", nrow(gene_region_map), "\n")
cat("Unique regions after foverlaps:", length(unique(gene_region_map$region)), "\n")

gene_region_map <- gene_region_map %>%
  group_by(region) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(CHROM, txStart) %>%
  mutate(region = as.numeric(factor(region, levels = unique(region))))

cat("Rows after n>1 filter:", nrow(gene_region_map), "\n")
cat("Unique regions after n>1 filter:", length(unique(gene_region_map$region)), "\n")

setDT(gene_region_map)

fwrite(gene_region_map, region_file, sep = "\t", quote = FALSE)
cat("Created region file with", nrow(gene_region_map), "genes\n")

region_gene_list_dir <- file.path(base_region_dir, "gene_lists")
dir.create(region_gene_list_dir, showWarnings = FALSE, recursive = TRUE)
cat("Gene list dir exists:", dir.exists(region_gene_list_dir), "\n")

cat("\nCreating gene list files per region...\n")
for (i in unique(gene_region_map$region)) {
  region_genes <- gene_region_map[region == i, GeneName]
  gene_list_file <- file.path(region_gene_list_dir, paste0("region", i, "_genes.txt"))
  writeLines(region_genes, gene_list_file)
  cat("Region", i, ":", length(region_genes), "genes\n")
}
cat("Created gene list files in", region_gene_list_dir, "\n")

options(scipen = 999)
regions <- unique(gene_region_map$region)
cat("\nStarting region folder creation for", length(regions), "regions...\n")

for (i in regions) {
  region <- gene_region_map[region == i]
  region_genes <- unique(region$GeneName)
  region_snps <- matched_sorted[GeneName %in% region_genes]

  region_dir <- file.path(base_region_dir, paste0("region", i))
  dir.create(region_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Creating folder:", region_dir, "| exists:", dir.exists(region_dir), "\n")

  snplist_file <- file.path(region_dir, "snplist.txt")
  if ("rsID" %in% names(region_snps)) {
    fwrite(region_snps[, .(rsID)], snplist_file, sep = "\t", quote = FALSE, col.names = FALSE)
  }

  snploc_file <- file.path(region_dir, "snp_loc.txt")
  fwrite(region_snps[, .(CHROM, POS)], snploc_file, sep = "\t", quote = FALSE)

  gwas_cols <- intersect(
    c("CHROM", "POS", "rsID", "REF", "ALT", "Beta_gwas", "Beta_sd_gwas", "PVAL_gwas", "Pvalue"),
    names(region_snps)
  )
  gwas_file <- file.path(region_dir, "GWAS.txt")
  fwrite(region_snps[, ..gwas_cols], gwas_file, sep = "\t", quote = FALSE)
}
cat("\nDone.\n")
