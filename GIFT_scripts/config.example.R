INPUT_BASE   <- "/path/to/GIFT/input"
GENO_BASE    <- "/path/to/geno/output"
WEIGHTS_BASE <- "/path/to/BGW_TPWAS/weights"
GEMMA_DIR    <- "/path/to/GEMMA/AD"
PRUNE_DIR    <- "/path/to/pruned/vcfs"

MATCHED_FILE     <- "/path/to/ALL_GENES_matched_combined_significant.txt"
GENE_COORDS_FILE <- "/path/to/ROSMAP_expr_TIGAR_format_WGS_IDs_b38_2023_chr1-22.tsv"
OUTPUT_SORTED    <- "/path/to/ALL_GENES_matched_sorted.txt"

EXPR_FILE  <- "/path/to/ROSMAP_expr_TIGAR_format_WGS_IDs_b38_2023_chr1-22.tsv"
ESTEP_MCMC <- "/path/to/BGW-xWAS-SS/bin/Estep_mcmc"

SLURM_PARTITIONS <- "yanglab,day-long-cpu,week-long-cpu,month-long-cpu"
CONDA_INIT       <- 'source "/path/to/miniconda3/etc/profile.d/conda.sh"'
CONDA_ENV        <- 'conda activate "/path/to/r_env/"'

N_RNASEQ         <- 931
PRUNE_REGION_IDS <- c(3, 9, 11, 13)

REGION_FILE <- file.path(INPUT_BASE, "region_from_matched.txt")
