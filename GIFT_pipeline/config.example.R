INPUT_BASE   <- "/path/to/GIFT/input"
GENO_BASE    <- "/path/to/geno/output"
WEIGHTS_BASE <- "/path/to/BGW/weights"
GEMMA_DIR    <- "/path/to/GEMMA"
PRUNE_DIR    <- "/path/to/pruned/vcfs"

MATCHED_FILE     <- "/path/to/ALL_GENES_matched_combined_significant.txt"
GENE_COORDS_FILE <- "/path/to/expression_TIGAR_format.tsv"
OUTPUT_SORTED    <- "/path/to/ALL_GENES_matched_sorted.txt"

EXPR_FILE  <- "/path/to/expression_TIGAR_format.tsv"
ESTEP_MCMC <- "/path/to/BGW-xWAS-SS/bin/Estep_mcmc"

SLURM_PARTITIONS <- "partition1,partition2"
CONDA_INIT       <- 'source "/path/to/miniconda3/etc/profile.d/conda.sh"'
CONDA_ENV        <- 'conda activate "/path/to/r_env/"'

N_SAMPLES        <- 931
PRUNE_REGION_IDS <- c()
DO_PRUNE         <- FALSE

REGION_FILE <- file.path(INPUT_BASE, "region_from_matched.txt")

XWAS_DIR     <- "/path/to/xwas/run"
LOG_DIR      <- file.path(XWAS_DIR, "logs")
MATCHED_DIR  <- file.path(XWAS_DIR, "matched_weights")

CIS_TRANS_FILTER <- ""
P_CUTOFF         <- 2.5e-6

ASSOC_DIR       <- "/path/to/assoc_results"
ASSOC_PREFIX    <- "CHR_"
ASSOC_SUFFIX    <- "_assoc.tsv"
ASSOC_GENE_COL  <- "GeneName"
ASSOC_PVAL_COL  <- "Pvalue"
ASSOC_STAT_COL  <- "ZStat"
ASSOC_STAT_TYPE <- "z"

WEIGHTS_FALLBACK <- ""

SS_DIRS   <- c("/path/to/summary_stats")
SS_SUFFIX <- ".SS.txt"

GENO_BASES  <- c("/path/to/geno/output")
GENO_SUFFIX <- ".geno.txt"
