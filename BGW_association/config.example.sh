# ============================================================================
# BGW association test configuration (shell / SLURM side)
#
# Copy this file to config.sh and edit:
#   cp config.example.sh config.sh
#
# config.sh is gitignored. It is sourced by slurm/run_association.sh and
# slurm/generate_geno_ss.sh.
# ============================================================================

# --- Repo / config ---------------------------------------------------------

# Absolute path to this module (the directory containing scripts/, slurm/).
PROJECT_DIR="/path/to/BGW_association"

# R config consumed by the association R scripts. The R scripts pick this up
# from the CONFIG_R environment variable, so export it here.
export CONFIG_R="${PROJECT_DIR}/config.R"

# --- Inputs ----------------------------------------------------------------

# Per-chromosome gene lists: <GENE_LIST_DIR>/Gene_list_<chr>.txt
GENE_LIST_DIR="/path/to/Gene_lists"

# BGW-xWAS weight directory (must match WEIGHT_DIR in config.R).
WEIGHT_DIR="/path/to/bgw/wkdir"

# --- Genotype / summary-stat generation step (generate_geno_ss.sh) ---------

# Per-gene (or per-region) VCFs used by Estep_mcmc -saveGeno.
# generate_geno_ss.sh expects ${TEMP_DIR}/<gene>.vcf.gz to exist.
TEMP_DIR="/path/to/temp_files"

# Where generated dosage files land (must match GENO_DIR in config.R).
GENO_OUT_DIR="/path/to/Test_Geno"

# Where generated summary-stat files land (must match SS_DIR in config.R).
SS_OUT_DIR="/path/to/Test_GWAS_SS"

# BGW-xWAS Estep_mcmc binary.
ESTEP_MCMC="/path/to/BGW-xWAS-SS/bin/Estep_mcmc"

# Placeholder trait/phenotype file required by Estep_mcmc -saveGeno.
TEST_PHENO="/path/to/trait.txt"

# --- Compute ---------------------------------------------------------------

CORES=8

# R module to load on your cluster. Leave empty ("") to skip `module load`.
R_MODULE="R"
