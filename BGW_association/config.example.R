# ============================================================================
# BGW association test configuration (R side)
#
# Copy this file to config.R and edit the paths to match your environment:
#   cp config.example.R config.R
#
# config.R is gitignored so your local, DUA-governed paths stay private.
# These variables are sourced by scripts/association_test.R and
# scripts/cis_association_test.R.
# ============================================================================

# Per-gene BGW-xWAS weight files, one subdirectory per gene:
#   <WEIGHT_DIR>/<gene>/<gene>_BGW_xQTL_weights.txt
WEIGHT_DIR <- "/path/to/bgw/wkdir"

# Per-gene genotype dosage files:
#   <GENO_DIR>/<gene>.geno.txt
GENO_DIR <- "/path/to/Test_Geno"

# Per-gene GWAS summary-statistic files:
#   <SS_DIR>/<gene>.SS.txt
SS_DIR <- "/path/to/Test_GWAS_SS"

# Output directory for per-chromosome association results (created if absent)
RESULT_DIR <- "/path/to/Results"
