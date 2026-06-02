#!/usr/bin/bash
#SBATCH --job-name=BGW-GENO
#SBATCH --partition=PARTITION_NAME        # EDIT: your cluster partition(s)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-22                       # one job per chromosome (1..22)
#SBATCH --time=24:00:00
#SBATCH --output=logs/BGW.Geno.%A_%a.out
#SBATCH --error=logs/BGW.Geno.%A_%a.err
#
# Generate per-gene genotype dosage (.geno.txt) and GWAS summary-stat
# (.SS.txt) files used as inputs to the association test.
#
# For each gene in the per-chromosome list it:
#   1. checks that a BGW-xWAS weight file exists and is non-empty,
#   2. validates the Trans column is integer-coded (0 = cis, 1 = trans),
#   3. runs Estep_mcmc -saveGeno on ${TEMP_DIR}/<gene>.vcf.gz,
#   4. moves the .geno.txt to GENO_OUT_DIR and .SS.txt to SS_OUT_DIR.
#
# Submit from the module root:
#   cd /path/to/BGW_association
#   mkdir -p logs
#   sbatch slurm/generate_geno_ss.sh
#
# Note: Estep_mcmc writes its outputs into the working directory; this script
# cd's into TEMP_DIR first, then moves the results. Adjust if your build
# writes elsewhere. All paths come from config.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

mkdir -p "${GENO_OUT_DIR}" "${SS_OUT_DIR}" "${TEMP_DIR}"
cd "${TEMP_DIR}"

gene_list="${GENE_LIST_DIR}/Gene_list_${SLURM_ARRAY_TASK_ID}.txt"

while read -r gene_name; do
  [ -z "${gene_name}" ] && continue
  echo "Processing gene ${gene_name} ..."

  weight_file="${WEIGHT_DIR}/${gene_name}/${gene_name}_BGW_xQTL_weights.txt"
  if [ ! -s "${weight_file}" ]; then
    echo "  weight file missing or empty, skipping"
    continue
  fi

  # Trans column (col 7) must be integer-coded (0 or 1).
  if ! awk 'NR>1 && $7!=0 && $7!=1 {exit 1}' "${weight_file}"; then
    echo "  non-integer Trans values, skipping"
    continue
  fi

  vcf="${TEMP_DIR}/${gene_name}.vcf.gz"
  if [ ! -s "${vcf}" ]; then
    echo "  VCF not found (${vcf}), skipping"
    continue
  fi

  echo "  generating genotype + summary-stat files ..."
  "${ESTEP_MCMC}" -vcf "${vcf}" -p "${TEST_PHENO}" -o "${gene_name}" \
    -GTfield GT -saveGeno -maf 0

  mv "${TEMP_DIR}/${gene_name}.geno.txt" "${GENO_OUT_DIR}/"
  mv "${TEMP_DIR}/${gene_name}.SS.txt"   "${SS_OUT_DIR}/"
done < "${gene_list}"

echo "Done."
