#!/usr/bin/bash
#SBATCH --job-name=BGW-ASSO
#SBATCH --partition=PARTITION_NAME        # EDIT: your cluster partition(s)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-22                       # one job per chromosome (1..22)
#SBATCH --time=24:00:00
#SBATCH --output=logs/BGW.Asso.%A_%a.out   # relative to your submit dir
#SBATCH --error=logs/BGW.Asso.%A_%a.err
#
# Run the BGW association test as a per-chromosome array job.
#
# Submit from the module root (so the logs/ paths above resolve), e.g.:
#   cd /path/to/BGW_association
#   mkdir -p logs
#   sbatch slurm/run_association.sh
#
# Edit the SBATCH --partition line for your cluster. All data paths come
# from config.sh (see config.example.sh).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

GENE_LIST_FILE="${GENE_LIST_DIR}/Gene_list_${SLURM_ARRAY_TASK_ID}.txt"

# Pick the test to run. Default is the full (cis + trans) test.
# Comment / uncomment to switch to the cis-only test.
R_SCRIPT="${SCRIPT_DIR}/../scripts/association_test.R"
# R_SCRIPT="${SCRIPT_DIR}/../scripts/cis_association_test.R"

if [ -n "${R_MODULE}" ]; then
  module load "${R_MODULE}"
fi

echo "Running association test on ${GENE_LIST_FILE} ..."
Rscript "${R_SCRIPT}" "${GENE_LIST_FILE}" "${SLURM_ARRAY_TASK_ID}" "${CORES}"
echo "Done."
