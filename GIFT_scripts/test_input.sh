#!/usr/bin/bash
#SBATCH --job-name=combine-matched
#SBATCH --partition=yanglab,day-long-cpu,week-long-cpu,month-long-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/combine_matched.%j.out
#SBATCH --error=logs/combine_matched.%j.err

# SLURM #SBATCH directives must be literal — edit above to change partitions/log paths.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if   [ -f "$SCRIPT_DIR/../config.sh" ]; then source "$SCRIPT_DIR/../config.sh"
elif [ -f "$SCRIPT_DIR/config.sh" ];    then source "$SCRIPT_DIR/config.sh"
else
  echo "ERROR: config.sh not found. Copy config.example.sh to config.sh and edit paths." >&2
  exit 1
fi

mkdir -p logs

echo "[$(date)] Starting combine job on node $(hostname)"

module load R
Rscript "$R_SCRIPT"

echo "[$(date)] Finished combine job"
