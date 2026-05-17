#!/usr/bin/bash
#SBATCH --job-name=xwas_gift_pipeline
#SBATCH --partition=yanglab,day-long-cpu,week-long-cpu,month-long-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/xwas_gift_pipeline.%j.out
#SBATCH --error=logs/xwas_gift_pipeline.%j.err

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if   [ -f "$SCRIPT_DIR/../config.sh" ]; then source "$SCRIPT_DIR/../config.sh"
elif [ -f "$SCRIPT_DIR/config.sh" ];    then source "$SCRIPT_DIR/config.sh"
else
  echo "ERROR: config.sh not found. Copy config.example.sh to config.sh and edit paths." >&2
  exit 1
fi

mkdir -p logs

PIPELINE_R="$SCRIPT_DIR/xwas_gift_pipeline.R"

START_STEP="${1:-1}"
END_STEP="${2:-5}"

echo "[$(date)] Running xWAS GIFT pipeline steps ${START_STEP} through ${END_STEP}"

module load R

for step in $(seq "$START_STEP" "$END_STEP"); do
  echo ""
  echo "[$(date)] === Step $step ==="
  Rscript "$PIPELINE_R" "$step"
  if [ $? -ne 0 ]; then
    echo "[$(date)] Step $step failed; aborting"
    exit 1
  fi
done

echo ""
echo "[$(date)] Pipeline complete. Submit GIFT jobs with:"
echo "  bash \$INPUT_BASE/GIFT_running/submit_all_regions.sh"
