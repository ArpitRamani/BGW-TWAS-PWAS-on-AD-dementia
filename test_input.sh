#!/usr/bin/bash
#SBATCH --job-name=combine-matched
#SBATCH --partition=yanglab,day-long-cpu,week-long-cpu,month-long-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/jyang51/YangLabData/aramani/GIFT2/logs/combine_matched.%j.out
#SBATCH --error=/home/jyang51/YangLabData/aramani/GIFT2/logs/combine_matched.%j.err

R_SCRIPT="/home/arama30/GIFT/retry_gift/GIFT_scripts/newGIFTGEN.R"

echo "[$(date)] Starting combine job on node $(hostname)"

module load R
Rscript "$R_SCRIPT"

echo "[$(date)] Finished combine job"
