#!/bin/bash
#SBATCH --job-name=forestMetrics_for_training
#SBATCH --partition=medium96s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=250G
#SBATCH --time=00:59:00
#SBATCH --array=1-35
#SBATCH --output=logs/forestMetrics_%A_%a.out
#SBATCH --error=logs/forestMetrics_%A_%a.err

module purge
module load gcc/14.2.0
module load udunits/2.2.28
module load gdal/3.10.0
module load r/4.5.2

mkdir -p logs

cd /user/pbasnet/u18501/ForestPulse_processing/Test_options/metrics_R

# --- Define all FORCE cube IDs here -----------------------------------------

FORCE_CUBE_IDS=(
  "X0061_Y0047" "X0061_Y0050" "X0062_Y0047" "X0062_Y0048" "X0062_Y0049"
  "X0062_Y0050" "X0062_Y0051" "X0063_Y0046" "X0063_Y0047" "X0063_Y0048"
  "X0063_Y0049" "X0063_Y0050" "X0063_Y0051" "X0063_Y0052" "X0064_Y0046"
  "X0064_Y0047" "X0064_Y0048" "X0064_Y0049" "X0064_Y0050" "X0064_Y0051"
  "X0064_Y0052" "X0065_Y0047" "X0065_Y0048" "X0065_Y0049" "X0065_Y0050"
  "X0065_Y0051" "X0066_Y0048" "X0066_Y0049" "X0066_Y0050" "X0066_Y0051"
  "X0067_Y0048" "X0067_Y0049" "X0067_Y0050" "X0068_Y0048" "X0068_Y0049"
)

FID=${FORCE_CUBE_IDS[$SLURM_ARRAY_TASK_ID - 1]}

# --- Paths defined here, passed to R via args --------------------------------
COPC_ROOT="/user/pbasnet/u18501/th_data/Thuringia_processed_tiles"
OUT_ROOT="/user/pbasnet/u18501/ForestPulse_processing/Test_options/metrics_R"
SPEC_ROOT="/mnt/vast-nhr/projects/nhr_ni_starter_24350/data/tree_species_fractions"

echo "[$(date)] Starting FORCE cube: $FID (task $SLURM_ARRAY_TASK_ID)"

Rscript --vanilla forestMetrics_for_training.R \
  --force_cube_id "$FID" \
  --copc_root     "$COPC_ROOT" \
  --out_root      "$OUT_ROOT" \
  --spec_root     "$SPEC_ROOT" \
  --workers       "$SLURM_CPUS_PER_TASK"

echo "[$(date)] Done: $FID"
