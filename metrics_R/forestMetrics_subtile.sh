#!/bin/bash
#SBATCH --job-name=forestMetrics_subtile
#SBATCH --partition=medium96s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=250G
#SBATCH --time=00:59:00
#SBATCH --array=1-3
#SBATCH --output=logs/forestMetrics_subtile_%A_%a.out
#SBATCH --error=logs/forestMetrics_subtile_%A_%a.err

module purge
module load gcc/14.2.0
module load udunits/2.2.28
module load gdal/3.10.0
module load r/4.5.2

cd /user/pbasnet/u18501/ForestPulse_processing/Test_options/metrics_R

FORCE_CUBE_IDS=(X0063_Y0049 X0063_Y0050 X0063_Y0051
)

FID=${FORCE_CUBE_IDS[$SLURM_ARRAY_TASK_ID-1]}

Rscript --vanilla forestMetrics_subtile.R \
  --force_cube_id "$FID" \
  --copc_root /user/pbasnet/u18501/th_data/Thuringia_processed_tiles \
  --out_root  /user/pbasnet/u18501/th_data/Thuringia_forestMetrics \
  --workers "$SLURM_CPUS_PER_TASK"
