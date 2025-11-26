#!/bin/bash
#SBATCH --job-name=preprocessing
#SBATCH --partition=standard96
#SBATCH --time=00:10:00
#SBATCH --array=0-99         #provide the array number based on the tiles
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

module load miniforge3
source activate lidar

INPUT_DIR="/user/pbasnet/u18501/jupyterhub-gwdg/Preprocessing_Trial/slurmTest/laz_tiles"
OUTPUT_DIR="/user/pbasnet/u18501/jupyterhub-gwdg/Preprocessing_Trial/slurmTest/process_tiles"
TEMPLATE="/user/pbasnet/u18501/jupyterhub-gwdg/Preprocessing_Trial/slurmTest/Preprocessing_PDAL.json"

mkdir -p "$OUTPUT_DIR"

# Collect all .laz files
LAS_FILES=("$INPUT_DIR"/*.laz)

# Pick file using SLURM index
INPUT_FILE="${LAS_FILES[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$INPUT_FILE" ]; then
    echo "ERROR: No input file for index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Output name
BASENAME=$(basename "$INPUT_FILE" .laz)
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}.copc.laz"

# Create temporary pipeline file
TEMP_PIPE=$(mktemp)

sed -e "s|%INPUT%|$INPUT_FILE|g" \
    -e "s|%OUTPUT%|$OUTPUT_FILE|g" \
    "$TEMPLATE" > "$TEMP_PIPE"

# Run PDAL pipeline
pdal pipeline "$TEMP_PIPE"

# Cleanup
rm "$TEMP_PIPE"
