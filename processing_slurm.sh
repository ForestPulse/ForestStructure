#!/bin/bash
#SBATCH --job-name=ALS_FORCE_ALIGN
#SBATCH --partition=medium96s
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=250G
#SBATCH -t 01:59:00
#SBATCH --array=0-14
#SBATCH -C ssd
#SBATCH -o logs/%x_%A_%a.out

set -euo pipefail

# ---------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------
FORCE_CUBE_ID="${FORCE_CUBE_ID:-}"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --force_cube_id) FORCE_CUBE_ID="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$FORCE_CUBE_ID" ]; then
    echo "ERROR: --force_cube_id is required."
    echo "Usage: sbatch --export=ALL,FORCE_CUBE_ID=X0064_Y0049 ./processing_slurm.sh"
    exit 1
fi

module load apptainer

# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------
CONTAINER="/$PROJECT/apptainer/pdal_container.sif"
HOST_HOME="/user/pbasnet/u18501"
HOST_DATA="/user/pbasnet/u18501/th_data"

# Host path to the raw Thuringia tiles directory (to be bound as /mnt/data)
HOST_RAW_TILES="/mnt/vast-nhr/projects/nhr_ni_starter_24350/data/shared/Thuringia_2020_2025_tiles"
# Container-internal path to the tiles (matches header cache)
INDIR="/mnt/data"

WORKDIR="${HOST_HOME}/ForestPulse_processing/Test_options/pdal_pc_level_reprojection"
FORCE_GPKG="${HOST_HOME}/ForestPulse_processing/DE_force-cube_25832.gpkg"
FORCE_LAYER="epsg25832"
HEADER_CACHE="${WORKDIR}/header_cache/header_cache_${FORCE_CUBE_ID}.pkl"
FINAL_OUTDIR_HOST="${HOST_DATA}/Thuringia_processed_tiles/copc_tiles_${FORCE_CUBE_ID}"

SRC_EPSG=25832
DST_EPSG=3035

# ---------------------------------------------------------------
# Tile batching
# ---------------------------------------------------------------
TOTAL_TILES=900
NUM_TASKS=15
BATCH_SIZE=$(( TOTAL_TILES / NUM_TASKS ))
START_ID=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
END_ID=$(( START_ID + BATCH_SIZE - 1 ))
if [ "$END_ID" -ge "$TOTAL_TILES" ]; then END_ID=$(( TOTAL_TILES - 1 )); fi

# ---------------------------------------------------------------
# SSD scratch — task-scoped to avoid collisions between concurrent array tasks
# ---------------------------------------------------------------
if   [ -n "${LOCAL_TMPDIR:-}"      ]; then FAST_TMP="$LOCAL_TMPDIR"
elif [ -n "${SHARED_SSD_TMPDIR:-}" ]; then FAST_TMP="$SHARED_SSD_TMPDIR"
else                                       FAST_TMP="$TMPDIR"
fi

TASK_TMP="${FAST_TMP}/task_${SLURM_ARRAY_TASK_ID}"
LOCAL_STATIC="${TASK_TMP}/static"
LOCAL_IN="${TASK_TMP}/input"
LOCAL_OUT="${TASK_TMP}/output"

mkdir -p "$LOCAL_STATIC" "$LOCAL_IN" "$LOCAL_OUT" "$FINAL_OUTDIR_HOST" logs

# Free SSD scratch on exit regardless of success/failure
trap 'echo "[cleanup] Removing ${TASK_TMP}"; rm -rf "${TASK_TMP}"' EXIT

# ---------------------------------------------------------------
# Stage static files to SSD
# ---------------------------------------------------------------
cp "$FORCE_GPKG"   "$LOCAL_STATIC/"
cp "$HEADER_CACHE" "$LOCAL_STATIC/"

# ---------------------------------------------------------------
# Extract only the .laz filenames needed for THIS batch slice
# ---------------------------------------------------------------
echo "Extracting required input files (tiles ${START_ID}-${END_ID}) for cube ${FORCE_CUBE_ID}..."

apptainer exec \
    --bind "${TASK_TMP}:${TASK_TMP}" \
    --bind "${HOST_HOME}:${HOST_HOME}" \
    "$CONTAINER" \
    python3 -c "
import pickle, os
cache_file = '${LOCAL_STATIC}/' + os.path.basename('${HEADER_CACHE}')
with open(cache_file, 'rb') as f:
    cache = pickle.load(f)
raw_headers = cache['tiles'] if isinstance(cache, dict) else cache

# Only files for this batch — avoids copying the full cube per task
batch = raw_headers[${START_ID}:${END_ID}+1]
needed = sorted(set(os.path.basename(h['filename']) for h in batch))
for fn in needed:
    print(fn)
" > "${TASK_TMP}/needed_files.txt"

FILE_COUNT=$(wc -l < "${TASK_TMP}/needed_files.txt")
echo "Copying ${FILE_COUNT} files to local SSD..."

# rsync runs on HOST — source is real host tiles directory, dest is LOCAL_IN
rsync -a --files-from="${TASK_TMP}/needed_files.txt" \
    "${HOST_RAW_TILES}/" \
    "${LOCAL_IN}/"

# ---------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------
echo "======================================================"
echo "Job Array ID : ${SLURM_ARRAY_TASK_ID}"
echo "Cube ID      : ${FORCE_CUBE_ID}"
echo "Tile range   : ${START_ID} – ${END_ID}"
echo "Input files  : ${FILE_COUNT}"
echo "Workers      : ${SLURM_CPUS_PER_TASK}"
echo "Node         : $(hostname)"
echo "Temp dir     : ${TASK_TMP}"
echo "Start time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================"

# ---------------------------------------------------------------
# Run reprojection
# HOST_RAW_TILES bound as /mnt/data → matches header cache paths exactly
# TASK_TMP bound as /local_ssd for SSD-staged input/output/static
# ---------------------------------------------------------------
apptainer exec \
    --bind "${HOST_RAW_TILES}:/mnt/data" \
    --bind "${HOST_HOME}:${HOST_HOME}" \
    --bind "${TASK_TMP}:/local_ssd" \
    "$CONTAINER" \
    python "${WORKDIR}/pdal_pc_level_reprojection.py" \
        --force_gpkg     "/local_ssd/static/$(basename "$FORCE_GPKG")" \
        --force_layer    "$FORCE_LAYER" \
        --force_cube_id  "$FORCE_CUBE_ID" \
        --input_laz_dir  "/mnt/data" \
        --output_laz_dir "/local_ssd/output" \
        --header_cache   "/local_ssd/static/$(basename "$HEADER_CACHE")" \
        --num_workers    "$SLURM_CPUS_PER_TASK" \
        --src_epsg       "$SRC_EPSG" \
        --dst_epsg       "$DST_EPSG" \
        --start_id       "$START_ID" \
        --end_id         "$END_ID"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "Processing successful. Syncing outputs to ${FINAL_OUTDIR_HOST}..."
    rsync -a "${LOCAL_OUT}/" "${FINAL_OUTDIR_HOST}/"
    echo "Done at $(date '+%Y-%m-%d %H:%M:%S')"
else
    echo "Job failed with exit code ${EXIT_CODE}."
fi

exit $EXIT_CODE
