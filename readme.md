# ForestStructure

ALS point cloud preprocessing pipeline for reprojection and retiling of Airborne Laser Scanning (ALS) data to align with the [FORCE Datacube](https://force-eo.readthedocs.io/) framework. Developed as part of the [ForestPulse](https://github.com/ForestPulse) project (Work Package 4) at HAWK Göttingen.

---

## Background

ALS data from German state survey administrations is typically delivered as **1 km × 1 km tiles** in **EPSG:25832** (UTM Zone 32N). The FORCE Datacube uses a different spatial scheme: **30 km × 30 km tiles** in **EPSG:3035** (Lambert Azimuthal Equal Area / ETRS89).

To spatially align ALS point clouds with Sentinel-2 imagery within FORCE, the data must be reprojected and retiled accordingly.

> A pilot dataset from Thüringia (2020–2025, ~16,945 tiles, ~33 pts/m²) was used for development and benchmarking. The workflow is designed to scale to all 480 FORCE cube tiles covering Germany.

---

## Terminology

| Term | Description |
|---|---|
| **Tile** | Input ALS LiDAR tile — 1 km × 1 km in EPSG:25832 |
| **FORCE cube tile** | FORCE Datacube tile — 30 km × 30 km in EPSG:3035 (480 tiles cover Germany) |
| **Subtile** | Output COPC-LAZ tile — 1 km × 1 km, axis-aligned within a FORCE cube tile |

Each FORCE cube tile contains **900 subtiles** (30 × 30 grid).

---

## Pipeline Overview

| Stage | Script | Description |
|---|---|---|
| **A** | `generate_header_cache.py` | Indexes intersecting ALS tiles for each FORCE cube tile |
| **B** | `pdal_pc_level_reprojection.py` | Reprojects, retiles, and preprocesses point clouds (EPSG:25832 → EPSG:3035) |

---

## Requirements

- [Apptainer](https://apptainer.org/) with `pdal_container.sif`
- FORCE cube tile index: `tiles-dev-25832.csv` (columns: `id`, `left`, `right`, `top`, `bottom`)
- FORCE cube GeoPackage: `DE_force-cube_25832.gpkg`
- Access to an HPC cluster (tested on GWDG HPC, `standard96` / `medium96s` partitions)

> **Note:** All paths in the commands below must be adapted to your user account and directory structure.

---

## Usage

### Interactive Apptainer Shell (optional)

Open an interactive shell within the container for testing:

```bash
apptainer shell \
  --bind /projects/extern/nhr/.../data:/mnt/data \
  ~/container/pdal_container.sif
```

> **Important:** Always resolve symlinks with `realpath` before binding paths on GWDG:
> ```bash
> realpath /projects/extern/nhr/nhr_ni/nhr_ni_starter/nhr_ni_starter_24350/dir.project/data/raw/th/Thuringia_2020_2025_tiles
> ```

---

## Stage A — Generate Header Cache

`generate_header_cache.py` scans input LAZ headers and identifies tiles that intersect the bounding extent of a given FORCE cube tile. This avoids reloading all input tiles on every run.

**Required inputs:**
- Path to input LAZ tiles
- `tiles-dev-25832.csv` — FORCE cube tile index in EPSG:25832
- A target FORCE cube ID (e.g., `X0061_Y0047`)

### Single tile

```bash
module load apptainer

apptainer exec \
  --bind /mnt/vast-nhr/projects/nhr_ni_starter_24350/data/shared/Thuringia_2020_2025_tiles:/mnt/data \
  --bind /user/pbasnet/u18501:/user/pbasnet/u18501 \
  /$PROJECT/apptainer/pdal_container.sif \
  python /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/generate_header_cache.py \
    --input /mnt/data \
    --force_csv /user/pbasnet/u18501/ForestPulse_processing/tiles-dev-25832.csv \
    --force_cube_id X0064_Y0050 \
    --workers 16 \
    --output /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/header_cache
```

To target a different FORCE cube tile, change only `--force_cube_id`.

### Batch processing (multiple tiles)

Save the following as `all_header_cache.sh` and run with `bash all_header_cache.sh`:

```bash
#!/bin/bash

CUBE_IDS=(
  "X0064_Y0050" "X0064_Y0051" "X0064_Y0052"
  "X0065_Y0047" "X0065_Y0048"
  "X0067_Y0049" "X0067_Y0050"
  "X0068_Y0048" "X0068_Y0049"
)

for ID in "${CUBE_IDS[@]}"; do
  echo "========================================"
  echo "Processing $ID..."
  echo "========================================"

  apptainer exec \
    --bind /mnt/vast-nhr/projects/nhr_ni_starter_24350/data/shared/Thuringia_2020_2025_tiles:/mnt/data \
    --bind /user/pbasnet/u18501:/user/pbasnet/u18501 \
    /$PROJECT/apptainer/pdal_container.sif \
    python /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/generate_header_cache.py \
      --input /mnt/data \
      --force_csv /user/pbasnet/u18501/ForestPulse_processing/tiles-dev-25832.csv \
      --force_cube_id "$ID" \
      --workers 16 \
      --output /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/header_cache
done

echo "All tasks completed!"
```

---

## Stage B — Reproject, Retile, and Preprocess

`pdal_pc_level_reprojection.py` reprojects point clouds from EPSG:25832 to EPSG:3035 and writes output as COPC-LAZ subtiles using PDAL.

**Required inputs:**
- `DE_force-cube_25832.gpkg` — FORCE cube GeoPackage
- Header cache from Stage A
- Input LAZ directory

### B.1 — Direct execution (login node, testing only)

> ⚠️ Not recommended for production. Use SLURM for full processing runs.

```bash
apptainer exec \
  --bind /mnt/vast-nhr/projects/nhr_ni_starter_24350/data/shared/Thuringia_2020_2025_tiles:/mnt/data \
  --bind /user/pbasnet/u18501:/user/pbasnet/u18501 \
  /$PROJECT/apptainer/pdal_container.sif \
  python /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/pdal_pc_level_reprojection.py \
    --force_gpkg /user/pbasnet/u18501/ForestPulse_processing/DE_force-cube_25832.gpkg \
    --force_layer epsg25832 \
    --force_cube_id X0065_Y0049 \
    --input_laz_dir /mnt/data \
    --output_laz_dir /user/pbasnet/u18501/test_tile_X0065_Y0049 \
    --header_cache /user/pbasnet/u18501/ForestPulse_processing/Test_options/pdal_pc_level_reprojection/header_cache/header_cache_X0065_Y0049.pkl \
    --num_workers 8 \
    --src_epsg 25832 \
    --dst_epsg 3035 \
    --start_id 0 \
    --end_id 899
```

**Expected output:**

FORCE cube X0065_Y0049 bounds (EPSG:3035): (4406026.36, 3074919.61, 4436026.36, 3104919.61)
Processing 900 subtiles with 8 workers.
Success (0, 0): .../las_4406026_3074920_x0065_y0049.copc.laz
Success (0, 1): .../las_4406026_3075920_x0065_y0049.copc.laz


**Output filename convention:**

las_<xmin><ymin><force_cube_id>.copc.laz


- `<xmin>` — minimum X of the subtile in EPSG:3035 (e.g., `4406026`)
- `<ymin>` — minimum Y of the subtile in EPSG:3035 (e.g., `3074920`)
- `<force_cube_id>` — FORCE cube tile ID (e.g., `X0065_Y0049`)

### B.2 — SLURM job submission (recommended)

For full production runs, submit via SLURM:

```bash
sbatch processing_slurm.sh --force_cube_id X0065_Y0049
```

**Benchmark (GWDG HPC, `standard96` partition):**
- Job arrays of 15, each processing 60 subtiles (30 in parallel)
- ~1 hour per FORCE cube tile
- ~1,700 core-hours per tile
- Estimated demand for 480 tiles (nationwide): ~816,000 core-hours

Log files are written to the working directory from which the SLURM job is submitted.

---

## Stage C — HPC Resource Monitoring

### During runtime

```bash
# Basic job statistics
sstat -a

# With horizontal scrolling
sstat -a | less -S
```

### After job completion

```bash
# Human-readable efficiency summary
module load py-reportseff
reportseff -u $USER

# Detailed Slurm accounting
sacct -j <job_id> \
  -o JobID,State,ExitCode,Elapsed,CPUTime,TotalCPU,NNodes,NCPUS,AllocCPUS,ReqMem,MaxRSS,Timelimit

# Export to CSV
sacct -j <job_id> \
  -o JobID,State,ExitCode,Elapsed,CPUTime,TotalCPU,NNodes,NCPUS,AllocCPUS,ReqMem,MaxRSS,Timelimit \
  --parsable2 > <job_id>_usage.csv
```

---

## Repository Structure

ForestStructure/
├── generate_header_cache.py # Stage A: header indexing
├── pdal_pc_level_reprojection.py # Stage B: reprojection and retiling
├── processing_slurm.sh # Stage B: SLURM job script
├── SLURM_Preprocessing.sh # Additional SLURM preprocessing script
├── Preprocessing_PDAL.json # PDAL pipeline configuration
├── generate_tiles.py # Tile generation utility
├── header_cache/ # Cached header outputs (Stage A)
├── force-cube/ # FORCE cube reference files
├── logs/ # Processing logs
└── README.md


---

## Data Source

Pilot ALS data for Thüringia (2020–2025) accessed from:
[Geoportal Thüringen — Höhendaten](https://geoportal.thueringen.de/gdi-th/download-offene-geodaten/download-hoehendaten)

- 16,945 tiles, 1 km × 1 km each
- Average point density: ~33 pts/m²
- Intersects 35 FORCE cube tiles

---

## License

To be defined by the ForestPulse project consortium.

---

## Contact

Developed at [HAWK Göttingen](https://www.hawk.de)/ForestPulse Work Package 4.

