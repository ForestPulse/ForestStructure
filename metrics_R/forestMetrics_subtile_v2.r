#!/usr/bin/env Rscript
# Unified metrics + PDM/acquisition dates + EarlyAccess V3 models.
# Sources: forestMetrics_subtile.r, date_layer.r, 260922_apply_new_models.r.
#
# Run inside an existing single-node Slurm allocation / R container:
# Rscript forestMetrics_subtile_v2.r \
#   --force_cube_id X0064_Y0050 \
#   --copc_root /path/to/copc_root \
#   --out_root /path/to/output_root \
#   --spec_root /path/to/tree_species_fractions \
#   --workers 18
# Optional: --mask_file /path/to/DE_mask.cog.tif
#
# Input: <copc_root>/copc_tiles_<FORCE_CUBE_ID>/*.copc.laz
# Output: <out_root>/metrics_<FORCE_CUBE_ID>/tiles_<cube_id_lower>/*_metrics.tif
# Mosaic: <out_root>/metrics_<FORCE_CUBE_ID>/metrics_<cube_id_lower>.tif + .vrt
#
# Keep the original 10 m FORCE grid, 1 km tiles, HeightAboveGround normalization,
# point/CHM functions, 0.5 m p2r(0.25) CHM, and date/PDM function definitions.
# Removed: most_frequent_day()/gpstime_day and the three old model equations.
# Added: normalized V3 wood density, ccover, PDM, four temporal bands, V3 models.
# One Float64 GeoTIFF preserves both decimal metrics and integer YYYYMMDD dates.
# Temporal/PDM NA values remain nodata, not zero. Old metric NA-to-zero handling
# is retained for metric/model bands only. No per-layer rounding is introduced.
#
# GPS assumptions are UNCHANGED from date_layer.r: Adjusted Standard GPS Time,
# +1e9 seconds, then -18 seconds for UTC. Input X/Y must already be EPSG:3035.
# PDM smoothing is tile-local (no neighboring-COPC halo), as in date_layer.r.
# Optional forest mask affects the V3 product bands only; dates/PDM stay intact.
# Without --mask_file, no forest mask is applied (the path cannot be inferred).
# The model script's Windows paths, regional SN/TH overlap exports, integer
# rounding, and retiling of precomputed vertical layers are not executed here:
# this pipeline calculates the full metric stack and PDM directly from COPC.
# Old-schema tile outputs are recomputed automatically. Complete new-schema
# outputs are reused. For changed inputs or mask settings, use a fresh out_root.

# Avoid nested OpenMP/BLAS/GDAL threads inside independent future workers.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", GDAL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(lidR)
  library(sf)
  library(terra)
  library(future)
  library(future.apply)
})
lidR::set_lidr_threads(1L)
options(lidR.raster.default = "terra")

# ---- Argument parsing -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args) ||
      startsWith(args[i + 1L], "--")) stop("Supply one value for ", flag)
  args[i + 1L]
}

FORCE_CUBE_ID <- get_arg("--force_cube_id")
COPC_ROOT     <- get_arg("--copc_root")
OUT_ROOT      <- get_arg("--out_root")
SPEC_ROOT     <- get_arg("--spec_root")
MASK_FILE     <- get_arg("--mask_file")
N_WORKERS     <- suppressWarnings(as.integer(get_arg(
  "--workers", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "18")
)))

if (is.null(FORCE_CUBE_ID)) stop("Must provide --force_cube_id")
if (is.null(COPC_ROOT)) stop("Must provide --copc_root")
if (is.null(OUT_ROOT)) stop("Must provide --out_root")
if (is.null(SPEC_ROOT)) stop("Must provide --spec_root")
if (is.na(N_WORKERS) || N_WORKERS < 1L) stop("--workers must be >= 1")
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (!is.na(slurm_cpus) && N_WORKERS > slurm_cpus) {
  stop("--workers exceeds SLURM_CPUS_PER_TASK=", slurm_cpus)
}

# ---- Paths and reference grid -----------------------------------------------
id_lower    <- tolower(FORCE_CUBE_ID)
COPC_DIR    <- file.path(COPC_ROOT, paste0("copc_tiles_", FORCE_CUBE_ID))
FINAL_DIR   <- file.path(OUT_ROOT, paste0("metrics_", FORCE_CUBE_ID))
OUT_DIR     <- file.path(FINAL_DIR, paste0("tiles_", id_lower))
vrt_file    <- file.path(FINAL_DIR, paste0("metrics_", id_lower, ".vrt"))
mosaic_file <- file.path(FINAL_DIR, paste0("metrics_", id_lower, ".tif"))
failed_log  <- file.path(FINAL_DIR, paste0("failed_tiles_", id_lower, ".txt"))
SPEC_RASTER <- file.path(SPEC_ROOT, paste0("tree_species_fractions_",
                                         FORCE_CUBE_ID, ".cog.tif"))
if (!dir.exists(COPC_DIR)) stop("COPC directory does not exist: ", COPC_DIR)
if (!dir.exists(SPEC_ROOT)) stop("Species root does not exist: ", SPEC_ROOT)

# Prefer the original exact filename; otherwise use date_layer.r's recursive
# lookup and require one unambiguous match.
if (!file.exists(SPEC_RASTER)) {
  spec_candidates <- list.files(
    SPEC_ROOT, pattern = paste0(FORCE_CUBE_ID, ".*\\.(tif|tiff)$"),
    full.names = TRUE, recursive = TRUE, ignore.case = TRUE
  )
  if (length(spec_candidates) != 1L) {
    stop("Expected exactly one matching species raster; found ",
         length(spec_candidates), " below ", SPEC_ROOT, "\n",
         paste(spec_candidates, collapse = "\n"))
  }
  SPEC_RASTER <- spec_candidates[1L]
}
SPEC_RASTER <- normalizePath(SPEC_RASTER, mustWork = TRUE)
force_ref <- terra::rast(SPEC_RASTER)
if (terra::nlyr(force_ref) != 14L) {
  stop("Species raster must contain the 14 bands in the original species order.")
}
if (!isTRUE(sf::st_crs(terra::crs(force_ref)) == sf::st_crs(3035))) {
  stop("Species/FORCE reference must be EPSG:3035; no reprojection is performed.")
}
FORCE_XMIN <- terra::ext(force_ref)$xmin
FORCE_YMIN <- terra::ext(force_ref)$ymin
FORCE_CRS  <- terra::crs(force_ref)
if (!is.null(MASK_FILE)) {
  MASK_FILE <- normalizePath(MASK_FILE, mustWork = TRUE)
  mask_ref <- terra::rast(MASK_FILE)
  if (terra::nlyr(mask_ref) != 1L ||
      !isTRUE(sf::st_crs(terra::crs(mask_ref)) == sf::st_crs(FORCE_CRS))) {
    stop("--mask_file must be a single-band raster in the FORCE CRS.")
  }
  rm(mask_ref)
}
rm(force_ref)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_DIR <- normalizePath(OUT_DIR, mustWork = TRUE)
message("FORCE_CUBE_ID : ", FORCE_CUBE_ID)
message("COPC_DIR      : ", COPC_DIR)
message("OUT_DIR       : ", OUT_DIR)
message("SPEC_RASTER   : ", SPEC_RASTER)
message("N_WORKERS     : ", N_WORKERS)
message("FORCE origin : xmin=", format(FORCE_XMIN, digits = 16),
        " ymin=", format(FORCE_YMIN, digits = 16))
message("Forest mask   : ", if (is.null(MASK_FILE)) "not applied" else MASK_FILE)

RES <- 10
TILE_SIZE <- 1000
GPS_ADJUSTMENT <- 1e9
GPS_UTC_OFFSET <- 18

# Acquisition dates: integer YYYYMMDD in UTC; span: integer calendar days.
# acq_flag (original order/precedence):
# 0 = cell has points but no finite GPS times; 1 = a single acquisition day;
# 4 = endpoint dates have different phenology classes (priority over span);
# otherwise 2 = <=14 days, 3 = 15-90 days, 5 = 91-180 days, 6 = >180 days.
# Phenology: leaf-off Oct-Apr; transition May/Sep; leaf-on Jun-Aug.
# Cells with no points retain nodata in all temporal bands.
# pdmodes3_smooth: 1/2 from the original second-largest 3x3 neighborhood rule.

OUTPUT_BANDS <- c(
  "topheight", "canopycover", "vci", "meanH", "maxH", "over02", "over05",
  "over10", "over15", "over20", "q95", "dens", "ccover", "pdmodes3_smooth",
  "acq_date_min", "acq_date_max", "acq_span_days", "acq_flag", "vol", "biom", "ba"
)
V3_BANDS <- c("topheight", "meanH", "maxH", "ccover", "dens", "vol", "biom", "ba")

# ---- Original point and CHM functions (unchanged) ----------------------------
forest_metrics <- function(hag, th = 2, zmax_cap = 60, by = 1) {
  total_pts <- length(hag)
  n_above2  <- sum(hag > th)
  hag       <- hag[hag >= 0]

  if (length(hag) < 2L)
    return(list(zp95 = NA_real_, zabove2 = NA_real_, vci = NA_real_))

  p999 <- quantile(hag, 0.999)
  hag  <- hag[hag <= min(p999, zmax_cap)]

  if (length(hag) < 2L)
    return(list(zp95 = NA_real_, zabove2 = NA_real_, vci = NA_real_))

  list(
    zp95    = as.numeric(quantile(hag, 0.95)),
    zabove2 = n_above2 / total_pts * 100,
    vci     = VCI(hag, zmax = zmax_cap, by = by)
  )
}

# ---- CHM metrics function ---------------------------------------------------
chm_metrics <- function(raster) {
  vec <- na.omit(as.vector(raster))
  if (length(vec) == 0)
    return(c(meanH = 0, maxH = 0, over02 = 0, over05 = 0, over10 = 0, over15 = 0, 
             over20 = 0, q95 = 0))
  c(meanH  = mean(vec),
    maxH   = max(vec),
    over02 = sum(vec > 2)  / length(vec),
    over05 = sum(vec > 5)  / length(vec),
    over10 = sum(vec > 10) / length(vec),
    over15 = sum(vec > 15) / length(vec),
    over20 = sum(vec > 20) / length(vec),
    q95 = unname(quantile(vec, 0.95)))
}

# ---- Original PDM and date functions (unchanged) ------------------------------
number_of_modes <- function(z, bandwidth = 3) {
  z <- z[is.finite(z)]

  if (length(z) < 3L || diff(range(z)) <= 0) {
    return(NA_real_)
  }

  tryCatch({
    dens <- stats::density(z, bw = bandwidth)

    maxima <- which(diff(sign(diff(dens$y))) == -2L) + 1L

    as.numeric(length(maxima))
  }, error = function(e) {
    NA_real_
  })
}

secondlargest <- function(x) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  if (length(x) == 1L) {
    return(min(max(x[1L], 1), 2))
  }

  x <- sort(x, decreasing = TRUE)

  # Binary result:
  # 1 = no locally supported bimodal profile
  # 2 = at least two cells in the 3x3 window have >= 2 modes
  min(max(x[2L], 1), 2)
}

gps_adjusted_to_utc <- function(gps_adjusted) {
  gps_epoch <- as.POSIXct(
    "1980-01-06 00:00:00",
    tz = "UTC"
  )

  gps_epoch + gps_adjusted + GPS_ADJUSTMENT - GPS_UTC_OFFSET
}

phenology_code <- function(date_value) {
  month_value <- as.integer(format(date_value, "%m"))

  if (month_value %in% c(6L, 7L, 8L)) {
    return(3L) # leaf-on
  }

  if (month_value %in% c(5L, 9L)) {
    return(2L) # transition
  }

  1L # leaf-off
}

time_metrics <- function(gpstime) {
  gpstime <- gpstime[is.finite(gpstime)]

  if (!length(gpstime)) {
    return(list(
      acq_date_min = NA_integer_,
      acq_date_max = NA_integer_,
      acq_span_days = NA_integer_,
      acq_flag = 0L
    ))
  }

  utc_min <- gps_adjusted_to_utc(min(gpstime))
  utc_max <- gps_adjusted_to_utc(max(gpstime))

  date_min <- as.Date(utc_min, tz = "UTC")
  date_max <- as.Date(utc_max, tz = "UTC")

  span_days <- as.integer(date_max - date_min)

  pheno_min <- phenology_code(date_min)
  pheno_max <- phenology_code(date_max)

  flag <- if (span_days == 0L) {
    1L
  } else if (pheno_min != pheno_max) {
    4L
  } else if (span_days <= 14L) {
    2L
  } else if (span_days <= 90L) {
    3L
  } else if (span_days <= 180L) {
    5L
  } else {
    6L
  }

  list(
    acq_date_min = as.integer(format(date_min, "%Y%m%d")),
    acq_date_max = as.integer(format(date_max, "%Y%m%d")),
    acq_span_days = span_days,
    acq_flag = as.integer(flag)
  )
}

# ---- Raster geometry and restart validation ---------------------------------
# Crop/extend whole cells only. Never interpolate temporal/categorical metrics.
fit_tile_grid <- function(x, template) {
  x <- terra::crop(x, terra::ext(template), snap = "near")
  x <- terra::extend(x, terra::ext(template), snap = "near")
  if (!terra::compareGeom(x, template, stopOnError = FALSE)) {
    stop("Metric raster does not match the 10 m FORCE tile grid.")
  }
  x
}

is_valid_output <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) ||
      file.info(path)$size <= 0) return(FALSE)
  tryCatch({
    r <- terra::rast(path)
    ok <- identical(names(r), OUTPUT_BANDS) &&
      terra::nlyr(r) == length(OUTPUT_BANDS) &&
      nrow(r) == TILE_SIZE / RES && ncol(r) == TILE_SIZE / RES &&
      all(abs(terra::res(r) - RES) < 1e-8) &&
      all(terra::datatype(r) == "FLT8S")
    if (!ok) return(FALSE)
    # Read all 100 x 100 cells to detect truncated/unreadable tile data.
    v <- terra::values(r, mat = TRUE)
    nrow(v) == nrow(r) * ncol(r) && ncol(v) == length(OUTPUT_BANDS)
  }, error = function(e) FALSE)
}

# ---- Per-tile worker --------------------------------------------------------
process_tile <- function(copc_file, out_dir) {
  lidR::set_lidr_threads(1L)
  options(lidR.raster.default = "terra")
  tile_id <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(copc_file)))
  out_file <- file.path(out_dir, paste0(tile_id, "_metrics.tif"))
  if (is_valid_output(out_file)) {
    message("[SKIP] ", basename(out_file))
    return(out_file)
  }
  if (file.exists(out_file)) message("[REPLACE] Old/incomplete: ", basename(out_file))

  tmp_copc <- tempfile(pattern = paste0(tile_id, "_"), fileext = ".copc.laz")
  # Write beside the final tile, then rename only when writing/validation finish.
  tmp_out <- tempfile(pattern = paste0(tile_id, "_"), tmpdir = out_dir,
                      fileext = ".tif")
  message("[START] ", basename(copc_file))
  tryCatch({
    if (!file.copy(copc_file, tmp_copc, overwrite = TRUE)) {
      stop("Copy to worker scratch failed.")
    }
    # Keep the original full read: HeightAboveGround is an extra-byte attribute.
    las <- lidR::readLAS(tmp_copc)
    if (is.null(las) || lidR::is.empty(las)) stop("LAS object is empty.")
    if (!all(c("HeightAboveGround", "gpstime") %in% names(las@data))) {
      stop("Input must contain 'HeightAboveGround' and 'gpstime'.")
    }
    sf::st_crs(las) <- sf::st_crs(3035)
    las$Z <- las$HeightAboveGround

    nx <- round((las@header$`Min X` - FORCE_XMIN) / TILE_SIZE)
    ny <- round((las@header$`Min Y` - FORCE_YMIN) / TILE_SIZE)
    tile_xmin <- FORCE_XMIN + nx * TILE_SIZE
    tile_ymin <- FORCE_YMIN + ny * TILE_SIZE
    tile_ext <- terra::ext(tile_xmin, tile_xmin + TILE_SIZE,
                           tile_ymin, tile_ymin + TILE_SIZE)
    template <- terra::rast(ext = tile_ext, resolution = RES, crs = FORCE_CRS)

    ## Original point metrics
    p.m <- lidR::pixel_metrics(
      las, ~forest_metrics(Z), res = RES, start = c(tile_xmin, tile_ymin)
    )
    terra::crs(p.m) <- FORCE_CRS
    p.m <- fit_tile_grid(p.m, template)
    names(p.m) <- c("topheight", "canopycover", "vci")
    p.m <- terra::subst(p.m, NA, 0)

    ## Updated V3 wood density: same 14 species/order, normalized by total fraction
    spec.crop <- terra::crop(terra::rast(SPEC_RASTER), tile_ext)
    names(spec.crop) <- c("Fi", "Ki", "Ta", "Dg", "La", "Bu", "Ei",
                          "Ah", "Bi", "Er", "Pp", "AL", "nn", "xx")
    spec.crop$dens <-
      (spec.crop$Fi * 0.38 + spec.crop$Ki * 0.44 + spec.crop$Ta * 0.42 +
       spec.crop$Dg * 0.44 + spec.crop$La * 0.49 + spec.crop$Bu * 0.60 +
       spec.crop$Ei * 0.60 + spec.crop$Ah * 0.51 + spec.crop$Bi * 0.55 +
       spec.crop$Er * 0.48 + spec.crop$Pp * 0.46 + spec.crop$AL * 0.48 +
       spec.crop$nn * 0.48 + spec.crop$xx * 0.48) /
      (spec.crop$Fi + spec.crop$Ki + spec.crop$Ta + spec.crop$Dg + spec.crop$La +
       spec.crop$Bu + spec.crop$Ei + spec.crop$Ah + spec.crop$Bi + spec.crop$Er +
       spec.crop$Pp + spec.crop$AL + spec.crop$nn + spec.crop$xx)
    spec.crop$dens[is.na(spec.crop$dens)] <- 0.48
    spec.aligned <- terra::resample(spec.crop$dens, p.m, method = "near")
    names(spec.aligned) <- "dens"

    ## Original 0.5 m CHM, aggregated to 10 m
    chm.template <- terra::disagg(p.m$topheight, 20)
    chm <- lidR::rasterize_canopy(las, p2r(0.25), res = chm.template)
    chm.agg <- terra::aggregate(chm, fact = 20, na.rm = TRUE,
                                fun = function(x, ...) chm_metrics(x))
    terra::crs(chm.agg) <- FORCE_CRS
    names(chm.agg) <- names(chm_metrics(1:4))
    chm.agg <- terra::resample(chm.agg, p.m, method = "near")
    chm.agg <- terra::subst(chm.agg, NA, 0)
    rm(chm, chm.template)

    ## PDM from date_layer.r; use the same normalized Z as all forest metrics
    las_above_ground <- lidR::filter_poi(las, Z >= 0.5)
    if (lidR::is.empty(las_above_ground)) {
      pdm <- terra::rast(template)
      terra::values(pdm) <- NA_real_
    } else {
      pdm <- lidR::pixel_metrics(
        las_above_ground,
        ~list(pdmodes3 = number_of_modes(Z, bandwidth = 3)),
        res = RES, start = c(tile_xmin, tile_ymin)
      )
      terra::crs(pdm) <- FORCE_CRS
      pdm <- fit_tile_grid(pdm, template)
    }
    names(pdm) <- "pdmodes3"
    pdm_smooth <- terra::focal(
      pdm, w = matrix(1, 3, 3), fun = secondlargest,
      na.policy = "omit", fillvalue = NA
    )
    names(pdm_smooth) <- "pdmodes3_smooth"
    rm(las_above_ground, pdm)
    gc()

    ## New temporal layers replace most_frequent_day()/gpstime_day.
    # Compute BEFORE removing las (the old script removed las too early).
    time_layers <- lidR::pixel_metrics(
      las, ~time_metrics(gpstime), res = RES, start = c(tile_xmin, tile_ymin)
    )
    terra::crs(time_layers) <- FORCE_CRS
    time_layers <- fit_tile_grid(time_layers, template)
    names(time_layers) <- c("acq_date_min", "acq_date_max", "acq_span_days", "acq_flag")
    rm(las)
    gc()
    if (!terra::compareGeom(pdm_smooth, time_layers, stopOnError = FALSE)) {
      stop("PDM and temporal layers are not geometrically identical.")
    }

    ## Preserve every original metric; add CHM-based ccover as a separate band
    m <- c(p.m, chm.agg, spec.aligned)
    m$ccover <- m$over02 * 100

    ## New V3 equations copied from 260922_apply_new_models.r
    m$vol <- 4.0543948 * m$meanH ^ 1.3322654 * m$dens ^ -0.8135485
    m$biom <- 6.1049792354 * m$meanH ^ 1.2452957 * m$dens ^ 0.1845087
    m$ba <- 2.3898955 * m$meanH ^ 0.7129061 * m$dens ^ -0.7446183
    m <- terra::subst(m, NA, 0)

    ## Optional original V3 forest-mask rule, only for V3 product bands
    if (!is.null(MASK_FILE)) {
      mask.crop <- terra::crop(terra::rast(MASK_FILE), tile_ext)
      mask.aligned <- terra::resample(mask.crop, template, method = "near")
      v3 <- m[[V3_BANDS]]
      v3[mask.aligned == 0] <- 0
      m[[V3_BANDS]] <- v3
    }
    m <- c(m, pdm_smooth, time_layers)
    m <- m[[OUTPUT_BANDS]]

    terra::writeRaster(
      m, tmp_out, overwrite = TRUE, NAflag = -9999, datatype = "FLT8S",
      gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256")
    )
    if (!is_valid_output(tmp_out)) stop("Written tile failed validation.")
    if (!file.rename(tmp_out, out_file)) stop("Could not finalize tile output.")
    # Avoid stale GDAL statistics from an older file at this same path.
    if (file.exists(paste0(out_file, ".aux.xml"))) unlink(paste0(out_file, ".aux.xml"))
    message("[DONE] ", basename(out_file))
    out_file
  }, error = function(e) {
    message("[ERROR] ", basename(copc_file), " : ", conditionMessage(e))
    NA_character_
  }, finally = {
    unlink(c(tmp_copc, tmp_out, paste0(tmp_out, ".aux.xml")))
  })
}

# ---- Discover and process tiles ---------------------------------------------
copc_files <- list.files(COPC_DIR, pattern = "\\.copc\\.laz$",
                         full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
if (!length(copc_files)) stop("No COPC tiles found in ", COPC_DIR)
copc_files <- normalizePath(copc_files, mustWork = TRUE)
message(length(copc_files), " COPC tiles found")
N_WORKERS <- min(N_WORKERS, length(copc_files))
future::plan(future::multisession, workers = N_WORKERS)
results <- tryCatch(
  future.apply::future_lapply(
    copc_files, process_tile, out_dir = OUT_DIR,
    # Explicitly export only functions/scalars/paths, never a live SpatRaster.
    future.globals = list(
      forest_metrics = forest_metrics, chm_metrics = chm_metrics,
      number_of_modes = number_of_modes, secondlargest = secondlargest,
      gps_adjusted_to_utc = gps_adjusted_to_utc, phenology_code = phenology_code,
      time_metrics = time_metrics, fit_tile_grid = fit_tile_grid,
      is_valid_output = is_valid_output, process_tile = process_tile,
      FORCE_XMIN = FORCE_XMIN, FORCE_YMIN = FORCE_YMIN, FORCE_CRS = FORCE_CRS,
      SPEC_RASTER = SPEC_RASTER, MASK_FILE = MASK_FILE,
      RES = RES, TILE_SIZE = TILE_SIZE, GPS_ADJUSTMENT = GPS_ADJUSTMENT,
      GPS_UTC_OFFSET = GPS_UTC_OFFSET, OUTPUT_BANDS = OUTPUT_BANDS,
      V3_BANDS = V3_BANDS
    ),
    future.packages = c("lidR", "terra", "sf"), future.seed = TRUE
  ),
  finally = future::plan(future::sequential)
)

# ---- Validate and mosaic ----------------------------------------------------
ok <- vapply(results, function(x) {
  length(x) == 1L && !is.na(x) && is_valid_output(x)
}, logical(1))
valid_files <- unlist(results[ok], use.names = FALSE)
failed_files <- copc_files[!ok]
message("Processed : ", length(valid_files), " / ", length(copc_files))

if (length(failed_files)) {
  writeLines(failed_files, failed_log)
  message(length(failed_files), " tile(s) failed; see ", failed_log)
} else if (file.exists(failed_log)) {
  unlink(failed_log)
}

if (!length(valid_files)) stop("No valid tile outputs; aborting mosaic.")

message("Building VRT from ", length(valid_files), " tiles ...")
terra::vrt(valid_files, filename = vrt_file, overwrite = TRUE, set_names = TRUE)
mosaic_vrt <- terra::rast(vrt_file)
names(mosaic_vrt) <- OUTPUT_BANDS

# Retain the original optional stats step, safely quoting the filename.
gdalinfo <- Sys.which("gdalinfo")
if (nzchar(gdalinfo)) {
  status <- system2(gdalinfo, c("-stats", shQuote(vrt_file)),
                    stdout = "/dev/null", stderr = "")
  if (status != 0L) warning("gdalinfo -stats returned status ", status)
}
# Keep Float64 here too; otherwise the mosaic would round YYYYMMDD dates.
mosaic_tmp <- tempfile(pattern = "metrics_mosaic_", tmpdir = FINAL_DIR, fileext = ".tif")
tryCatch({
  terra::writeRaster(
    mosaic_vrt, mosaic_tmp, overwrite = TRUE, NAflag = -9999, datatype = "FLT8S",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
             "BIGTIFF=YES")
  )
  check <- terra::rast(mosaic_tmp)
  if (!identical(names(check), OUTPUT_BANDS) ||
      !all(terra::datatype(check) == "FLT8S")) stop("Mosaic failed validation.")
  rm(check)
  if (!file.rename(mosaic_tmp, mosaic_file)) stop("Could not finalize mosaic.")
  if (file.exists(paste0(mosaic_file, ".aux.xml"))) unlink(paste0(mosaic_file, ".aux.xml"))
}, finally = unlink(c(mosaic_tmp, paste0(mosaic_tmp, ".aux.xml"))))
message("Done -> ", mosaic_file)
