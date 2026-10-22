#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lidR)
  library(terra)
  library(sf)
  library(future)
  library(future.apply)
})

# =============================================================================
# Usage:
#
# Rscript date_layer.r \
#   --force_cube_id X0064_Y0050 \
#   --copc_root /path/to/copc_root \
#   --out_root /path/to/output_root \
#   --spec_root /path/to/tree_species_fractions \
#   --workers 10
#
# Input COPC directory:
#   <copc_root>/copc_tiles_<FORCE_CUBE_ID>/*.copc.laz
#
# The script searches recursively under <spec_root> for exactly one GeoTIFF
# whose filename contains the current FORCE cube ID.
#
# Per-tile output:
#   <out_root>/metrics_<FORCE_CUBE_ID>/tiles_<cube_id_lower>/*_pdm_time.tif
#
# Output bands:
#   1. pdmodes3_smooth
#   2. acq_date_min
#   3. acq_date_max
#   4. acq_span_days
#   5. acq_flag
# =============================================================================

# ---- Arguments ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)

  if (!length(i)) {
    return(default)
  }

  if (i == length(args)) {
    stop("No value supplied for ", flag)
  }

  args[i + 1L]
}

FORCE_CUBE_ID <- get_arg("--force_cube_id")
COPC_ROOT     <- get_arg("--copc_root")
OUT_ROOT      <- get_arg("--out_root")
SPEC_ROOT     <- get_arg("--spec_root")
N_WORKERS     <- as.integer(get_arg("--workers", "18"))

if (is.null(FORCE_CUBE_ID)) stop("Missing --force_cube_id")
if (is.null(COPC_ROOT)) stop("Missing --copc_root")
if (is.null(OUT_ROOT)) stop("Missing --out_root")
if (is.null(SPEC_ROOT)) stop("Missing --spec_root")

if (is.na(N_WORKERS) || N_WORKERS < 1L) {
  stop("--workers must be an integer >= 1")
}

message("FORCE_CUBE_ID: ", FORCE_CUBE_ID)
message("COPC_ROOT:     ", COPC_ROOT)
message("OUT_ROOT:      ", OUT_ROOT)
message("SPEC_ROOT:     ", SPEC_ROOT)
message("N_WORKERS:     ", N_WORKERS)

# ---- Paths -------------------------------------------------------------------

id_lower  <- tolower(FORCE_CUBE_ID)
COPC_DIR  <- file.path(COPC_ROOT, paste0("copc_tiles_", FORCE_CUBE_ID))
FINAL_DIR <- file.path(OUT_ROOT, paste0("metrics_", FORCE_CUBE_ID))
OUT_DIR   <- file.path(FINAL_DIR, paste0("tiles_", id_lower))

FAILED_LOG <- file.path(
  FINAL_DIR,
  paste0("failed_tiles_", id_lower, "_pdm_time.txt")
)

dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(COPC_DIR)) {
  stop("COPC directory does not exist: ", COPC_DIR)
}

if (!dir.exists(SPEC_ROOT)) {
  stop("Species-raster root directory does not exist: ", SPEC_ROOT)
}

# ---- Find matching species raster -------------------------------------------

# Searches SPEC_ROOT recursively and requires exactly one matching .tif/.tiff.
# Example valid filenames:
#   tree_species_fractions_X0064_Y0050.tif
#   species_x0064_y0050.tif
#   X0064_Y0050_species_fractions.tiff

spec_candidates <- list.files(
  path = SPEC_ROOT,
  pattern = paste0(FORCE_CUBE_ID, ".*\\.(tif|tiff)$"),
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (!length(spec_candidates)) {
  stop(
    "No species raster matching FORCE cube ID '", FORCE_CUBE_ID,
    "' was found below:\n",
    SPEC_ROOT
  )
}

if (length(spec_candidates) > 1L) {
  stop(
    "More than one species raster matches FORCE cube ID '",
    FORCE_CUBE_ID,
    "'. Refusing to choose one arbitrarily:\n",
    paste(spec_candidates, collapse = "\n")
  )
}

SPEC_RASTER <- spec_candidates[1L]

message("Species/FORCE raster selected: ", SPEC_RASTER)

# ---- FORCE cube grid from species raster ------------------------------------

force_ref <- terra::rast(SPEC_RASTER)

FORCE_XMIN <- terra::ext(force_ref)$xmin
FORCE_YMIN <- terra::ext(force_ref)$ymin

FORCE_CRS <- terra::crs(force_ref)

if (is.na(FORCE_CRS) || FORCE_CRS == "") {
  stop("The selected species/FORCE reference raster has no CRS: ", SPEC_RASTER)
}

message(
  "FORCE origin: xmin=", format(FORCE_XMIN, digits = 16),
  " ymin=", format(FORCE_YMIN, digits = 16)
)

message("FORCE CRS: ", FORCE_CRS)

# ---- Settings ----------------------------------------------------------------

RES       <- 10
TILE_SIZE <- 1000

# GPS time assumption:
# Values around 3.9e8–4.3e8 indicate Adjusted Standard GPS Time:
# standard GPS seconds since 1980-01-06 minus 1e9.
GPS_ADJUSTMENT <- 1e9
GPS_UTC_OFFSET <- 18

# acq_flag:
# 0 = points exist in the raster cell but all GPS-time values are invalid
# 1 = one acquisition day
# 2 = <= 14 days and same phenology class
# 3 = 15-90 days and same phenology class
# 4 = minimum and maximum dates have different phenology classes
# 5 = 91-180 days and same phenology class
# 6 = > 180 days
#
# Phenology:
# leaf-off   = October-April
# transition = May, September
# leaf-on    = June-August

# ---- Metric functions --------------------------------------------------------

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

# ---- Output validation -------------------------------------------------------

is_valid_output <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    return(FALSE)
  }

  tryCatch({
    r <- terra::rast(path)

    terra::nlyr(r) == 5L &&
      identical(
        names(r),
        c(
          "pdmodes3_smooth",
          "acq_date_min",
          "acq_date_max",
          "acq_span_days",
          "acq_flag"
        )
      )
  }, error = function(e) {
    FALSE
  })
}

# ---- Tile processing ---------------------------------------------------------

process_tile <- function(
    copc_file,
    out_dir,
    force_xmin,
    force_ymin,
    force_crs
) {
  tile_id <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(copc_file))
  )

  out_file <- file.path(
    out_dir,
    paste0(tile_id, "_pdm_time.tif")
  )

  if (is_valid_output(out_file)) {
    message("[SKIP] ", basename(out_file))
    return(out_file)
  }

  if (file.exists(out_file)) {
    message("[REPLACE] Invalid/incomplete file: ", basename(out_file))
    file.remove(out_file)
  }

  message("[START] ", basename(copc_file))

  tryCatch({
    # Reads X, Y, Z, GPS time, and classification only.
    las <- lidR::readLAS(
      copc_file,
      select = "xyztc"
    )

    if (lidR::is.empty(las)) {
      stop("LAS object is empty.")
    }

    if (!"gpstime" %in% names(las@data)) {
      stop("No 'gpstime' attribute was loaded.")
    }

    # Assign the known FORCE CRS label.
    # This does not transform X/Y coordinates.
    sf::st_crs(las) <- sf::st_crs(force_crs)

    # ---- FORCE-aligned 1 km tile extent -------------------------------------

    # The species raster defines the authoritative FORCE origin.
    # Each COPC tile is snapped to a 1 km grid derived from that origin.
    nx <- round((las@header$`Min X` - force_xmin) / TILE_SIZE)
    ny <- round((las@header$`Min Y` - force_ymin) / TILE_SIZE)

    tile_xmin <- force_xmin + nx * TILE_SIZE
    tile_ymin <- force_ymin + ny * TILE_SIZE

    tile_ext <- terra::ext(
      tile_xmin,
      tile_xmin + TILE_SIZE,
      tile_ymin,
      tile_ymin + TILE_SIZE
    )

    message(
      "[GRID] ", tile_id,
      " xmin=", format(tile_xmin, digits = 16),
      " ymin=", format(tile_ymin, digits = 16)
    )

    # ---- PDM modes -----------------------------------------------------------

    # Assumes Z is already normalized height / HeightAboveGround.
    las_above_ground <- lidR::filter_poi(las, Z >= 0.5)

    if (lidR::is.empty(las_above_ground)) {
      pdm <- terra::rast(
        ext = tile_ext,
        resolution = RES,
        crs = force_crs
      )

      terra::values(pdm) <- NA_real_
      names(pdm) <- "pdmodes3"

    } else {
      pdm <- lidR::pixel_metrics(
        las_above_ground,
        ~list(
          pdmodes3 = number_of_modes(Z, bandwidth = 3)
        ),
        res = RES,
        start = c(tile_xmin, tile_ymin)
      )

      terra::crs(pdm) <- force_crs

      # Enforce the exact 1 km FORCE-aligned tile extent.
      pdm <- terra::crop(
        pdm,
        tile_ext,
        snap = "out"
      )

      names(pdm) <- "pdmodes3"
    }

    # 3x3 cells at 10 m resolution = 30 m x 30 m neighborhood.
    pdm_smooth <- terra::focal(
      pdm,
      w = matrix(1, 3, 3),
      fun = secondlargest,
      na.policy = "omit",
      fillvalue = NA
    )

    names(pdm_smooth) <- "pdmodes3_smooth"

    # ---- Per-pixel temporal information -------------------------------------

    time_layers <- lidR::pixel_metrics(
      las,
      ~time_metrics(gpstime),
      res = RES,
      start = c(tile_xmin, tile_ymin)
    )

    terra::crs(time_layers) <- force_crs

    # Enforce the same FORCE-derived 1 km extent.
    time_layers <- terra::crop(
      time_layers,
      tile_ext,
      snap = "out"
    )

    names(time_layers) <- c(
      "acq_date_min",
      "acq_date_max",
      "acq_span_days",
      "acq_flag"
    )

    # Both metric products must share precisely the same grid.
    # Do not resample: mismatch means an alignment problem.
    if (!terra::compareGeom(
      pdm_smooth,
      time_layers,
      stopOnError = FALSE
    )) {
      stop(
        "PDM and temporal layers are not geometrically identical. ",
        "Check FORCE origin and pixel_metrics start coordinates."
      )
    }

    out <- c(pdm_smooth, time_layers)

    # GeoTIFF requires one output datatype for this five-band file.
    # INT4S safely stores dates in YYYYMMDD format and -9999 nodata.
    terra::writeRaster(
      out,
      out_file,
      overwrite = TRUE,
      NAflag = -9999,
      datatype = "INT4S",
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES",
        "BLOCKXSIZE=256",
        "BLOCKYSIZE=256"
      )
    )

    message("[DONE]  ", basename(out_file))

    out_file

  }, error = function(e) {
    message(
      "[ERROR] ",
      basename(copc_file),
      ": ",
      conditionMessage(e)
    )

    NA_character_
  })
}

# ---- Discover input COPC tiles ----------------------------------------------

copc_files <- list.files(
  COPC_DIR,
  pattern = "\\.copc\\.laz$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (!length(copc_files)) {
  stop("No '*.copc.laz' files found in: ", COPC_DIR)
}

message("COPC tiles found: ", length(copc_files))

# ---- Parallel processing -----------------------------------------------------

# Avoid nested OpenMP/BLAS threading inside parallel R workers.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

future::plan(
  future::multisession,
  workers = N_WORKERS
)

results <- future.apply::future_lapply(
  copc_files,
  process_tile,
  out_dir = OUT_DIR,
  force_xmin = FORCE_XMIN,
  force_ymin = FORCE_YMIN,
  force_crs = FORCE_CRS,
  future.packages = c("lidR", "terra", "sf"),
  future.seed = TRUE
)

future::plan(future::sequential)

# ---- Processing summary ------------------------------------------------------

valid_files <- Filter(
  function(x) !is.na(x) && is_valid_output(x),
  results
)

input_stems <- tools::file_path_sans_ext(
  tools::file_path_sans_ext(basename(copc_files))
)

output_stems <- tools::file_path_sans_ext(
  gsub("_pdm_time\\.tif$", "", basename(unlist(valid_files)))
)

failed_files <- copc_files[
  !(input_stems %in% output_stems)
]

message("Completed: ", length(valid_files), " / ", length(copc_files))

if (length(failed_files)) {
  writeLines(failed_files, FAILED_LOG)
  message("Failed tiles written to: ", FAILED_LOG)
} else if (file.exists(FAILED_LOG)) {
  file.remove(FAILED_LOG)
}
