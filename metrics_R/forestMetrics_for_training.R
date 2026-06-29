#!/usr/bin/env Rscript

library(lidR)
library(sf)
library(terra)
library(future)
library(future.apply)

# ---- Argument parsing -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(paste("No value for", flag))
  args[i + 1]
}

FORCE_CUBE_ID <- get_arg("--force_cube_id")
COPC_ROOT     <- get_arg("--copc_root")
OUT_ROOT      <- get_arg("--out_root")
SPEC_ROOT     <- get_arg("--spec_root")
N_WORKERS     <- as.integer(get_arg("--workers", "18"))

if (is.null(FORCE_CUBE_ID)) stop("Must provide --force_cube_id")
if (is.null(COPC_ROOT))     stop("Must provide --copc_root")
if (is.null(OUT_ROOT))      stop("Must provide --out_root")
if (is.null(SPEC_ROOT))     stop("Must provide --spec_root")

message("FORCE_CUBE_ID : ", FORCE_CUBE_ID)
message("COPC_ROOT     : ", COPC_ROOT)
message("OUT_ROOT      : ", OUT_ROOT)
message("SPEC_ROOT     : ", SPEC_ROOT)
message("N_WORKERS     : ", N_WORKERS)

# ---- Paths ------------------------------------------------------------------
id_lower    <- tolower(FORCE_CUBE_ID)
COPC_DIR    <- file.path(COPC_ROOT, paste0("copc_tiles_",   FORCE_CUBE_ID))
FINAL_DIR   <- file.path(OUT_ROOT,  paste0("metrics_",      FORCE_CUBE_ID))
OUT_DIR     <- file.path(FINAL_DIR, paste0("tiles_",        id_lower))
vrt_file    <- file.path(FINAL_DIR, paste0("metrics_",      id_lower, ".vrt"))
mosaic_file <- file.path(FINAL_DIR, paste0("metrics_",      id_lower, ".tif"))
failed_log  <- file.path(FINAL_DIR, paste0("failed_tiles_", id_lower, ".txt"))
SPEC_RASTER <- file.path(SPEC_ROOT, paste0("tree_species_fractions_",
                                           FORCE_CUBE_ID, ".cog.tif"))

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- FORCE cube grid from species raster ------------------------------------
force_ref  <- terra::rast(SPEC_RASTER)
FORCE_XMIN <- terra::ext(force_ref)$xmin
FORCE_YMIN <- terra::ext(force_ref)$ymin

message("FORCE origin: xmin=", format(FORCE_XMIN, digits = 16),
        " ymin=", format(FORCE_YMIN, digits = 16))

# ---- Per-tile worker --------------------------------------------------------
process_tile <- function(copc_file, out_dir) {
  tile_id  <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(copc_file)))
  out_file <- file.path(out_dir, paste0(tile_id, "_metrics.tif"))
  
  if (file.exists(out_file)) return(out_file)
  
  tmp_copc <- file.path(tempdir(), basename(copc_file))
  
  tryCatch({
    file.copy(copc_file, tmp_copc, overwrite = TRUE)
    
    las <- lidR::readLAS(tmp_copc)
    sf::st_crs(las) <- sf::st_crs(3035)
    las$Z <- las$HeightAboveGround
    
    ## FORCE-aligned tile extent
    nx <- round((las@header$`Min X` - FORCE_XMIN) / 1000)
    ny <- round((las@header$`Min Y` - FORCE_YMIN) / 1000)
    
    tile_xmin <- FORCE_XMIN + nx * 1000
    tile_ymin <- FORCE_YMIN + ny * 1000
    tile_ext  <- terra::ext(
      tile_xmin, tile_xmin + 1000,
      tile_ymin, tile_ymin + 1000
    )
    ################Bis hier schon überarbeitet
    
    ## species data as raster template
    spec.tile <- terra::rast(SPEC_RASTER)
    spec.crop <- terra::crop(spec.tile, tile_ext)
    raster.template <- terra::disagg(spec.crop[[1]], 10)
    
    # create canopy height model
    chm <- lidR::rasterize_canopy(las, res = raster.template)
    
    # aggregate
    chm.mean10 <- terra::aggregate(chm, fact  =10, na.rm = TRUE, fun = mean)
    
    chm.max05 <- terra::aggregate(chm, fact = 5, na.rm = TRUE, fun = max)
    twostep <- terra::aggregate(chm.max05, fact = 2, na.rm = TRUE, fun = mean)
    
    ext(chm.mean10) == ext(raster.template)
    
    rm(las); gc()
    
    # combine 10m-mean and two-step aggregation
    metrics <- c(chm.mean10, twostep)
    names(metrics) <-c("mean10", "max5mean10")
    terra::crs(metrics) <- "EPSG:3035"
    
    
    m <- terra::subst(m, NA, 0)
    
    terra::writeRaster(
      m, out_file, overwrite = TRUE,
      NAflag = -9999,
      gdal = c("COMPRESS=LZW", "TILED=YES",
               "BLOCKXSIZE=256", "BLOCKYSIZE=256")
    )
    
    out_file
    
  }, error = function(e) {
    message("[ERROR] ", basename(copc_file), " : ", conditionMessage(e))
    NA_character_
  }, finally = {
    if (file.exists(tmp_copc)) unlink(tmp_copc)
  })
}

# ---- 1. Discover tiles ------------------------------------------------------
copc_files <- list.files(
  COPC_DIR, pattern = "\\.copc\\.laz$",
  full.names = TRUE, recursive = FALSE
)
message(length(copc_files), " COPC tiles found in ", COPC_DIR)

if (length(copc_files) == 0L) stop("No COPC tiles found in ", COPC_DIR)

# ---- 2. Parallel processing -------------------------------------------------
plan(multisession, workers = N_WORKERS)

results <- future_lapply(
  copc_files,
  process_tile,
  out_dir        = OUT_DIR,
  future.globals = list(
    forest_metrics = forest_metrics,
    chm_metrics    = chm_metrics,
    process_tile   = process_tile,
    FORCE_CUBE_ID  = FORCE_CUBE_ID,
    SPEC_RASTER    = SPEC_RASTER,
    FORCE_XMIN     = FORCE_XMIN,
    FORCE_YMIN     = FORCE_YMIN
  ),
  future.packages = c("lidR", "terra"),
  future.seed     = TRUE
)

plan(sequential)

# ---- 3. Validate results ----------------------------------------------------
valid_files <- Filter(function(x) !is.na(x) && file.exists(x), results)

processed_stems <- gsub("_metrics\\.tif$", "", basename(unlist(valid_files)))
input_stems     <- tools::file_path_sans_ext(
  tools::file_path_sans_ext(basename(copc_files)))
failed_files    <- copc_files[!input_stems %in% processed_stems]

message("Processed : ", length(valid_files), " / ", length(copc_files))

if (length(failed_files) > 0) {
  message("Failed : ", length(failed_files), " (see ", failed_log, ")")
  writeLines(failed_files, failed_log)
} else {
  message("All tiles completed successfully")
  if (file.exists(failed_log)) file.remove(failed_log)
}

if (length(valid_files) == 0L) stop("No valid tile outputs; aborting mosaic.")

# ---- 4. Mosaic --------------------------------------------------------------
message("Building VRT from ", length(valid_files), " tiles ...")
terra::vrt(unlist(valid_files), filename = vrt_file,
           overwrite = TRUE, set_names = TRUE)
system(paste("gdalinfo -stats", vrt_file, "> /dev/null"))
mosaic_vrt <- terra::rast(vrt_file)

terra::writeRaster(
  mosaic_vrt, mosaic_file, overwrite = TRUE,
  NAflag = -9999,
  gdal = c("COMPRESS=LZW", "TILED=YES",
           "BLOCKXSIZE=512", "BLOCKYSIZE=512",
           "BIGTIFF=YES")
)

message("Done -> ", mosaic_file)