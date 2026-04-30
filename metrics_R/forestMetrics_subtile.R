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

FORCE_CUBE_ID <- "X0063_Y0050"  # get_arg("--force_cube_id") # modified FOR LOCAL TESTING
COPC_ROOT     <- get_arg("--copc_root", ".")
OUT_ROOT      <- get_arg("--out_root",  ".")
N_WORKERS     <- as.integer(get_arg("--workers", "32"))  # use many cores on full node

if (is.null(FORCE_CUBE_ID)) {
  stop("Must provide --force_cube_id, e.g. --force_cube_id X0063_Y0048")
}

id_lower  <- tolower(FORCE_CUBE_ID)

COPC_DIR  <- file.path(COPC_ROOT, paste0("copc_tiles_", FORCE_CUBE_ID))
FINAL_DIR <- file.path(OUT_ROOT, paste0("metrics_",    FORCE_CUBE_ID))
OUT_DIR   <- file.path(FINAL_DIR, paste0("tiles_",     id_lower))
vrt_file  <- file.path(FINAL_DIR, paste0("metrics_",   id_lower, ".vrt"))
mosaic_file <- file.path(FINAL_DIR, paste0("metrics_", id_lower, ".tif"))
failed_log  <- file.path(FINAL_DIR, paste0("failed_tiles_", id_lower, ".txt"))

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Point metrics function -------------------------------------------------------
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
chm_metrics <- function(raster){
  
  vec <- na.omit(as.vector(raster))
  
  res <- c(meanH = mean(vec),
           maxH = max(vec),
           over02 = sum(vec>2)/length(vec),
           over05 = sum(vec>5)/length(vec),
           over20 = sum(vec>20)/length(vec))
  return(res)
}

# ---- Per-tile worker (with tempdir + parallel) ------------------------------
process_tile <- function(copc_file, out_dir) {
  tile_id  <- tools::file_path_sans_ext(
                tools::file_path_sans_ext(basename(copc_file)))
  out_file <- file.path(out_dir, paste0(tile_id, "_metrics.tif"))

  if (file.exists(out_file)) return(out_file)

  tmp_copc <- file.path(tempdir(), basename(copc_file))

  tryCatch({
    # copy COPC to fast local storage
    file.copy(copc_file, tmp_copc, overwrite = TRUE)
    
    ## read LAS
    ctg <- lidR::readLAScatalog(tmp_copc)
    lidR::opt_chunk_size(ctg)   <- 0
    lidR::opt_chunk_buffer(ctg) <- 0
    lidR::opt_progress(ctg)     <- FALSE

    xmin <- ctg@data$Min.X
    ymin <- ctg@data$Min.Y
    
    ## read species data
    spec.tile.name <- paste0("mnt/vast-nhr/projects/nhr_ni_starter_24350/data/", 
                             "tree_species_fractions/tree_species_fractions_",
                             FORCE_CUBE_ID,
                             ".cog.tif")
    spec.tile <- terra::rast(spec.tile.name)
    spec.crop <- terra::crop(spec.tile, terra::ext(p.m), snap = "in")
    # calculate density layer from species composition
    spec.crop$dens <-spec.crop$tree_species_fractions.cog_1 * 0.0038 +
      spec.crop$tree_species_fractions.cog_2 * 0.0044 +
      spec.crop$tree_species_fractions.cog_3 * 0.0042 +
      spec.crop$tree_species_fractions.cog_4 * 0.0044 +
      spec.crop$tree_species_fractions.cog_5 * 0.0048 +
      spec.crop$tree_species_fractions.cog_6 * 0.006 +
      spec.crop$tree_species_fractions.cog_7 * 0.0061 +
      spec.crop$tree_species_fractions.cog_8 * 0.0052 +
      spec.crop$tree_species_fractions.cog_9 * 0.0055 +
      spec.crop$tree_species_fractions.cog_10 * 0.0045 +
      spec.crop$tree_species_fractions.cog_11 * 0.0043 +
      spec.crop$tree_species_fractions.cog_12 * 0.005 +
      spec.crop$tree_species_fractions.cog_13 * 0.005 +
      spec.crop$tree_species_fractions.cog_14 * 0.005
    # all missing or erroneous values will be set to 0.5
    spec.crop[is.na(spec.crop$dens) | 
                spec.crop$dens < 0.38 | 
                spec.crop$dens > 0.61 ] <- 0.5
    

    ## pixel metrics
    p.m <- lidR::pixel_metrics(
      ctg,
      ~forest_metrics(HeightAboveGround),
      res   = 10,
      start = c(xmin, ymin)
    )
    p.m <- terra::crop(p.m, terra::ext(ctg))
    names(p.m) <- c("topheight", "canopycover", "vci")

    ## chm metrics
    chm.resolution.template <- terra::disagg(p.m$topheight, 20)
    # create chm
    chm <- lidR::rasterize_canopy(ctg, p2r(0.25), res = chm.resolution.template)
    # aggregate to 10m pixels
    chm.agg <- terra::aggregate(chm, fact = 20, fun = "chm_metrics")
    names(chm.agg) <- names(chm_metrics(c(1,2,3,4)))

    # merge into one raster
    spec.crop.aligned <- terra::resample(spec.crop, p.m) # This should not be necessary
    m <- c(p.m, chm.agg, spec.crop.aligned$dens)
    
    
    # calculate target variables
    m$vol <- 3.4838 * m$meanH^1.3921 * m$dens^-0.76431
    m$biom <- 4.988 * m$meanH^1.343 * m$dens^0.3047
    m$ba <- 2.1713 * m$meanH^0.75521 * m$dens^-0.71128

    
    # write
    terra::writeRaster(
      m, out_file, overwrite = TRUE,
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

# ---- 2. Parallel processing over tiles -------------------------------------
plan(multisession, workers = N_WORKERS)

results <- future_lapply(
  copc_files,
  process_tile,
  out_dir         = OUT_DIR,
  future.globals  = list(forest_metrics = forest_metrics,
                         chm_metrics    = chm_metrics,
                         process_tile   = process_tile),
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

if (length(valid_files) == 0L) {
  stop("No valid tile outputs; aborting mosaic.")
}

# ---- 4. Mosaic --------------------------------------------------------------
message("Building VRT from ", length(valid_files), " tiles ...")
terra::vrt(unlist(valid_files), filename = vrt_file,
           overwrite = TRUE, set_names = TRUE)
system(paste("gdalinfo -stats", vrt_file, "> /dev/null"))
mosaic_vrt <- terra::rast(vrt_file)
# names already set from first tile; keep as-is
terra::writeRaster(
  mosaic_vrt, mosaic_file, overwrite = TRUE,
  gdal = c("COMPRESS=LZW", "TILED=YES",
           "BLOCKXSIZE=512", "BLOCKYSIZE=512",
           "BIGTIFF=YES")
)

message("Done -> ", mosaic_file)
