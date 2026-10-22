################################################################################
# Calculate modeled targetvariables from Felix's preliminary germanywide heights
################################################################################
# Timon Miesner, 22.10.26
#
# Models calculated from BWI data with data-driven outlier removal
# and without spatial balancing for sampling density

require(terra)
require(sf)

# define data paths
HEIGHTS_VRT <- "O:/a7bze/TimonMiesner/ForestPulse_EarlyAccess_Products_Oct26/top_canopy_height_2023/height_2023.vrt"
SPECIES_VRT <- "O:/a7bze/TimonMiesner/germany_TreeSpecies.vrt"
TREED_MASK <- "C:/Users/miesner/Data/FP-Earlyaccess/DE_mask.cog.tif"
FORCE_GRID <- "C:/Users/miesner/Data/FP-Earlyaccess/datacube_template_rasters/data/datacube-grid_DEU_10km.gpkg"
WOOD_DENS_CSV <- "C:/Users/miesner/Data/Global Wood Density Database v2/density_extracted.csv"
OUTPUT_DIR <- "~/v4_germanywide/2023/"


# read germany-wide data
force_grid_sf <- sf::read_sf(FORCE_GRID)
heights_de <- rast(HEIGHTS_VRT)
spec_de <- rast(SPECIES_VRT)
mask_de <- rast(TREED_MASK)
densities <- read.csv2(WOOD_DENS_CSV, sep  = "\t")
default_dens <- densities[densities$spec_group == "Default",]$dens_new

# loop over force tiles
for(tile_i in 1:nrow(force_grid_sf)){
  
  tile <- force_grid_sf[tile_i, ]
  
  # crop everything to current tile
  heights_tile <- terra::crop(heights_de, tile)
  spec_tile <- terra::crop(spec_de, tile)
  mask_tile <- terra::crop(mask_de, tile)
  
  masked_rast <- heights_tile[[1]] * mask_tile[[1]]
  masked_rast[masked_rast[[1]] < 0] <- 0
  names(masked_rast) <- c("height")
  
  # calculate wood density 
  spec_tile$dens <- (spec_tile[[1]] * densities[1, "dens_new"]+
                       spec_tile[[2]] * densities[2, "dens_new"]+
                       spec_tile[[3]] * densities[3, "dens_new"]+
                       spec_tile[[4]] * densities[4, "dens_new"]+
                       spec_tile[[5]] * densities[5, "dens_new"]+
                       spec_tile[[6]] * densities[6, "dens_new"]+
                       spec_tile[[7]] * densities[7, "dens_new"]+
                       spec_tile[[8]] * densities[8, "dens_new"]+
                       spec_tile[[9]] * densities[9, "dens_new"]+
                       spec_tile[[10]] * densities[10, "dens_new"]+
                       spec_tile[[11]] * densities[11, "dens_new"]+
                       spec_tile[[12]] * default_dens +
                       spec_tile[[13]] * default_dens +
                       spec_tile[[14]] * default_dens 
                     ) / (
                       spec_tile[[1]] + spec_tile[[2]] + spec_tile[[3]] +
                         spec_tile[[4]] + spec_tile[[5]] + spec_tile[[6]] + 
                         spec_tile[[7]] + spec_tile[[8]] + spec_tile[[9]] + 
                         spec_tile[[10]] + spec_tile[[12]] + spec_tile[[13]] + 
                         spec_tile[[13]] + spec_tile[[14]]) 
  spec_tile[is.na(spec_tile$dens)] <- default_dens
  
  # calculate modeled variables
  masked_rast$vol <- 2.3564686 * masked_rast[[1]] ^ 1.4394929 * spec_tile$dens ^ -0.7341489
  masked_rast$biom <- 3.7250261 * masked_rast[[1]] ^ 1.3447925 * spec_tile$dens ^ 0.2913108
  masked_rast$ba <- 2.1813681 * masked_rast[[1]] ^ 0.6868182 * spec_tile$dens ^ -0.7883417
  
  # export layer by layer
  for(layer in names(masked_rast)){
    outfilename <- paste0(OUTPUT_DIR,  "ForStr_v04_", layer, "_", tile$Tile_ID, 
                          ".tif")
    data.format <- ifelse(layer %in% c("vol", "biom"), "INT2U", "INT1U")
    mylayer <- round(masked_rast[layer])
    writeRaster(mylayer, outfilename, datatype = data.format)
  }
}


