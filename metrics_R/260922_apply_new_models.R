################################################################################
# This script calculates modeled forest structure variables for EarlyAccess V3
####
# Since the previous early access products were found to have a bias, 
# it was necessary to calculate the modeled variables again.
# Now, models were used which were optimized on data 
# from Thüringen, Rheinland-Pfalz and Sachsen together.
# Wood density literature values were slightly altered compared to earlyer versions. 
# Input data was generated with forestMetrics_subtile.R
# File paths are all for Timons local machine.


require(terra)

th_data <- list.files("C:/Users/miesner/Data/FP-Earlyaccess/Thuringia_forestMetrics_eA/V2/", 
                      full.names =T, pattern = ".tif")
rp_data <- list.files("C:/Users/miesner/Data/FP-Earlyaccess/Rheinland_metrics/", 
                      full.names =T, pattern = ".tif")
sn_data <- list.files("C:/USers/miesner/Data/FP-Earlyaccess/sn_ht/", 
                      full.names=  T, pattern = ".tif")
rp_data <- rp_data[substr(rp_data, nchar(rp_data)-7, nchar(rp_data)) != ".aux.xml"]

spec <- rast("C:/Users/miesner/Data/FP-Earlyaccess/tree_species_fractions.cog.tif")
mask <- rast("C:/Users/miesner/Data/FP-Earlyaccess/DE_mask.cog.tif")


## For Thürungen and RlP ##
for(filename in c(th_data, rp_data))try({

  metr <- rast(filename)
  
  metr$ccover <- metr$over02*100
  
  metr.base <- subset(metr, c("topheight", "meanH", "maxH", "ccover"))
  
  # calculate wood density
  spec.crop <- terra::crop(spec, ext(metr.base))
  names(spec.crop) <- c("Fi", "Ki", "Ta", "Dg", "La", "Bu", "Ei", "Ah", "Bi", "Er", "Pp", "AL", "nn", "xx")
  spec.crop$dens <-
    (spec.crop$Fi * 0.38 + spec.crop$Ki  * 0.44 + spec.crop$Ta  * 0.42 + 
    spec.crop$Dg  * 0.44 + spec.crop$La  * 0.49 + spec.crop$Bu  * 0.60 +
    spec.crop$Ei  * 0.60 + spec.crop$Ah  * 0.51 + spec.crop$Bi  * 0.55 + 
    spec.crop$Er  * 0.48 + spec.crop$Pp  * 0.46 + spec.crop$AL  * 0.48  +
    spec.crop$nn  * 0.48 + spec.crop$xx  * 0.48) /
    (spec.crop$Fi + spec.crop$Ki + spec.crop$Ta + spec.crop$Dg + spec.crop$La + 
     spec.crop$Bu + spec.crop$Ei + spec.crop$Ah + spec.crop$Bi + spec.crop$Er + 
     spec.crop$Pp + spec.crop$AL + spec.crop$nn + spec.crop$xx)
  # replace NA with 48
  spec.crop$dens[is.na(spec.crop$dens)] <- 0.48
  metr.new <- c(metr.base, spec.crop$dens)
  
  # calculate target variables
  
  metr.new$vol <- 4.0543948 * metr.new$meanH ^ 1.3322654 * metr.new$dens ^ -0.8135485
  metr.new$biom <- 6.1049792354 * metr.new$meanH ^ 1.2452957 * metr.new$dens ^ 0.1845087
  metr.new$ba <- 2.3898955 * metr.new$meanH ^  0.7129061 * metr.new$dens ^ -0.7446183
  
  # mask with forest mask
  mask.crop <- terra::crop(mask, ext(metr.new))
  metr.new[ mask.crop==0 ] <- 0
  
  # export layer by layer
  forcetile.suffix <- substr(filename, nchar(filename)- 14, nchar(filename))
  print(forcetile.suffix)
  
  for(layer in names(metr.new)){
    outfilename <- paste0("C:/Users/miesner/Data/FP-Earlyaccess/v3/ForStr_v03_", 
                       layer, "_", forcetile.suffix)
    data.format <- ifelse(layer %in% c("vol", "biom"), "INT2U", "INT1U")
    mylayer <- round(metr.new[layer])
    writeRaster(mylayer, outfilename, datatype = data.format)
  }
})

## For SN ##

for(filename in sn_data)try({
  forcetile.suffix <- substr(filename, nchar(filename)- 14, nchar(filename))
  print(forcetile.suffix)
  
  if(file.exists(paste0("C:/Users/miesner/Data/FP-Earlyaccess/v3/ForStr_v03_ba_",
                        forcetile.suffix))){
    next
  }
  
  
  metr <- rast(filename)
  
  metr$ccover <- metr$cc.chm.p95*100
  
  metr.base <- subset(metr, c("mean.chm.p95", "ccover"))
  names(metr.base)[1] <- "meanH"
  
  # calculate wood density
  spec.crop <- terra::crop(spec, ext(metr.base))
  names(spec.crop) <- c("Fi", "Ki", "Ta", "Dg", "La", "Bu", "Ei", "Ah", "Bi", "Er", "Pp", "AL", "nn", "xx")
  spec.crop$dens <-
    (spec.crop$Fi * 0.38 + spec.crop$Ki  * 0.44 + spec.crop$Ta  * 0.42 + 
       spec.crop$Dg  * 0.44 + spec.crop$La  * 0.49 + spec.crop$Bu  * 0.60 +
       spec.crop$Ei  * 0.60 + spec.crop$Ah  * 0.51 + spec.crop$Bi  * 0.55 + 
       spec.crop$Er  * 0.48 + spec.crop$Pp  * 0.46 + spec.crop$AL  * 0.48  +
       spec.crop$nn  * 0.48 + spec.crop$xx  * 0.48) /
    (spec.crop$Fi + spec.crop$Ki + spec.crop$Ta + spec.crop$Dg + spec.crop$La + 
       spec.crop$Bu + spec.crop$Ei + spec.crop$Ah + spec.crop$Bi + spec.crop$Er + 
       spec.crop$Pp + spec.crop$AL + spec.crop$nn + spec.crop$xx)
  
  metr.new <- c(metr.base, spec.crop$dens)
  
  # calculate target variables
  
  metr.new$vol <- 4.0543948 * metr.new$meanH ^ 1.3322654 * metr.new$dens ^ -0.8135485
  metr.new$biom <- 6.1049792354 * metr.new$meanH ^ 1.2452957 * metr.new$dens ^ 0.1845087
  metr.new$ba <- 2.3898955 * metr.new$meanH ^  0.7129061 * metr.new$dens ^ -0.7446183
  
  
  # export layer by layer
  for(layer in names(metr.new)){
    outfilename <- paste0("C:/Users/miesner/Data/FP-Earlyaccess/v3/ForStr_v03_", 
                          layer, "_", forcetile.suffix)
    data.format <- ifelse(layer %in% c("vol", "biom"), "INT2U", "INT1U")
    mylayer <- round(metr.new[layer])
    # if a file of the name already exists (because it is in the overlap between SN and TH), merge the two
    if(file.exists(outfilename)){
      existing.rast <- rast(outfilename)
      merge(mylayer, existing.rast, filename = outfilename, overwrite = TRUE)
    }else{
      writeRaster(mylayer, outfilename, datatype = data.format)
    }
  }
})


## process vertical-layering layer ##

require(sf)
force <- read_sf("C:/USers/miesner/Data/FP-Earlyaccess/datacube_template_rasters/data/datacube-grid_DEU_10km.gpkg")
lyr.th <- rast("C:/Users/miesner/Data/FP-Earlyaccess/Layer_th_cog.tif")
lyr.rp <- rast("C:/USers/miesner/Data/FP-Earlyaccess/Layer_Date_rp.cog.tif")
mask.th <- crop(mask, lyr.th)
lyr.th[mask.th == 0 & !is.na(lyr.th$Layer)] <- 0
mask.rp <- crop(mask, lyr.rp)
lyr.rp[mask.rp == 0 & !is.na(lyr.rp$Layer)] <- 0

for(i in 1:nrow(force)){
  tile <- force[i,]
  outfilename <- paste0("C:/Users/miesner/Data/FP-Earlyaccess/v3/ForStr_v03_",
                        "vertclyr_x00", tile$Tile_X, "_y00", tile$Tile_Y, ".tif")
  try({
    small.lyr <- crop(lyr.th, tile)$Layer
    if(length(unique(small.lyr))>0){ # if the layer contains meaningful content
      writeRaster(small.lyr, outfilename, datatype = "INT1U")
    }
  })
  try({
    small.lyr <- crop(lyr.rp, tile)$Layer
    if(length(unique(small.lyr))>0){ # if the layer contains meaningful content
      writeRaster(small.lyr, outfilename, datatype = "INT1U")
    }
  })
}
