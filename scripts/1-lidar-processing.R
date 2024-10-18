################################################################################
# Code to calculate stand structure metrics from LiDAR data
  # Canopy height
  # Canopy cover %
  # Understory density (NRD)

# Amanda Zak
# July 2024
################################################################################
library(lidR)
library(terra)

res_chm <- 16.4 # 5 meter resolution for canopy height model
res_ud <- 196.85 # 60 meter for understory density

# Each study area will be processed separately for more manageable processing times
# Resulting .tif rasters from this code will be mosaiced in arcpy 

# List of study areas
studyareas <- c("14_51","33","46","44","45","39","49","50","52")

# Set study area number for following code
sa <- studyareas[2]

###############################################################################

# 1. Load catalog

saCat <- eval(parse(text = paste0("readLAScatalog('LiDAR tifs/area_", sa, "/tiles/',
                                  filter = '-drop_class 7 -drop_withheld')")))

# check catalog
las_check(saCat)
plot(saCat)



# 2. Normalize ground

# specify output location
opt_output_files(saCat) <- eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/normalized/{ORIGINALFILENAME}_normalized'")))

# normalize
normCat <- normalize_height(saCat,knnidw())



# 3. Filter high and low noise points by re-initializing the catalog

# using a 180-foot height cut-off because the tallest tree in PA is 180 feet tall
filtCat <- eval(parse(text = paste0("readLAScatalog('LiDAR tifs/area_", sa, "/normalized/',
                                  filter = '-drop_z_below 0 -drop_z_above 180')")))



# 4. Canopy height model

# specify output location
opt_output_files(filtCat) <- eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/chm/{ORIGINALFILENAME}_chm'")))

# calculate
chmCat <- rasterize_canopy(filtCat, res = res_chm, algorithm = dsmtin()) # 5 meter pixel size

# check for errors
plot(chmCat)
max(chmCat)

# write as raster
writeRaster(chmCat, eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/chm/",sa,"_chm.tif'"))), overwrite = TRUE)

# scale chm raster up from 5-meter pixel to 60-meter pixel

# read in the completed canopy height tif
chmTif <- eval(parse(text = paste0("rast('LiDAR tifs/area_",sa,"/chm/",sa,"_chm.tif')")))

# aggregate from pixel size of 5 meters to 60 meters (factor of 12)
chmAgg <- aggregate(chmTif,12)
plot(chmAgg)

# save as raster
writeRaster(chmAgg, eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/chm/",sa,"_chm_60.tif'"))), overwrite = TRUE)



# 5. Canopy cover

# set up reclass matrix
h <- 10 # minimum canopy height in feet
m <- c(-100000,h,0,
       h+0.00001,100000,1) # reclass matrix
mat <- matrix(m, ncol=3,byrow=T)

# reclassify pixels as 0 (below 10 feet) or 1 (above 10 feet)
chm_classCat <- classify(chmCat, mat)

# calculate canopy cover at larger pixel size
can_covCat <- aggregate(chm_classCat,fact=12)
can_covCat <- subst(can_covCat,NA,0)

# check for errors
plot(can_covCat)

# save as raster
writeRaster(can_covCat, eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/canopycover/",sa,"_cancov_60.tif'"))), overwrite = TRUE)



# 6. Understory density

# create function
calc_nrd <- function(Z,i,j) {
  length(which(Z >= i & Z < j))/length(which(Z < j))
}
i <- 0.82 # lower understory bound - 0.25 m
j <- 6.56 # upper understory bound - 2 m

# specify output location
opt_output_files(filtCat) <- eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/ud/{ORIGINALFILENAME}_ud'")))

# calculate
udCat <- pixel_metrics(filtCat, ~calc_nrd(Z,i,j), res = res_ud, overwrite = TRUE)

# check for errors
length(which(is.na(udCat[])))
plot(udCat)

# write function to fill holes by taking the mean of its neighboring cells
# (requires more than one neighboring cell)
requireNeighbors <- function(x){
  if (length(na.omit(x)) < 2) { # if fewer than 2 neighbors with values
    return(NA) # keep NA
  } else { # otheriwse, return the mean
    return(mean(na.omit(x)))
  }
}

# fill holes
udCat_filled <- focal(udCat, w = matrix(c(1,1,1,1,0,1,1,1,1),nrow=3), fun = requireNeighbors, na.policy = "only")
plot(udCat_filled)
i <- 1
while (length(which(is.na(udCat_filled[]))) > 0) { # as long as some cells are NA, continue filling holes
  udCat_filled <- focal(udCat_filled, w = matrix(c(1,1,1,1,0,1,1,1,1),nrow=3), fun = requireNeighbors, na.policy = "only")
  i <- i+1 # track how many iterations required to fill all holes
}
print(i)

# check for errors again
length(which(is.na(udCat_filled[])))
plot(udCat_filled)

# write raster
writeRaster(udCat_filled, eval(parse(text = paste0("'LiDAR tifs/area_",sa,"/ud/",sa,"_ud.tif'"))), overwrite = TRUE)

