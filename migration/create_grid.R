#!/bin/env Rscript 

.libPaths(c("/public/groups/meyerlab/esteiger/R_libs", .libPaths()))

library(dggridR)
library(sf)
library(terra)

# this doesn't appear to be necessary... uses all the mem anyhow
# memory.limit(size=1000000)   
# print("increased the memory limit to 1000gb)


# ---- Create Global Grid at Resolution 14 ----
dgg_res14 <- dgconstruct(res = 14, projection = "ISEA", aperture = 4, topology = "TRIANGLE")
print("Constructed a discrete global grid system (dggs) object")

grid_res14 <- dgearthgrid(dgg_res14)
print("Returned the coordinates constituting the boundary of cells for the entire Earth")
rm(dgg_res14)

st_write(grid_res14, "res14_dgg.shp", delete_layer = TRUE)
print("saved to external shape file")
rm(grid_res14)

# ---- Crop to Study Area ----
#bbox_extent_14 <- ext(-120.7, -120.1, 34.4, 34.9)
#grid_res14_cropped <- st_crop(grid_res14, bbox_extent_14)
#rm(grid_res14)

# Export cropped grid
#st_write(grid_res14_cropped, "study_area_res14_dgg.shp", delete_layer = TRUE)
#rm(grid_res14_cropped)

# ---- Overwrite with Resolution 17 Grid ----
dgg_res17 <- dgconstruct(res = 17, projection = "ISEA", aperture = 4, topology = "TRIANGLE")
print("Constructed a discrete global grid system (dggs) object")

grid_res17 <- dgearthgrid(dgg_res17)
print("Returned the coordinates constituting the boundary of cells for the entire Earth")
rm(dgg_res17)

st_write(grid_res17, "res17_dgg.shp", delete_layer = TRUE) 
print("saved to external shape file")

# Crop to second bounding box
#bbox_extent_17 <- ext(-120.55, -120.50, 34.56, 34.605)
#grid_res17_cropped <- st_crop(grid_res17, bbox_extent_17)
#rm(grid_res17)

# Export cropped grid
#st_write(grid_res17_cropped, "strauss_area_res17_dgg.shp", delete_layer = TRUE) 
#rm(grid_res17_cropped)

###################
# Attempt where grid created and saved at same time
# Don't like this bc I can't tell whether failing at grid production or export
# print("starting grid construction")
# dgearthgrid(dgconstruct(res=14, projection = "ISEA", aperture = 4, topology = "TRIANGLE"),savegrid = "world_triangle_res14.shp")
# print("Saved grid resolution 14 to shapefile... ")

# print("starting grid construction")
# dgearthgrid(dgconstruct(res=17, projection = "ISEA", aperture = 4, topology = "TRIANGLE"),savegrid = "world_triangle_res17.shp")
# print("Saved grid resolution 17 to shapefile... ")

