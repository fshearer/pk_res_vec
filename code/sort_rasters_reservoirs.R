# organise covariates for analysis

# list the covariates of interest
covars <- c('data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2001.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2002.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2003.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2004.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2005.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2006.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2007.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2008.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2009.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2010.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2011.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2012.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2001.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2002.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2003.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2004.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2005.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2006.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2007.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2008.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2009.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2010.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2011.fractional.5km.tif',
            'data/raw/raster/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2012.fractional.5km.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class07_Open_Shrublands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class08_Woody_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class09_Savannas.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class10_Grasslands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class11_Permanent_Wetlands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class12_Croplands.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2001.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2002.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2003.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2004.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2005.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2006.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2007.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2008.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2009.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2010.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2011.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/landcover 2001-12/2012.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            'data/raw/raster/covariates/MODIS satellite data/tasseled cap brightness/TCB_Mean_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/tasseled cap brightness/TCB_SD_5km_Mean.tif',
            'data/raw/raster/covariates/SEAsia_pop5k/SEAsia_pop5k.tif',
            'data/raw/raster/covariates/MODIS satellite data/enhanced vegetation index/EVI_Mean_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/enhanced vegetation index/EVI_SD_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/land surface temp/LST_Day_Mean_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/land surface temp/LST_Day_SD_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/tasselled cap wetness/TCW_Mean_5km_Mean.tif',
            'data/raw/raster/covariates/MODIS satellite data/tasselled cap wetness/TCW_SD_5km_Mean.tif',
            'data/raw/raster/covariates/SRTM elevation/mod_dem_5k_CLEAN.flt')


# give them readable names
names <- c('forest_intact_2001',
           'forest_intact_2002',
           'forest_intact_2003',
           'forest_intact_2004',
           'forest_intact_2005',
           'forest_intact_2006',
           'forest_intact_2007',
           'forest_intact_2008',
           'forest_intact_2009',
           'forest_intact_2010',
           'forest_intact_2011',
           'forest_intact_2012',
           'forest_disturbed_2001',
           'forest_disturbed_2002',
           'forest_disturbed_2003',
           'forest_disturbed_2004',
           'forest_disturbed_2005',
           'forest_disturbed_2006',
           'forest_disturbed_2007',
           'forest_disturbed_2008',
           'forest_disturbed_2009',
           'forest_disturbed_2010',
           'forest_disturbed_2011',
           'forest_disturbed_2012',
           'open_shrublands_2001',
           'open_shrublands_2002',
           'open_shrublands_2003',
           'open_shrublands_2004',
           'open_shrublands_2005',
           'open_shrublands_2006',
           'open_shrublands_2007',
           'open_shrublands_2008',
           'open_shrublands_2009',
           'open_shrublands_2010',
           'open_shrublands_2011',
           'open_shrublands_2012',
           'woody_savannas_2001',
           'woody_savannas_2002',
           'woody_savannas_2003',
           'woody_savannas_2004',
           'woody_savannas_2005',
           'woody_savannas_2006',
           'woody_savannas_2007',
           'woody_savannas_2008',
           'woody_savannas_2009',
           'woody_savannas_2010',
           'woody_savannas_2011',
           'woody_savannas_2012',
           'savannas_2001',
           'savannas_2002',
           'savannas_2003',
           'savannas_2004',
           'savannas_2005',
           'savannas_2006',
           'savannas_2007',
           'savannas_2008',
           'savannas_2009',
           'savannas_2010',
           'savannas_2011',
           'savannas_2012',
           'grasslands_2001',
           'grasslands_2002',
           'grasslands_2003',
           'grasslands_2004',
           'grasslands_2005',
           'grasslands_2006',
           'grasslands_2007',
           'grasslands_2008',
           'grasslands_2009',
           'grasslands_2010',
           'grasslands_2011',
           'grasslands_2012',
           'permanent_wetlands_2001',
           'permanent_wetlands_2002',
           'permanent_wetlands_2003',
           'permanent_wetlands_2004',
           'permanent_wetlands_2005',
           'permanent_wetlands_2006',
           'permanent_wetlands_2007',
           'permanent_wetlands_2008',
           'permanent_wetlands_2009',
           'permanent_wetlands_2010',
           'permanent_wetlands_2011',
           'permanent_wetlands_2012',
           'croplands_2001',
           'croplands_2002',
           'croplands_2003',
           'croplands_2004',
           'croplands_2005',
           'croplands_2006',
           'croplands_2007',
           'croplands_2008',
           'croplands_2009',
           'croplands_2010',
           'croplands_2011',
           'croplands_2012',
           'cropland_natural_vegetation_mosaic_2001',
           'cropland_natural_vegetation_mosaic_2002',
           'cropland_natural_vegetation_mosaic_2003',
           'cropland_natural_vegetation_mosaic_2004',
           'cropland_natural_vegetation_mosaic_2005',
           'cropland_natural_vegetation_mosaic_2006',
           'cropland_natural_vegetation_mosaic_2007',
           'cropland_natural_vegetation_mosaic_2008',
           'cropland_natural_vegetation_mosaic_2009',
           'cropland_natural_vegetation_mosaic_2010',
           'cropland_natural_vegetation_mosaic_2011',
           'cropland_natural_vegetation_mosaic_2012',
           'TCB_mean', 
           'TCB_SD',
           'SEAsia_pop', 
           'EVI_mean',
           'EVI_SD',
           'LSTday_mean',
           'LSTday_SD',
           'TCW_mean',
           'TCW_SD',
           'SRTM_elevation')

# loop through opening links to the rasters
covs <- lapply(covars, raster)

# get all extents
extents <- t(sapply(covs, function (x) as.vector(extent(x))))

# get the smallest extent squares of all layers
ext <- extent(c(max(extents[, 1]),
                min(extents[, 2]),
                max(extents[, 3]),
                min(extents[, 4])))

# crop all layers by this
covs <- lapply(covs, crop, ext)

# stack the rasters
covs <- stack(covs) 

# loop through covariate layers, find all NAs, find cov with most NAs, and mask the whole stack by this layer
#max_nas <-0
#nas_layer <-0

#for (i in 1:nlayers(covs)){
#  nas <- length(which(is.na(getValues(covs[[i]]))))
#  if (nas > max_nas) {
#    max_nas <- nas
#    nas_layer <- i
#  }
#}

#covs <- mask(covs, covs[[nas_layer]])

# loop through covs masking one layer by every other layer to create a master mask
master_mask <- covs[[1]]

for (i in 1:nlayers(covs)){
  
  master_mask <- mask(master_mask, covs[[i]])
  
}

# mask all covariates by master mask
covs <- mask(covs, master_mask)

cbind(names(covs), names)

warning('check the names match up!')

# give them nicer names
names(covs) <- names

# create subsets of temporal, non-temporal and the most contemporary covariates for running monkey and vector models 
covs_current <- subset(covs, c(12, 24, 36, 48, 60, 72, 84, 96, 108:118))
covs_temporal <- subset(covs, 1:108)
covs_nontemporal <- subset(covs, 109:118)

# output them as multiband .grd files
writeRaster(covs_current,
            file = 'data/clean/raster/raster',
            overwrite = TRUE)

writeRaster(covs_temporal,
            file = 'data/clean/raster/raster_temporal',
            overwrite = TRUE)

writeRaster(covs_nontemporal,
            file = 'data/clean/raster/raster_nontemporal',
            overwrite = TRUE)

