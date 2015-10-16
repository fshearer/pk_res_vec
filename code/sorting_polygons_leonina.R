# prepare leonina occurrence data for modelling

# load data

# raw occurrence data
dat <- read.csv('data/raw/occurrence/model leonina data JAN2015.csv',
                      stringsAsFactors = FALSE)

# use occurrence data only
dat <- subset(dat, dat$Presence==1)

# create data subset without polygons
dat_new <- subset(dat, dat$Geometry_type!='polygon')

# subset the data into polygons
poly1_dat <- subset(dat, dat$Admin_level=='leonina_ras1')

poly2_dat <- subset(dat, dat$Admin_level=='leonina_ras2')

poly3_dat <- subset(dat, dat$Admin_level=='admin2')

# load matching rasters for polygons
poly1 <- brick('data/raw/raster/polygon_rasters/leonina_ras1.tif')

poly2 <- brick('data/raw/raster/polygon_rasters/leonina_ras2.tif')

admin2 <- brick('data/raw/raster/polygon_rasters/Admin2_Raster.flt')

# load a 5km template raster of the study area extent
template <- raster('data/clean/raster/raster.grd')

# crop admin2 raster to this extent
poly3 <- crop(admin2, template)

# ~~~~~~~~~~~
# sample random points from within polygons on a 5km grid
for (polygon in c('poly1_dat', 'poly2_dat', 'poly3_dat')){
  
  # get data subset for polygon
  poly_dat <- get(polygon)
  
  for (i in 1:nrow(poly_dat)){
    
    # get polygon code from dataset
    
    if (polygon %in% c('poly1_dat', 'poly2_dat')) {
      polycode <- poly_dat$Polygon_code[i]
      
    } else {
      polycode <- poly_dat$Gaul_code[i]
    }  
    
    # get raster for this polygon
    substring <- substr(polygon, 0, 5)
    raster <- get(substring)
    
    # obtain index for pixels in raster with this identifier 
    idx <- which(getValues(raster)==polycode)
    
    # get coordinates of the cells there
    pts <- xyFromCell(raster, idx)
    
    # the number of points
    n <- nrow(pts)
    
    # if its larger than 100, pick 100 of them at random
    if (n>100){
      pts<- pts[sample(1:100, replace=FALSE),]
      n <- 100
    }
    
    # pull out the info for that record
    info_mat <- poly_dat[rep(i, nrow(pts)),]
    
    # stick the coordinates in
    info_mat[,c('Latitude', 'Longitude')] <- pts[,2:1]
    
    # append it to new_dat
    dat_new <- rbind(dat_new, info_mat)
}
}

# remove the polygon columns
dat_new <- subset(dat_new, select= -c(Admin_level, Polygon_code, Gaul_code))


# if the resulting dataset looks fairly clean, write it to disk
if (!(any(is.na(dat_new)))) {
  # output the resulting table
  write.csv(dat_new,
          file = 'data/clean/occurrence/polygon_leonina_data.csv',
          row.names = FALSE)
}  



