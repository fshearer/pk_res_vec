### function for sampling data containing polygon, sampling from the unique_ID column

subsamplePolys <- function (data, ...) {
  # given a presence-background dataset, with multiple rows for some of the
  # occurrence records, subset it so that there's only one randomly selected
  # point from each polygon and then take a bootstrap it using `subsample`.
  # Dots argument is passed to subsample.
  
  # get the point data
  point_dat <- subset(data, data$Geometry_type!='polygon') 
  
  # subset to get polygon data only
  poly_dat <- subset(data, data$Geometry_type=='polygon')
  
  if (nrow(poly_dat) > 0) {
    
    # get the different IDs
    u <- unique(poly_dat$Unique_ID)
  
    # remove the NA values
    u <- u[!is.na(u)]
  
    # loop through, picking an index for each based on the number available
    data_idx <- sapply(u,
                     function (identifier, data) {
                       idx <- which(data$Unique_ID == identifier)
                       sample(idx, 1)
                     },
                     data)
  
    # get the subsetted dataset
    dat <- data[data_idx, ]
    
    # append the subsetted polygon data to the point data
    point_dat <- rbind(dat, point_dat)
  
  }

  # randomly subsample the dataset
  ans <- subsample(point_dat,
                   n = nrow(point_dat),
                   ...)
                    
  # remove the polygon columns
  #ans <- subset(ans, select= -c(Geometry_type, Polygon_code, Gaul_code))
  
  return (ans)
}




