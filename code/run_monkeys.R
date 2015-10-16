# run monkey niche models

# set the RNG seed
set.seed(1)

# whether to run with GPs (rather than bootstrapped BRTs)
useGPs <- FALSE

# set output path
outpath <- 'output/monkeys/'

# load the data

# raster covariates
covs_current <- brick('data/clean/raster/raster.grd')
covs_temporal <- brick('data/clean/raster/raster_temporal.grd')
covs_nontemporal <- brick('data/clean/raster/raster_nontemporal.grd')

# rename the temporal layers (remove year)
names(covs_current)[names(covs_current)=='forest_intact_2012']<-'forest_intact'
names(covs_current)[names(covs_current)=='forest_disturbed_2012']<-'forest_disturbed'
names(covs_current)[names(covs_current)=='open_shrublands_2012']<-'open_shrublands'
names(covs_current)[names(covs_current)=='woody_savannas_2012']<-'woody_savannas'
names(covs_current)[names(covs_current)=='savannas_2012']<-'savannas'
names(covs_current)[names(covs_current)=='grasslands_2012']<-'grasslands'
names(covs_current)[names(covs_current)=='permanent_wetlands_2012']<-'permanent_wetlands'
names(covs_current)[names(covs_current)=='croplands_2012']<-'croplands'
names(covs_current)[names(covs_current)=='cropland_natural_vegetation_mosaic_2012']<-'cropland_natural_vegetation_mosaic'

# load background dataset
mammal_bias <- read.csv('data/clean/occurrence/mammals-bias.csv')

# change column names to match monkey species data sets
colnames(mammal_bias)[colnames(mammal_bias)=='decimalLatitude']<-'Latitude'
colnames(mammal_bias)[colnames(mammal_bias)=='decimalLongitude']<-'Longitude'
colnames(mammal_bias)[colnames(mammal_bias)=='year']<-'Year'

# add columns for later to deal with polygon data in leonina occurrence dataset
mammal_bias$Geometry_type <- 'point'
mammal_bias$Unique_ID <- NA

#colnames(mammal_bias)[colnames(mammal_bias)=='gbifID']<-'Survey_id'

# occurrence data
fascicularis <- read.csv('data/clean/occurrence/model fascicularis data JAN2015.csv')
nemestrina <- read.csv('data/clean/occurrence/model nemestrina data JAN2015.csv')
leonina <- read.csv('data/clean/occurrence/polygon_leonina_data.csv')

# remove extra columns in fascicularis dataset 
fascicularis <- subset(fascicularis, select= -c(9:19))

# add columns for later to deal with polygon data in leonina occurrence dataset
fascicularis$Unique_ID <- NA 
nemestrina$Unique_ID <- NA 

# loop through these species fitting models
for (species in c('fascicularis', 'nemestrina', 'leonina')) {
  
  # get occcurrence data for this species
  occ <- get(species)
  
  # only use presence points 
  occ <- subset(occ, occ$Presence==1)
  
  # set output path for this species
  outpath_species <- paste0(outpath,
                            species,
                            '/')
  
  # create the directory if it doesn't already exist
  if (!file.exists(outpath_species)) {
    dir.create(outpath_species)
  }
  
  # combine the occurrence and background records
  dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                     occ[, c('Unique_ID', 'Longitude','Latitude', 'Year', 'Geometry_type')]),
               cbind(PA = rep(0, nrow(mammal_bias)),
                     mammal_bias[ ,c('Unique_ID', 'Longitude', 'Latitude', 'Year','Geometry_type')]))
  
  
  # create an empty dataframe
  dat_covs <- extract(covs_current, dat[, c('Longitude', 'Latitude')])
  
  dat_covs[] <- NA
  
  
  # loop through extracting covs for each year
  for (year in 2001:2012) {
    
    # truncate years to the ones we have covariates for
    years <- dat$Year
    years <- pmax(years, 2001)
    years <- pmin(years, 2012)
    
    # index for this year's data
    idx <- which(years == year)
    
    # covs for this year
    
    # get index for temporal covariates for relevant year  
    
    idx2 <- which (substr(names(covs_temporal),
                          nchar(names(covs_temporal)) - 3,
                          nchar(names(covs_temporal)))==year)
    
    covs_year <- covs_temporal[[idx2]]
    names(covs_year)[names(covs_year)==paste0('forest_intact_', year)]<-'forest_intact'
    names(covs_year)[names(covs_year)==paste0('forest_disturbed_', year)]<-'forest_disturbed'
    names(covs_year)[names(covs_year)==paste0('open_shrublands_', year)]<-'open_shrublands'
    names(covs_year)[names(covs_year)==paste0('woody_savannas_', year)]<-'woody_savannas'
    names(covs_year)[names(covs_year)==paste0('savannas_', year)]<-'savannas'
    names(covs_year)[names(covs_year)==paste0('grasslands_', year)]<-'grasslands'
    names(covs_year)[names(covs_year)==paste0('permanent_wetlands_', year)]<-'permanent_wetlands'
    names(covs_year)[names(covs_year)==paste0('croplands_', year)]<-'croplands'
    names(covs_year)[names(covs_year)==paste0('cropland_natural_vegetation_mosaic_', year)]<-'cropland_natural_vegetation_mosaic'
    
    # add nontemporal covariates
    covs_year <- addLayer(covs_year, covs_nontemporal)
    
    # extract data
    covs_year_extract <- extract(covs_year, dat[idx, c('Longitude', 'Latitude')])
    
    # check they're all there
    stopifnot(all(colnames(dat_covs) %in% colnames(covs_year_extract)))
    
    # match up the column names so they're in the right order
    match <- match(colnames(dat_covs), colnames(covs_year_extract))
    
    # extract covariates for all points
    dat_covs[idx, ] <- covs_year_extract[, match]
    
  } 
  
  # combine covariates with the other info
  dat_all <- cbind(dat, dat_covs)
  
  # remove any records with missing covariates (outside mask)
  
  # which rows in dat_covs contain an NA
  outside_idx <- attr(na.omit(dat_covs), 'na.action')
  
  stopifnot(is.null(outside_idx))
  
  if (!is.null(outside_idx)) {
    # remove these from dat_all
    dat_all <- dat_all[-outside_idx, ]
  }
  
  
  # if we're running Gaussian process models
  if (useGPs) {
    
    # get the point data
    # a single row for each polygon will be added later
    dat_final <- subset(dat_all, dat_all$Geometry_type!='polygon')
    
    
    # need to calculate the error matrix for covariates
    # (matrix of standard deviations for covariates for each datapoint)
    
    # create matrix for standard deviations of covariate values
    # set standard deviation to 0 for point data
    # a single row for each polygon will be added later
    e <- matrix(data=0, nrow=nrow(dat_final), ncol=ncol(dat_covs))
    
    # need to account for polygons (not replicate each as 100s of points)
    
    # get the polygon data
    dat_poly <- subset(dat_all, dat_all$Geometry_type=='polygon')
    
    # get the different polygon codes
    u <- unique(dat_poly$Polygon_code)
    
    # remove the NA values
    u <- u[!is.na(u)]
      
    # loop through obtaining the column means for each polygon and add a single row to dat_final for each
    for (poly in u){
      
      #poly <- u[21] 
      # get the indexes of all polygons with a particular poly code
      idx <- which(dat_poly$Polygon_code==poly)
      # subset the data into covariate columns only
      sub <- dat_poly[idx, colnames(dat_covs)] 
      # calculate the the column means
      means <- colMeans(sub) 
      # create a data frame of single row using the first index
      newrow <- dat_poly[idx[1],] 
      
      # loop through the covariate values replacing with means 
      for (cov in colnames(dat_covs)){
        newrow[,cov]<-means[cov]
      }
      
      # add new row to the point data
      dat_final <- rbind(dat_final, newrow)
      
      if (length(idx) > 1){
        # convert subsetted data frame to a matrix
        mat <- as.matrix(sub)
        # get the standard deviation for each column (covariate)
        stdev <- apply (mat, 2, FUN=sd) 
        
      } else {
        stdev <- rep(0, ncol(dat_covs))
      }
              
      # add standard deviations to the e matrix
      e <- rbind(e, stdev)           
      
    }
    
    # get the different Gaul codes
    u2 <- unique(dat_poly$Gaul_code)
    
    # remove the NA values
    u2 <- u2[!is.na(u2)]
    
    # loop through obtaining the column means for each polygon
    for (poly in u2){
      
      # get the indexes of all polygons with a particular poly code
      idx <- which(dat_poly$Gaul_code==poly)
      # subset the data into covariate columns only
      sub <- dat_poly[idx, colnames(dat_covs)] 
      # calculate the column means
      means <- colMeans(sub)
      # create a data frame of single row using the first index    
      newrow <- dat_poly[idx[1],]   
      
      # loop through the covariate values replacing with means 
      for (cov in colnames(dat_covs)){
        newrow[,cov]<-means[cov]
      }
      
      # add new row to the point data
      dat_final <- rbind(dat_final, newrow)
      
      if (length(idx) > 1){
        # convert subsetted data frame to a matrix
        mat <- as.matrix(sub)
        # get the standard deviation for each column (covariate)
        stdev <- apply (mat, 2, FUN=sd) 
        
      } else {
        stdev <- rep(0, ncol(dat_covs))
      }
      
      # add standard deviations to the e matrix
      e <- rbind(e, stdev)     
      
    }
    
    # make a vector denoting presence and absence records
    y <- dat_final$PA
    
    # make a dataframe of covariates
    x <- data.frame(dat_final[, colnames(dat_covs)])
    
    # add weights as in the BRT models
    wts <- ifelse(dat_final$PA == 1, 1, sum(dat_final$PA) / sum(1 - dat_final$PA))
    
    # fit the model
    m <- graf(y = y,
              x = x,
              error = e,
              weights = wts,
              opt.l = TRUE)  # whether to optimise the lengthscales
    
    # make predictions
    preds <- predict(m, covs_current)
    
    # save prediction raster
    writeRaster(preds,
                file = paste0(outpath_species,
                              species,
                              '_GP'),
                format = 'GTiff',
                overwrite = TRUE)
    
    # plot the rasters
    png(paste0(outpath_species,
               'prediction_GP.png'),
        width = 2000,
        height = 2000,
        pointsize = 30)
    
    par(oma = rep(0, 4),
        mar = c(0, 0, 0, 2))
    
    plot(preds[[1]],
         axes = FALSE,
         box = FALSE)
    
    points(dat[dat$PA == 1, 3:4],
           pch = 16,
           cex = 1,
           col = 'blue')
    
    dev.off()
    
    # plot the effect curves
    png(paste0(outpath_species,
               'effects_GP.png'),
        width = 2000,
        height = 2500,
        pointsize = 30)
    
    par(mfrow = n2mfrow(ncol(dat_covs)),
        mar = c(4, 4, 0.5, 1))
    plot(m, peak = TRUE)
    
    dev.off()
    
    # print the lengthscales
    lengthscales <- m$ls
    names(lengthscales) <- colnames(dat_covs)
    write.csv(lengthscales,
              file = paste0(outpath_species,
                            'GP_lengthscales.csv'))
    
  } else {
    # if we're running bootstrapped BRT models
    
    ncpu <- 60
    nboot <- ncpu*2
    
    # create list of random permutations of dat_all
    # for leonina dataset, sampling one cell from each polygon for each iteration,
    # then bootstrapping as usual (minimum 5 pres/5 abs)
    data_list <- replicate(nboot,
                           subsamplePolys(dat_all,
                                          minimum = c(5, 5),
                                          replace = TRUE),
                           simplify = FALSE)
    
    # initialize the cluster
    sfInit(parallel = TRUE, cpus = ncpu)
    sfLibrary(seegSDM)
    
   # gbm.x = 7:ncol(data_list[[1]]) - specifies number of columns in which to look for covariates, based on the first row in data_list
    
    model_list <- sfLapply(data_list,
                           runBRT,
                           gbm.x = 7:ncol(data_list[[1]]),
                           gbm.y = 1,
                           n.folds = 5,
                           pred.raster = covs_current,
                           gbm.coords = 3:4,
                           wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)))
    
    # get cv statistics in parallel
    stat_lis <- sfLapply(model_list, getStats)
    
    # summarise all the ensembles
    preds <- stack(lapply(model_list, '[[', 4))
    # same as stack(lapply(model_list, function(x) x[[4]])) 
    
    # summarise the predictions in parallel
    preds_sry <- combinePreds(preds, parallel = TRUE)
    
    # stop the cluster
    sfStop()
    
    # convert the stats list into a matrix using the do.call function
    stats <- do.call("rbind", stat_lis)
    
    # save them
    write.csv(stats,
              paste0(outpath_species, 
                     'stats.csv'))
    
    names(preds_sry) <- c('mean',
                          'median',
                          'lowerCI',
                          'upperCI')
    
    # save the prediction summary
    writeRaster(preds_sry,
                file = paste0(outpath_species,
                              species),
                format = 'GTiff',
                overwrite = TRUE)
    
    # save the relative influence scores
    relinf <- getRelInf(model_list)
    
    write.csv(relinf,
              file = paste0(outpath_species,
                            'relative_influence.csv'))
    
    # plot and the marginal effect curves
    png(paste0(outpath_species,
               'effects.png'),
        width = 2000,
        height = 2500,
        pointsize = 30)
    
    par(mfrow = n2mfrow(nlayers(covs_current)))
    
    effects <- getEffectPlots(model_list,
                              plot = TRUE)
    
    dev.off()
    
    # plot the risk map
    png(paste0(outpath_species,
               'prediction_mean.png'),
        width = 2000,
        height = 2000,
        pointsize = 30)
    
    par(oma = rep(0, 4),
        mar = c(0, 0, 0, 2))
    
    plot(preds_sry[[1]],
         axes = FALSE,
         box = FALSE)
    
    points(dat[dat$PA == 1, 3:4],
           pch = 16,
           cex = 1,
           col = 'blue')
    
    dev.off()
    
  }
  
}
