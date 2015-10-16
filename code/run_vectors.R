# run vector niche models

# set the RNG seed
set.seed(1)

# set output path
outpath <- 'output/vectors/'

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
anopheles_bias <- read.csv('data/clean/occurrence/anopheles_bias.csv')

# change column names to match vector species data sets
colnames(anopheles_bias)[colnames(anopheles_bias)=='latitude']<-'Latitude'
colnames(anopheles_bias)[colnames(anopheles_bias)=='longitude']<-'Longitude'
colnames(anopheles_bias)[colnames(anopheles_bias)=='year']<-'Year'

#colnames(anopheles_bias)[colnames(anopheles_bias)=='id']<-'Unique_id'

# check all points fall within the extent
# get extent of covs_current
ext <- extent(covs_current)

# find index of all points falling outside the extent
outside_ext_idx <- which((anopheles_bias$Latitude < ext[3]) 
                         | (anopheles_bias$Latitude > ext[4]) 
                         | (anopheles_bias$Longitude < ext[1]) 
                         | (anopheles_bias$Longitude > ext[2]))

# remove these from dataset
anopheles_bias <- anopheles_bias[-outside_ext_idx, ]

# some anopheles_bias data points have missing values (NA) for'Year'
# make an index of data points with and without NA for year
NA_year_idx <- which(is.na(anopheles_bias$Year))
year_idx <- which(!(is.na(anopheles_bias$Year)))
x <- anopheles_bias[year_idx, 'Year'] 

# sample from distribution of years within data set
years <- sample(x, size=length(NA_year_idx), replace=TRUE)

# replace NAs with year by sampling from distribution of years within data set
for (i in 1:length(NA_year_idx)) {
  anopheles_bias[NA_year_idx[[i]], 'Year']<- years [[i]]
}

# create data subsets for data that has been identified to species or complex level
bias_species <- subset(anopheles_bias, anopheles_bias$species_identified=='TRUE')

# occurrence data
vector_group<- read.csv('data/clean/occurrence/model vector data JAN2015.csv',
                        stringsAsFactors =FALSE)

#vector_group$Unique_id <- 999

# create data subsets for the two complexes and each species
dirus_complex<-subset(vector_group, vector_group$Complex=='Dirus')
leucosphyrus_complex<-subset(vector_group, vector_group$Complex=='Leucosphyrus')

# create data subsets for each species
dirus <- subset(vector_group, vector_group$Species=='Anopheles dirus')
#cracens <- subset(vector_group, vector_group$Species=='Anopheles cracens')
#balabacensis <- subset(vector_group, vector_group$Species=='Anopheles balabacensis')
#introlatus <-subset(vector_group, vector_group$Species=='Anopheles introlatus')
#latens <-subset(vector_group, vector_group$Species=='Anopheles latens ')

# loop through these species fitting models
for (species in c('dirus', 'dirus_complex', 'leucosphyrus_complex', 'vector_group')) {
  
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
  
  
  if (species %in% c('dirus_complex', 'leucosphyrus_complex', 'vector_group')) {
    bias <- anopheles_bias
    
  } else {
    bias <- bias_species
  }
  
  # combine the occurrence and background records
  dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                     occ[, c('Longitude', 'Latitude', 'Year')]),
               cbind(PA = rep(0, nrow(bias)),
                     bias[ ,c('Longitude', 'Latitude', 'Year')]))
  
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
    
    idx2 <- which (substr(names(covs_temporal), nchar(names(covs_temporal)) - 3, nchar(names(covs_temporal)))==year)
    
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
  #dat_all <- na.omit(dat_all)
  
  # which rows in dat_covs contain an NA
  outside_idx <- attr(na.omit(dat_covs), 'na.action')
  
  #stopifnot(is.null(outside_idx))
  
  if (!is.null(outside_idx)) {
    # remove these from dat_all
    dat_all <- dat_all[-outside_idx, ] # this step removes two points from anopheles_bias dataset which fall within the extent but off land
  }
  
  
  ncpu <- 20
  nboot <- ncpu
  
  # get random bootstraps of the data (minimum 5 pres/5 abs)
  data_list <- replicate(nboot,
                         subsample(dat_all,
                                   nrow(dat_all),
                                   replace = TRUE,
                                   minimum = c(5, 5)),
                         simplify = FALSE)
  
  # initialize the cluster
  sfInit(parallel = TRUE, cpus = ncpu)
  sfLibrary(seegSDM)
  
  model_list <- sfLapply(data_list,
                         runBRT,
                         gbm.x = 5:ncol(dat_all),
                         gbm.y = 1,
                         n.folds = 5,
                         pred.raster = covs_current,
                         gbm.coords = 2:3,
                         wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)))
  
  # get cv statistics in parallel
  stat_lis <- sfLapply(model_list, getStats)
  
  # summarise all the ensembles
  preds <- stack(lapply(model_list, '[[', 4))
  
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
  
  points(dat[dat$PA == 1, 2:3],
         pch = 16,
         cex = 1,
         col = 'blue')
  
  dev.off()
  
}
