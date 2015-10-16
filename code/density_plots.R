# make density plots of land cover types for occurrence maps 

# set output path
outpath <- 'output/monkeys/'

# load model outputs
fascicularis <- brick('output/monkeys/fascicularis/fascicularis.tif')
leonina <- brick('output/monkeys/leonina/leonina.tif')
nemestrina <- brick('output/monkeys/nemestrina/nemestrina.tif') 

# load raster covariates
covs_current <- brick('data/clean/raster/raster.grd')

# add urban and build up 5km percentage covariate
# load raster
urban <- brick('data/raw/raster/covariates/MODIS satellite data/landcover 2012 - each year available - includes wetlands/Proportional/13_Urban_And_Built_Up_5km_Percentage.tif')

# crop urban raster to study extent
urban <- crop(urban, covs_current)

# add urban layer to covs_current raster stack
covs_current <- addLayer(covs_current, urban)

# remove non-land cover covariates 
nonlc_covs <- c('TCB_mean', 
            'TCB_SD',
            'SEAsia_pop', 
            'EVI_mean',
            'EVI_SD',
            'LSTday_mean',
            'LSTday_SD',
            'TCW_mean',
            'TCW_SD',
            'SRTM_elevation')

covs_current <- covs_current[[which(!(names(covs_current)%in% nonlc_covs))]] 

names(covs_current)
# layers with scales of 0-100

smallscale <- which(maxValue(covs_current) == 1) #c('forest_intact_2012', 'forest_disturbed_2012')
for(i in smallscale) covs_current[[i]] <- covs_current[[i]] * 100


# rename the urban layer
names(covs_current)[names(covs_current)=='X13_Urban_And_Built_Up_5km_Percentage']<-'Urban_and_built-up'

######## monkeys

# loop through generating a density plot for each covariate for each monkey species
for (species in c('fascicularis', 'leonina', 'nemestrina')) {
  
  # species <- 'fascicularis'
  
  output <- get(species)
  #plot(output[[1]])
  
  #set presence threshold 
  output <- output[[1]] > 0.75
  #plot (output)
  
  # set all 0 values to NA
  output[output==0] <- NA
  
  
  # set output path for this species
  outpath_species <- paste0(outpath,
                            species,
                            '/')
  
  png(paste0(outpath_species,
             'density_plots.png'),
      width = 2400,
      height = 3000,
      pointsize = 30,)
  
  par(mfrow = n2mfrow(nlayers(covs_current)))
    
  for (i in 1:nlayers(covs_current)) {
    
    i <- 1
    
    ras <- covs_current[[i]]
    #plot(ras)
    
    # mask covariate layer with output 
    ras2 <- mask(ras, output)
    
    # remove any NAs
    ras3 <- na.omit(getValues(ras2))
    
    # all cell values
    all_vals <- na.omit(getValues(ras))


    # remove 0 cover values
    all_vals <- all_vals[all_vals > 0]
    ras3 <- ras3[ras3 > 0]

    nb <- seq(from=0, to=100, by=1)
    
    # get histograms
    all <- hist(all_vals, breaks = nb, plot=TRUE)
    pres <- hist(ras3, breaks = nb, plot=TRUE)

    # scale to 0
    all$counts <- all$counts / sum(all$counts)
    pres$counts <- pres$counts / sum(pres$counts)

    ratio <- pres
    ratio$counts <- ratio$counts / all$counts

    newvals <- rep(ratio$mids, round(ratio$counts * 10000))

    # make density estimate
    d <- density(newvals, from = 10, to = 90, bw=6)
    
    # make data frame using density values so the plot can be coloured
    poly_d <- data.frame(x=d$x, y=d$y)
    toprow <- c(poly_d[1,1],0)
    bottomrow <- c(poly_d[nrow(poly_d), 1], 0)
    
    # add bottomrow to bottom of dataframe
    poly_d <- rbind(poly_d, bottomrow)
    # add toprow to top of dataframe
    poly_d <- rbind(toprow, poly_d)
    
    # plotting window
    plot(d,
         main = paste0(names(covs_current[[i]]), '_', species),
         xlab='% cover',
         xlim=c(10,90),
         ylim=c(0, 0.04),
         type='n')
    
    # plot density points
    polygon(poly_d, col=grey(0.4), border='grey')
    
  }
  
  dev.off()
}

########## vectors

# set output path
outpath <- 'output/vectors/'

# load model outputs
vector_group <- brick('output/vectors/vector_group/vector_group.tif')
dirus_complex <- brick('output/vectors/dirus_complex/dirus_complex.tif')
leucosphyrus_complex <- brick('output/vectors/leucosphyrus_complex/leucosphyrus_complex.tif')  
dirus <- brick('output/vectors/dirus/dirus.tif')

# loop through plotting a histogram for each covariate for each vector species
for (species in c('vector_group', 'dirus_complex', 'leucosphyrus_complex', 'dirus')) {
    
  output <- get(species)
  #plot(output[[1]])
  
  #set presence threshold 
  output <- output[[1]] > 0.75
  #plot (output)
  
  # set all 0 values to NA
  output[output==0] <- NA
  
  
  # set output path for this species
  outpath_species <- paste0(outpath,
                            species,
                            '/')
  
  png(paste0(outpath_species,
             'density_plots.png'),
      width = 2400,
      height = 3000,
      pointsize = 30,)
  
  par(mfrow = n2mfrow(nlayers(covs_current)))
  
  for (i in 1:nlayers(covs_current)) {
    
    #i <- 1
    
    ras <- covs_current[[i]]
    #plot(ras)
    
    # mask covariate layer with output 
    ras2 <- mask(ras, output)
    
    
    ras3 <- na.omit(getValues(ras2))
    
    # all cell values
    all_vals <- na.omit(getValues(ras))
    
    
    # remove 0 cover values
    all_vals <- all_vals[all_vals > 0]
    ras3 <- ras3[ras3 > 0]
    
    # set number of bins
    nb <- seq(from=0, to=100, by=1)
    
    # get histograms
    all <- hist(all_vals, breaks = nb, plot=FALSE)
    pres <- hist(ras3, breaks = nb, plot=FALSE)
    
    # scale to 0
    all$counts <- all$counts / sum(all$counts)
    pres$counts <- pres$counts / sum(pres$counts)
    
    ratio <- pres
    ratio$counts <- ratio$counts / all$counts
    
    newvals <- rep(ratio$mids, round(ratio$counts * 10000))
    
    # make density estimate
    d <- density(newvals, from = 10, to = 90, bw=6)
    
    # make data frame using density values so the plot can be coloured
    poly_d <- data.frame(x=d$x, y=d$y)
    toprow <- c(poly_d[1,1],0)
    bottomrow <- c(poly_d[nrow(poly_d), 1], 0)
    
    # add bottomrow to bottom of dataframe
    poly_d <- rbind(poly_d, bottomrow)
    # add toprow to top of dataframe
    poly_d <- rbind(toprow, poly_d)
    
    # plotting window
    plot(d,
         main = paste0(names(covs_current[[i]]), '_', species),
         xlab='% cover',
         xlim=c(10,90),
         ylim=c(0, 0.045),
         type='n')
    
    # plot density points
    polygon(poly_d, col=grey(0.4), border='grey')
    
  }
  
  dev.off()
}




