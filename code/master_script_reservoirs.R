# script to run all analysis

# clear workspace
rm(list = ls())

# ~~~~~~~~~~~~
# load packages
library(seegSDM)
library(raster)
library(snowfall)
library(GRaF)

#~~~~~~~~~~~~~
# load functions file
source('code/functions_reservoirs.R')

# ~~~~~~~~~~~~
# prepare model inputs

# organise covariate layers
source('code/sort_rasters_reservoirs.R')

# organise occurrence data containing polygons 
source('code/sorting_polygons_leonina.R')

# ~~~~~~~~~~~~
# run models, calculate stats and plot outputs
source('code/run_monkeys.R')
source('code/run_vectors.R')

# ~~~~~~~~~~~~
# generate covariate density plots
source('code/density_plots.R')