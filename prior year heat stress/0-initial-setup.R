##################################################################
##       Prior Year Heat Stress Analysis (Initial Setup)        ##
##                   By: Migdonio A. González                   ##
##                   Edited by: Sean Connolly                   ##
##################################################################
# Loading libraries.
library(ggplot2)
library(gridExtra)
library(geoR)
library(MASS)
library(nlme) 
library(boot) 
library(modelr) 
library(dplyr) 
library(Matrix) 

# Spatial Eigenvector Libraries.
library(spdep) 
library(spatialreg) 
library(ncf) 
library(adespatial) 

# Library to show ggplots. 
library(gridExtra)

# Loading the spatial weighting matrix selection function.
source("analysis-functions.R")

# Loading the data.
load("datasets/bleach_datasets.RData")

# Adding a column to the 2020 dataset that contains the max between DHW_5km_max_16 & 17. 
bleach_2020 <- bleach_2020 %>% mutate(DHW_5km_max_1617 = pmax(DHW_5km_max_16, DHW_5km_max_17))

# Creating coordinates matrix.
coordinates_2020 <- as.matrix(bleach_2020[,c("X","Y")])
