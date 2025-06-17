##################################################################
##             DHW * year Analysis (Initial Setup)              ##
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

# Loading the spatial weighting matrix selection function.
source("analysis-functions.R")

# Loading the data.
load("datasets/bleach_datasets.RData")
load("datasets/bleach_all.RData")

# Coordinates Matrices.
coordinates_1998 <- as.matrix(bleach_1998[,c("X","Y")])
coordinates_2002 <- as.matrix(bleach_2002[,c("X","Y")])
coordinates_2016 <- as.matrix(bleach_2016[,c("X","Y")])
coordinates_2017 <- as.matrix(bleach_2017[,c("X","Y")])
coordinates_2020 <- as.matrix(bleach_2020[,c("X","Y")])

# For the DHW*year analysis (year subsets).
subset_1998 <- which(bleach_all$year == "1998")
subset_2002 <- which(bleach_all$year == "2002")
subset_2016 <- which(bleach_all$year == "2016")
subset_2017 <- which(bleach_all$year == "2017")
subset_2020 <- which(bleach_all$year == "2020")