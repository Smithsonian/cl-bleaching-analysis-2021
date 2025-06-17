##################################################################
##            DHW * region Analysis (Initial Setup)             ##
##                   By: Migdonio A. González                   ##
##                   Edited by: Sean Connolly                   ##
##################################################################
# Loading libraries.
library(ggplot2)
library(grid)
library(gridExtra)
library(geoR)
library(MASS)
library(nlme) 
library(boot) 
library(modelr) 
library(dplyr) 

# Spatial Eigenvector Libraries.
library(spdep) 
library(spatialreg) 
library(ncf) 
library(adespatial) 

# Loading the spatial weighting matrix selection function.
source("analysis-functions.R")

# Loading the data.
load("datasets/bleach_datasets.RData")

# Coordinates Matrices.
coordinates_1998 <- as.matrix(bleach_1998[,c("X","Y")])
coordinates_2002 <- as.matrix(bleach_2002[,c("X","Y")])
coordinates_2016 <- as.matrix(bleach_2016[,c("X","Y")])
coordinates_2017 <- as.matrix(bleach_2017[,c("X","Y")])
coordinates_2020 <- as.matrix(bleach_2020[,c("X","Y")])