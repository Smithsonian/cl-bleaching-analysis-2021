#########################################################################
##  Prior year heat stress analysis (2002,2016,2017) (Initial Setup)   ##
##                       By: Migdonio A. González                      ##
##                       Edited by: Sean Connolly                      ##
#########################################################################
# Loading libraries.
library(ggplot2)
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

# Library to show ggplots. 
library(gridExtra)

# Library to read XLSX.
library(readxl)

# Loading the functions necessary to run all scripts in this analysis.
source("analysis-functions.R")

# Loading the data.
load("datasets/bleach_datasets.RData")

# Coordinates Matrices.
coordinates_2002 <- as.matrix(bleach_2002[,c("X","Y")])
coordinates_2016 <- as.matrix(bleach_2016[,c("X","Y")])
coordinates_2017 <- as.matrix(bleach_2017[,c("X","Y")])