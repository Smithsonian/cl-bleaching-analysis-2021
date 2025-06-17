#################################################################
##                 2002 Analysis - DHW * year                  ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
### Top 3 candidates.
candidate.number <- 81 # 3 EV
# candidate.number <- 41 # 4 EV 
# candidate.number <- 61 # 4 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- "DHW_25km_max"
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column,"+. "))

# Extracting the spatial structure candidate for this year.
listw.2002 <- candidates_2002[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_2002$vector.index[selected_2002$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.2002 <- bleach_2002[, c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.2002 <- mem(listw = listw.2002, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.2002 <- glm(formula1, 
              family = binomial, 
              data = bleach.reduced.2002)

# Calculating the residuals of this model.
residuals.2002 <- residuals(m.2002, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.2002 <- mems.2002[,mem.index]
print(ncol(eigenvectors.2002))
colnames(eigenvectors.2002) <- paste0("X",1:ncol(eigenvectors.2002))
bleach.reduced.2002 <- cbind(bleach.reduced.2002, eigenvectors.2002)

# Fitting GLM model with SAC accounted for (for comparison).
mS.2002 <- glm(formula2, 
               family = binomial, 
               data = bleach.reduced.2002)

# Calculating the residuals for this model.
residuals.S.2002 <- residuals(mS.2002, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.2002 <- create_variogram_object_RA("2002", 
                                          bleach_2002, 
                                          residuals.2002)

variogS.2002 <- create_variogram_object_RA("2002", 
                                           bleach_2002, 
                                           residuals.S.2002)

par(mfrow = c(2,1))
plot(variog.2002, main = "2002 - No SAC")
abline(h = variog.2002$var.mark)
plot(variogS.2002, main = "2002 - SAC")
abline(h = variogS.2002$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.2002 <- correlog(coordinates_2002[,1],
#                           coordinates_2002[,2],
#                           increment = 5,
#                           residuals.2002)
# 
# correlogS.2002 <- correlog(coordinates_2002[,1],
#                            coordinates_2002[,2],
#                            increment = 5,
#                            residuals.S.2002)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.2002, main = "2002 - No SAC")
# abline(h = 0)
# plot(correlogS.2002, main = "2002 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.2002 <- generate_prediction_dataframe(bleach.reduced.2002, 
                                              dhw.column,
                                              nevs = ncol(eigenvectors.2002))

# Calculating predictions.
pr.2002 <- generate_predictions(m.2002, 
                                newdata.2002$no.ev,
                                dhw.column)
pr.S.2002 <- generate_predictions(mS.2002, 
                                  newdata.2002$ev,
                                  dhw.column)

# Creating plot objects.
plot.2002 <- create_plot_object("2002 - No Eigenvectors", 
                                dhw.column,
                                x.label = " ",
                                pr.2002)
plot.S.2002 <- create_plot_object("2002 - Eigenvectors", 
                                  dhw.column,
                                  x.label = " ",
                                  pr.S.2002)

# Plotting.
grid.arrange(plot.2002, plot.S.2002, nrow = 2)
