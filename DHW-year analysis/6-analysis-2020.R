#################################################################
##                 2020 Analysis - DHW * year                  ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
### Top 3 candidates.
candidate.number <- 121 # 6 EV
# candidate.number <- 201 # 9 EV 
# candidate.number <- 141 # 9 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- "DHW_5km_max"
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column,"+. "))

# Extracting the spatial structure candidate for this year.
listw.2020 <- candidates_2020[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_2020$vector.index[selected_2020$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.2020 <- bleach_2020[, c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.2020 <- mem(listw = listw.2020, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.2020 <- glm(formula1, 
              family = binomial, 
              data = bleach.reduced.2020)

# Calculating the residuals of this model.
residuals.2020 <- residuals(m.2020, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.2020 <- mems.2020[,mem.index]
print(ncol(eigenvectors.2020))
colnames(eigenvectors.2020) <- paste0("X",1:ncol(eigenvectors.2020))
bleach.reduced.2020 <- cbind(bleach.reduced.2020, eigenvectors.2020)

# Fitting GLM model with SAC accounted for (for comparison).
mS.2020 <- glm(formula2, 
               family = binomial, 
               data = bleach.reduced.2020)

# Calculating the residuals for this model.
residuals.S.2020 <- residuals(mS.2020, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.2020 <- create_variogram_object_RA("2020", 
                                          bleach_2020, 
                                          residuals.2020)

variogS.2020 <- create_variogram_object_RA("2020", 
                                           bleach_2020, 
                                           residuals.S.2020)

par(mfrow = c(2,1))
plot(variog.2020, main = "2020 - No SAC")
abline(h = variog.2020$var.mark)
plot(variogS.2020, main = "2020 - SAC")
abline(h = variogS.2020$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.2020 <- correlog(coordinates_2020[,1],
#                           coordinates_2020[,2],
#                           increment = 5,
#                           residuals.2020)
# 
# correlogS.2020 <- correlog(coordinates_2020[,1],
#                            coordinates_2020[,2],
#                            increment = 5,
#                            residuals.S.2020)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.2020, main = "2020 - No SAC")
# abline(h = 0)
# plot(correlogS.2020, main = "2020 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.2020 <- generate_prediction_dataframe(bleach.reduced.2020, 
                                              dhw.column,
                                              nevs = ncol(eigenvectors.2020))

# Calculating predictions.
pr.2020 <- generate_predictions(m.2020, 
                                newdata.2020$no.ev,
                                dhw.column)
pr.S.2020 <- generate_predictions(mS.2020, 
                                  newdata.2020$ev,
                                  dhw.column)

# Creating plot objects.
plot.2020 <- create_plot_object("2020 - No Eigenvectors", 
                                dhw.column,
                                x.label = " ",
                                pr.2020)
plot.S.2020 <- create_plot_object("2020 - Eigenvectors", 
                                  dhw.column,
                                  x.label = " ",
                                  pr.S.2020)

# Plotting.
grid.arrange(plot.2020, plot.S.2020, nrow = 2)
