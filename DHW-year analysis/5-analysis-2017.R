#################################################################
##                 2017 Analysis - DHW * year                  ##
##                   By: Migdonio A. Gonz�lez                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
### Top 3 candidates.
candidate.number <- 21 # 4 EV
# candidate.number <- 41 # 5 EV 
# candidate.number <- 61 # 5 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- "DHW_5km_ytd"
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column,"+. "))

# Extracting the spatial structure candidate for this year.
listw.2017 <- candidates_2017[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_2017$vector.index[selected_2017$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.2017 <- bleach_2017[, c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.2017 <- mem(listw = listw.2017, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.2017 <- glm(formula1, 
              family = binomial, 
              data = bleach.reduced.2017)

# Calculating the residuals of this model.
residuals.2017 <- residuals(m.2017, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.2017 <- mems.2017[,mem.index]
print(ncol(eigenvectors.2017))
colnames(eigenvectors.2017) <- paste0("X",1:ncol(eigenvectors.2017))
bleach.reduced.2017 <- cbind(bleach.reduced.2017, eigenvectors.2017)

# Fitting GLM model with SAC accounted for (for comparison).
mS.2017 <- glm(formula2, 
               family = binomial, 
               data = bleach.reduced.2017)

# Calculating the residuals for this model.
residuals.S.2017 <- residuals(mS.2017, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.2017 <- create_variogram_object_RA("2017", 
                                          bleach_2017, 
                                          residuals.2017)

variogS.2017 <- create_variogram_object_RA("2017", 
                                           bleach_2017, 
                                           residuals.S.2017)

par(mfrow = c(2,1))
plot(variog.2017, main = "2017 - No SAC")
abline(h = variog.2017$var.mark)
plot(variogS.2017, main = "2017 - SAC")
abline(h = variogS.2017$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.2017 <- correlog(coordinates_2017[,1],
#                           coordinates_2017[,2],
#                           increment = 5,
#                           residuals.2017)
# 
# correlogS.2017 <- correlog(coordinates_2017[,1],
#                            coordinates_2017[,2],
#                            increment = 5,
#                            residuals.S.2017)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.2017, main = "2017 - No SAC")
# abline(h = 0)
# plot(correlogS.2017, main = "2017 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.2017 <- generate_prediction_dataframe(bleach.reduced.2017, 
                                              dhw.column,
                                              nevs = ncol(eigenvectors.2017))

# Calculating predictions.
pr.2017 <- generate_predictions(m.2017, 
                                newdata.2017$no.ev,
                                dhw.column)
pr.S.2017 <- generate_predictions(mS.2017, 
                                  newdata.2017$ev,
                                  dhw.column)

# Creating plot objects.
plot.2017 <- create_plot_object("2017 - No Eigenvectors", 
                                dhw.column,
                                x.label = " ",
                                pr.2017)
plot.S.2017 <- create_plot_object("2017 - Eigenvectors", 
                                  dhw.column,
                                  x.label = " ",
                                  pr.S.2017)

# Plotting.
grid.arrange(plot.2017, plot.S.2017, nrow = 2)
