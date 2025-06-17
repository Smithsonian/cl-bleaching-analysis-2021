####################################################################
##  2020 Analysis - DHW(2020) * DHWprior(2016) - DHWprior(2016)   ##
##                    By: Migdonio A. González                    ##
##                     Edited by: Sean Connolly                   ##
####################################################################
# Top 3 candidates.
candidate.number <- 121 # 7 EV
# candidate.number <- 201 # 8 EV
# candidate.number <- 221 # 8 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- c("DHW_5km_max","DHW_5km_max_16")
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," - ",dhw.column[2]))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," +. "," - ",dhw.column[2]))

# Extracting the spatial structure candidate for this year.
listw.DHWM16 <- candidates_DHWM[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_DHWM16$vector.index[selected_DHWM16$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.DHWM16 <- bleach_2020[,c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.DHWM16 <- mem(listw = listw.DHWM16, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.DHWM16 <- glm(formula1, 
                family = binomial, 
                data = bleach.reduced.DHWM16)

# Calculating the residuals of this model.
residuals.DHWM16 <- residuals(m.DHWM16, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.DHWM16 <- mems.DHWM16[,mem.index]
print(ncol(eigenvectors.DHWM16))
colnames(eigenvectors.DHWM16) <- paste0("X",1:ncol(eigenvectors.DHWM16))
bleach.reduced.DHWM16 <- cbind(bleach.reduced.DHWM16, eigenvectors.DHWM16)

# Fitting GLM model with SAC accounted for (for comparison).
mS.DHWM16 <- glm(formula2, 
                 family = binomial, 
                 data = bleach.reduced.DHWM16)

# Calculating the residuals for this model.
residuals.S.DHWM16 <- residuals(mS.DHWM16, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.DHWM16 <- create_variogram_object_RA("DHW_5km_max_2016", 
                                            bleach_2020, 
                                            residuals.DHWM16)

variogS.DHWM16 <- create_variogram_object_RA("DHW_5km_max_2016", 
                                             bleach_2020, 
                                             residuals.S.DHWM16)

par(mfrow = c(2,1))
plot(variog.DHWM16, main = "DHW Max 2016 - No SAC")
abline(h = variog.DHWM16$var.mark)
plot(variogS.DHWM16, main = "DHW Max 2016 - SAC")
abline(h = variogS.DHWM16$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.DHWM16 <- correlog(coordinates_2020[,1],
#                           coordinates_2020[,2],
#                           increment = 5,
#                           residuals.DHWM16)
# 
# correlogS.DHWM16 <- correlog(coordinates_2020[,1],
#                            coordinates_2020[,2],
#                            increment = 5,
#                            residuals.S.DHWM16)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.DHWM16, main = "DHWM16 - No SAC")
# abline(h = 0)
# plot(correlogS.DHWM16, main = "DHWM16 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.DHWM16 <- generate_prediction_dataframe(dhw.column,
                                                bleach.reduced.DHWM16, 
                                                nevs = ncol(eigenvectors.DHWM16))

# Calculating predictions.
pr.DHWM16 <- generate_predictions(dhw.column, m.DHWM16, newdata.DHWM16$no.ev)
pr.S.DHWM16 <- generate_predictions(dhw.column, mS.DHWM16, newdata.DHWM16$ev)

# Creating plot objects.
plot.DHWM16 <- create_plot_object(dhw.column, 
                                  " ", 
                                  "DHW in 2016",
                                  " ",
                                  pr.DHWM16)
plot.S.DHWM16 <- create_plot_object(dhw.column,
                                    " ", 
                                    "DHW in 2016",
                                    " ",
                                    pr.S.DHWM16)

# Plotting.
grid.arrange(plot.DHWM16, plot.S.DHWM16, nrow = 2)