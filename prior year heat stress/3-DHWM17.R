####################################################################
##  2020 Analysis - DHW(2020) * DHWprior(2017) - DHWprior(2017)   ##
##                    By: Migdonio A. González                    ##
##                     Edited by: Sean Connolly                   ##
####################################################################
# Top 3 candidates.
candidate.number <- 121 # 7 EV
# candidate.number <- 141 # 7 EV
# candidate.number <- 201 # 9 EV 

# Setting up variables that will be needed throughout the script. 
dhw.column <- c("DHW_5km_max","DHW_5km_max_17")
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," - ",dhw.column[2]))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," +. "," - ",dhw.column[2]))

# Extracting the spatial structure candidate for this year.
listw.DHWM17 <- candidates_DHWM[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_DHWM17$vector.index[selected_DHWM17$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.DHWM17 <- bleach_2020[,c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.DHWM17 <- mem(listw = listw.DHWM17, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.DHWM17 <- glm(formula1, 
                family = binomial, 
                data = bleach.reduced.DHWM17)

# Calculating the residuals of this model.
residuals.DHWM17 <- residuals(m.DHWM17, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.DHWM17 <- mems.DHWM17[,mem.index]
print(ncol(eigenvectors.DHWM17))
colnames(eigenvectors.DHWM17) <- paste0("X",1:ncol(eigenvectors.DHWM17))
bleach.reduced.DHWM17 <- cbind(bleach.reduced.DHWM17, eigenvectors.DHWM17)

# Fitting GLM model with SAC accounted for (for comparison).
mS.DHWM17 <- glm(formula2, 
                 family = binomial, 
                 data = bleach.reduced.DHWM17)

# Calculating the residuals for this model.
residuals.S.DHWM17 <- residuals(mS.DHWM17, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.DHWM17 <- create_variogram_object_RA("DHW_5km_max_2017", 
                                            bleach_2020, 
                                            residuals.DHWM17)

variogS.DHWM17 <- create_variogram_object_RA("DHW_5km_max_2017", 
                                             bleach_2020, 
                                             residuals.S.DHWM17)

par(mfrow = c(2,1))
plot(variog.DHWM17, main = "DHW Max 2017 - No SAC")
abline(h = variog.DHWM17$var.mark)
plot(variogS.DHWM17, main = "DHW Max 2017 - SAC")
abline(h = variogS.DHWM17$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.DHWM17 <- correlog(coordinates_2020[,1],
#                           coordinates_2020[,2],
#                           increment = 5,
#                           residuals.DHWM17)
# 
# correlogS.DHWM17 <- correlog(coordinates_2020[,1],
#                            coordinates_2020[,2],
#                            increment = 5,
#                            residuals.S.DHWM17)

# Loading correlogram objects.
# load("correlograms/correlogramsDHWM17.RData")

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.DHWM17, main = "DHWM17 - No SAC")
# abline(h = 0)
# plot(correlogS.DHWM17, main = "DHWM17 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.DHWM17 <- generate_prediction_dataframe(dhw.column,
                                                bleach.reduced.DHWM17, 
                                                nevs = ncol(eigenvectors.DHWM17))

# Calculating predictions.
pr.DHWM17 <- generate_predictions(dhw.column, m.DHWM17, newdata.DHWM17$no.ev)
pr.S.DHWM17 <- generate_predictions(dhw.column, mS.DHWM17, newdata.DHWM17$ev)

# Creating plot objects.
plot.DHWM17 <- create_plot_object(dhw.column, 
                                  " ", 
                                  "DHW in 2017",
                                  " ",
                                  pr.DHWM17)
plot.S.DHWM17 <- create_plot_object(dhw.column,
                                    " ", 
                                    "DHW in 2017",
                                    " ",
                                    pr.S.DHWM17)

# Plotting.
grid.arrange(plot.DHWM17, plot.S.DHWM17, nrow = 2)