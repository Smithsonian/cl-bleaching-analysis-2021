####################################################################################
##  2020 Analysis - DHW(2020) * DHWpriormax(2016,2017) - DHWpriormax(2016,2017)   ##
##                            By: Migdonio A. González                            ##
##                             Edited by: Sean Connolly                           ##
####################################################################################
# Top 3 candidates.
candidate.number <- 121 # 7 EV
# candidate.number <- 201 # 8 EV
# candidate.number <- 141 # 8 EV 

# Setting up variables that will be needed throughout the script. 
dhw.column <- c("DHW_5km_max","DHW_5km_max_1617")
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," - ",dhw.column[2]))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column[1],"*",dhw.column[2]," +. "," - ",dhw.column[2]))

# Extracting the spatial structure candidate for this year.
listw.DHWM1617 <- candidates_DHWM[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_DHWM1617$vector.index[selected_DHWM1617$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.DHWM1617 <- bleach_2020[,c("bin.score",dhw.column)]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.DHWM1617 <- mem(listw = listw.DHWM1617, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.DHWM1617 <- glm(formula1, 
                  family = binomial, 
                  data = bleach.reduced.DHWM1617)

# Calculating the residuals of this model.
residuals.DHWM1617 <- residuals(m.DHWM1617, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.DHWM1617 <- mems.DHWM1617[,mem.index]
print(ncol(eigenvectors.DHWM1617))
colnames(eigenvectors.DHWM1617) <- paste0("X",1:ncol(eigenvectors.DHWM1617))
bleach.reduced.DHWM1617 <- cbind(bleach.reduced.DHWM1617, eigenvectors.DHWM1617)

# Fitting GLM model with SAC accounted for (for comparison).
mS.DHWM1617 <- glm(formula2, 
                   family = binomial, 
                   data = bleach.reduced.DHWM1617)

# Calculating the residuals for this model.
residuals.S.DHWM1617 <- residuals(mS.DHWM1617, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.DHWM1617 <- create_variogram_object_RA("DHW_5km_max_2016/2017", 
                                              bleach_2020, 
                                              residuals.DHWM1617)

variogS.DHWM1617 <- create_variogram_object_RA("DHW_5km_max_2016/2017", 
                                               bleach_2020, 
                                               residuals.S.DHWM1617)

par(mfrow = c(2,1))
plot(variog.DHWM1617, main = "DHW Max 2016/2017 - No SAC")
abline(h = variog.DHWM1617$var.mark)
plot(variogS.DHWM1617, main = "DHW Max 2016/2017 - SAC")
abline(h = variogS.DHWM1617$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.DHWM1617 <- correlog(coordinates_2020[,1],
#                           coordinates_2020[,2],
#                           increment = 5,
#                           residuals.DHWM1617)
# 
# correlogS.DHWM1617 <- correlog(coordinates_2020[,1],
#                            coordinates_2020[,2],
#                            increment = 5,
#                            residuals.S.DHWM1617)

# Loading correlogram objects.
# load("correlograms/correlogramsDHWM1617.RData")

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.DHWM1617, main = "DHWM1617 - No SAC")
# abline(h = 0)
# plot(correlogS.DHWM1617, main = "DHWM1617 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.DHWM1617 <- generate_prediction_dataframe(dhw.column,
                                                  bleach.reduced.DHWM1617, 
                                                  nevs = ncol(eigenvectors.DHWM1617))

# Calculating predictions.
pr.DHWM1617 <- generate_predictions(dhw.column, m.DHWM1617, newdata.DHWM1617$no.ev)
pr.S.DHWM1617 <- generate_predictions(dhw.column, mS.DHWM1617, newdata.DHWM1617$ev)

# Creating plot objects.
plot.DHWM1617 <- create_plot_object(dhw.column, 
                                    " ", 
                                    "Max DHW from 2016/2017",
                                    "DHW (°C-week)",
                                    pr.DHWM1617)
plot.S.DHWM1617 <- create_plot_object(dhw.column,
                                      " ", 
                                      "Max DHW from 2016/2017",
                                      "DHW (°C-week)",
                                      pr.S.DHWM1617)

# Plotting.
grid.arrange(plot.DHWM1617, plot.S.DHWM1617, nrow = 2)