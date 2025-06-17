#################################################################
##                1998 Analysis - DHW * region                 ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
### Top 3 candidates.
candidate.number <- 41 # 5 EV
# candidate.number <- 61 # 5 EV
# candidate.number <- 81 # 5 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- "DHW_25km_max"
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column,"*region - region"))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column,"*region +. ","-region"))

# Extracting the spatial structure candidate for this year.
listw.1998 <- candidates_1998[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_1998$vector.index[selected_1998$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.1998 <- bleach_1998[,c("bin.score",dhw.column,"region")]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.1998 <- mem(listw = listw.1998, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.1998 <- glm(formula1, 
              family = binomial, 
              data = bleach.reduced.1998)

# Calculating the residuals of this model.
residuals.1998 <- residuals(m.1998, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.1998 <- mems.1998[,mem.index]
colnames(eigenvectors.1998) <- paste0("X",1:ncol(eigenvectors.1998))
bleach.reduced.1998 <- cbind(bleach.reduced.1998, eigenvectors.1998)

# Fitting GLM model with SAC accounted for (for comparison).
mS.1998 <- glm(formula2, 
               family = binomial, 
               data = bleach.reduced.1998)

# Calculating the residuals for this model.
residuals.S.1998 <- residuals(mS.1998, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.1998 <- create_variogram_object_RA("1998", 
                                          bleach_1998, 
                                          residuals.1998)

variogS.1998 <- create_variogram_object_RA("1998", 
                                           bleach_1998, 
                                           residuals.S.1998)

par(mfrow = c(2,1))
plot(variog.1998, main = "1998 - No SAC")
abline(h = variog.1998$var.mark)
plot(variogS.1998, main = "1998 - SAC")
abline(h = variogS.1998$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.1998 <- correlog(coordinates_1998[,1],
#                           coordinates_1998[,2],
#                           increment = 5,
#                           residuals.1998)
# 
# correlogS.1998 <- correlog(coordinates_1998[,1],
#                            coordinates_1998[,2],
#                            increment = 5,
#                            residuals.S.1998)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.1998, main = "1998 - No SAC")
# abline(h = 0)
# plot(correlogS.1998, main = "1998 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.1998 <- generate_prediction_dataframe(bleach.reduced.1998, 
                                              dhw.column,
                                              nevs = ncol(eigenvectors.1998))

# Calculating predictions.
pr.1998 <- generate_predictions(m.1998, newdata.1998$no.ev, dhw.column)
pr.1998$region <- factor(pr.1998$region, levels = c("North","Central","South"))

pr.S.1998 <- generate_predictions(mS.1998, newdata.1998$ev, dhw.column)
pr.S.1998$region <- factor(pr.S.1998$region, levels = c("North","Central","South"))

# Creating plot objects.
plot.1998 <- create_plot_object("1998", dhw.column, " ", " ", pr.1998)
plot.S.1998 <- create_plot_object("1998", dhw.column, " ", " ", pr.S.1998)

# Plotting.
grid.arrange(plot.1998, plot.S.1998, nrow = 2)
