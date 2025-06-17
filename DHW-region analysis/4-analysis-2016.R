#################################################################
##                2016 Analysis - DHW * region                 ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
### Top 3 candidates.
candidate.number <- 1 # 5 EV
# candidate.number <- 61 # 5 EV
# candidate.number <- 41 # 5 EV

# Setting up variables that will be needed throughout the script. 
dhw.column <- "DHW_5km_ytd"
formula1 <- as.formula(paste0("bin.score ~ ",dhw.column,"*region - region"))
formula2 <- as.formula(paste0("bin.score ~ ",dhw.column,"*region +. ","-region"))

# Extracting the spatial structure candidate for this year.
listw.2016 <- candidates_2016[[candidate.number]]

# Eigenvector indeces.
mem.index <- as.integer(unlist(strsplit(selected_2016$vector.index[selected_2016$candidate == candidate.number], " ")))

# Reduced dataframe.
bleach.reduced.2016 <- bleach_2016[,c("bin.score",dhw.column,"region")]

##################################################################
##         Fitting models with and without eigenvectors         ##
##################################################################
# Calculating eigenvectors for the model. 
mems.2016 <- mem(listw = listw.2016, MEM.autocor = "positive")

# Fitting GLM model with no SAC accounted for (for comparison).
m.2016 <- glm(formula1, 
              family = binomial, 
              data = bleach.reduced.2016)

# Calculating the residuals of this model.
residuals.2016 <- residuals(m.2016, type = "pearson")

# Adding eigenvectors to the dataframe. 
eigenvectors.2016 <- mems.2016[,mem.index]
print(ncol(eigenvectors.2016))
colnames(eigenvectors.2016) <- paste0("X",1:ncol(eigenvectors.2016))
bleach.reduced.2016 <- cbind(bleach.reduced.2016, eigenvectors.2016)

# Fitting GLM model with SAC accounted for (for comparison).
mS.2016 <- glm(formula2, 
               family = binomial, 
               data = bleach.reduced.2016)

# Calculating the residuals for this model.
residuals.S.2016 <- residuals(mS.2016, type = "pearson")

##################################################################
##               Creating and plotting variograms               ##
##################################################################
# Calculating and plotting both variograms (before and after accounting for SAC).
variog.2016 <- create_variogram_object_RA("2016", 
                                          bleach_2016, 
                                          residuals.2016)

variogS.2016 <- create_variogram_object_RA("2016", 
                                           bleach_2016, 
                                           residuals.S.2016)

par(mfrow = c(2,1))
plot(variog.2016, main = "2016 - No SAC")
abline(h = variog.2016$var.mark)
plot(variogS.2016, main = "2016 - SAC")
abline(h = variogS.2016$var.mark)

par(mfrow = c(1,1))

##################################################################
##              Creating and plotting correlograms              ##
##################################################################
# Calculating correlogram objects.
# correlog.2016 <- correlog(coordinates_2016[,1],
#                           coordinates_2016[,2],
#                           increment = 5,
#                           residuals.2016)
# 
# correlogS.2016 <- correlog(coordinates_2016[,1],
#                            coordinates_2016[,2],
#                            increment = 5,
#                            residuals.S.2016)

# Plotting the correlograms.
# par(mfrow = c(2,1))
# plot(correlog.2016, main = "2016 - No SAC")
# abline(h = 0)
# plot(correlogS.2016, main = "2016 - SAC")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##              Calculating predictions - Response              ##
##################################################################
# Creating newdata dataframe.
newdata.2016 <- generate_prediction_dataframe(bleach.reduced.2016, 
                                              dhw.column,
                                              nevs = ncol(eigenvectors.2016))

# Calculating predictions.
pr.2016 <- generate_predictions(m.2016, newdata.2016$no.ev, dhw.column)
pr.2016$region <- factor(pr.2016$region, levels = c("North","Central","South"))

pr.S.2016 <- generate_predictions(mS.2016, newdata.2016$ev, dhw.column)
pr.S.2016$region <- factor(pr.S.2016$region, levels = c("North","Central","South"))

# Creating plot objects.
plot.2016 <- create_plot_object("2016", dhw.column, " ", " ", pr.2016)
plot.S.2016 <- create_plot_object("2016", dhw.column, " ", " ", pr.S.2016)

# Plotting.
grid.arrange(plot.2016, plot.S.2016, nrow = 2)
