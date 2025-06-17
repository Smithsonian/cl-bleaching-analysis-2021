##################################################################
##          Including all years in model - DHW * year           ##
##                   By: Migdonio A. González                   ##
##                   Edited by: Sean Connolly                   ##
##################################################################
# Library to plot. 
library(gridExtra)

# Reducing the bleach dataset to relevant columns. 
bleach_red <- bleach_all[,c("bin.score","DHW","year")]

# Adding zeroes to each eigenvector.
eigenvector.1998 <- create_mem_dataframe(eigenvectors.1998, subset_1998)
eigenvector.2002 <- create_mem_dataframe(eigenvectors.2002, subset_2002)
eigenvector.2016 <- create_mem_dataframe(eigenvectors.2016, subset_2016)
eigenvector.2017 <- create_mem_dataframe(eigenvectors.2017, subset_2017)
eigenvector.2020 <- create_mem_dataframe(eigenvectors.2020, subset_2020)

# Matrix of eigenvectors. 
ev.matrix <- cbind(eigenvector.1998,
                   eigenvector.2002,
                   eigenvector.2016,
                   eigenvector.2017,
                   eigenvector.2020)

colnames(ev.matrix) <- paste0("X",seq(1,ncol(ev.matrix)))

# Adding matrix of eigenvectors to the bleach dataset. 
bleach_red <- cbind(bleach_red, ev.matrix)

# Fitting the model with no eigenvectors. 
m_all <- glm(bin.score ~ DHW * year,
             family = binomial,
             data = bleach_red)

# Extracting the residuals. 
residuals.all <- residuals(m_all, type = "pearson")

# Fitting the model with eigenvectors. 
mS_all <- glm(bin.score ~ DHW * year +.,
              family = binomial,
              data = bleach_red)

# Extracting the residuals.
residuals.S.all <- residuals(mS_all, type = "pearson")

#################################################################
##                     Plotting Variograms                     ##
#################################################################
# Variograms without accounting for SAC. 
variog.1998 <- create_variogram_object_RA("1998",
                                          bleach_all[subset_1998,],
                                          residuals.all[subset_1998])
variog.2002 <- create_variogram_object_RA("2002",
                                          bleach_all[subset_2002,],
                                          residuals.all[subset_2002])
variog.2016 <- create_variogram_object_RA("2016",
                                          bleach_all[subset_2016,],
                                          residuals.all[subset_2016])
variog.2017 <- create_variogram_object_RA("2017",
                                          bleach_all[subset_2017,],
                                          residuals.all[subset_2017])
variog.2020 <- create_variogram_object_RA("2020",
                                          bleach_all[subset_2020,],
                                          residuals.all[subset_2020])

# Variograms accounting for SAC. 
variogS.1998 <- create_variogram_object_RA("1998",
                                           bleach_all[subset_1998,],
                                           residuals.S.all[subset_1998])
variogS.2002 <- create_variogram_object_RA("2002",
                                           bleach_all[subset_2002,],
                                           residuals.S.all[subset_2002])
variogS.2016 <- create_variogram_object_RA("2016",
                                           bleach_all[subset_2016,],
                                           residuals.S.all[subset_2016])
variogS.2017 <- create_variogram_object_RA("2017",
                                           bleach_all[subset_2017,],
                                           residuals.S.all[subset_2017])
variogS.2020 <- create_variogram_object_RA("2020",
                                           bleach_all[subset_2020,],
                                           residuals.S.all[subset_2020])

# Plotting variograms for comparison.
par(mfrow = c(2,1))
plot(variog.1998, main = "1998")
abline(h = variog.1998$var.mark)
plot(variogS.1998, main = "1998")
abline(h = variogS.1998$var.mark)

plot(variog.2002, main = "2002")
abline(h = variog.2002$var.mark)
plot(variogS.2002, main = "2002")
abline(h = variogS.2002$var.mark)

plot(variog.2016, main = "2016")
abline(h = variog.2016$var.mark)
plot(variogS.2016, main = "2016")
abline(h = variogS.2016$var.mark)

plot(variog.2017, main = "2017")
abline(h = variog.2017$var.mark)
plot(variogS.2017, main = "2017")
abline(h = variogS.2017$var.mark)

plot(variog.2020, main = "2020")
abline(h = variog.2020$var.mark)
plot(variogS.2020, main = "2020")
abline(h = variogS.2020$var.mark)

par(mfrow = c(1,1))

#################################################################
##                    Plotting Correlograms                    ##
#################################################################
# Correlograms without accounting for SAC.
# correlog.1998 <- correlog(bleach_all[subset_1998,c("X")],
#                           bleach_all[subset_1998,c("Y")],
#                           residuals.all[subset_1998],
#                           increment = 5)
# correlog.2002 <- correlog(bleach_all[subset_2002,c("X")],
#                           bleach_all[subset_2002,c("Y")],
#                           residuals.all[subset_2002],
#                           increment = 5)
# correlog.2016 <- correlog(bleach_all[subset_2016,c("X")],
#                           bleach_all[subset_2016,c("Y")],
#                           residuals.all[subset_2016],
#                           increment = 5)
# correlog.2017 <- correlog(bleach_all[subset_2017,c("X")],
#                           bleach_all[subset_2017,c("Y")],
#                           residuals.all[subset_2017],
#                           increment = 5)
# correlog.2020 <- correlog(bleach_all[subset_2020,c("X")],
#                           bleach_all[subset_2020,c("Y")],
#                           residuals.all[subset_2020],
#                           increment = 5)
# 
# # Correlograms accounting for SAC.
# correlogS.1998 <- correlog(bleach_all[subset_1998,c("X")],
#                            bleach_all[subset_1998,c("Y")],
#                            residuals.S.all[subset_1998],
#                            increment = 5)
# correlogS.2002 <- correlog(bleach_all[subset_2002,c("X")],
#                            bleach_all[subset_2002,c("Y")],
#                            residuals.S.all[subset_2002],
#                            increment = 5)
# correlogS.2016 <- correlog(bleach_all[subset_2016,c("X")],
#                            bleach_all[subset_2016,c("Y")],
#                            residuals.S.all[subset_2016],
#                            increment = 5)
# correlogS.2017 <- correlog(bleach_all[subset_2017,c("X")],
#                            bleach_all[subset_2017,c("Y")],
#                            residuals.S.all[subset_2017],
#                            increment = 5)
# correlogS.2020 <- correlog(bleach_all[subset_2020,c("X")],
#                            bleach_all[subset_2020,c("Y")],
#                            residuals.S.all[subset_2020],
#                            increment = 5)


# # Plotting correlograms for comparison. 
# par(mfrow = c(2,1))
# plot(correlog.1998, main = "1998")
# abline(h = 0)
# plot(correlogS.1998, main = "1998")
# abline(h = 0)
# 
# plot(correlog.2002, main = "2002")
# abline(h = 0)
# plot(correlogS.2002, main = "2002")
# abline(h = 0)
# 
# plot(correlog.2016, main = "2016")
# abline(h = 0)
# plot(correlogS.2016, main = "2016")
# abline(h = 0)
# 
# plot(correlog.2017, main = "2017")
# abline(h = 0)
# plot(correlogS.2017, main = "2017")
# abline(h = 0)
# 
# plot(correlog.2020, main = "2020")
# abline(h = 0)
# plot(correlogS.2020, main = "2020")
# abline(h = 0)
# 
# par(mfrow = c(1,1))

##################################################################
##                     Moran's I Comparison                     ##
##################################################################
# Moran's I without accounting for SAC. 
# moran.1998 <- moran.mc(residuals.all[subset_1998],
#                        listw.1998,
#                        999)
# moran.ci.1998 <- calculate_moran_ci(moran.1998)
# 
# moran.2002 <- moran.mc(residuals.all[subset_2002],
#                        listw.2002,
#                        999)
# moran.ci.2002 <- calculate_moran_ci(moran.2002)
# 
# moran.2016 <- moran.mc(residuals.all[subset_2016],
#                        listw.2016,
#                        999)
# moran.ci.2016 <- calculate_moran_ci(moran.2016)
# 
# moran.2017 <- moran.mc(residuals.all[subset_2017],
#                        listw.2017,
#                        999)
# moran.ci.2017 <- calculate_moran_ci(moran.2017)
# 
# moran.2020 <- moran.mc(residuals.all[subset_2020],
#                        listw.2020,
#                        999)
# moran.ci.2020 <- calculate_moran_ci(moran.2020)
# 
# # Moran's I accounting for SAC.
# moran.S.1998 <- moran.mc(residuals.S.all[subset_1998],
#                          listw.1998,
#                          999)
# moran.S.ci.1998 <- calculate_moran_ci(moran.S.1998)
# 
# moran.S.2002 <- moran.mc(residuals.S.all[subset_2002],
#                          listw.2002,
#                          999)
# moran.S.ci.2002 <- calculate_moran_ci(moran.S.2002)
# 
# moran.S.2016 <- moran.mc(residuals.S.all[subset_2016],
#                          listw.2016,
#                          999)
# moran.S.ci.2016 <- calculate_moran_ci(moran.S.2016)
# 
# moran.S.2017 <- moran.mc(residuals.S.all[subset_2017],
#                          listw.2017,
#                          999)
# moran.S.ci.2017 <- calculate_moran_ci(moran.S.2017)
# 
# moran.S.2020 <- moran.mc(residuals.S.all[subset_2020],
#                          listw.2020,
#                          999)
# moran.S.ci.2020 <- calculate_moran_ci(moran.S.2020)
# 
# # Plotting comparisons. 
# par(mfrow = c(2,1))
# plot(moran.1998, main = "1998")
# abline(v = c(moran.ci.1998$lower, moran.ci.1998$upper), col = "red", lty = 2)
# plot(moran.S.1998, main = "1998")
# abline(v = c(moran.S.ci.1998$lower, moran.S.ci.1998$upper), col = "red", lty = 2)
# 
# plot(moran.2002, main = "2002")
# abline(v = c(moran.ci.2002$lower, moran.ci.2002$upper), col = "red", lty = 2)
# plot(moran.S.2002, main = "2002")
# abline(v = c(moran.S.ci.2002$lower, moran.S.ci.2002$upper), col = "red", lty = 2)
# 
# plot(moran.2016, main = "2016")
# abline(v = c(moran.ci.2016$lower, moran.ci.2016$upper), col = "red", lty = 2)
# plot(moran.S.2016, main = "2016")
# abline(v = c(moran.S.ci.2016$lower, moran.S.ci.2016$upper), col = "red", lty = 2)
# 
# plot(moran.2017, main = "2017")
# abline(v = c(moran.ci.2017$lower, moran.ci.2017$upper), col = "red", lty = 2)
# plot(moran.S.2017, main = "2017")
# abline(v = c(moran.S.ci.2017$lower, moran.S.ci.2017$upper), col = "red", lty = 2)
# 
# plot(moran.2020, main = "2020")
# abline(v = c(moran.ci.2020$lower, moran.ci.2020$upper), col = "red", lty = 2)
# plot(moran.S.2020, main = "2020")
# abline(v = c(moran.S.ci.2020$lower, moran.S.ci.2020$upper), col = "red", lty = 2)
# 
# par(mfrow = c(1,1))

#################################################################
##                   Calculating Predictions                   ##
#################################################################
# "New data" tables. 
newdata <- generate_prediction_dataframe_full(bleach_all,
                                              nevs = ncol(ev.matrix))

# Predictions for model without accounting for SAC. 
pr.no.sac <- generate_predictions_full(m_all, newdata$no.ev)
plot.no.sac <- create_plot_object_full(" ", " ", pr.no.sac)

# Predictions for model accounting for SAC.
pr.sac <- generate_predictions_full(mS_all, newdata$ev)
plot.sac <- create_plot_object_full(" ", " ", pr.sac)

# Showing the plots. 
grid.arrange(plot.no.sac, plot.sac, nrow = 2)