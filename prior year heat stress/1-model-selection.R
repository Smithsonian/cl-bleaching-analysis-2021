#################################################################
##       Spatial Structure Candidates Selection Process        ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################

#### DHWM16 analysis. 
# Generating candidates (same candidates for all models).
candidates_DHWM <- listw.candidates(coordinates_2020,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = seq(100,320,20),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_DHWM16 <- select_weighting_matrix(bin.score ~ DHW_5km_max * DHW_5km_max_16 - DHW_5km_max_16,
                                           bin.score ~ DHW_5km_max * DHW_5km_max_16 +. - DHW_5km_max_16,
                                           "DHW_5km_max_16",
                                           bleach_2020,
                                           candidates_DHWM)


#### DHWM17 analysis. 
# Candidate selection based on the ME function.
selected_DHWM17 <- select_weighting_matrix(bin.score ~ DHW_5km_max * DHW_5km_max_17 - DHW_5km_max_17,
                                           bin.score ~ DHW_5km_max * DHW_5km_max_17 +. - DHW_5km_max_17,
                                           "DHW_5km_max_17",
                                           bleach_2020,
                                           candidates_DHWM)


#### DHWM1617 analysis. 
# Candidate selection based on the ME function.
selected_DHWM1617 <- select_weighting_matrix(bin.score ~ DHW_5km_max * DHW_5km_max_1617 - DHW_5km_max_1617,
                                             bin.score ~ DHW_5km_max * DHW_5km_max_1617 +. - DHW_5km_max_1617,
                                             "DHW_5km_max_1617",
                                             bleach_2020,
                                             candidates_DHWM)


# Adding a column that has the absolute value of the Moran's I statistic.
selected_DHWM16 <- selected_DHWM16 %>% mutate(moran_abs = abs(moran))
selected_DHWM17 <- selected_DHWM17 %>% mutate(moran_abs = abs(moran))
selected_DHWM1617 <- selected_DHWM1617 %>% mutate(moran_abs = abs(moran))

#################################################################
##            Saving OR loading the selected tables            ##
#################################################################
# Saving model selection tables.
save(selected_DHWM16, file = "selected/selected16.RData")
save(selected_DHWM17, file = "selected/selected17.RData")
save(selected_DHWM1617, file = "selected/selected1617.RData")

# Loading model selection tables (generated in this script).
load("selected/selected16.RData")
load("selected/selected17.RData")
load("selected/selected1617.RData")