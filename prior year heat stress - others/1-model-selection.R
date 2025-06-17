#################################################################
##       Spatial Structure Candidates Selection Process        ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
#### 2002 analysis. 
# Generating candidates.
candidates_2002 <- listw.candidates(coordinates_2002,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = c(50,75,100,125,150),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_2002 <- select_weighting_matrix(bin.score ~ DHW_25km_max * DHW_25km_max_98 - DHW_25km_max_98,
                                         bin.score ~ DHW_25km_max * DHW_25km_max_98 +. - DHW_25km_max_98,
                                         c("DHW_25km_max","DHW_25km_max_98"),
                                         bleach_2002,
                                         candidates_2002)

#### 2016 analysis. 
# Generating candidates.
candidates_2016 <- listw.candidates(coordinates_2016,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = c(50,75,100,125,150),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_2016 <- select_weighting_matrix(bin.score ~ DHW_5km_ytd * DHW_25km_max_02 - DHW_25km_max_02,
                                         bin.score ~ DHW_5km_ytd * DHW_25km_max_02 +. - DHW_25km_max_02,
                                         c("DHW_5km_ytd","DHW_25km_max_02"),
                                         bleach_2016,
                                         candidates_2016)

#### 2017 analysis. 
# Generating candidates.
candidates_2017 <- listw.candidates(coordinates_2017,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = c(100,120,140,160),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_2017 <- select_weighting_matrix(bin.score ~ DHW_5km_ytd * DHW_5km_max_16 - DHW_5km_max_16,
                                         bin.score ~ DHW_5km_ytd * DHW_5km_max_16 +. - DHW_5km_max_16,
                                         c("DHW_5km_ytd","DHW_5km_max_16"),
                                         bleach_2017,
                                         candidates_2017)

# Adding a column that has the absolute value of the Moran's I statistic.
selected_2002 <- selected_2002 %>% mutate(moran_abs = abs(moran))
selected_2016 <- selected_2016 %>% mutate(moran_abs = abs(moran))
selected_2017 <- selected_2017 %>% mutate(moran_abs = abs(moran))

#################################################################
##            Saving OR loading the selected tables            ##
#################################################################
# Saving model selection tables.
save(selected_2002, file = "selected/selected2002.RData")
save(selected_2016, file = "selected/selected2016.RData")
save(selected_2017, file = "selected/selected2017.RData")

# Loading model selection tables (generated in this script).
load("selected/selected2002.RData")
load("selected/selected2016.RData")
load("selected/selected2017.RData")