#################################################################
##       Spatial Structure Candidates Selection Process        ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
#### 1998 analysis. 
# Generating candidates.
candidates_1998 <- listw.candidates(coordinates_1998,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = c(50,75,100,125,150),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_1998 <- select_weighting_matrix(bleach_1998,
                                         "DHW_25km_max",
                                         bin.score ~ DHW_25km_max,
                                         bin.score ~ DHW_25km_max +.,
                                         candidates_1998)

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
selected_2002 <- select_weighting_matrix(bleach_2002,
                                         "DHW_25km_max",
                                         bin.score ~ DHW_25km_max,
                                         bin.score ~ DHW_25km_max +.,
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
selected_2016 <- select_weighting_matrix(bleach_2016,
                                         "DHW_5km_ytd",
                                         bin.score ~ DHW_5km_ytd,
                                         bin.score ~ DHW_5km_ytd +.,
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
selected_2017 <- select_weighting_matrix(bleach_2017,
                                         "DHW_5km_ytd",
                                         bin.score ~ DHW_5km_ytd,
                                         bin.score ~ DHW_5km_ytd +.,
                                         candidates_2017)

#### 2020 analysis.
# Generating candidates.
candidates_2020 <- listw.candidates(coordinates_2020,
                                    style = "W",
                                    nb = "dnear",
                                    d1 = 0,
                                    d2 = seq(100,320,20),
                                    weights = "fup",
                                    y_fup = seq(0.25,5,0.25))

# Candidate selection based on the ME function.
selected_2020 <- select_weighting_matrix(bleach_2020,
                                         "DHW_5km_max",
                                         bin.score ~ DHW_5km_max,
                                         bin.score ~ DHW_5km_max +.,
                                         candidates_2020)

# Adding a column that has the absolute value of the Moran's I statistic.
selected_1998 <- selected_1998 %>% mutate(moran_abs = abs(moran))
selected_2002 <- selected_2002 %>% mutate(moran_abs = abs(moran))
selected_2016 <- selected_2016 %>% mutate(moran_abs = abs(moran))
selected_2017 <- selected_2017 %>% mutate(moran_abs = abs(moran))
selected_2020 <- selected_2020 %>% mutate(moran_abs = abs(moran))

#################################################################
##            Saving OR loading the selected tables            ##
#################################################################
# Saving model selection tables.
save(selected_1998, file = "selected/selected1998.RData")
save(selected_2002, file = "selected/selected2002.RData")
save(selected_2016, file = "selected/selected2016.RData")
save(selected_2017, file = "selected/selected2017.RData")
save(selected_2020, file = "selected/selected2020.RData")

# Loading model selection tables (generated in this script).
load("selected/selected1998.RData")
load("selected/selected2002.RData")
load("selected/selected2016.RData")
load("selected/selected2017.RData")
load("selected/selected2020.RData")
