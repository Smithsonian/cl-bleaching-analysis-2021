#################################################################
##                Confusion Matrix Calculation                 ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
# Loading caret package for confusion matrix function.
library(caret)
load("datasets/bleach_datasets.RData")

##---------------------
##  Non-spatial Model  
##---------------------
### 1998
bleach_1998 <- bleach_1998[, c("bin.score","DHW_25km_max")]
prt.1998 <- predict(m.1998, newdata = bleach_1998, type = "response")
cm.1998 <- confusionMatrix(data = as.factor(as.numeric(prt.1998 > 0.5)), 
                           reference = as.factor(bleach_1998$bin.score))

### 2002
bleach_2002 <- bleach_2002[, c("bin.score","DHW_25km_max")]
prt.2002 <- predict(m.2002, newdata = bleach_2002, type = "response")
cm.2002 <- confusionMatrix(data = as.factor(as.numeric(prt.2002 > 0.5)), 
                           reference = as.factor(bleach_2002$bin.score))

### 2016
bleach_2016 <- bleach_2016[, c("bin.score","DHW_5km_ytd")]
prt.2016 <- predict(m.2016, newdata = bleach_2016, type = "response")
cm.2016 <- confusionMatrix(data = as.factor(as.numeric(prt.2016 > 0.5)), 
                           reference = as.factor(bleach_2016$bin.score))

### 2017
bleach_2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd")]
prt.2017 <- predict(m.2017, newdata = bleach_2017, type = "response")
cm.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                           reference = as.factor(bleach_2017$bin.score))

### 2020
bleach_2020 <- bleach_2020[, c("bin.score","DHW_5km_max")]
prt.2020 <- predict(m.2020, newdata = bleach_2020, type = "response")
cm.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                           reference = as.factor(bleach_2020$bin.score))

load("datasets/bleach_datasets.RData")

##---------------------
##  Zero Eigenvectors  
##---------------------
### 1998
bleach_1998 <- bleach_1998[, c("bin.score","DHW_25km_max")]
bleach_1998 <- cbind(bleach_1998, eigenvectors.1998)
bleach_1998[,3:ncol(bleach_1998)] <- 0
prt.1998 <- predict(mS.1998, newdata = bleach_1998, type = "response")
cmZS.1998 <- confusionMatrix(data = as.factor(as.numeric(prt.1998 > 0.5)), 
                             reference = as.factor(bleach_1998$bin.score))

### 2002
bleach_2002 <- bleach_2002[, c("bin.score","DHW_25km_max")]
bleach_2002 <- cbind(bleach_2002, eigenvectors.2002)
bleach_2002[,3:ncol(bleach_2002)] <- 0
prt.2002 <- predict(mS.2002, newdata = bleach_2002, type = "response")
cmZS.2002 <- confusionMatrix(data = as.factor(as.numeric(prt.2002 > 0.5)), 
                             reference = as.factor(bleach_2002$bin.score))

### 2016
bleach_2016 <- bleach_2016[, c("bin.score","DHW_5km_ytd")]
bleach_2016 <- cbind(bleach_2016, eigenvectors.2016)
bleach_2016[,3:ncol(bleach_2016)] <- 0
prt.2016 <- predict(mS.2016, newdata = bleach_2016, type = "response")
cmZS.2016 <- confusionMatrix(data = as.factor(as.numeric(prt.2016 > 0.5)), 
                             reference = as.factor(bleach_2016$bin.score))

### 2017
bleach_2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd")]
bleach_2017 <- cbind(bleach_2017, eigenvectors.2017)
bleach_2017[,3:ncol(bleach_2017)] <- 0
prt.2017 <- predict(mS.2017, newdata = bleach_2017, type = "response")
cmZS.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                             reference = as.factor(bleach_2017$bin.score))

### 2020
bleach_2020 <- bleach_2020[, c("bin.score","DHW_5km_max")]
bleach_2020 <- cbind(bleach_2020, eigenvectors.2020)
bleach_2020[,3:ncol(bleach_2020)] <- 0
prt.2020 <- predict(mS.2020, newdata = bleach_2020, type = "response")
cmZS.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                             reference = as.factor(bleach_2020$bin.score))

load("datasets/bleach_datasets.RData")

##-------------------------
##  Non-zero Eigenvectors  
##-------------------------
### 1998
bleach_1998 <- bleach_1998[, c("bin.score","DHW_25km_max")]
bleach_1998 <- cbind(bleach_1998, eigenvectors.1998)
prt.1998 <- predict(mS.1998, newdata = bleach_1998, type = "response")
cmS.1998 <- confusionMatrix(data = as.factor(as.numeric(prt.1998 > 0.5)), 
                            reference = as.factor(bleach_1998$bin.score))

### 2002
bleach_2002 <- bleach_2002[, c("bin.score","DHW_25km_max")]
bleach_2002 <- cbind(bleach_2002, eigenvectors.2002)
prt.2002 <- predict(mS.2002, newdata = bleach_2002, type = "response")
cmS.2002 <- confusionMatrix(data = as.factor(as.numeric(prt.2002 > 0.5)), 
                            reference = as.factor(bleach_2002$bin.score))

### 2016
bleach_2016 <- bleach_2016[, c("bin.score","DHW_5km_ytd")]
bleach_2016 <- cbind(bleach_2016, eigenvectors.2016)
prt.2016 <- predict(mS.2016, newdata = bleach_2016, type = "response")
cmS.2016 <- confusionMatrix(data = as.factor(as.numeric(prt.2016 > 0.5)), 
                            reference = as.factor(bleach_2016$bin.score))

### 2017
bleach_2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd")]
bleach_2017 <- cbind(bleach_2017, eigenvectors.2017)
prt.2017 <- predict(mS.2017, newdata = bleach_2017, type = "response")
cmS.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                            reference = as.factor(bleach_2017$bin.score))

### 2020
bleach_2020 <- bleach_2020[, c("bin.score","DHW_5km_max")]
bleach_2020 <- cbind(bleach_2020, eigenvectors.2020)
prt.2020 <- predict(mS.2020, newdata = bleach_2020, type = "response")
cmS.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                            reference = as.factor(bleach_2020$bin.score))
