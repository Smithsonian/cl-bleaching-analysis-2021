#################################################################
##                Confusion Matrix Calculation                 ##
##                   By: Migdonio A. González                  ##
##                   Edited by: Sean Connolly                  ##
#################################################################
# Loading caret package for confusion matrix calculation.
library(caret)

##---------------------
##  Non-spatial Model  
##---------------------
bleach.new.2020 <- bleach_2020[, c("bin.score","DHW_5km_max","DHW_5km_max_17")]
prt.2020 <- predict(m.DHWM17, newdata = bleach.new.2020, type = "response")
cm.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                           reference = as.factor(bleach.new.2020$bin.score))

##---------------------
##  Zero Eigenvectors  
##---------------------
bleach.new.2020 <- bleach_2020[, c("bin.score","DHW_5km_max","DHW_5km_max_17")]
bleach.new.2020 <- cbind(bleach.new.2020, eigenvectors.DHWM17)
bleach.new.2020[,4:ncol(bleach.new.2020)] <- 0
prt.2020 <- predict(mS.DHWM17, newdata = bleach.new.2020, type = "response")
cmZS.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                             reference = as.factor(bleach.new.2020$bin.score))


##-------------------------
##  Non-zero Eigenvectors  
##-------------------------
### 2017
bleach.new.2020 <- bleach_2020[, c("bin.score","DHW_5km_max","DHW_5km_max_17")]
bleach.new.2020 <- cbind(bleach.new.2020, eigenvectors.DHWM17)
prt.2020 <- predict(mS.DHWM17, newdata = bleach.new.2020, type = "response")
cmS.2020 <- confusionMatrix(data = as.factor(as.numeric(prt.2020 > 0.5)), 
                            reference = as.factor(bleach.new.2020$bin.score))