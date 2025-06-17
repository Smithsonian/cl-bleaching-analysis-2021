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
bleach.new.2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd","DHW_5km_max_16")]
prt.2017 <- predict(m.2017, newdata = bleach.new.2017, type = "response")
cm.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                           reference = as.factor(bleach.new.2017$bin.score))

##---------------------
##  Zero Eigenvectors  
##---------------------
bleach.new.2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd","DHW_5km_max_16")]
bleach.new.2017 <- cbind(bleach.new.2017, eigenvectors.2017)
bleach.new.2017[,4:ncol(bleach.new.2017)] <- 0
prt.2017 <- predict(mS.2017, newdata = bleach.new.2017, type = "response")
cmZS.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                             reference = as.factor(bleach.new.2017$bin.score))

##-------------------------
##  Non-zero Eigenvectors  
##-------------------------
### 2017
bleach.new.2017 <- bleach_2017[, c("bin.score","DHW_5km_ytd","DHW_5km_max_16")]
bleach.new.2017 <- cbind(bleach.new.2017, eigenvectors.2017)
prt.2017 <- predict(mS.2017, newdata = bleach.new.2017, type = "response")
cmS.2017 <- confusionMatrix(data = as.factor(as.numeric(prt.2017 > 0.5)), 
                            reference = as.factor(bleach.new.2017$bin.score))