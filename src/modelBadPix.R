library(tidyverse)
library(e1071)

dataset <- read_csv('OSBS/badPixModel/modelToGetBad.csv')
dat <- cbind(seq(1, nrow(dataset)), dataset)
colnames(dat)[1]<- "pix_id"
dat$badgood <- factor(dat$badgood)
sp <- split(dat, list(dat$badgood, dat$pixel_crownID))
samples <- lapply(sp, function(x) x[sample(1:nrow(x), nrow(x) * 0.7, FALSE),])
calibration <- do.call(rbind, samples)
validation <- anti_join(dat, calibration, by="pix_id")

sp <- split(calibration, list(calibration$badgood, calibration$pixel_crownID))
samples <- lapply(sp, function(x) x[sample(1:nrow(x), nrow(x) * 0.7, FALSE),])
train <- do.call(rbind, samples)
test <- anti_join(calibration,train, by="pix_id")

svm.train <- train[,-c(1:2)]
svm.train <- Filter(function(x)!all(is.na(x)), svm.train)

svm.test <- test[,-c(1:2)]
svm.test <- Filter(function(x)!all(is.na(x)), svm.test)

#svm_tune <- tune(svm, badgood ~ ., data = svm.train, kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
svm_tune <- tune(svm, train.y = svm.train$badgood, train.x = svm.train[-(colnames(svm.train) == "badgood")], kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

print(svm_tune)
best_mod <- svm_tune$best.model
best_mod_pred <- predict(best_mod, svm.train[-(colnames(svm.train) == "badgood")]) 

error_best_mod <- as.integer(svm.train$badgood) - as.integer(best_mod_pred) 

# this value can be different on your computer
# because the tune method randomly shuffles the data
best_mod_RMSE <- sqrt(mean(error_best_mod^2)) # 1.290738

pred_test <- predict(best_mod, svm.test[-(colnames(svm.test) == "badgood")]) 
error_best_mod.test <- as.integer(svm.test$badgood) - as.integer(pred_test) 
best_mod_RMSE <- sqrt(mean(error_best_mod.test^2)) # 1.290738
summary(error_best_mod.test)

svm.val <- validation[,-c(1:2)]
svm.val <- Filter(function(x)!all(is.na(x)), svm.val)
pred_val <- predict(best_mod, svm.val[-(colnames(svm.val) == "badgood")]) 
error_val <- as.integer(svm.val$badgood) - as.integer(pred_val) 
best_mod_RMSE <- sqrt(mean(error_val^2)) # 1.290738
summary(error_val)
table(error_val)

