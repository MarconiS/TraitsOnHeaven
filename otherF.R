library (leaps)
for (j in 1:5) {
  train.data = data.frame(Y=Y[j], X)
  colnames(train.data)[1] <- "Y"
  head(train.data)
  regfit.full <- regsubsets(Y~., data=train.data, nvmax = 50, method="forward")
}
set.seed (1)

library (leaps)

folds = cut(seq(1,nrow(aug.spectra)), breaks=10, labels=FALSE)
data = data.frame(Y, X)
for (i in 1:10) {
  test_index = which(folds==i, arr.ind=TRUE)
  data.test = data[test_index, ]
  data.train = data[-test_index, ]
  for (j in 1:5) {
    train.data = data.frame(Y=Y[j], X)
    regfit.full <- regsubsets(Y~., data=train.data, nvmax = 50, method="forward")
    
  }
}