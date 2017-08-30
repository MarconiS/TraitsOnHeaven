library(dplyr)
library(data.table)
ensemble.bagging = list()
yhat.bagging = list()
y_train.bagging = list()
data_train.bagging = list()
i=6

fname <- Sys.glob(paste(in.dir,'FinalSet/chosenSets/onePix1Crown_', names[i], '*',sep="" ))
wholeCrowns = pd.DataFrame.from_csv('./OSBS/CrownPix_norm.csv')
cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
nCrowns <- dim(cr.Traits)[1]


token <- 0
ridge.aug <- matrix(nrow = length(eval.set$cr.id -1), ncol = length(fname)+1)
ridge.aug[,1] <- train.data[,names(aug.spectra) %in% names[i]]
for(f in fname){
  token <- token +1
  aug.spectra = read.csv(f)
  aug.spectra$X=aug.spectra$X.1 = aug.spectra$X.2 = aug.spectra$X.3 = aug.spectra$X.4 = aug.spectra$X.5 = 
    aug.spectra$X.6 = aug.spectra$X.7 = aug.spectra$X.8 = aug.spectra$X.9 = aug.spectra$X.10 = NULL
  
  aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
  X<- aug.spectra[grepl("band", names(aug.spectra))]
  X=X[, colSums(is.na(X)) == 0]
  Y <- aug.spectra[,names(aug.spectra) %in% names]
  X_corr = cor(as.matrix(X), Y)
  aug.X <- data.frame(aug.spectra$name, Y, X)
  # Subset data into cal/val by site
  eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
  train.data <- eval.set$train
  test.data <- eval.set$test
  
  # Run calibration PLSR analysis to select optimal number of components
  pls.mod.train <- pls.cal(train.data, 25,nm = names, j, norm = F)
  #calculate number of components given min test PRESS or RMSEP
  optim.ncomps <- opt.comps(pls.mod.train, Y, j)
  pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,optim.ncomps,j, norm = F)
  out.data <- res.out(pred.val.data, train.data,nm = names, test.data, j)

  #bag <- substr(x, nchar(x)-n+1, nchar(x))
  ensemble.bagging[[token]] <- pls.mod.train
  
  #data for bagging
  foo <- predict.pls(pls.mod.train, train.data, nm = names,optim.ncomps,j, norm = F)
  foo
  #data_train.bagging[[token]] <- cbind(eval.set$cr.id, train.data[,names(aug.spectra) %in% names[i]], as.vector(unlist(predict.pls(pls.mod.train, train.data, nm = names,optim.ncomps,j, norm = F))))
  ridge.aug[,token+1] <-as.vector(unlist(predict.pls(pls.mod.train, train.data, nm = names,optim.ncomps,j, norm = F)))[1:42]
  yhat.bagging[[token]] <- pred.val.data
  #y_train.bagging[[token]] <- predict.pls(pls.mod.train, train.data, nm = names,optim.ncomps,j, norm = F)
  #data_train.bagging[[token]] <- train.data[eval.set$cr.id,names(aug.spectra) %in% names[i]]
}
colnames(ridge.aug) <- c("Y", paste("X", seq(1,length(fname)), sep=""))

lambdas <- 10^seq(3, -2, by = -.1)
cv_fit <- cv.glmnet(ridge.aug[,-1], ridge.aug[,1], alpha = 0, lambda = lambdas)
plot(cv_fit)
opt_lambda <- cv_fit$lambda.min
fit <- cv_fit$glmnet.fit
summary(fit)

