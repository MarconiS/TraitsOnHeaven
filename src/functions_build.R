#-------imp.spectra-------------------------------------------------------------------------------------------

imp.spectra <- function(f, pwd)
{
  # Import dry spectra dataset
  aug.data <- read.table(paste(pwd,f,sep=""), header=TRUE,sep=",")
  return(aug.data)
}
#-------cut.set-------------------------------------------------------------------------------------------

cut.set<-function(aug.X,out.dir, c.id, prop = 0.7){
  species <- unique(aug.X$aug.spectra.name)
  train.data <- 0
  test.data <- 0
  j <- 1
  cr.id <- NULL
  for (i in as.character(species)){
    set.seed(1)
    #subset by species
    temp.data <- aug.X[which(aug.X$aug.spectra.name==i),]
    rows <- sample(1:nrow(temp.data),ceiling(prop*nrow(temp.data)))
    foo <- c.id[which(aug.X$aug.spectra.name==i)]
    cr.id <- c(cr.id, foo[rows])
    cal.data = droplevels(temp.data[rows,])
    val.data = droplevels(temp.data[-rows,])

    if(j==1){
      train.data <- cal.data
      test.data <- val.data
    } else {
      train.data <- rbind(train.data,cal.data)
      test.data <- rbind(test.data,val.data)
    }

    j <- j+1
  }
  return(list(train=train.data, test=test.data, cr.id=cr.id))
}

# pls.cal ---------------------------------------------------
pls.cal <- function(train.data, numcomps = 15, nm,j,  normalz = F, lab.train = NULL){

  spectra <- as.matrix(train.data[grepl("band", names(train.data))])
  traits <- as.matrix(train.data[,names(train.data) %in% nm])

  pls.summ <- list()
  spectra_log_dif_snv <- spectra
  #spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  if(normalz){spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))}

  leaf.trait <- traits[,j]
  if(is.character(leaf.trait)){
    #clean
    spectra_log_dif_snv=spectra_log_dif_snv[, colSums(is.na(spectra_log_dif_snv)) == 0]
    tmp.y<-matrix(as.numeric(lab.train$id,ncol=1)) # make numeric matrix
    train.PLS = data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    #tmp.pls<-plsda(y = factor(tmp.y), x = spectra_log_dif_snv,  ncomp = comps,  probMethod = "softmax", type = "class")
    tmp.pls<- plsr(Y ~ X, ncomp=numcomps,validation="LOO", probMethod = "softmax", trace=TRUE, method = "oscorespls", data = train.PLS, type = "class")
  }else{
    train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = plsr(Y ~ X,scale=F, ncomp=numcomps,validation="LOO", trace=TRUE, method = "oscorespls", probMethod = "softmax", data = train.PLS)
  }
  pls.summ[[nm[j]]] <- tmp.pls

  return(pls.summ)
}

#-------opt.comps-------------------------------------------------------------------------------------------
opt.comps <- function(pls.mod.train, Y, j, Class = F){
  ncomps <- NA
  tmp.pls = pls.mod.train[[1]]
  if(!Class){
    ncomps <- which(tmp.pls$validation$PRESS==min(tmp.pls$validation$PRESS[3:length(tmp.pls$validation$PRESS)]))
  }
  return(c(ncomps))
}
#-----predict.pls------------------------------------------------------------------------------------------
predict.pls <- function(pls.mod.train, test.data, optim.ncomps=NULL,j, nm, norm = F, lab.test = NULL){
  spectra <- as.matrix(test.data[grepl("band", names(test.data))])
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  out <- list()
  spectra_log_dif_snv <- spectra
  if(norm){  spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))}
  leaf.trait <- traits[,j]
  if(is.character(leaf.trait)){
    spectra_log_dif_snv=spectra_log_dif_snv[, colSums(is.na(spectra_log_dif_snv)) == 0]

    tmp.y<-matrix(as.numeric(lab.test$id, ncol=1)) # make numeric matrix
    test.PLS <- data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    pred<-predict(pls.mod.train$name,newdata=test.PLS, comps = optim.ncomps)
    accu = (sum(round(pred)==lab.test$id)/length(lab.test$id))
    print(accu)
    #out.data = data.frame(lab.test$name,pred.val.data$,rep(NA, length(lab.test$name)), -9999, -9999)
    temp = data.frame(traits[,j],round(pred), rep(accu,length(traits[,j])))
    names(temp) = c("obs","pred", "R2")

    out[[nm[j]]] <- temp

  }else{
    test.PLS <- data.frame(X=I(spectra_log_dif_snv))
    out[[nm[j]]] <- as.vector(predict(pls.mod.train[[1]], newdata = test.PLS, ncomp=optim.ncomps, type="response"))
  }
  return(out)
}
#-----res.out------------------------------------------------------------------------------------------
res.out <- function(pred.val.data, train.data, test.data, j, nm, labelsSp = NULL)
{
  out <- list()
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  # Build output dataset
  # PLSR Summary statistics
  if(is.character(traits[,j])){
    pred.data <- round(as.vector(pred.val.data))
    correct <- labelsSp$id == as.numeric(pred.data)
    accuracy <- sum(as.numeric(correct))/length(correct)
    temp = data.frame(as.numeric(factor(traits[,j])),as.numeric(pred.data),correct, accuracy)
    names(temp) = c("obs","pred","res", "R2")
  }else{
    pred.data <- as.vector(eval(parse(text = paste('pred.val.data$',nm[j],sep=""))))
    res <- traits[,j]- pred.data
    MSE.test <- mean(res^2)
    RMSE.test <- sqrt(MSE.test)
    R2 <- 1- length(pred.data) * MSE.test / sum((traits[,j]- mean(traits[,j]))^2)
    ### Output val dataset
    temp = data.frame(traits[,j],pred.data,res, R2, MSE.test)
    names(temp) = c("obs","pred","res", "R2", "MSE")
  }

  out[[names(train.data[,names(train.data) %in% names])[j]]] <- temp

  return(temp)
}

#-----scale to 0:1------------------------------------------------------------------------------------------
scalar1 <- function(x) {x / sqrt(sum(x^2))}
#-----decimalplaces------------------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#-----pls_glm--------------------------------------------------------------------------------------$

pls_glm <- function(j = "P_pct", loops = NULL,cores = 2, nrmlz = F, out.dir = NULL, in.dir = NULL, nsites = 2){
  library(foreach)
  library(doParallel)

  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())

  ll = as.character(1:loops)
  mod.out <- foreach(laps = ll, .verbose = T, .multicombine = TRUE) %dopar% {

    library('plsRglm')
    library(readr)
    library(gnm)

    source(paste('../../src/functions_build.R', sep=""))
    source(paste('../../src/src_build.R', sep=""))
    out.dir = paste(getwd(), "/outputs/", sep="")
    in.dir = paste(getwd(), "/inputs/", sep="")

    cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
    nCrowns <- dim(cr.Traits)[1]
    aug.spectra <- read.csv(paste(in.dir, 'Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    aug.spectra$X= NULL
    aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
    X<- aug.spectra[grepl("band", names(aug.spectra))]
    X=X[, colSums(is.na(X)) == 0]
    Y <- as.data.frame(aug.spectra[,names(aug.spectra) %in% j])
    Y <- (round(Y,3))
    if(dim(Y)[2]==1){names(Y) = j}
    aug.X <- data.frame(aug.spectra$name, Y, X)

    # Subset data into cal/val by site
    eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
    train.data <- eval.set$train
    test.data <- eval.set$test
    colnames(train.data)[1] <- "name"
    colnames(test.data)[1] <- "name"
    # write_csv(train.data, paste(in.dir, "datSets/train_for_baseline_", laps,".csv", sep=""))
    # write_csv(test.data, paste(in.dir, "datSets/test_for_baseline_", laps,".csv", sep=""))

    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    Y <- as.vector(train.data[,names(train.data) %in% j])
    #K=nrow(Y)
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }

    train.PLS<- cv.plsRglm(dataY=log(Y),dataX=X, nt=15,NK=1, K=10,modele="pls-glm-family",family=gaussian())
    out <- list()
    #out["summary"] <- summary(train.PLS)
    press <- unlist(kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))

    #chosen to be built
    mod <- plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]),modele="pls-glm-gaussian")
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }

    Y.test <- as.vector(test.data[,names(test.data) %in% j])
    out[["pred"]] <- predict(mod,newdata=X.tst,type="response",comps=as.integer(out["ncomp"]), se.fit = T)
    out[["mod"]] <- mod
    out[["score"]] <- gnm(out$pred$fit ~ (Y.test), family = gaussian())
    out["aic"]<- AIC(out$score)
    #r2 <- summary(lm(pred$fit~ Y.test))
    out
  }
  #print(paste(j, "done"))
  #save(mod.out, file = paste("pls_glm_", j,  sep = ""))
  stopCluster(cl)
  return(mod.out)
}

#rm(list=ls(all=TRUE))   # clear workspace
get_stats <- function()
  library("tsensembler")
  library(readr)
  library(plsRglm)
  library(dplyr)
  NeonSite = "ALL"
  wd = "/ufrc/ewhite/s.marconi/Marconi2018/ModelBuild/"
  setwd(wd)
  source('../src/functions_build.R')
  source('../src/src_build.R')
  source( '../src/functions_infer.R')
  source('../src/src_infer.R')
  #source(paste(getwd(), '/src/pls_par.R', sep=""))


  path = paste(wd, NeonSite, sep="/")
  setwd(path)
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")

  nm = "C"
  j="C_pct"
  get100=T
  method = 'press'
  normz = T
  LOG=F
  glm = T
  if(get100){
    loops = 1000
    if(glm){
      load(file = paste('./pls_glm_',j, sep="" ))
    }else{
      load(file = paste('./pls_',j, sep="" ))
    }
    mod.aic <- rep(NA, loops)
    mod.r2 <- rep(NA, loops)
    weights = rep(0,loops)
    for(bb in 1: length(mod.aic)){
      #mod.aic[bb] <-foo[[bb]]$score$aic
      mod.r2[bb] <- 1 - sum((foo[[bb]]$pred$fit - (foo[[bb]]$score$data$Y.test))^2) /
        sum((foo[[bb]]$score$data$Y.test - mean(foo[[bb]]$score$data$Y.test))^2)
    }
    mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 15)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 15)])
    #mask <- cbind(which(mod.aic %iwen% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])
    #mask <- mask[order(mask[,2]),]
    plsglm <- foo[mask[,1]]
    save(plsglm, file = paste("plsglm_", j,  sep = ""))
  }else{
    load(file=paste("./plsglm_", j,  sep = ""))

  }
  mod.r2=rep(0,length(plsglm))
  mod.aic=rep(0,length(plsglm))

  for(bb in 1: length(plsglm)){
    #mod.aic[bb] <-plsglm[[bb]]$score$aic #plsglm[[bb]]$aic
    mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
      sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
  }
  mod.aic
  mod.r2

  if(method=='aic'){
    delta.aic <- mod.aic - min(mod.aic)
    weights <- softmax(-0.5*delta.aic)
  }else if(method =='press'){
    delta.aic <- (max(mod.r2)-mod.r2)
    weights <- softmax(-0.5*delta.aic)

  }
  #
  #weights <- softmax(-0.5*delta.aic)
  #weights[mask[,1]] <- softmax(0.5*delta.press)
  weights
  #weights <- softmax(0.5*delta.press)

  #out of bag
  test.data.y <- read.csv(paste(in.dir, "Spectra/CrownTraits_outBag.csv", sep=""))
  test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
  if(LOG){test.data.y[names(test.data.y)==j] <- log(test.data.y[names(test.data.y)==j])}
  test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag_new.csv", sep=""))

  ndvi <- (test.data.x$band_90- test.data.x$band_58)/(test.data.x$band_58 + test.data.x$band_90)
  nir860 <- (test.data.x$band_96 + test.data.x$band_97)/2
  test.data.x[which(ndvi < 0.7 | nir860 < .3),]=NA
  #allData[which(nir860 < .3),] = NA
  test.data.x <- test.data.x[complete.cases(test.data.x), ]
  pName<- test.data.x$pixel_crownID

  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  test.data.x[,bad_Bands]=NA

  test.data.x <- test.data.x[grepl("band", names(test.data.x))]

  test.data.x[test.data.x>1]=NA
  #test.data.x$pixel_crownID <- pName
  specMat=as.matrix(test.data.x[,-c(1:2)])

  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  test.data.x[,-(1:2)] <- normMat

  crownID <- as.data.frame(pName)
  spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
  if(normz){
    spectra_log_dif_snv[spectra_log_dif_snv==0] <- 0.0000001
    if(NeonSite =="ALL"){
      spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv[,-(1:2)]))),differences=1, lag=3))
      spectra_log_dif_snv<- cbind(test.data.x[,1:2], spectra_log_dif_snv)
    }else{
      spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv))),differences=1, lag=3))
    }
  }
  #test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
  test.PLS = (as.matrix(spectra_log_dif_snv))

  rm(output)
  rm(output.sd)
  output.daic =  output.sd.daic <- crownID
  output.daic$yhat = output.sd.daic$yhat = 0
  w.daic <- rep(0, length(weights))
  w.daic <- weights
  out <- list()
  pred.val.data <- list()
  for(bb in 1:length(w.daic)){
    pls.mod.train <- plsglm[[bb]]$mod
    optim.ncomps <- plsglm[[bb]]$ncomp
    if(glm){
      pred.val.data <- predict(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, type='response',  se.fit = T)
      output.sd.daic$yhat <-  output.sd.daic$yhat + as.vector(pred.val.data$se.fit) * weights[bb]
    }else{
      pred.val.data$fit <- predict(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, type='response', se.fit = T)
    }
    output.daic$yhat <- output.daic$yhat + pred.val.data$fit * w.daic[bb]
    #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
    if(!exists("output")){
      output <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
      if(glm){
        output.sd <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))
      }
    }else{
      output <- rbind(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))))
      if(glm){
        output.sd <- rbind(output.sd, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))))
      }
    }
  }
  colnames(output) <- c("pixel_crownID", "modelID", "yhat")
  output.daic <- as.data.frame(as.matrix(output.daic))
  colnames(output.daic) <- c("pixel_crownID", "yhat")
  if(glm){
    colnames(output.sd) <- c("pixel_crownID", "modelID", "yhat_unc")
    colnames(output.sd.daic) <- c("pixel_crownID", "yhat_unc")

    output.sd <- output.sd[!is.infinite(output.sd$yhat), ]
    output.sd <- output.sd[!is.nan(output.sd$yhat), ]
    output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]
  }

  output <- output[!is.infinite(output$yhat), ]
  output <- output[!is.nan(output$yhat), ]
  output.daic <- output.daic[!is.nan(output.daic$yhat), ]


  #output$yhat <- output$yhat * weights
  pixel.mat <- inner_join(output, test.data.y, by = "pixel_crownID")
  pixel.daic.mat <- inner_join(output.daic, test.data.y, by = "pixel_crownID")
  colnames(pixel.daic.mat)[3] <- "y"
  colnames(pixel.mat)[4] <- "y"
  pixel.based.r2 <-  1 - sum((pixel.mat$yhat - (pixel.mat$y))^2) / sum((pixel.mat$y - mean(pixel.mat$y))^2)
  pixel.daic.r2 <-  1 - sum((pixel.daic.mat$yhat - (pixel.daic.mat$y))^2) / sum((pixel.daic.mat$y - mean(pixel.daic.mat$y))^2)

  #pixel.based <- cor(pixel.mat$yhat, pixel.mat$y)^2
  crown.based <- pixel.mat %>%
    group_by(pixel_crownID) %>%
    summarise(yhat = mean(yhat))
  crown.based <- inner_join(crown.based, test.data.y, by = "pixel_crownID")
  colnames(crown.based) <- c("pixel_crownID", "yhat", "y")
  crown.based.r2 <-  1 - sum((crown.based$yhat - (crown.based$y))^2) / sum((crown.based$y - mean(crown.based$y))^2)
  #crown.based.r2 <- cor(crown.based$yhat, crown.based$y)^2

  crown.based.daic <- pixel.daic.mat %>%
    group_by(pixel_crownID) %>%
    summarise(yhat = mean(yhat))
  crown.based.daic <- inner_join(crown.based.daic, test.data.y, by = "pixel_crownID")
  colnames(crown.based.daic) <- c("pixel_crownID", "yhat", "y")
  crown.daic.r2 <-  1 - sum((crown.based.daic$yhat - (crown.based.daic$y))^2) / sum((crown.based.daic$y - mean(crown.based.daic$y))^2)


  cor(crown.based.daic$yhat,crown.based.daic$y)^2

  crown.daic.r2
  crown.based.r2
  pixel.based.r2
  pixel.daic.r2

  print("sd")
  mean(output.sd.daic$yhat_unc, na.rm = T)


  out.sd <- output.sd.daic %>%
    group_by(pixel_crownID) %>%
    summarise(yhat = mean(yhat_unc, na.rm = T))

  mean(out.sd$yhat)

  print("rmse")
  rmse(pixel.daic.mat$yhat, pixel.daic.mat$y)
  rmse(crown.based.daic$yhat, crown.based.daic$y)



  mlim = c(min(min(pixel.daic.mat$yhat), min(pixel.daic.mat$y)))
  Mlim = c(max(min(pixel.daic.mat$yhat), max(pixel.daic.mat$y)))
  plot(pixel.daic.mat$yhat, pixel.daic.mat$y, xlim = c(mlim,Mlim),ylim = c(mlim,Mlim))
  par(new=TRUE)
  plot(crown.based.daic$yhat, crown.based.daic$y,xlim = c(mlim,Mlim),ylim = c(mlim,Mlim), col="red")
  abline(0, 1)

  write.csv(crown.based.daic, paste(out.dir, nm, "_crown_res.csv", sep=""))
  write.csv(pixel.daic.mat, paste(out.dir, nm, "_pix_res.csv", sep=""))
}
