###--- normalize
normalize <- function(CrownIDS, NeonSite, rescale = T, in.dir = NULL, out.dir = NULL){
	# function used to reduce peripheral illumination noise. The idea is to equalize
	# the spectral reflectance of each pixel by 
	# credits to Adytia Singh
	
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  #for space reason, better keep NEON data to be integer: rescale used to get to the real reflectance [0:1]
  if(rescale){allData[,-c(1,2)] <- allData[,-c(1,2)] /10000}
  # Set bad bands to zero
  ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
  nir860 <- (allData$band_96 + allData$band_97)/2 
  allData[which(ndvi < 0.7 | nir860 < .3),]=NA
  allData <- allData[complete.cases(allData), ]
  allData[,bad_Bands]=NA
  #check if there are at least 1+ pixels per crown. If not, remove that crown from the analysis
  pixPerCrown <- allData %>%
    group_by(pixel_crownID) %>%
    summarise(n = n())
  colnames(CrownIDS)[colnames(CrownIDS) %in% NeonSite] <- "pixel_crownID"
  badCrowns <- anti_join(CrownIDS, pixPerCrown, by = "pixel_crownID")$pixel_crownID
  badCrowns <- na.omit(badCrowns)
  #remove any reflectance bigger than 1
  pixel_crownID <- allData[,colnames(allData) %in% "pixel_crownID"]
  #allData <- allData[,-1]
  allData <-allData[grepl("band", names(allData))]
  allData[allData>1]=NA
  
  allData <- cbind(pixel_crownID,allData)
  
  # Find unique crowns (only for plotting purposes)
  unqCrown=unique(allData$pixel_crownID)
  # Extract spectra into matrix
  specMat=as.matrix(allData[,all_Bands])
  # Vector normalize spectra
  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  
  # Write vector normalized spectra back into dataframe
  normDF=allData
  normDF[,all_Bands]=normMat
  
  # Write vector normalized spectra to CSV
  write.csv(normDF, "./inputs/Spectra/CrownPix_norm.csv",row.names=FALSE)
  return(badCrowns)
}



### -- pixPerm

pixPerm <-function(loops, unqCrown, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                    path = NULL){
	# randomply extract 1 pixel per crown in order to develop a single calibration-validation set. The quality of such combination of pixels will be used to weight importance of instances (pixels of a bag) to predict the response variable
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  # Read data
  allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
  allData=imp.spectra("Spectra/fCrownPix_norm.csv",in.dir)  
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  # Find unique crowns (only for plotting purposes)
  bootDat <- allData[1:length(unqCrown), ]
  bootPix <- allData[1:length(unqCrown), 1:2]
  bootDat[] <- NA
  bootPix[] <- NA  
  for (laps in 1:loops){
    tk = 1
    for(i in unqCrown){
      if(sum(allData$pixel_crownID==i)==1){
        bootDat[tk,] <- allData[which(allData$pixel_crownID==i),]
        bootPix[tk,] <- c(i,1 + which(allData$pixel_crownID==i))
      }else{
        set.seed(laps) # todo: change laps * odd number?
        bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
        set.seed(laps)
        bootPix[tk,] <- c(i,1 + sample(which(allData$pixel_crownID==i), 1))
      }
      tk = tk +1
      
    }
    print(laps)
    names(bootPix) <- c("pixel_crownID", "ChosenPix")
    write.csv(bootDat, paste('inputs/Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    write.csv(bootPix, paste('inputs/Permutations/onePix1Position_', laps, '.csv', sep = ''))
  }
}

### --- pls_glm
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
    out
  }
  stopCluster(cl)
  return(mod.out)
}



### --- 
pPLS <- function(j = NULL, loops = 1000,nrmlz = F,numcomps=20, cores =2,nsites = 2,  out.dir = NULL, in.dir = NULL){
  library(foreach)
  library(doParallel)
  
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  ll = as.character(1:loops)
  mod.out <- foreach(laps = ll, .verbose = T, .multicombine = TRUE) %dopar% {
    
    library('pls')
    library(readr)
    library(gnm)
    
    source(paste('../../src/functions_build.R', sep=""))
    source(paste('../../src/src_build.R', sep=""))
    out.dir = paste(getwd(), "/outputs/", sep="")
    in.dir = paste(getwd(), "/inputs/", sep="")
    
    out <- list()
    
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
    
    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    Y <- as.vector(train.data[,names(train.data) %in% j])
    #K=nrow(Y)
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }
    train.PLS = data.frame(Y = I(Y), X=I(X))
    pls.mod.train = plsr(Y ~ X,scale=F, ncomp=numcomps,validation="LOO", 
                         trace=TRUE, method = "oscorespls", probMethod = "softmax", 
                         data = train.PLS)     
    
    
    optim.ncomps <-which(pls.mod.train$validation$PRESS==min(pls.mod.train$validation$PRESS[3:length(pls.mod.train$validation$PRESS)]))
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }
    
    Y.test <- as.vector(test.data[,names(test.data) %in% j])
    test.PLS <- data.frame(X=I(X.tst))
    pred.val.data <- as.vector(predict(pls.mod.train, newdata = test.PLS, 
                                       ncomp=optim.ncomps, type="response"))
    
    out[["pred"]]$fit <- pred.val.data 
    out[["mod"]] <- pls.mod.train
    out[["ncomp"]] <- optim.ncomps
    out[["score"]]$data$Y.test <- Y.test
    out[["score"]]$score$aic <- 1 - sum((pred.val.data - (Y.test))^2) /
      sum((Y.test - mean(Y.test))^2)
    out
  } 
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

