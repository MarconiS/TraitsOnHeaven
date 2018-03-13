#-----pixPerm------------------------------------------------------------------------------------------
pixPerm <- function(loops, unqCrown, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                    path = NULL){
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  # Read data
  allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
  allData=imp.spectra("Spectra/CrownPix_norm.csv",in.dir)  
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
        rraand <- laps *runif(1, 1, 10^6)
        set.seed(rraand) # todo: change laps * odd number?
        bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
        set.seed(rraand)
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

#-----PLS------------------------------------------------------------------------------------------
PLS <- function(names = NULL, loops = 1000,norm = F, out.dir = NULL, in.dir = NULL){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_1.csv', sep = ''), in.dir)
  for (j in 1:length(names)) {
    print(j)
    if(names[j] %in% c("C_pct", "P_pct")){
      norm = T
    }else{
      norm = F
    }
    #matrix to store performances in
    mod.out = vector("list", loops)
    mod.stats = vector("list", loops)
    mod.comps = rep(NA, loops)
    for(laps in 1:loops) {
      aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_', laps, '.csv', sep = ''), in.dir)
      aug.spectra$X= NULL
      aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
      X<- aug.spectra[grepl("band", names(aug.spectra))]
      X=X[, colSums(is.na(X)) == 0]
      Y <- aug.spectra[,names(aug.spectra) %in% names]
      if(is.null(dim(Y))){names(Y) = names}
      aug.X <- data.frame(aug.spectra$name, Y, X)
      # Subset data into cal/val by site
      eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
      train.data <- eval.set$train
      test.data <- eval.set$test
      colnames(train.data)[1] <- "name"
      colnames(test.data)[1] <- "name"
      write_csv(train.data, "train_for_baseline.csv")
      write_csv(test.data, "test_for_baseline.csv")
      
      print(paste(laps, names(Y[j])))
      # Run calibration PLSR analysis to select optimal number of components
      pls.mod.train <- pls.cal(train.data, 15,nm = names, j, norm = norm)
      #calculate number of components given min test PRESS or RMSEP
      
      optim.ncomps <- opt.comps(pls.mod.train, Y, j, Class = F)
      
      pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,optim.ncomps,j, norm = norm)
      out.data <- res.out(pred.val.data, train.data,nm = names, test.data, j)
      
      mod.out[[laps]] <- pls.mod.train 
      mod.stats[[laps]] <- out.data
      mod.comps[laps] <- optim.ncomps
      setwd(out.dir)
    } 
    save(mod.out, file = paste("models_out_", names[j],  sep = ""))
    save(mod.stats,  file = paste("models_stats_", names[j],  sep = ""))
    save(mod.comps, file = paste("models_comps_", names[j],  sep = ""))
  }
}

#-----pPls------------------------------------------------------------------------------------------
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
    nsites = 2
    print(j)
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
  stopCluster(cl)
  return(mod.out)
}

#-----normalize------------------------------------------------------------------------------------------
normalize<-function(CrownIDS, NeonSite, rescale = T, in.dir = NULL, out.dir = NULL){
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  if(rescale){allData[,-c(1,2)] <- allData[,-c(1,2)] /10000}
  # Set bad bands to zero
  ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
  nir860 <- (allData$band_96 + allData$band_97)/2 
  allData[which(ndvi < 0.7 | nir860 < .3),]=NA
  #allData[which(nir860 < .3),] = NA
  allData <- allData[complete.cases(allData), ]
  allData[,bad_Bands]=NA
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
  #equalizer=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  #mu_train = apply(specMat, FUN =mean, MARGIN = 2, na.rm=T)
  #sd_train =  apply(specMat, FUN =sd, MARGIN = 2, na.rm=T)
  #normMat=matrix(data=rep(equalizer,ncol(specMat)),ncol=ncol(specMat))
  #for(i in 1:length(sd_train)){
  #  normMat[,i]= (specMat[,i] - mu_train[i])/(sd_train[i])
  #}
  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  #normMat=(specMat - mu_train)/t(sd_train)
  
  # Write vector normalized spectra back into dataframe
  normDF=allData
  normDF[,all_Bands]=normMat
  
  # Write vector normalized spectra to CSV
  write.csv(normDF, "./inputs/Spectra/CrownPix_norm.csv",row.names=FALSE)
  #write.csv(cbind(mu_train, sd_train), "./outputs/Equal_norm_vectors.csv",row.names=FALSE)
  
  return(badCrowns)
}
#-----perform_summary------------------------------------------------------------------------------------------
perform_summary <- function(names=NULL,out.name = NULL, in.dir = NULL, out.dir = NULL, out.of.bag = T, normlz = F, weighted = T){
  for(j in names){
    load(file = paste(out.dir, 'models_comps_',j, sep="" ))
    load(file = paste(out.dir, 'models_out_',j, sep="" ))
    load(file = paste(out.dir, 'models_stats_',j, sep="" ))
    mod.r2 <- rep(NA, length(mod.stats))
    if(j == "name"){
      for(bb in 1: length(mod.r2)){
        mod.r2[bb] <- mean(mod.stats[[bb]]$name$R2)
      }
      mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 100)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 100)])
      mask <- mask[order(mask[,2]),]
      pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
      pred.weighted <- rep(0, length(mod.stats[[1]]$name$pred))
      norm.R2 <- scalar1(mask[,2])
      multiplier <- 10^2 #decimalplaces(min(norm.R2))
      weights <- floor(norm.R2 * multiplier)
      
      test.data <- read_csv("test_classification.csv")
      test.data <- test.data[colnames(test.data) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]
      
      lab.sp <- read.csv("spLabels.csv", stringsAsFactors = F)
      
    }else{
      for(bb in 1: length(mod.r2)){
        mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
      }
      mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 100)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 100)])
      mask <- mask[order(mask[,2]),]
      pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
      if(out.of.bag){
        test.data.y <- read_csv(paste(in.dir, "Spectra/CrownTraits_outBag.csv", sep=""))
        test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
        test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag.csv", sep=""))
        crownID <- as.data.frame(test.data.x$pixel_crownID)
        test.data.x <- test.data.x[grepl("band", names(test.data.x))]
        spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
        if(normlz){  spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv))),differences=1, lag=3))}    
        
      }else{
        predictions <- sapply(mod.stats, "[[", 2)
      }
    }
    
    if(j == "name"){
      rm(pred.weighted, pred.weight, pred, predictions)
      tkn <- 0
      for(bb in (mask[,1])){
        tkn <- tkn + 1
        if(exists("pred.weighted")){
          predictions <- cbind(predictions, matrix(mod.stats[[bb]]$name$pred, nrow(test.data)))
          pred.weighted <- cbind(pred.weighted, matrix(mod.stats[[bb]]$name$pred, nrow(test.data), as.integer(mask[tkn,2])))
        }else{
          predictions <- matrix(mod.stats[[bb]]$name$pred, nrow(test.data))
          pred.weighted <-  matrix(mod.stats[[bb]]$name$pred, nrow(test.data), as.integer(mask[tkn,2]))
        }
      }
      pred=rep(NA, length(mod.stats[[1]]$name$pred))
      pred.weight=rep(NA, length(mod.stats[[1]]$name$pred))
      for(ii in 1: length(mod.stats[[bb]]$name$pred)) {
        temp.freq <- table(predictions[ii,])
        temp.freq.weight <- table(pred.weighted[ii,])
        pred[ii] <- names(temp.freq)[which(temp.freq == max(temp.freq))]
        pred.weight[ii] <- names(temp.freq.weight)[which(temp.freq.weight == max(temp.freq.weight))]
      }
      pred <- as.data.frame(cbind(test.data$pixel_crownID, pred))
      pred.weight <- as.data.frame(cbind(test.data$pixel_crownID, pred.weight))
      colnames(pred) = c("pixel_crownID","spID")
      colnames(pred.weight) = c("pixel_crownID","spID_w")
    }else{
      if(out.of.bag){
        test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
        rm(output)
        for(jj in mask[,1]){
          pls.mod.train <- mod.out[[jj]]
          optim.ncomps <- mod.comps[jj]
          pred.val.data <- predict(pls.mod.train[[j]], newdata = test.PLS, ncomp=optim.ncomps, type='response')
          
          #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
          if(!exists("output")){
            output <- cbind(crownID, rep(jj, dim(test.data.x)[1]), as.vector(pred.val.data))
          }else{
            output <- rbind(output, cbind(crownID, rep(jj, dim(test.data.x)[1]), as.vector(pred.val.data)))
          }
        }
      }else{
        #????????????????
        pred=rep(NA, length(mod.stats[[1]]$pred))
        pred.weight=rep(NA, length(mod.stats[[1]]$pred))
        if(!exists("pred.weighted")){
          pred.weighted <- mod.stats[[bb]]$pred * mod.r2[bb]
          pred <- mod.stats[[bb]]$pred
        }else{
          pred.weighted <- pred.weighted + mod.stats[[bb]]$pred * mod.r2[bb]
          pred <- pred + mod.stats[[bb]]$pred
        }
      }
    }
    
    if(j == "name"){
      rm(evalFinal)
      evalFinal <- merge(lab.sp, test.data, by = "name")
      evalFinal <-unique(evalFinal)
      evalFinal <- merge(evalFinal, pred.weight, by = "pixel_crownID")
      evalFinal <- merge(evalFinal, pred, by = "pixel_crownID")
      correct <- evalFinal$spID_train == evalFinal$spID
      accuracy <- sum(as.numeric(correct))/length(correct)
      correct.weigth <- evalFinal$spID_train == evalFinal$spID_w
      accuracy.weight <- sum(as.numeric(correct.weigth))/length(correct.weigth)
      print(accuracy)
      print(accuracy.weight)
      
      if(!exists("out")){
        out <- c(j, accuracy, accuracy.weight)
      }else{
        out <- rbind(out, c(j, accuracy, accuracy.weight))
      }
    }else{
      if(out.of.bag){
        colnames(output) <- c("pixel_crownID", "modelID", "yhat")
        if(weighted){
          foo <- mask[,2]/sum(mask[,2])
          weights <- NULL
          for(temp in 1:length(foo)){
            weights<- c(weights, rep(foo[temp], dim(test.data.x)[1]))
          }
        }else{
          weights <- rep(1, length(output$yhat))
        }
        #output$yhat <- output$yhat * weights
        pixel.mat <- inner_join(output, test.data.y, by = "pixel_crownID")
        pixel.based <- cor(pixel.mat$yhat, eval(parse(text = paste('pixel.mat$',j,sep=""))))^2
        crown.based <- pixel.mat %>%
          group_by(pixel_crownID) %>%
          summarise(yhat = mean(yhat))
        crown.based <- inner_join(crown.based, test.data.y, by = "pixel_crownID")
        colnames(crown.based) <- c("pixel_crownID", "yhat", "y")
        #crown.based <- cor(crown.based$yhat, crown.based$y)^2
        #plot(lm(crown.based$yhat ~ crown.based$y))
        par(mfrow=c(2,1))
        plot(crown.based$yhat ~ crown.based$y, main = "1:1 crown object leaf N (%)", cex = 2, xlab = "OOB observed N (%)",ylab = "OOB predicted N (%)", xlim = c(0.5, 2.5), ylim = c(0.5, 2.5), pch=21, bg="blue")
        abline(0,1, col="red")
        plot(pixel.mat$yhat ~ eval(parse(text = paste('pixel.mat$',j,sep=""))), main = "1:1 pixel leaf N (%)", cex = 2, xlab = "OOB observed N (%)",ylab = "OOB predicted N (%)", xlim = c(0.5, 2.5), ylim = c(0.5, 2.5), pch=21, bg="blue")
        abline(0,1, col="red")
        print(paste(pixel.based, "|", crown.based))
        if(!exists("out")){
          out <- c(j, pixel.based, crown.based)
        }else{
          out <- rbind(out, c(j, pixel.based, crown.based))
        }
      }else{
        predictions <- predictions / length(mod.stats[[1]]$pred)
        pred.weighted <- pred.weighted /sum(mod.r2[mask])
        if(!exists("out")){
          out <- c(j, cor(pred, mod.stats[[1]]$obs)^2, cor(pred.weighted, mod.stats[[1]]$obs)^2)
        }else{
          out <- rbind(out, c(j, cor(pred, mod.stats[[1]]$obs)^2, cor(pred.weighted, mod.stats[[1]]$obs)^2))
        }
      }
    }
    write_csv(data.frame(pred.weighted), paste(out.dir, "predicted_",j, ".csv", sep=""))
  }     
  write_csv(data.frame(out), paste(out.dir, out.name, sep=""))
  return(out)
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

plot_upper_cor <- function(cormat){
  
  #this is when uncertainty associated to species is 0
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
}
