#-----pixPerm------------------------------------------------------------------------------------------
pixPerm <- function(rounds, loops, unqCrown, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                    path = ("/Users/sergiomarconi/Projects/OSBS")){
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
      set.seed(laps + (rounds -1) * loops) # todo: change laps * odd number?
      bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
      set.seed(laps + (rounds -1) * loops)
      bootPix[tk,] <- c(i,1 + sample(which(allData$pixel_crownID==i), 1))
      tk = tk +1
    }
    names(bootPix) <- c("pixel_crownID", "ChosenPix")
    write.csv(bootDat, paste('inputs/Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    write.csv(bootPix, paste('inputs/Permutations/onePix1Position_', laps, '.csv', sep = ''))
  }
}

#-----PLS------------------------------------------------------------------------------------------
PLS <- function(names = NULL, loops = 1000, out.dir = paste(getwd(), "/outputs/", sep=""), in.dir = paste(getwd(), "/inputs/", sep="")){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_1.csv', sep = ''), in.dir)
  for (j in 1:length(names)) {
    print(j)
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
      
      print(paste(laps, names(Y[j])))
      # Run calibration PLSR analysis to select optimal number of components
      pls.mod.train <- pls.cal(train.data, 15,nm = names, j, norm = F)
      #calculate number of components given min test PRESS or RMSEP
      
      optim.ncomps <- opt.comps(pls.mod.train, Y, j)
      
      pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,optim.ncomps,j, norm = F)
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
#-----PLS_DA------------------------------------------------------------------------------------------
PLS_DA <- function(loops = 1000, names = c("name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                   out.dir = paste(getwd(), "/outputs/", sep=""), in.dir = paste(getwd(), "/inputs/", sep="")){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/crownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  names = "name"
  aug.spectra <- imp.spectra(paste('Spectra/CrownPix_norm.csv', sep = ''), in.dir)
  
  #matrix to store performances in
  mod.out = vector("list", loops)
  mod.stats = vector("list", loops)
  mod.comps = rep(NA, loops)
  aug.spectra$X= NULL
  aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
  X<- aug.spectra[grepl("band", names(aug.spectra))]
  X=X[, colSums(is.na(X)) == 0]
  Y <- aug.spectra[,names(aug.spectra) %in% names]
  names(Y) <- "name"
  aug.X <- data.frame(aug.spectra$name, Y, X)
  # Subset data into cal/val by site
  prova <- data.frame(unique(cbind(aug.spectra$pixel_crownID, aug.spectra$name)))
  colnames(prova) <- c("pixel_crownID", "species")
  proportions <- 0.7
  set.seed(14)
  out <- prova %>% 
    group_by(species) %>%
    filter(pixel_crownID %in% sample(pixel_crownID, ceiling(proportions*length(pixel_crownID))))
  
  pixID <- data.frame(as.integer(levels(out$pixel_crownID)[out$pixel_crownID]))
  names(pixID) <- "pixel_crownID"
  train.data <- inner_join(pixID, aug.spectra, by = "pixel_crownID")
  test.data <- anti_join(aug.spectra, by = "pixel_crownID", pixID)
  
  print(paste(laps, names(Y[j])))
  # Run calibration PLSR analysis to select optimal number of components
  pls.mod.train <- pls.cal(train.data, 15, nm = names, j=1, norm = F)
  #calculate number of components given min test PRESS or RMSEP
  
  optim.ncomps <- opt.comps(pls.mod.train, Y, j)
  
  pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,optim.ncomps,j, norm = F)
  out.data <- res.out(pred.val.data, train.data,nm = names, test.data, j)
  
  mod.out[[laps]] <- pls.mod.train 
  mod.stats[[laps]] <- out.data
  mod.comps[laps] <- optim.ncomps
  setwd(out.dir)
  save(mod.out, file = paste("models_out_", names[j],  sep = ""))
  save(mod.stats,  file = paste("models_stats_", names[j],  sep = ""))
  save(mod.comps, file = paste("models_comps_", names[j],  sep = ""))
  
}
#-----normalize------------------------------------------------------------------------------------------
normalize<-function(CrownIDS){
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  #allData[,-c(1,2)] <- allData[,-c(1,2)] /10000
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
  allData <- allData[,-1]
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
#-----perform_summary------------------------------------------------------------------------------------------
perform_summary <- function(names=NULL,out.name = NULL){
  
  for(j in names){
    load(file = paste(out.dir, 'models_comps_',j, sep="" ))
    load(file = paste(out.dir, 'models_out_',j, sep="" ))
    load(file = paste(out.dir, 'models_stats_',j, sep="" ))
    mod.r2 <- rep(NA, length(mod.stats))
    
    for(bb in 1: length(mod.r2)){
      mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
    }
    mask <- which(mod.r2 %in% tail(sort(mod.r2), 100) & mod.r2 > 0.0)
    pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
    if(j == "name"){
      norm.R2 <- scalar1(mod.r2[mask])
      multiplier <- 10^4 #decimalplaces(min(norm.R2))
      weights <- floor(norm.R2 * multiplier)
    }  
    predictions <- sapply(mod.stats, "[[", 2)
    #rm(pred.weighted, pred.weight, pred)
    if(j == "name"){
      tkn <- 0
      for(bb in (mask)){
        tkn <- tkn + 1
        if(exists("pred.weighted")){
          pred.weighted <- cbind(pred.weighted, matrix(predictions[, bb], nrow(predictions), as.integer(weights[tkn])))
        }else{
          pred.weighted <- matrix(predictions[, bb], nrow(predictions), as.integer(weights[tkn]))
        }
      }
    }
    for(bb in mask){
      if(j == "name"){
        pred=rep(NA, length(mod.stats[[1]]$pred))
        pred.weight=rep(NA, length(mod.stats[[1]]$pred))
        
        for(ii in 1: length(mod.stats[[1]]$pred)) {
          temp.freq <- table(predictions[ii,])
          temp.freq.weight <- table(pred.weighted[ii,])
          pred[ii] <- names(temp.freq)[which(temp.freq == max(temp.freq))]
          pred.weight[ii] <- names(temp.freq.weight)[which(temp.freq.weight == max(temp.freq.weight))]
        }
      }else{
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
      correct <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred)
      accuracy <- sum(as.numeric(correct))/length(correct)
      correct.weigth <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred.weight)
      accuracy.weight <- sum(as.numeric(correct.weigth))/length(correct.weigth)
      if(!exists("out")){
        out <- c(j, accuracy, accuracy.weight)
      }else{
        out <- rbind(out, c(j, accuracy, accuracy.weight))
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
  write_csv(data.frame(out), paste(out.dir, out.name, sep=""))
}
