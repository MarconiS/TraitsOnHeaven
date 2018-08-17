#-----pixPerm------------------------------------------------------------------------------------------
pixPerm <- function(in_dir, out_dir, loops){
  
  # Read data
  allData=read_csv(paste(in_dir, "Spectra/CrownPix_norm.csv", sep="/"))  
  unqCrown = unique(allData$pixel_crownID)
  bootDat <- allData[1:length(unqCrown), ]
  bootDat[] <- NA
  for (laps in 1:loops){
    tk = 1
    for(i in unqCrown){
      if(sum(allData$pixel_crownID==i)==1){
        bootDat[tk,] <- allData[which(allData$pixel_crownID==i),]
      }else{
        rraand <- laps *runif(1, 1, 10^6)
        set.seed(rraand) # todo: change laps * odd number?
        bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
      }
      tk = tk +1
    }
    print(laps)
    write_csv(bootDat, paste(in_dir, '/Permutations/onePix1Crown_', laps, '.csv', sep = ''))
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
  #allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  # Find bad bands
  #all_Bands=as.character(allBand$BandName)
  #bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  if(rescale){allData[,-c(1:3)] <- allData[,-c(1:3)] /10000}
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
  
  #remove any reflectance bigger than 1
  pixel_crownID <- allData[,colnames(allData) %in% "pixel_crownID"]
  #allData <- allData[,-1]
  allData <-allData[grepl("band", names(allData))]
  allData[allData>1]=NA
  
  allData <- cbind(pixel_crownID,allData)
  
  # Find unique crowns (only for plotting purposes)
  unqCrown=unique(allData$pixel_crownID)
  # Extract spectra into matrix
  specMat=as.matrix(allData[-c(1:3)])
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
  normMat=data.frame(specMat/normMat)
  colnames(normMat) = colnames(allData[-c(1:3)])
  #normMat=(specMat - mu_train)/t(sd_train)
  
  # Write vector normalized spectra back into dataframe
  #normDF=data.frame(allData
  #normDF[,all_Bands]=normMat
  
  # Write vector normalized spectra to CSV
  write.csv(cbind(allData[c(1:3)],normMat),  "./inputs/Spectra/CrownPix_norm.csv",row.names=FALSE)
  #write.csv(cbind(mu_train, sd_train), "./outputs/Equal_norm_vectors.csv",row.names=FALSE)
  
  return(unqCrown)
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

pls_glm <- function(nm = "N_pct",cores = 2, wd , loops, nrmlz){
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
    source(paste(wd, '/src/src_build.R', sep=""))
    out_dir = paste(wd, '/outputs/', sep="")
    in_dir = paste(wd, '/inputs/', sep="")
    
    cr.Traits <- read_csv(paste(in_dir, "Spectra/CrownTraits.csv",sep=""))
    aug.spectra <- read_csv(paste(in_dir, 'Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    aug.spectra <- merge(cr.Traits, aug.spectra, by = "pixel_crownID")
    tmp_features<- aug.spectra[grepl("band", names(aug.spectra))]
    tmp_features <- tmp_features[-(which(colnames(tmp_features) %in% c("band_1", "band_370")))]
    tmp_variables <- aug.spectra[names(aug.spectra) %in% nm]
    tmp_variables <- (round(tmp_variables,3))
    aug.X <- data.frame(aug.spectra$name, tmp_features, tmp_variables)
    
    # Subset data into cal/val by site
    eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
    train.data <- eval.set$train
    test.data <- eval.set$test
    colnames(train.data)[1] <- "name"
    colnames(test.data)[1] <- "name"
    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    Y <- as.vector(train.data[,names(train.data) %in% nm])
    #K=nrow(Y)
    nsites <- length(unique(aug.spectra$site))
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }
    
    train.PLS<- cv.plsRglm(dataY=(Y),dataX=X, nt=15,NK=1, K=10,modele="pls-glm-family",family=gaussian(), verbose = F)
    out <- list()
    #out["summary"] <- summary(train.PLS)
    press <- unlist(kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))
    
    #chosen to be built
    mod <- plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]),modele="pls-glm-gaussian", verbose = F)
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }
    
    Y.test <- as.vector(test.data[,names(test.data) %in% nm])
    out[["pred"]] <- predict(mod,newdata=X.tst,type="response",comps=as.integer(out["ncomp"]))
    out[["mod"]] <- mod
    out[["score"]] <- gnm(out$pred ~ (Y.test), family = gaussian())
    out["aic"]<- AIC(out$score)
    out
  }
  stopCluster(cl)
  return(mod.out)
}

predict.withsd <- function(object,newdata,comps=object$computed_nt,type=c("link", "response", "terms", "scores", "class", "probs"),se.fit=FALSE, wt = NULL, dispersion = NULL,methodNA="adaptative",verbose=TRUE,...)
{
  if (!inherits(object, "plsRglmmodel"))
    stop("Primary argument much be a plsRglmmodel object")
  if(missing(type)){type="link"}
  if (!(type %in% c("link",  "response", "terms", "scores", "class", "probs")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="link"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="response"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="terms"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="class"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}
    }
    if (type=="probs"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}
    }
    if (type=="scores"){return(object$tt[,1:comps])}
  } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$call$formula)){
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("subset", "weights"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$data <- newdata
      mf$formula <- object$call$formula[-2]
      mf$drop.unused.levels <- TRUE
      mf$na.action <- na.pass
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      mt <- attr(mf, "terms")
      #    attr(mt,"intercept")<-0L
      newdata.frame <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)[,-1]
      else matrix(, nrnd, 0L)
      weights <- as.vector(model.weights(mf))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
    newdataNA <- !is.na(newdata)
    newdata.scaledwotNA <- as.matrix(newdata.scaled)
    newdata.scaledwotNA[!newdataNA] <- 0
    ttpredY <- NULL
    if (methodNA=="adaptative") {
      for(ii in 1:nrnd){
        if (all(newdataNA[ii,])){
          ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps)))
        }
        else {
          if(verbose){cat("Missing value in row ",ii,".\n")}
          ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
        }}}
    if(methodNA=="missingdata") {
      if(verbose){cat("Prediction as if missing values in every row.\n")}
      for (ii in 1:nrnd) {
        ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
      }
    }
    colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY)
    if("glm" %in% class(object$FinalModel)){
      if (type=="link"){return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit,dispersion = dispersion,...))
      }
      #if (type=="response"){return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit,dispersion = dispersion,...))
      #}
      if (type=="terms"){return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit,dispersion = dispersion,...))
      }
      if (type=="response"){return(predict.lm(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit,interval = 'prediction',weights = wt, dispersion = dispersion,...))
      }
    }
    if("polr" %in% class(object$FinalModel)){
      if (type=="class"){return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}
      if (type=="probs"){return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}
    }
    if (type=="scores"){colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="");
    return(ttpredY[,1:comps,drop=FALSE])}
  }
}
