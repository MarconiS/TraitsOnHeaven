PLSandLasso <- function( rounds = 10, loops = 400, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                         out.dir = paste(getwd(), "/outputs/", sep=""), in.dir = paste(getwd(), "/inputs/", sep="")){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  
  for (j in 1:6) {
    R = matrix(NA, ncol=4, nrow=loops)
    
    for(laps in 1:loops) {
      aug.spectra <- imp.spectra(paste('Bootstrap_',rounds,'/onePix1Crown_',names[j], laps, '.csv', sep = ''), in.dir)
      aug.spectra$X=aug.spectra$X.1 = aug.spectra$X.2 = aug.spectra$X.3 = aug.spectra$X.4 = aug.spectra$X.5 = 
        aug.spectra$X.6 = aug.spectra$X.7 = aug.spectra$X.8 = aug.spectra$X.9 = aug.spectra$X.10 = NULL
      
      aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
      X<- aug.spectra[grepl("band", names(aug.spectra))]
      
      
      X=X[, colSums(is.na(X)) == 0]
      # filter <- apply(X, 2, max)
      # if(any(filter > 1)){
      #   X <- X[-which(filter > 1)]
      # }
      # X=X[, colSums(is.na(X)) == 0]
      Y <- aug.spectra[,names(aug.spectra) %in% names]
      X_corr = cor(as.matrix(X), Y)
      aug.X <- data.frame(aug.spectra$name, Y, X)
      # Subset data into cal/val by site
      eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
      train.data <- eval.set$train
      test.data <- eval.set$test
      
      print(paste(laps, names(Y[j])))
      # Run calibration PLSR analysis to select optimal number of components
      pls.mod.train <- pls.cal(train.data, 25,nm = names, j, norm = F)
      #calculate number of components given min test PRESS or RMSEP
      optim.ncomps <- opt.comps(pls.mod.train, Y, j)
      pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,optim.ncomps,j, norm = F)
      out.data <- res.out(pred.val.data, train.data,nm = names, test.data, j)
      
      #-- run lasso -------------------------------------------------------------------------------------#
      lasso.reg <- lasso.asd(test.data,train.data, j, nm = names,  norm = F)
      print.coeffs(pls.mod.train,optim.ncomps, rounds, lasso.reg, Y, j, out.dir) 
      
      
      # #----run GP ------#
      # library("GPfit")
      # GPmodel = GP_fit(train.data[,8:368], train.data[,j+1] )
      # 
      setwd(out.dir)
      
      #I apology for the lack of organization in the plotting scheme
      R[laps,1:2] <- c(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$r2",  sep = ""))),
                       mean(eval(parse(text = paste("out.data$", names(Y[j]), "$R2", sep = "")))))
      R[laps,3:4] <- c(sqrt(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$mse",  sep = "")))),
                       sqrt(mean(eval(parse(text = paste("out.data$", names(Y[j]), "$MSE", sep = ""))))))
    } 
    assign(paste("R.", names(Y[j]),  sep = ""), R)
  }
  
  R2est <- as.data.frame(cbind(R.LMA_g.m2[,1:2],R.d13C[,1:2], R.d15N[,1:2], R.C_pct[,1:2], R.N_pct[,1:2],R.P_pct[,1:2]))
  colnames(R2est) <- c("LMA Lasso", "LMA PLS","δ13C Lasso","δ13C PLS","δ15N Lasso", "δ15N PLS","C Lasso", "C PLS","N Lasso","N PLS", "P Lasso", ".P_pct")
  
  MSEest <- as.data.frame(cbind(R.LMA_g.m2[,3:4], R.d13C[,3:4], R.d15N[,3:4],  R.C_pct[,3:4], R.N_pct[,3:4],R.P_pct[,3:4]))
  colnames(MSEest) <- c("LMA_g.m2", ".LMA_g.m2","d13C",".d13C","d15N", ".d15N","C_pct", ".C_pct","N_pct",".N_pct", "P_pct", ".P_pct")
  write.csv(R2est, file = paste(names[j],'Boost', rounds,'R2.csv',sep="_"))
  write.csv(MSEest, file = paste(names[j],'Boost', rounds,'MSE.csv',sep="_"))
}


normalize<-function(){
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  #allData[,-c(1,2)] <- allData[,-c(1,2)] /10000
  # Set bad bands to zero
  ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
  allData[which(ndvi < 0.7),]=NA
  allData <- allData[complete.cases(allData), ]
  allData[,bad_Bands]=NA
  
  #remove any reflectance bigger than 1
  pixel_crownID <- allData[,2]
  allData <- allData[,-c(1,2)]
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
}



rPix <-function(rounds, loops, unqCrown, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                path = ("/Users/sergiomarconi/Projects/OSBS")){
  
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  # Read data
  allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
  
  for (j in 1:6){
    if(rounds > 1){
      allData=read.csv(paste(in.dir, "Spectra", "/CrownPix_norm_",rounds-1,names[j], ".csv", sep=""))
      all_Bands=as.character(allBand$BandName)
      bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
    }else{
      allData=imp.spectra("Spectra/CrownPix_norm.csv",in.dir)  
      all_Bands=as.character(allBand$BandName)
      bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
    }
    # Set bad bands to zero
    allData[,bad_Bands]=NA
    #allData=allData[, colSums(is.na(allData)) == 0]
    
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
      write.csv(bootDat, paste('inputs/Bootstrap_', rounds,'/onePix1Crown_',names[j], laps, '.csv', sep = ''))
      write.csv(bootPix, paste('inputs/Bootstrap_',rounds,'/onePix1Position_',names[j], laps, '.csv', sep = ''))
    }
  }
}



DiscardAndRerun <- function(names, rounds, nsim, nCrowns, lass.included = F){
  library(readr)
  for (j in 1:6){
    pls_coef <- read_csv(paste(out.dir,'pls_coeffs_',rounds,names[j],'.csv',sep=""),col_names = FALSE)
    if(lass.included = T){
      lasso_coef <- read_csv(paste(out.dir,'las_coeffs_',rounds,names[j],'.csv',sep=""), col_names = FALSE)
    }
    # pls_coef <- read.csv(paste(out.dir,'pls_coeffs_',names[j],'.csv',sep=""),header = FALSE)
    # lasso_coef <- read.csv(paste(out.dir,'las_coeffs_',names[j],'.csv',sep=""), header = FALSE)
    nEntries = 100
    press <- data.frame(cbind(seq(1,nsim),pls_coef$X1))
    if(lass.included = T){
      r2.las <-data.frame(cbind(seq(1,nsim),lasso_coef$X1))
      p.ls <- r2.las[order(r2.las$X2),]
      p.ls <- p.ls[1:100,1]
      nEntries = 200
    }
    p <- press[order(-press$X2),]
    head(p)
    p <- p[1:100,1]
    
    pixels <- data.frame(matrix(NA, ncol=nCrowns, nrow = nEntries))
    tk = 0
    for(k in p){
      tk = tk + 1
      #import ith permutation associated with the horrible correlation
      tmp <- imp.spectra(paste('Bootstrap_',rounds,'/onePix1Position_',names[j], k, '.csv', sep = ''), in.dir)
      pixels[tk,] <- tmp$ChosenPix 
    }
    if(lass.included = T){
      for(k in p.ls){
        tk = tk + 1
        #import ith permutation associated with the horrible correlation
        #check why 1
        tmp <- imp.spectra(paste('Bootstrap_',rounds,'/onePix1Position_',names[j], k, '.csv', sep = ''), in.dir)
        #tmp <- imp.spectra(paste('Bootstrap_1/onePix1Position_', names[j], k, '.csv', sep = ''), in.dir)
        pixels[tk,] <- tmp$ChosenPix #it is called band_1 just because I didn't change the name when created 1pix1crown, but is the actual pixel  chosen from the ith crown
      }
    }
    discarded = kept = NULL
    for(k in 1:nCrowns){
      if(dim(table(pixels[,k])) >4){
        rank.p <- as.data.frame(table(pixels[,k]))
        #order frequency of bad pixels. Cut off pixels summing up to 85%?
        rank.p <- rank.p[order(-rank.p$Freq),]
        tk = qtl = 0
        while(qtl < 40){
          tk <- tk +1
          qtl <- rank.p$Freq[tk] + qtl
        }
        discarded <- c(discarded, as.character(rank.p$Var1[1:tk-1]))
        kept <- c(kept, as.character(rank.p$Var1[tk:dim(rank.p)[1]]))
      }
    }
    write.csv(discarded, paste(in.dir, 'Bootstrap_2/badPix_', names[j],'.csv', sep = ''))
    write.csv(kept, paste(in.dir, 'Bootstrap_2/goodPix_', names[j],'.csv', sep = ''))
    
    allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
    allData=imp.spectra("Spectra/CrownPix_norm.csv",in.dir)  
    # Find bad bands
    all_Bands=as.character(allBand$BandName)
    bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
    
    # Subset bad bands to zero > here you want to delete the bad pixel from CrownPix.csv [now maintaning the good one]
    allData <- allData[-as.integer(discarded),]
    write.csv(allData, paste(in.dir, 'Spectra/CrownPix_norm_',rounds, names[j],'.csv', sep = ''))
  }
}


## GPr
GProute<- function(train, test, name, kernel){
}
