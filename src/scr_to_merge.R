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


PLS <- function( rounds = 10, loops = 1000, names = NULL,
                         out.dir = paste(getwd(), "/outputs/", sep=""), in.dir = paste(getwd(), "/inputs/", sep="")){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_1.csv', sep = ''), in.dir)
  for (j in 1) {
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
      
      #I apology for the lack of organization in the plotting scheme
      # R[laps,1:2] <- c(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$r2",  sep = ""))),
      #                  mean(eval(parse(text = paste("out.data$", names(Y[j]), "$R2", sep = "")))))
      # R[laps,3:4] <- c(sqrt(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$mse",  sep = "")))),
      #                  sqrt(mean(eval(parse(text = paste("out.data$", names(Y[j]), "$MSE", sep = ""))))))
    } 
    save(mod.out, file = paste("models_out_", names[j],  sep = ""))
    save(mod.stats,  file = paste("models_stats_", names[j],  sep = ""))
    save(mod.comps, file = paste("models_comps_", names[j],  sep = ""))
  }
}
