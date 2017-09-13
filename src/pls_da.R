PLS_DA <- function( rounds = 10, loops = 1000, names = c("name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
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
  
  out <- prova %>% 
    group_by(species) %>%
    filter(pixel_crownID %in% sample(pixel_crownID, ceiling(0.7*length(pixel_crownID))))
  
  pixID <- data.frame(as.integer(levels(out$pixel_crownID)[out$pixel_crownID]))
  names(pixID) <- "pixel_crownID"
  train.data <- inner_join(pixID, aug.spectra, by = "pixel_crownID")
  test.data <- anti_join(aug.spectra, by = "pixel_crownID", pixID)
  
  #  write_csv(train.data, paste(in.dir, "Spectra/train_species.csv", sep =""))
  #  write_csv(test.data, paste(in.dir, "Spectra/test_species.csv", sep =""))
  
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
  
  #I apology for the lack of organization in the plotting scheme
  # R[laps,1:2] <- c(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$r2",  sep = ""))),
  #                  mean(eval(parse(text = paste("out.data$", names(Y[j]), "$R2", sep = "")))))
  # R[laps,3:4] <- c(sqrt(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$mse",  sep = "")))),
  #                  sqrt(mean(eval(parse(text = paste("out.data$", names(Y[j]), "$MSE", sep = ""))))))
  save(mod.out, file = paste("models_out_", names[j],  sep = ""))
  save(mod.stats,  file = paste("models_stats_", names[j],  sep = ""))
  save(mod.comps, file = paste("models_comps_", names[j],  sep = ""))
}
