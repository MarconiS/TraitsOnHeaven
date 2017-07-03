DiscardandRerun <- function (names, round, loops){
  library(readr)
  for (j in 1:6){
    pls_coef <- read_csv(paste(out.dir,'pls_coeffs_',names[j],'.csv',sep=""), col_names = FALSE)
    lasso_coef <- read_csv(paste(out.dir,'las_coeffs_',names[j],'.csv',sep=""), col_names = FALSE)
    
    press <- data.frame(cbind(seq(1,loops),pls_coef$X1))
    r2.las <-data.frame(cbind(seq(1,loops),lasso_coef$X1))
    p <- press[order(-press$X2),]
    p.ls <- r2.las[order(r2.las$X2),]
    head(p)
    p <- p[1:100,1]
    p.ls <- p.ls[1:100,1]
    
    pixels <- data.frame(matrix(NA, ncol=nCrowns, nrow = 200))
    tk = 0
    for(k in p){
      tk = tk + 1
      #import ith permutation associated with the horrible correlation
      tmp <- imp.spectra(paste('Bootstrap_',round,'/onePix1Position_', k, '.csv', sep = ''), in.dir)
      pixels[tk,] <- tmp[,3]
    }
    for(k in p.ls){
      tk = tk + 1
      #import ith permutation associated with the horrible correlation
      tmp <- imp.spectra(paste('Bootstrap_1/onePix1Position_', k, '.csv', sep = ''), in.dir)
      pixels[tk,] <- tmp[,3]
    }
    
    discarded = kept = NULL
    for(k in 1:nCrowns){
      flush.console()
      rank.p <- as.data.frame(table(pixels[,k]))
      #order frequency of bad pixels. Cut off pixels summing up to 85%?
      rank.p <- rank.p[order(-rank.p$Freq),]
      tk = qtl = 0
      while(qtl < 90){
        tk <- tk +1
        qtl <- rank.p$Freq[tk] + qtl
      }
      discarded <- c(discarded, as.character(rank.p$Var1[1:tk-1]))
      kept <- c(kept, as.character(rank.p$Var1[tk:dim(rank.p)[1]]))
    }
    write.csv(discarded, paste(in.dir, 'Bootstrap_2/badPix_', names[j],'.csv', sep = ''))
    write.csv(kept, paste(in.dir, 'Bootstrap_2/goodPix_', names[j],'.csv', sep = ''))
    
    allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
    allData=imp.spectra("Spectra/CrownPix.csv",in.dir)  
    # Find bad bands
    all_Bands=as.character(allBand$BandName)
    bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
    
    # Subset bad bands to zero
    allData <- allData[kept,]
    write.csv(allData, paste(in.dir, 'Bootstrap_2/secondCrownRound_',names[j],'.csv', sep = ''))
  }
}