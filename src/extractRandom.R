rPix <-function(names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"), 
                path = ("/Users/sergiomarconi/Documents/PhD/Projects/MappingTraits"), round){
  
  setwd(path)
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  # Read data
  allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
  allData=imp.spectra("Spectra/CrownPix.csv",in.dir)  
  
  for (j in 1:6){
    allData=read.csv(paste(in.dir, "Bootstrap_", round, "/secondCrownRound_",names[j], ".csv", sep=""))
    all_Bands=as.character(allBand$BandName)
    bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
    
    # Set bad bands to zero
    allData[,bad_Bands]=NA
    allData=allData[, colSums(is.na(allData)) == 0]
    
    # Find unique crowns (only for plotting purposes)
    unqCrown=unique(allData$pixel_crownID)
    bootDat = bootPix <- allData[1:length(unqCrown), ]
    bootDat[] <- NA
    bootPix[] <- NA
    for (laps in 1:400){
      tk = 1
      for(i in unqCrown){
        bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
        bootPix[tk,] <- sample(which(allData$pixel_crownID==i), 1)
        tk = tk +1
      }
      write.csv(bootDat, paste('inputs/Bootstrap_', round,'/onePix1Crown_',names[j], laps, '.csv', sep = ''))
      write.csv(bootPix, paste('inputs/Bootstrap_',round,'/onePix1Position_',names[j], laps, '.csv', sep = ''))
    }
  }
# }else{
#   # Find bad bands
#   all_Bands=as.character(allBand$BandName)
#   bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
#   
#   # Set bad bands to zero
#   allData[,bad_Bands]=NA
#   allData=allData[, colSums(is.na(allData)) == 0]
#   
#   # Find unique crowns (only for plotting purposes)
#   unqCrown=unique(allData$pixel_crownID)
#   bootDat = bootPix <- allData[1:length(unqCrown), ]
#   bootDat[] <- NA
#   bootPix[] <- NA
#   for (laps in 1:1000){
#     tk = 1
#     for(i in unqCrown){
#       bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
#       bootPix[tk,] <- sample(which(allData$pixel_crownID==i), 1)
#       tk = tk +1
#     }
#     write.csv(bootDat, paste('onePix1Crown_', laps, '.csv', sep = ''))
#     write.csv(bootPix, paste('onePix1Position_', laps, '.csv', sep = ''))
#   }
# }
}