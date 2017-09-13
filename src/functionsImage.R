normalizeImg<-function(hsp){
  allData = hsp
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  #allData[,-c(1,2)] <- allData[,-c(1,2)] /10000
  # Set bad bands to zero
  ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
  nir860 <- (allData$band_96 + allData$band_97)/2 
  allData[which(ndvi < 0.7 | nir860 < .3),]=NA
  #allData <- allData[complete.cases(allData), ]
  allData[,colnames(allData) %in% bad_Bands]=NA
  
  #remove any reflectance bigger than 1
  #pixel_crownID <- allData[,2]
  #allData <- allData[,-c(1,2)]
  allData[allData>1]=NA
  #allData <- cbind(pixel_crownID,allData)
  
  # Find unique crowns (only for plotting purposes)
  #unqCrown=unique(allData$pixel_crownID)
  # Extract spectra into matrix
  #specMat=as.matrix(allData[,all_Bands])
  specMat = allData[,colnames(allData) %in% all_Bands]
  
  # Vector normalize spectra
  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  
  # Write vector normalized spectra back into dataframe
  #normDF=allData
  #normDF[,all_Bands]=normMat
  normDF = as.data.frame(normMat)
  # Write vector normalized spectra to CSV
  return(normDF)
}
