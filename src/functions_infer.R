#-----normalizeImg------------------------------------------------------------------------------------------

normalizeImg<-function(hsp, scaled = T){
  allData = hsp
  # Read data
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  if(scaled = F) {allData[,-c(1,2)] <- allData[,-c(1,2)] /10000}
  # Set bad bands to zero
  ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
  nir860 <- (allData$band_96 + allData$band_97)/2 
  allData[which(ndvi < 0.7 | nir860 < .3),]=NA
  allData[,colnames(allData) %in% bad_Bands]=NA
  
  #remove any reflectance bigger than 1
  allData[allData>1]=NA
  specMat = allData[,colnames(allData) %in% all_Bands]
  
  # Vector normalize spectra
  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  
  normDF = as.data.frame(normMat)
  # Write vector normalized spectra to CSV
  return(normDF)
}
#-----read.centroids------------------------------------------------------------------------------------------

read.centroids <- function(file){
  plot.centr <- readOGR("./inputs/Geofiles/vector", file, stringsAsFactors = FALSE)
  return(plot.centr)
}

#-----loadIMG------------------------------------------------------------------------------------------

loadIMG <- function(path, proj, img.lb = 1){
  path.raster <- paste(in.dir, "Geofiles/RSdata/", sep="")
  #chm.sample <- raster(paste(path.raster, "LiDAR/", 'plot_', img.lb, "_chm.tif", sep = ""))
  hs.sample <- brick(paste(path.raster, "hs/", img.lb, "_hyper.tif", sep = ""))
  return(list(hsp = hs.sample))
}


