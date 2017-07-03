
# Read data
allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("./inputs/Spectra/CrownPix.csv")
#allData=read.csv("/Users/sergiomarconi/Documents/Classes/Dropbox/FinalProject_Marconi/FinalProject_Marconi/Inputs/foliar_spectra_full.csv")
# Find bad bands
all_Bands=as.character(allBand$BandName)
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])

# Set bad bands to zero
allData[,bad_Bands]=NA

# Find unique crowns (only for plotting purposes)
unqCrown=unique(allData$pixel_crownID)

# Extract spectra into matrix
specMat=as.matrix(allData[,all_Bands])
#matplot(t(specMat),type="l")

# Vector normalize spectra
normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
normMat=specMat/normMat
matplot(t(normMat),type="l")

# Write vector normalized spectra back into dataframe
normDF=allData
normDF[,all_Bands]=normMat

# Write vector normalized spectra to CSV
write.csv(normDF, "./inputs/Spectra/CrownPix_norm.csv",row.names=FALSE)

# Write comparison plots out
#pdf(file="./Outputs/44_CrownPix.pdf",width=11,height=6.5,paper="USr")
crown = 44
{
  print(crown)
  flush.console()
  par(mfrow=c(1,2))
  wv <- seq(385, 2510, 5)
  foo <- t(specMat[which(allData$pixel_crownID==crown),all_Bands])
  plot(wv, foo[,1], type="l",lty=1,col="red",ylab="Reflectance",main='Original spectra',xlab="wavelength",cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3, ylim=c(0,0.5))
  for(iii in 2:20){
    par(new = TRUE)
    plot(wv, foo[,iii],axes=F, type="l", xlab='', lty=1,col="red",ylab='',  ylim=c(0,0.5))
  }
  
  foo.norm <- t(normMat[which(allData$pixel_crownID==crown),all_Bands])
  plot(wv, foo.norm[,1], type="l",lty=1,col="blue",ylab="SVN(reflectance)",main='Normalized spectra',xlab="wavelength",cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3, ylim=c(0,0.12))
  for(iii in 2:20){
    par(new = TRUE)
    plot(wv, foo.norm[,iii], type="l", axes=F,xlab='', lty=1,col="blue",ylab='',ylim=c(0,0.12))
  }
}
#dev.off()

