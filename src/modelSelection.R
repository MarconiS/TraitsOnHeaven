#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#' 
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}

#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
#'              the PRESS statistic.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}


extractPerm <- function(sz=1000){
  # Read data
  allBand=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/neon_aop_bands.csv")
  allData=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/CrownPix_norm.csv")
  
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  
  # Set bad bands to zero
  allData[,bad_Bands]=NA
  allData=allData[, colSums(is.na(allData)) == 0]
  
  # Find unique crowns (only for plotting purposes)
  unqCrown=unique(allData$pixel_crownID)
  bootDat <- allData[1:length(unqCrown), ]
  bootDat[] <- NA
  for (laps in 1:sz){
    tk = 1
    for(i in unqCrown){
      bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
      tk = tk +1
    }
    write.csv(bootDat, paste('Inputs/BootNorm/onePix1Crown_', laps, '.csv', sep = ''))
  }
}



normAditya <- function(){
  # Read data
  allBand=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/neon_aop_bands.csv")
  allData=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  # Set bad bands to zero
  allData[,bad_Bands]=NA
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
  write.csv(normDF, "/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/CrownPix_norm.csv",row.names=FALSE)
}