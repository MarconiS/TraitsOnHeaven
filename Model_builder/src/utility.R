# utility functions
na.rm <- function(dat){
  na_rm <- rowSums(is.na(dat[grepl("band", names(dat))])) == 0
  dat=dat[na_rm, ]
}

no.val.rm <- function(dat, no_value = 10000){
  dat <- dat[!rowSums(dat[grepl("band", names(dat))] >no_value),]
  return(dat)
}

optical.filter<- function(dat){
  ndvi <- (dat$band_90- dat$band_58)/(dat$band_58 + dat$band_90) <0.7
  nir860 <- (dat$band_96 + dat$band_97)/20000 < 0.3
  naval = as.logical(ndvi | nir860)
  dat[naval,] = NA
  return(dat)
}

hiper_features <- function(dat, normalization = "singh", augment = "none"){
  dat <- dat[grepl("band", names(dat))] %>%
    as.matrix
  
  if(normalization == "singh"){
    normMat <- sqrt(apply(dat^2,FUN=sum,MAR=1, na.rm=TRUE)) 
    normMat <- matrix(data=rep(normMat,ncol(dat)),ncol=ncol(dat)) 
    normMat=dat/normMat
  }
  
}

mod.r2 <- function(pred, obs){
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
} 

category_to_dummy <- function(cat){
  #dummy variable example
  foo=list()
  foo$band_ <- cat
  mat <- model.matrix( ~ band_ - 1, data=foo)
  return(mat)
}
