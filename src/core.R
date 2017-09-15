rm(list=ls(all=TRUE))   # clear workspace


main <- function(NeonSite = "OSBS", loops=1000, rebuild = T,
                 wd = '~/Documents/Projects/TraitsOnHeaven/', 
                 names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")){
  
  #---------------- Load required libraries ---------------------------------------------------------#
  library (boot)
  library (glmnet)
  library (leaps)
  library(raster)
  library(readr)
  library(reshape2)
  library(rjson)
  library(rgdal)
  library(itcSegment)
  require(MASS)
  library(pls)
  library(sm)
  library(tidyverse)
  
  setwd(wd)
  
  source(paste(getwd(), '/src/functions_build.R', sep=""))
  source(paste(getwd(), '/src/src_build.R', sep=""))
  source(paste(getwd(), '/src/functions_infer.R', sep=""))
  source(paste(getwd(), '/src/src_infer.R', sep=""))
  
  path = paste(wd, NeonSite, sep="/")
  setwd(path)
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  CrownIDS <- read_csv(paste(path, "inputs", "Spectra", 'CrownIDS.csv', sep="/"))
  cid <- eval(parse(text = paste("CrownIDS", NeonSite, sep="$")))
  cid <- cid[!is.na(cid)]
  nCrowns = length(cid) #85
  if(!file.exists(paste(in.dir, "/Spectra/CrownPix_norm.csv", sep=""))){
    emptyCrowns <- normalize(CrownIDS, NeonSite, in.dir = in.dir, out.dir = out.dir)
    if(!is_empty(emptyCrowns)){warnings(paste("there are crowns which are empty. run whoIsEmpty() to get which one"))}
    pixPerm(loops= 1000, unqCrown = cid, names, path=path)
    print("pixPerm ok")
  }
  if(rebuild == T){
    PLS(loops=1000, names, in.dir = in.dir, out.dir = out.dir)
    print("PLS ok")
    setwd(path)
    PLS_DA(loops=1000, names = "name", in.dir = in.dir, out.dir = out.dir)
    setwd(path)
    print("PLS_DA ok")
  }
  performance <- perform_summary(names, "baggedTraits.csv", in.dir = in.dir, out.dir = out.dir)
  print(performance)
  getPointCloud(NeonString = "2014_OSBS_", in.dir = in.dir, out.dir = out.dir)
  print("pointCloud ok")
  getSpatialRegression(names = "name", in.dir = in.dir, out.dir = out.dir)
  print("getSpatialRegression class ok")
  getSpatialRegression(names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"), in.dir = in.dir, out.dir = out.dir)
  print("getSpatialRegression regr ok")
  setwd(path)
  getTreeRegression(in.dir = in.dir, out.dir = out.dir)
}

main(NeonSite = "OSBS", rebuild = T)
