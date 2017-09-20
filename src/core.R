#rm(list=ls(all=TRUE))   # clear workspace


main <- function(NeonSite = "OSBS", year = 2014, epsg = NULL, loops=1000, rebuild = T, rescale = T,
                 wd = '~/Documents/Projects/TraitsOnHeaven/', tile = 80, 
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
    emptyCrowns <- normalize(CrownIDS, NeonSite,rescale = rescale, in.dir = in.dir, out.dir = out.dir)
    if(!is_empty(emptyCrowns)){warnings(paste("there are crowns which are empty. run whoIsEmpty() to get which one"))}
    pixPerm(loops= loops, unqCrown = cid, names, path=path)
    print("pixPerm ok")
  } 
  if(rebuild == T){
    PLS(loops=loops, names, in.dir = in.dir, out.dir = out.dir)
    print("PLS ok")
    setwd(path)
    PLS_DA(loops=loops, names = "name", in.dir = in.dir, out.dir = out.dir)
    setwd(path)
    print("PLS_DA ok")
  }
  #performance <- perform_summary(names, "baggedTraits.csv", in.dir = in.dir, out.dir = out.dir)
  #print(performance)
  print(paste(year, "_", NeonSite,"_", sep=""))
#  getPointCloud(NeonString = paste(year, "_", NeonSite,"_", sep=""), in.dir = in.dir, out.dir = out.dir)
  print("pointCloud ok")
  getSpatialRegression(NeonSite = NeonSite, names = "name", in.dir = in.dir, out.dir = out.dir,
                       tile = tile, epsg = epsg, proj = paste("+init=epsg:",epsg, sep=""))
  print("getSpatialRegression class ok")
  getSpatialRegression(NeonSite = NeonSite, names = names, in.dir = in.dir, out.dir = out.dir,
                       tile = tile, epsg = epsg, proj = paste("+init=epsg:",epsg, sep=""))
  print("getSpatialRegression regr ok")
  setwd(path)
  getTreeRegression(in.dir = in.dir, out.dir = out.dir)
}

#main(NeonSite = "TALL", year = 2015, epsg = 2154, rebuild = F)
main(NeonSite = "OSBS", year = 2014, epsg = 32617, rebuild = F)
