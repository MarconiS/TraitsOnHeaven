rm(list=ls(all=TRUE))   # clear workspace
setwd('/Users/sergiomarconi/Documents/Projects/TraitsOnHeaven/')


main <- function(NeonSite = "OSBS", 
                 loops=1000,
                 rebuild = T,
                 wd = '~/Documents/Projects/TraitsOnHeaven/', 
                 out.dir = paste(getwd(), "/outputs/", sep=""),
                 in.dir = paste(getwd(), "/inputs/", sep=""),  
                 names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")){
  
  #---------------- Load required libraries ---------------------------------------------------------#
  library (boot)
  library (glmnet)
  library (leaps)
  library(raster)
  library(readr)
  library(reshape2)
  library(rgdal)
  library(itcSegment)
  require(MASS)
  library(pls)
  library(sm)
  library(tidyverse)
  
  
  source(paste(getwd(), '/src/functions_build.R', sep=""))
  source(paste(getwd(), '/src/scr_build.R', sep=""))
  source(paste(getwd(), '/src/functions_infer.R', sep=""))
  source(paste(getwd(), '/src/scr_infer.R', sep=""))
  
  path = paste(wd, NeonSite, sep="/")
  setwd(path)
  rebuild = F
  CrownIDS <- read_csv(paste(path, "inputs", "Spectra", 'CrownIDS.csv', sep="/"))
  cid <- eval(parse(text = paste("CrownIDS", NeonSite, sep="$")))
  cid <- cid[!is.na(cid)]
  nCrowns = length(cid) #85
  if(!file.exists(paste(in.dir, "/Spectra/CrownPix_norm.csv", sep=""))){
    emptyCrowns <- normalize(CrownIDS)
    if(!is_empty(emptyCrowns)){warnings(paste("there are crowns which are empty. run whoIsEmpty() to get which one"))}
    pixPerm(rounds, 1000, unqCrown = cid, names, path=path)
  }
  if(rebuild = T){
    PLS(rounds,1000, names)
    setwd(path)
    PLS_DA(rounds,1000, names = "name")
    setwd(path)
  }
  performance <- perform_summary(names, "baggedTraits.csv")
  getSpatialRegression()
  getTreeRegression()
}

main()