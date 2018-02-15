args = commandArgs(trailingOnly=TRUE)

buildCore <- function(NeonSite = "ALL", nm = NULL,  year = NULL, epsg = NULL, loops=1000, wd = "/ufrc/ewhite/s.marconi/Marconi2018/", 
rebuild = T, rescale = F,nrmlz = F, tile = 80, spatial = T, deriv.reg = F, deriv.class = T,perf = T,
                 nm = c("LMA_g.m2", "C_pct","N_pct", "P_pct"),method = 'pls',cores = 32){
  
  #---------------- Load required libraries ---------------------------------------------------------#
  library(plsRglm)
  library(gnm)
  library (boot)
  library (glmnet)
  library (leaps)
  library(raster)
  library(readr)
  library(reshape2)
  library(rjson)
  library(rgdal)
  require(MASS)
  library(pls)
  library(sm)
  library(tidyverse)
  library(parallel)
  
  
  setwd(wd)
  
  source(paste(getwd(), '/src/build_functions.R', sep=""))
  path = paste(wd, "ModelBuild", NeonSite, sep="/")
  setwd(path)
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  CrownIDS <- read_csv(paste('./inputs/Spectra/', 'CrownIDS.csv', sep="/"))
  cid <- eval(parse(text = paste("CrownIDS", NeonSite, sep="$")))
  cid <- cid[!is.na(cid)]
  nCrowns = length(cid) #85
  if(!file.exists(paste(in.dir, "/Spectra/CrownPix_norm.csv", sep=""))){
    emptyCrowns <- normalize(CrownIDS, NeonSite,rescale = rescale, in.dir = in.dir, out.dir = out.dir)
    if(!is_empty(emptyCrowns)){warnings(paste("there are crowns which are empty. run whoIsEmpty() to get which one"))}
    pixPerm(loops= loops, unqCrown = cid, names, path=path)
  } 
  if(rebuild == T){
	  if(method=='plsglm'){
        foo <- pls_glm(j = nm, loops = 1000, cores = cores, nrmlz = nrmlz)
        save(foo, file = paste("./pls_glm_", j,  sep = ""))
      }else if(method == 'pls'){
        foo <- pPLS(j = nm, loops = 1000, cores = cores, nrmlz = nrmlz)
        save(foo, file = paste("./pls_", j,  sep = ""))
      }
    }
  }
}

#run the whole pipeline
buildCore(NeonSite = args[1], nm = args[2],  cores = args[3], loops = 1000,  rebuild = T, spatial = F, nrmlz = F, method = "pls")
