main_build <- function(NeonSite = "ALL", year = NULL, epsg = NULL, loops=1000, rebuild = T, rescale = F,nrmlz = F,
                 wd = "/ufrc/ewhite/s.marconi/Marconi2018/", tile = 80, spatial = T, deriv.reg = F, deriv.class = T,perf = T,
                 nm = "N_pct",method = 'pls',cores = 32){

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

  source(paste(getwd(), '/src/functions_build.R', sep=""))
  source(paste(getwd(), '/src/src_build.R', sep=""))
  #source(paste(getwd(), '/src/functions_infer.R', sep=""))
  #source(paste(getwd(), '/src/src_infer.R', sep=""))
  #source(paste(getwd(), '/src/pls_par.R', sep=""))


  path = paste(wd, "ModelBuild", NeonSite, sep="/")
  setwd(path)
  out.dir = paste(getwd(), "/outputs/", sep="")
  in.dir = paste(getwd(), "/inputs/", sep="")
  CrownIDS <- read_csv(paste('./inputs/Spectra/', 'CrownIDS.csv', sep="/"))
  cid <- eval(parse(text = paste("CrownIDS", NeonSite, sep="$")))
  cid <- cid[!is.na(cid)]
  nCrowns = length(cid) #85
  #if(!file.exists(paste(in.dir, "/Spectra/CrownPix_norm.csv", sep=""))){
    emptyCrowns <- normalize(CrownIDS, NeonSite,rescale = rescale, in.dir = in.dir, out.dir = out.dir)
    if(!is_empty(emptyCrowns)){
      warnings(paste("there are crowns which are empty. run whoIsEmpty() to get which one"))
      whoIsEmpty()
      }
    pixPerm(loops= loops, unqCrown = cid, names, path=path)
  #}
  if(rebuild == T){
    foo <- pls_glm(j = nm, loops = 1000, cores = cores, nrmlz = T)
    save(foo, file = paste("/ufrc/ewhite/s.marconi/Marconi2018/ModelBuild/",NeonSite,  "/pls_glm_", nm,  sep = ""))
  }
  get_stats()
}
#c("LMA_g.m2", "C_pct","N_pct", "P_pct")
#run the whole pipeline
main_build(NeonSite = "ALL",loops = 1000,cores = 63, nm = "N_pct", rebuild = T, spatial = F, nrmlz = T)
main_build(NeonSite = "ALL",loops = 1000,cores = 63, nm = "P_pct", rebuild = T, spatial = F, nrmlz = T)
main_build(NeonSite = "ALL",loops = 1000,cores = 63, nm = "LMA_g.m2", rebuild = T, spatial = F, nrmlz = T)

#main(NeonSite = "OSBS",loops = 1000,cores = 63, nm = c("C_pct"), rebuild = T, spatial = F, nrmlz = T)
#main(NeonSite = "TALL",loops = 1000,cores = 63, nm = c("C_pct"), rebuild = T, spatial = F, nrmlz = T)
