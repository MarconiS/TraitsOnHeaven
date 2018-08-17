#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)
print(args)

build_model <- function(loops=1000, cores = NULL, nrmlz = T,
                        wd = "/ufrc/ewhite/s.marconi/Chapter1/Model_builder", nm = "P_pct"){
  
  #---------------- Load required libraries ---------------------------------------------------------#
  library(plsRglm)
  library(gnm)
  library (boot)
  library (glmnet)
  library(raster)
  library(readr)
  library(reshape2)
  library(rgdal)
  require(MASS)
  source(paste(wd, '/src/src_build.R', sep=""))
  
  
  out_dir = paste(wd, "/outputs/", sep="")
  in_dir = paste(wd, "/inputs/", sep="")
  if(length(list.files(paste(in_dir, 'Permutations', sep="/"), pattern=".csv"))<loops){
    pixPerm(in_dir, out_dir, loops= loops)
  }
  foo <- pls_glm(nm = nm, wd = wd, loops = loops, cores = cores, nrmlz = nrmlz)
  save(foo, file = paste(out_dir,  "/pls_glm_nodiff", nm,  sep = ""))
}
build_model(loops = 1000,cores = 64, nm = "P_pct", nrmlz = T)
