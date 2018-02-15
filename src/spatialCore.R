#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

spatialCore <- function(pt ="./spatialPhase/hpsOSBS/", f = "NEON_D03_OSBS_DP1_20140507_150703_reflectance.h5", 
                         nm = c("LMA_g.m2", "N_pct", "P_pct"), cores = 16, itc_method = "silva", scaled = F){
  setwd("/ufrc/ewhite/s.marconi/Marconi2018/")
  .libPaths(c("/ufrc/ewhite/s.marconi/Marconi2018/RstudioLibs", 
              "/home/s.marconi/R/x86_64-pc-linux-gnu-library/3.4",
              "/apps/gcc/6.3.0/R/3.4.0/lib64/R/library"))
  
  print(.libPaths())
  library(rhdf5)
  library(doParallel)
  source(paste(getwd(), '/src/sp_functions.R', sep=""))
	
  
  spInfo <- h5read(paste(pt,f,sep=""),"/OSBS/Reflectance/Metadata")
  pathSize <- dim(spInfo$Ancillary_Imagery$Aspect)
  epsg <- spInfo$Coordinate_System$`EPSG Code`
  mapInfo <- as.integer(unlist(strsplit(spInfo$Coordinate_System$Map_Info, ",")))
  xmin <- mapInfo[4]
  ymax <- mapInfo[5]
  xmax <- xmin + pathSize[1]
  ymin <- ymax - pathSize[2]
  nx = as.integer(xmax/1000) - as.integer(xmin/1000)
  ny = as.integer(ymax/1000) - as.integer(ymin/1000)
  rm(spInfo)
  #get subtile
  
  tile.args <- expand.grid(0:nx,0:ny)
  tile.args <- data.frame(tile.args, xmin, xmax, ymin, ymax, epsg, pt, f, scaled, itc_method, stringsAsFactors = F)
  colnames(tile.args) <- c("xx","yy","xmin", "xmax","ymin","ymax", "epsg", "path", "f", "scaled", "itc_method")

  registerDoSEQ()
  cl <- makeCluster(cores, outfile = "")
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  registerDoParallel(cl)
  parApply(cl, tile.args, 1, traitsFromITC, nm = nm)
  stopCluster(cl)
}

spatialCore(f = args[1], cores = 15)
sessionInfo()