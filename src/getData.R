#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)

traitsFromITC <- function(tile.args, nm = NULL, old=F, ntop = 20){
  library(raster)
  source(paste(getwd(), '/src/function_spatial.R', sep=""))
  library(rgdal)
  library("tsensembler")
  library(plsRglm)
  library(maptools)
  library(data.table)
  rasterOptions(tmpdir="./tmp")
  tile.args <- data.frame((tile.args), stringsAsFactors = F)

  tile.args <- transform(tile.args,
                         xmin = as.integer(xmin),xmax = as.integer(xmax), scaled = as.logical(scaled),
                         ymin = as.integer(ymin),ymax = as.integer(ymax), epsg = as.integer(epsg))
  #create layer of delineated crowns
  NeonSite <- unlist(strsplit(tile.args$f, "_"))[3]
  f_nm_1 = as.integer(tile.args$xmin)
  f_nm_2 = as.integer(tile.args$ymin)
  #python much faster: produce the ith hiperspectral tile
  system2("python", args=(sprintf('"%1$s" "%2$s" "%3$s"',"./src/stripeToRaster.py",
                                  tile.args$f, tile.args$path)))
  #produce the itc layer
  crownIT(las_id = paste(f_nm_1, f_nm_2, sep ="_"), epsg=tile.args$epsg,
          NeonSite = NeonSite, method = 'silva', max_cr_factor = 0.6, exclusion = 0.5, mv = 7, minh = 7)
  print("ITC saved")
}

# pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/"
# f = "NEON_D03_OSBS_DP3_408000_3280000_reflectance.h5"
# pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018"
# itc_method = "silva"
# scaled = F
# NeonSite = "OSBS"
getITC_and_raster <- function(pt =NULL, pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018",
                         f = NULL, nm = NULL, NeonSite = "NULL", cores = 1, itc_method = "silva", scaled = F){
  setwd(pathSpatial)
  .libPaths(c("/usr/lib/R/site-library","/usr/local/lib/R/site-library","/usr/lib/R/library",
              "/home/s.marconi/R/x86_64-pc-linux-gnu-library/3.4",
              "/apps/gcc/6.3.0/R/3.4.0/lib64/R/library"))


  library(rhdf5)
  pt = paste(pt, NeonSite, "/Reflectance/", sep = "")

  spInfo <- h5read(paste(pt,f,sep=""),paste("/", NeonSite, "/Reflectance/Metadata", sep=""))
  pathSize <- dim(spInfo$Ancillary_Imagery$Aspect)
  epsg <- spInfo$Coordinate_System$`EPSG Code`
  mapInfo <- as.integer(unlist(strsplit(spInfo$Coordinate_System$Map_Info, ",")))

  xmin <- mapInfo[4]
  ymax <- mapInfo[5]
  xmax <- xmin + pathSize[1]
  ymin <- ymax - pathSize[2]
  rm(spInfo)

  tile.args <- data.frame(xmin, xmax, ymin, ymax, epsg, pt, f, scaled, itc_method, stringsAsFactors = F)
  colnames(tile.args) <- c("xmin", "xmax","ymin","ymax", "epsg", "path", "f", "scaled", "itc_method")
  traitsFromITC(tile.args, nm = nm, old = !newmeta)
}
getITC_and_raster(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = args[8], nm = "LMA_g.m2")
#getITC_and_raster(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = "TALL", nm = "N_pct")
