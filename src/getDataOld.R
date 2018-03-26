#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)

myFun <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE))
  a <- paste0(a, sprintf("%04d", sample(9999, 1, TRUE)), sample(LETTERS, n, TRUE))
}


traitsFromITC <- function(tile.args, nm = NULL, old=F, ntop = 20){
  library(raster)
  source(paste(getwd(), '/src/function_spatial.R', sep=""))
  library(rgdal)
  library(maptools)
  rasterOptions(tmpdir="./tmp")
  tile.args <- data.frame(t(tile.args), stringsAsFactors = F)

  tile.args <- transform(tile.args, xx = as.integer(xx),yy = as.integer(yy),
                         xmin = as.integer(xmin),xmax = as.integer(xmax), scaled = as.logical(scaled),
                         ymin = as.integer(ymin),ymax = as.integer(ymax), epsg = as.integer(epsg))
  #create layer of delineated crowns
  NeonSite <- tile.args$site
  if(NeonSite=="OSBS"){
    ras_dir = "/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/D03/OSBS/L4/Rasters"
    itc_dir = "/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/D03/OSBS/L4/ITCs"
    Domain <- "D03"
  }else if (NeonSite=="TALL"){
    ras_dir = "/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/D08/TALL/L4/Rasters"
    itc_dir = "/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/D08/TALL/L4/ITCs"
    Domain <- "D08"

  }
  whichMeta <- "./src/stripeToRasterOld.py"
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE))
  tileID <- paste0(a, sprintf("%04d", sample(9999, 1, TRUE)), sample(LETTERS, 1, TRUE))

  system2("python", args=(sprintf('"%1$s" "%2$s" "%3$s" "%4$i" "%5$i" "%6$s" "%7$s" "%8$s"', whichMeta,
                                  tile.args$f, tile.args$path_h5, tile.args$xx, tile.args$yy, as.character(tile.args$epsg), ras_dir, tileID)))


  f_nm_1 = as.integer(tile.args$xmin/1000)*1000 +1000*as.integer(tile.args$xx)
  f_nm_2 = as.integer(tile.args$ymin/1000)*1000 +1000*(as.integer(tile.args$yy))
  itc_path <-paste("/ufrc/ewhite/s.marconi/NeonData/2015_Campaign", Domain, NeonSite, "L1/DiscreteLidar/Classified_point_cloud/",sep="/")
  itc <- list.files(itc_path, pattern = paste(f_nm_1, f_nm_2, sep ="_"))[1] # paste(f_nm_1, f_nm_2, sep ="_"))

  #produce the itc layer
  crownIT(las_id = paste(f_nm_1, f_nm_2, sep ="_"), pt = itc_path,  tileID = tileID, out_path = itc_dir, epsg=tile.args$epsg,
          NeonSite = NeonSite, method = 'silva', max_cr_factor = 0.6, exclusion = 0.5, mv = 7, minh = 7)
  print("ITC saved")
}

# pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/"
# f= "NIS1_20140507_143910_atmcor.h5"
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

  print(NeonSite)
  library(rhdf5)
  library(doParallel)
  if(NeonSite=="OSBS"){
    pt_h5 = paste(pt, "D03/OSBS/L1/Spectrometer/Reflectance/", sep = "")
    spInfo <- h5read(paste(pt_h5,f,sep=""), "/map info")
    mapInfo <- as.integer(unlist(strsplit(spInfo, ",")))
    pathSize <- dim(h5read(paste(pt_h5,f,sep=""), "/Aspect"))[1:2]
    #there was no epsg in the former version??? check!
    epsg <- "32617"
  }else if(NeonSite=="TALL"){
    pt_h5 = paste(pt, "D08/TALL/L1/Spectrometer/Reflectance/", sep = "")
    spInfo <- h5read(paste(pt_h5,f,sep=""), "/map info")
    mapInfo <- as.integer(unlist(strsplit(spInfo, ",")))
    pathSize <- dim(h5read(paste(pt_h5,f,sep=""), "/Aspect"))[1:2]
    #there was no epsg in the former version??? check!
    epsg <- "32616"
  }


  xmin <- mapInfo[4]
  ymax <- mapInfo[5]
  xmax <- xmin + pathSize[1]
  ymin <- ymax - pathSize[2]
  nx = as.integer(xmax/1000) - as.integer(xmin/1000)
  ny = as.integer(ymax/1000) - as.integer(ymin/1000)
  rm(spInfo)
  tile.args <- expand.grid(0:nx,0:ny)
  tile.args <- data.frame(tile.args, xmin, xmax, ymin, ymax, epsg, NeonSite, pt_h5, f, scaled, itc_method, stringsAsFactors = F)
  colnames(tile.args) <- c("xx","yy","xmin", "xmax","ymin","ymax", "epsg", "site", "path_h5", "f", "scaled", "itc_method")

  cores <-32
  registerDoSEQ()
  cl <- makeCluster(cores, outfile = "")
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  registerDoParallel(cl)
  parApply(cl, tile.args, 1, traitsFromITC, nm = nm, old = T)
  stopCluster(cl)
}
getITC_and_raster(f = args[7], pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/", NeonSite = "TALL", nm = "LMA_g.m2")
#getITC_and_raster(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = "TALL", nm = "N_pct")
