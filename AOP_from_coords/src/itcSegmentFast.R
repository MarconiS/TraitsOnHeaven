crownITC <- function(pt = NULL,wd = NULL, pattern, cores = 2,pybin = "/home/s.marconi/.conda/envs/quetzal3/bin", epsg=NULL, chm_f = NULL){
  library(foreach)
  library(doParallel)
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  results <- foreach(i = pattern) %dopar% {
    #for(i in pattern){
    library(raster)
    library(lidR)
    source(paste(wd, "src/polygonize.R", sep=""))
    #tryCatch({
      if(!length(list.files(paste(wd, "/outputs/itcShp/", sep=""), pattern=i))>0){
        f = list.files(pt, pattern = i)
        
        las = readLAS(paste(pt, f, sep=""))
          
        # normalization
        lasnormalize(las, method = "knnidw", k = 10L)
        
        # compute a canopy image
        chm = grid_canopy(las, res = 0.5, subcircle = 0.2, na.fill = "knnidw", k = 4, p=1)
        chm = as.raster(chm)
        kernel = matrix(1,3,3)
        chm = raster::focal(chm, w = kernel, fun = mean)
        chm = raster::focal(chm, w = kernel, fun = mean)
        
        writeRaster(chm, paste(chm_f, i, "_chm.tif",sep=""), format="GTiff", overwrite=TRUE)
        #silva 2016
        ttops = tree_detection(chm, 5, 2)
        crowns <-lastrees_silva(las, chm, ttops, max_cr_factor = 0.6, exclusion = 0.5, 
                                mv = 7, minh = 7, extra = T)
        Sys.setenv(PATH = paste(pybin, Sys.getenv("PATH"),sep=":"))
        
        x <- gdal_polygonizeR(crowns, pypath = paste(wd, "src/polygonize/", sep=""))
        proj4string(x) <-  CRS(paste("+init=epsg:", epsg, sep=""))
        writeOGR(x, dsn=paste(wd,"outputs/itcShp/", sep=""), paste(i, "silva", sep="_"), overwrite_layer = T, check_exists = T, driver="ESRI Shapefile")
        print(i)
      }else{
        print("tile exists")
      }
   # },error=function(e){})
  }
  stopCluster(cl)
  return(results)
}

