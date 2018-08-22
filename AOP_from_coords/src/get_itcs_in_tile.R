get_itcs_in_tile <- function(hps_fi = NULL, centroids = NULL, f_path = NULL, NeonSites = NULL, buffer = 2){
  library(rhdf5)
  h5name <- paste(f_path,hps_fi,sep = "/")
  #mapInfo <- h5read(h5name,"map info")
  #mapInfo<-unlist(strsplit(mapInfo, ","))
  
  #new metadata ['Metadata']['Coordinate_System']['Map_Info']
  mapInfo <- h5read(h5name,paste(NeonSites, "/Reflectance/Metadata", sep=""))
  mapInfo <-unlist(strsplit(mapInfo$Coordinate_System$Map_Info, ","))
  #grab the utm coordinates of the lower left corner
  xMin<-as.numeric(mapInfo[4])
  yMax<-as.numeric(mapInfo[5])
  reflInfo <- h5ls(h5name)
  #reflInfo <- unlist(strsplit(reflInfo$dim[3], split = "x"))
  reflInfo <- unlist(strsplit(reflInfo$dim[7], split = "x"))
  
  nRows <- as.integer(reflInfo[2])
  nCols <- as.integer(reflInfo[1])
  xMax <- xMin + nCols
  yMin <- yMax - nRows
  #be sure you canclip a 60 by 60 pic
  mask.x <- (centroids$easting < xMax-buffer) & (centroids$easting > buffer + xMin)
  mask.y <- (centroids$northing < yMax-buffer) & (centroids$northing > yMin+buffer)
  itcextract <- centroids[mask.x & mask.y,]
  h5closeAll()
  return(itcextract)
}
