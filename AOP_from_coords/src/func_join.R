file_path <- "/Users/sergiomarconi/Documents/Chapter2/JointTraitsDist/Extracted_by_NEON/outputs/buffer/Spectra/"
flist <- list.files(file_path)
#flist <- flist[-grep("OSBS", flist)]

#init database
dataset = readr::read_csv(paste(file_path,flist[1], sep="/"))
idinfo <- unlist(strsplit(flist[1], split = "_"))
dataset <- cbind(gsub('.{4}$', '', idinfo[4]) , idinfo[1], idinfo[2], idinfo[3], dataset)
colnames(dataset) <- c("individualID", "flightpath", "band_site", "band_species", "band_chm", paste("band", 1:368, sep="_"))
for(i in flist[-1]){
  tryCatch({
    tmp <- readr::read_csv(paste(file_path,i, sep="/"))
    idinfo <- unlist(strsplit(i, split = "_"))
    tmp <- cbind(gsub('.{4}$', '', idinfo[4]) , idinfo[1], idinfo[2], idinfo[3], tmp)
    colnames(tmp) <- colnames(dataset)
    
    dataset <- rbind(dataset, tmp)
  },error=function(e){print(e)})
}
readr::write_csv(dataset, "/Users/sergiomarconi/Documents/Chapter2/JointTraitsDist/buffered_features.csv")



get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32617, 32616, 32616, 32612, 32613, 32617) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", "JORN", "OSBS")
  return(dictionary[colnames(dictionary)==utm])
}

get_buffer_itc <- function(){
  library(raster)
  require(rgeos)
  require(rgdal)
  library(tidyverse)
  library(doParallel)
  source("./AOP_from_coords/src/get_crown_area.R")
  
  # place a circle around the crown
  file_path <- "./AOP_from_coords/outputs/itcTiff/"
  wd = "./AOP_from_coords/"
  flist <- list.files(file_path)
  #flist <- flist[-grep("OSBS", flist)]
  
  # estimate crown allometry from dbh, where needed
  tos_data <- crowne_size_data <- readr::read_csv("./TOS_retriever/out/field_data.csv") %>%
    dplyr::select("individualID","taxonID", "siteID","nlcdClass","elevation",
                  "UTM_E","UTM_N","stemDiameter","height","baseCrownHeight","maxCrownDiameter","ninetyCrownDiameter") %>%
    unique %>%
    get_crown_area
  
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  #results <- foreach(z = clean_hps, .combine = 'cbind', .verbose = T) %:%
  results <- foreach(i = 1:length(flist), .verbose = T) %dopar% {
    library(raster)
    require(rgeos)
    require(rgdal)
    library(tidyverse)
    clip<- raster::stack(paste(file_path, flist[i], sep="/"))
    centroid <- data.frame(individualID = gsub('.{4}$', '', flist[i]) , 
                           UTM_E = extent(clip)[1] + 20, 
                           UTM_N = extent(clip)[3] + 20, stringsAsFactors = F)
    coordinates(centroid) <- ~UTM_E+UTM_N   
    NeonSites <- unlist(strsplit(gsub('.{4}$', '', flist[i]), split = "_"))
    proj4string(centroid)=CRS(paste("+init=epsg:", get_epsg_from_utm(NeonSites[2]), sep=""))
    proj4string(clip) <- CRS(paste("+init=epsg:", get_epsg_from_utm(NeonSites[2]), sep=""))
    
    #create buffer
    max_crd <- tos_data[tos_data$individualID == NeonSites[4],"maxCrownDiameter"]
    if(is.na(max_crd) || length(unlist(max_crd))==0){
      max_crd <- tos_data[tos_data$taxonID==NeonSites[3], "maxCrownDiameter"] %>%
        unlist %>%
        mean(na.rm = T)
      if(is.nan(max_crd)){
        max_crd = 8
      }
      #max_crd = mean(tos_data[tos_data$taxonID==NeonSites[3], "maxCrownDiameter"], na.rm =T)
    }
    buf = gBuffer(centroid,width=max_crd/2)
    #extract data
    spectra_from_round <- raster::extract(clip, buf)
    #create outputs
    writeOGR(SpatialPolygonsDataFrame(buf, data.frame(t(NeonSites)),  match.ID = F), 
             paste(wd,'/outputs/buffer/itcID/', sep = ""), gsub('.{4}$', '', flist[i]),  
             driver="ESRI Shapefile", overwrite_layer = T)
    readr::write_csv(as.data.frame(spectra_from_round), paste(wd,'/outputs/buffer/Spectra/',gsub('.{4}$', '', flist[i]),".csv", sep=""))
  }
  stopCluster(cl)
}  