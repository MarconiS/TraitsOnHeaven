foreach( ii = list_of_tiffs) %dopar% {
  tryCatch({
    library(raster)
    library(rgdal)
    library(rgeos)
    library(tidyverse)
    library(sf)
    # token = token +1
    img <- brick(paste("./AOP_from_coords/outputs/itcTiff/", ii, sep=""))
    itc_name = strsplit(ii, split = "_")[[1]][4]
    itc_name = gsub('.{4}$', '', itc_name)
    itc_db <- dataset[dataset$individualID == itc_name,]
    if(!is_empty(itc_db$individualID)){
      epsg <- get_epsg_from_utm(itc_db$siteID)
      center <- data.frame(extent(img)[1]+20, extent(img)[3]+20)
      colnames(center) <- c("umt_e", "umt_n")
      radius <- itc_db$maxCrownDiameter / 2
      # dat_sf <- st_as_sf(center, coords = c("umt_e", "umt_n"), crs = epsg) 
      # dat_circles <- st_buffer(dat_sf, dist = radius)
      # #extract data
      crs(img) <- paste("+init=epsg:", epsg, sep="")
      centroid_spdf = SpatialPointsDataFrame(center, proj4string=CRS(paste("+init=epsg:", epsg, sep="")), itc_db)
      itc_spectra <- raster::extract(x = img, y = centroid_spdf, buffer=min(20,radius), df = T)
      entryID = gsub('.{4}$', '', ii)
      write_csv(as.data.frame(itc_spectra), paste(wd,'outputs/buffer_spectra/',entryID, ".csv", sep=""))
    }}, error(cond){print(ii)})
  }
  stopCluster(cl)

  
  cr_per_path <- matrix(NA, nrow = length(hps_f), ncol = 4)
  for(ii in 1:length(hps_f)){
    cr_per_path[ii,] <- tryCatch(
      {get_itcs_in_tile2015(hps_f[ii], centroids=centroids, NeonSites = NeonSites, f_path=f_path)},
      error=function(cond){
        message(paste("hiperspectral flight path corrupted or incoplete:", ii))
        return(0)
      })
  }
  
  
  
  
  
  
  
  
  