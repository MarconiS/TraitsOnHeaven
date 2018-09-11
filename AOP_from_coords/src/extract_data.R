#32611
itcExtract <- function(x, f, chm, itc.f, epsg, token, wd,pybin = "/home/s.marconi/.conda/envs/quetzal3/bin",  buffer = 20){
  library(raster)
  library(rgeos)
  library(rgdal)
  library(readr)
  #256196.6 4108745
  clip.xmin <- max((x$easting) - buffer, as.integer(x$easting/1000)*1000)
  clip.ymin <- max((x$northing) - buffer, as.integer(x$northing/1000)*1000)
  clip.xmax <- min(clip.xmin + 2*buffer, as.integer(x$easting/1000+1)*1000)
  clip.ymax <- min(clip.ymin + 2*buffer, as.integer(x$northing/1000+1)*1000)
  
  #256166.6 4108775
  entryID <- paste(token, x$siteID, x$taxonID,  x$individualID, sep="_")
  #launch python "/Users/sergiomarconi/anaconda3/bin/"
  Sys.setenv(PATH = paste(pybin, Sys.getenv("PATH"),sep=":"))
  system2("python3", args=(sprintf('"%1$s" "%2$s" "%3$s" "%4$s" "%5$f" "%6$f" "%7$f" "%8$f" "%9$s" "%10$s"',paste(wd, "src/extractCrown.py", sep=""), 
                                   f, chm, entryID, clip.xmin, clip.xmax, clip.ymin, clip.ymax, as.character(epsg), wd)))
  clip<- stack(paste(wd, "/outputs/itcTiff/", entryID, '.tif', sep=""))
  #if(max(clip) == min(clip)){
  
  itc <- paste(as.integer(clip.xmin/1000)*1000, as.integer(clip.ymax/1000)*1000, 'silva', sep="_")
  itc_shp <-readOGR(paste(wd, "outputs/itcShp/", sep=""), itc, stringsAsFactors = FALSE)
  proj4string(itc_shp) <-  CRS(paste("+init=epsg:", epsg, sep=""))
  proj4string(clip) <-  CRS(paste("+init=epsg:", epsg, sep=""))
  
  #which cell is that poin in
  centDat <- data.frame(long = clip.xmin+buffer, lat = clip.ymin+buffer, name = entryID)
  coordinates(centDat) <- ~ long + lat
  proj4string(centDat) <- CRS(paste("+init=epsg:", epsg, sep=""))
  
  datExtract <- over(centDat, itc_shp)
  if(is.na(datExtract)){
    print("not really")
    dist <- gDistance(centDat, itc_shp,byid=TRUE)
    datDist <-which(dist==min(dist))
    area = 0
    for(pp  in 1:length(datDist)){
      #check which polygon at the same distance has the greatest area
      foo <- itc_shp[datDist[pp],]
      afoo <-  unlist(lapply(foo@polygons, function(x) slot(x, "area")))
      if(afoo>=area){
        area = max(area,afoo)
        datExtract <- foo
      }
    }
    itc_spectra <- extract(clip,datExtract)
    write_csv(as.data.frame(itc_spectra), paste(wd,'outputs/Spectra/',entryID, ".csv", sep=""))
  }else{
    #extract infor from that clump
    datExtract<- itc_shp[which(itc_shp$DN==datExtract$DN),]
    for(pp  in 1:length(datExtract)){
      #check which polygon at the same distance has the greatest area
      foo <- datExtract[pp,]
      afoo <- over(centDat, foo)
      if(!is.na(afoo)){
        ffoo <- foo
      }
    }
    datExtract<-ffoo
    itc_spectra <- extract(clip, datExtract)
  }
  writeOGR(datExtract, paste(wd,'outputs/itcID/', sep = ""), entryID,  driver="ESRI Shapefile", overwrite_layer = T)
  write_csv(as.data.frame(itc_spectra), paste(wd,'outputs/Spectra/',entryID, ".csv", sep=""))
  
}



maskTiff <- function(x, epsg, token){
  library(raster)
  library(rgeos)
  library(rgdal)
  library(readr)
  #256196.6 4108745
  clip.xmin <- (x$easting) - buffer
  clip.ymin <- (x$northing) - buffer
  clip.xmax <- clip.xmin + 2*buffer
  clip.ymax <- clip.ymin + 2*buffer
  entryID <- paste(token, x$siteID, x$taxaID,  x$treeID, sep="_")
  clip<- stack(paste(entryID, '.tif', sep=""))
  itc <- paste(as.integer(clip.xmax/1000)*1000, as.integer(clip.ymax/1000)*1000, 'silva', sep="_")
  itc_shp <-readOGR("./itcShp/", itc, stringsAsFactors = FALSE)
  proj4string(itc_shp) <-  CRS(paste("+init=epsg:", epsg, sep=""))
  proj4string(clip) <-  CRS(paste("+init=epsg:", epsg, sep=""))
  centDat <- data.frame(long = clip.xmin+buffer, lat = clip.ymin+buffer, name = entryID)
  coordinates(centDat) <- ~ long + lat
  proj4string(centDat) <- CRS(paste("+init=epsg:", epsg, sep=""))
  datExtract <- over(centDat, itc_shp)
  
  if(is.na(datExtract)){
    print("not really")
    dist <- gDistance(centDat, itc_shp,byid=TRUE)
    datDist <-which(dist==min(dist))
    area = 0
    for(pp  in 1:length(datDist)){
      #check which polygon at the same distance has the greatest area
      foo <- itc_shp[datDist[pp],]
      afoo <-  unlist(lapply(foo@polygons, function(x) slot(x, "area")))
      if(afoo>=area){
        area = max(area,afoo)
        datExtract <- foo
      }
    }
  }else{
    #extract infor from that clump
    datExtract<- itc_shp[which(itc_shp$DN==datExtract$DN),]
    for(pp  in 1:length(datExtract)){
      #check which polygon at the same distance has the greatest area
      foo <- datExtract[pp,]
      afoo <- over(centDat, foo)
      if(!is.na(afoo)){
        ffoo <- foo
      }
    }
    datExtract<-ffoo
  }
  msk <- mask(clip[[1]], datExtract)
  onlyIn <- msk * clip
  writeRaster(onlyIn, paste('./output/masked/', entryID, '.tif',sep=""), overwrite = T)
  
}

convert_stei <- function(dat){
  library(rgdal)
  what_to_keep <- which(dat$easting > 270000)
  tmp_keep <- dat[what_to_keep,]
  dat$utmZone <- "15N"
  dat$siteID <- "CHEQ"
  coordinates(dat) <- c("easting", "northing")
  proj4string(dat) <- CRS("+init=epsg:32616") # WGS 84
  dat <- spTransform(dat, CRS("+init=epsg:32615"))
  new_dat <- cbind(dat@data, dat@coords)
  new_dat[what_to_keep,] <- tmp_keep
  return(new_dat)
}

