
loadITC <- function(in.dir, proj){
  path.vect <- paste(in.dir, "Geofiles/vector/", sep="")
  tree.tg <- readOGR(path.vect, "Final Tagged Trees", stringsAsFactors=F)
  tree.untg <- readOGR(path.vect, "Final Untagged Trees", stringsAsFactors=F)
  plot.cnt <- read.centroids("OSBS_diversity_plot_centroids")
  utm.tg <- spTransform(tree.tg, CRS(proj))
  utm.untg <- spTransform(tree.untg, CRS(proj))
  return(list(tg = utm.tg, untg  = tree.untg))
}

loadIMG <- function(path, proj, img.lb = 1){
  path.raster <- paste(in.dir, "Geofiles/RSdata/", sep="")
  #chm.sample <- raster(paste(path.raster, "LiDAR/", 'plot_', img.lb, "_chm.tif", sep = ""))
  hs.sample <- brick(paste(path.raster, "hs/", img.lb, "_hyper.tif", sep = ""))
  return(list(hsp = hs.sample))
}

spectra.profile <-function(f, xi, yi){
  #first call required libraries
  library(rhdf5)
  library(plyr)
  library(ggplot2)
  
  pre.path <- getwd()
  setwd("/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D3/OSBS/2014/OSBS_L1/OSBS_Spectrometer/Reflectance/")
  spInfo <- h5readAttributes(f,"map info")
  
  #read in the wavelength information from the HDF5 file
  wavelengths<- h5read(f,"wavelength")
  #extract Some Spectra from a single pixel
  aPixel<- h5read(f,"Reflectance",index=list(xi,yi,NULL))
  #reshape the data and turn into dataframe
  b <- adply(aPixel,c(3))
  #create clean data frame
  aPixeldf <- b[2]
  #add wavelength data to matrix
  aPixeldf$Wavelength <- wavelengths
  reflInfo <- h5readAttributes(f,"Reflectance")
  #add scaled data column to DF
  aPixeldf$scaled <- (aPixeldf$V1/reflInfo$`Scale Factor`)
  aPixeldf$scaled[aPixeldf$scaled>1] <- NA
  #make nice column names
  names(aPixeldf) <- c('Reflectance','Wavelength','ScaledReflectance')
  
  #plot spectral profile
  sp.plot <- plot(x=aPixeldf$Wavelength, 
                  y=aPixeldf$ScaledReflectance,
                  xlab="Wavelength (nm)",
                  ylab="Reflectance")
  
  return(list(prof.frame = aPixeldf, prof.plot = sp.plot))
}

flight.path <- function(plot.ext, crs.plot){
  #first call required libraries
  library(rhdf5)
  library(plyr)
  library(rgeos)
  library(neonAOP)
  
  #walk trough the h5 to find the right image
  files <- read_csv("inputs/OSBS_h5.csv", col_types = cols(NIS1_20140507_143910_atmcor.h5 = col_character()))
  pre.path <- getwd()
  #TODO this becomes an argument as soon as you go for more than one site
  setwd("/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D3/OSBS/2014/OSBS_L1/OSBS_Spectrometer/Reflectance/")
  #a polygon object with several polygons inside. refer to them by giving the ith number
  #TODO: for now I'll focus on plot OSBS_001. Change in future to work with all plots in batch 
  plot.ID <- polys.df@polygons[[27]]@ID
  plot.ext <- polys.df@polygons[[27]]@Polygons[[1]]@coords
  plot.ext <- cbind(plot.ext[,2], plot.ext[,1])
  p <- Polygon(plot.ext)
  p <- Polygons(list(p), 1)
  plot.ext <- SpatialPolygons(list(p))
  proj4string(plot.ext) <- CRS(crs.plot)  
  
  j = 0
  is.included <- rep(0,dim(files)[1])
  for(i in t(files)){
    j <- j + 1
    i <- t(files)[5]
    h5ls(i,all=T) 
    
    #calculate ith object extent
    dims <- get_data_dims(i)
    
    ras.ext <- get.extent(i, dims)
    #see if the jth plot is contained in the ith flight
    ei <- as(extent(rasExt), "SpatialPolygons")
    proj4string(ei) <- CRS(crs.plot)  
    
    #check who contains the plot
    if (gContainsProperly(plot.ext, ei)) {
      is.included[j] <- 1
    }
    print(is.included[j])
  }
  fl.paths <- files[as.logical(is.included)]
  setwd(pre.path)
}
get.plot.extent <- function(plots, buffersize){
  centroids <- coordinates(plots)
  xPlus <- centroids[,2]+buffersize
  yPlus <- centroids[,1]+buffersize
  xMinus <- centroids[,2]-buffersize
  yMinus <- centroids[,1]-buffersize  
  ext.pts <- cbind(xMinus, yMinus, xMinus, yPlus, xPlus, yPlus, xPlus, yMinus, xMinus, yMinus)
  ID = plots$Plot_ID
  
  #credits: http://stackoverflow.com/questions/26620373/spatialpolygons-creating-a-set-of-polygons-in-r-from-coordinates
  polys <- SpatialPolygons(mapply(function(poly, id) {
    xy <- matrix(poly, ncol=2, byrow=TRUE)
    Polygons(list(Polygon(xy)), ID=id)
  }, split(ext.pts, row(ext.pts)), ID))
  
  # Create SPDF
  polys.df <- SpatialPolygonsDataFrame(polys, data.frame(id=ID, row.names=ID))
  plot(polys.df, col=rainbow(50, alpha=0.5))
  return(polys.df)
}

get.extent <- function(f, dims){
  h5.ext <- h5read(i, "map info")
  #Extract each element of the map info information 
  #so we can extract the lower left hand corner coordinates.
  h5.ext<-unlist(strsplit(h5.ext, ","))
  res <- as.numeric(c(h5.ext[2], h5.ext[3]))
  xMin <- as.numeric(h5.ext[4]) 
  yMax <- as.numeric(h5.ext[5])
  xMax <- (xMin + (dims[1])*res[1])
  yMin <- (yMax - (dims[2])*res[2]) 
  
  #h5 extent
  rasExt <- extent(xMin,xMax,yMin,yMax)
  return(rasExt)
}

read.centroids <- function(file){
  plot.centr <- readOGR("./inputs/Geofiles/vector", file, stringsAsFactors = FALSE)
  return(plot.centr)
}

filter.bad <- function(img, bb, img.lb){
  img <- dropLayer(img, bb)
  writeRaster(img, paste('./outputs/filtered/', img.lb, "_bb350_2512.tif", sep= ""), format = 'GTiff')
  return(brick(paste('./outputs/filtered/', img.lb, "_bb350_2512.tif", sep= "")))
}

itc.plot <- function(trees, rasters, just.mnf, hybrid = NULL, f=FALSE){
  
  # setwd("/Users/sergiomarconi/Documents/PhD/Projects/scalingHiper/Outputs/")
  #pdf(paste('./outputs/crowns/', f,'_itc.pdf',sep=""),height=12,width=8)
  trees$untg <- spTransform(trees$untg, CRS("+init=epsg:32617"))
  trees$tg <- spTransform(trees$tg, CRS("+init=epsg:32617"))
  just.mnf <- spTransform(just.mnf, CRS("+init=epsg:32617"))
  
  eval.extent <-  union(trees$tg, trees$untg)
  eval.extent <- eval.extent[just.mnf,]
  tst.hsp <- just.mnf[eval.extent,]
  rasters <- crop(rasters, tst.hsp)
  
  if(!f){
    plot(rasters, col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
  }else{
    plotRGB(brick(rasters), r=1,g=2,b=3, stretch="Lin")
  }
  plot(tst.hsp,axes=T, border="blue", add=TRUE, lwd = 2)
  plot(eval.extent,axes=T, border="red3", add=TRUE, lwd = 2)
  if(!is.null(hybrid)){
    hybrid <- spTransform(hybrid, CRS("+init=epsg:32617"))
    hybrid <- hybrid[eval.extent,]
    plot(hybrid,axes=T, border="green", add=TRUE, lwd = 1.5)
  }
  #dev.off()
}


create.geoJSON <- function(tmp.poly, out){
  tmp <- (tmp.poly@Polygons[[1]]@coords)
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c("x", "y")
  coordinates(tmp) = c("x", "y")
  tmp <- Polygon(tmp)
  tmp <- Polygons(list(tmp), 1)
  tmp <- SpatialPolygons(list(tmp))
  proj4string(tmp) <- proj
  tmp <- SpatialPolygonsDataFrame(tmp,data)
  writeOGR(obj = tmp, out,  layer = 'plot_001', driver="GeoJSON")
}

myCluster <-function(Xt,nC, k, ext, proj){
  require(cluster)
  require(raster)
  X <- as.array(Xt)
  X <- X[,,1:nC]
  x.size <- dim(X)
  dim(X) <- c(x.size[1] * x.size[2], x.size[3])
  clust.X <- clara(X, k) 
  r <- as.matrix(clust.X$clustering)
  dim(r) <- c(x.size[1], x.size[2])
  r <-raster(r)
  extent(r) <- ext
  crs(r) <- proj
  return(r)
}

myHCluster <-function(Xt,nC, k, ext, proj){
  Xt <- hyp.clust
  require(cluster)
  require(raster)
  X <- as.array(Xt)
  X <- X[,,1:nC]
  x.size <- dim(X)
  dim(X) <- c(x.size[1] * x.size[2], x.size[3])
  clust.X <- hclust(dist(X))
  ct <- cutree(clust.X, k)
  r <- as.matrix(ct)
  dim(r) <- c(x.size[1], x.size[2])
  r <-raster(r)
  extent(r) <- ext
  crs(r) <- proj
  return(r)
}

itcHhps <- function(hyp.clust, nC, lasITC, par=T, fast=T, verbose = F){
  library(dynamicTreeCut)
  library(cluster)
  for( i in 1:length(lasITC@polygons)){
    p.id <- hyp.clust
    if(verbose){
      plot(lasITC)
      plot(SpatialPolygons(list(lasITC@polygons[[i]])), border = 'red',add=T)
    }
    p.id <- r.mask <- crop(hyp.clust, (SpatialPolygons(list(lasITC@polygons[[i]]))))
    if(verbose){
      plot(p.id[[1]])
    }
    ext <- extent(p.id$pca_plot1.1)
    proj <- p.id@crs
    p.id <- as.array(p.id)
    p.id <- p.id[,,1:nC]
    p.dim <- dim(p.id)
    dim(p.id) <- c(p.dim[1] * p.dim[2], p.dim[3])
    # Ward Hierarchical Clustering
    d <- dist(p.id, method = "manhattan") # distance matrix
    #rdist.w.na
    #d <- rdist.w.na(as.numeric(p.id), as.numeric(p.id)) # distance matrix
    
    fit <- hclust(d, method="ward.D") 
    
    #use k-nn given number of k from hierarchical clustering p-val
    #groups <- clara(p.id, k=get.k(p.id, par, fast))
    #groups <- as.matrix(groups$clustering)
    
    #use the same clustering given by hierarchical trees
    #groups <- cutree(fit, k=get.k(p.id, par, fast)) # cut tree into n clusters given a = 0.95
    d <- daisy(p.id, metric="manhattan")
    groups <- cutreeDynamic(fit, minClusterSize=10, distM = as.matrix(d))
    
    dim(groups) <- c(p.dim[1], p.dim[2])
    dim(p.id) <- c(p.dim[1], p.dim[2], p.dim[3])
    groups[p.id[,,1] ==-1L] <- NA
    if(max(groups, na.rm = TRUE)==0){groups <-groups + 1}
    #plot(raster(groups))
    r <- raster(groups)
    extent(r) <- ext
    crs(r) <- proj
    mask <- rasterize(SpatialPolygons(list(lasITC@polygons[[i]])), r.mask)
    #plot(mask)
    extent(mask) <- ext
    r <- mask(r, mask)
    r2 <- gapfill(r, mask, ext)
    # plot(r2)
    # r2 <- gapfill(r2, mask, ext)
    
    if(length(r2)>0){
      r <- rasterToPolygons(r2, dissolve = T, n = 16)
      if(i ==1){
        poly.list <- r
      }else{
        poly.list <- rbind(poly.list, r)
      }
    }
  }
  return(poly.list)
}

get.k <- function(p.id, par, fast){
  library(pvclust)
  if(fast){
    fit <- pvclust(t(p.id), method.hclust="ward.D",
                   method.dist="manhattan", parallel = par, r = 1)
    pp <-pvpick(fit, pv="bp", alpha =.2)
  }else{
    fit <- pvclust(t(p.id), method.hclust="ward.D",
                   method.dist="manhattan", parallel = par, r = c(0.4,0.8))
    pp <-pvpick(fit, alpha =.95)
  }
  return(length(pp$clusters))
}




gapfill <- function(r,mask, ext){
  # http://stackoverflow.com/questions/24465627/clump-raster-values-depending-on-class-attribute
  # extend raster, otherwise left and right edges are 'touching'
  r <- extend(r, c(1,1))
  mask <- extend(mask, c(1,1))
  # get al unique class values in the raster
  clVal <- unique(r)
  # remove '0' (background)
  clVal <- clVal[!clVal==0]
  # create a 1-value raster, to be filled in with NA's
  r.NA <- setValues(raster(r), 1)
  #plot(r.NA)
  # set background values to NA
  r.NA <- mask(r, mask, updatevalue = 0)
  r.NA[r==0]<- NA
  #plot(r.NA)
  # loop over all unique class values
  for (ii in clVal) {
    # create & fill in class raster
    r.class <- setValues(raster(r), NA)
    r.class[r == ii]<- 1
    # clump class raster
    clp <- clump(r.class, directions = 4)
    # calculate frequency of each clump/patch
    cl.freq <- as.data.frame(freq(clp))
    # store clump ID's with frequency 1
    rmID <- cl.freq$value[which(cl.freq$count <= 4)]
    # assign NA to all clumps whose ID's have frequency 1
    r.NA[clp %in% rmID] <- 0
  } 
  
  # multiply original raster by the NA raster
  r2 <- r * r.NA
  #r2 <- crop(r2, ext)
  
  # Create focal weights matrix.
  nbrhd <- matrix(c(0,1,0,1,1,1,0,1,0), ncol = 3, nrow = 3)
  r3 <- focal(x=r2, w=nbrhd, fun = fill.na)
  r3 <- r3 + r2
  r3 <- crop(r3, ext)
  #re-remove background 0 
  # r.NA <- setValues(raster(r3), 1)
  # r.NA <- mask(r3, r, updatevalue = 0)
  # r.NA[r==0]<- NA
  # plot(r.NA)
  return(r3)
}
# idea from http://stackoverflow.com/questions/27906111/r-focal-raster-package-how-to-apply-filter-to-subsets-of-background-data
fill.na <- function(x) {
  center <- x[ceiling(length(x)/2)]
  if(!is.na(center)){
    if( center==0 ) {
      fr <-x
      fr <- fr[fr!= 0]
      fr <- fr[!is.na(fr)]
      fr <- as.integer(names(which(table(fr) == max(table(fr)))[1]))
      if(length(fr)>0)
        return((fr))
    }
  }
}

# http://stackoverflow.com/questions/22785475/does-a-function-like-dist-rdist-exist-which-handles-nas
rdist.w.na <- function(X,Y)
{
  library(pdist)
  
  if (!is.matrix(X)) 
    X = as.matrix(X)
  if (!is.matrix(Y)) 
    Y = as.matrix(Y)
  distances <- matrix(pdist(X,Y)@dist, ncol=nrow(X), byrow = TRUE)
  #count NAs
  na.count <- sapply(1:nrow(X),function(i){rowSums(is.na(Y) | is.na(X[i,]))})
  #scaling to number of cols
  distances * sqrt(ncol(X)/(ncol(X) - na.count))
}
