#--- traitsFromITC
traitsFromITC <- function(args, nm = c("LMA_g.m2", "N_pct")){
	#main function to get traits from remote sensing data, applying the models developed in phase 1
  library(raster)
  source(paste(getwd(), '/src/function_spatial.R', sep=""))
  library(rgdal)
  library("tsensembler")
  library(plsRglm)
  library(maptools)
  library(data.table)
  rasterOptions(tmpdir="./tmp")
  args <- data.frame(t(args), stringsAsFactors = F)
  args <- transform(args, xx = as.integer(xx),yy = as.integer(yy),
                    xmin = as.integer(xmin),xmax = as.integer(xmax), scaled = as.logical(scaled),
                    ymin = as.integer(ymin),ymax = as.integer(ymax), epsg = as.integer(epsg))
  (print(args))
  #python much faster: produce the ith hiperspectral tile
  system2("python", args=(sprintf('"%1$s" "%2$s" "%3$s" "%4$i" "%5$i"',"./src/stripeToRaster.py", 
                                  args$f, args$path, args$xx, args$yy)))
  f_nm_1 = max(as.integer(args$xmin/1000)*1000 +1000*as.integer(args$xx), as.integer(args$xmin))
  f_nm_2 = min(as.integer(args$ymin/1000)*1000 +1000*(as.integer(args$yy)+1), as.integer(args$ymax))
  hsp.str <- list.files("./tmp/", pattern = paste(f_nm_1, f_nm_2, sep =""))
  hsp <- brick(paste("./tmp/",hsp.str, sep=""))
  NeonSite <- unlist(strsplit(args$f, "_"))[3]
  f_nm_1 = as.integer(args$xmin/1000)*1000 +1000*as.integer(args$xx)
  f_nm_2 = as.integer(args$ymin/1000)*1000 +1000*(as.integer(args$yy))
  # #tryCatch({ 
  itc_path <-paste("./spatialPhase/segmentation", NeonSite, sep="/")
  itc <- list.files(itc_path, pattern = paste(f_nm_1, f_nm_2, sep ="_"))[1] # paste(f_nm_1, f_nm_2, sep ="_"))
  if(!is.na(itc) && !is.nan(itc) && length(itc)!=0){
    itc <- paste(c(head(unlist(strsplit(itc, "_")),-1), args$itc_method), collapse="_")
    warning(print(itc))
    #apply method on tif file
    
    x.size <- dim(hsp)
    extnt <- extent(hsp)
    hsp <- as.array(hsp)
    dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
    hsp <- as.data.table(hsp)
    if(args$scaled == F) {hsp <- hsp /10000}
    # Set bad bands to zero
    ndvi <- (hsp$band_90- hsp$band_58)/(hsp$band_58 + hsp$band_90)
    nir860 <- (hsp$band_96 + hsp$band_97)/2
    hsp[which(ndvi < 0.7 | nir860 < .3),]=NA
    hsp[hsp <0 |hsp>1]=NA
    # Vector normalize spectra
    normMat=sqrt(apply(hsp^2,FUN=sum,MAR=1, na.rm=TRUE))
    normMat=matrix(data=rep(normMat,ncol(hsp)),ncol=ncol(hsp))
    normMat=hsp/normMat
    hsp = as.data.table(normMat)
    rm(ndvi, nir860, normMat)
    
    goodPix.pos <- which(!is.na(hsp[,1]))
    hsp=as.matrix(hsp)
    site_bands <- matrix(0, dim(hsp)[1],2)
    #change in the future
    if(NeonSite=="OSBS"){
      site_bands[,1]<-1
    }else{
      site_bands[,2]<-1
    }
    hsp<- cbind(site_bands, hsp)
    
    #end of normalization of the spectra
    
    sp_path <- './spatialPhase'
    itc_shp <-readOGR(paste(sp_path, "segmentation", NeonSite, sep="/"), itc, stringsAsFactors = FALSE)
    proj4string(itc_shp) <-  CRS(paste("+init=epsg:", args$epsg, sep=""))
    warning(print("for cicle"))
    for(j in nm){
      if(j == "P_pct"){
        warning(print(j))
        
        hsp[hsp ==0] <- 0.000001
        hsp<- hsp[, -c(1:2)]
        hsp <-t(diff(t(log(hsp)),differences=1, lag=3))
        hsp<- cbind(site_bands, hsp)
      }
      hsp <- hsp[complete.cases(hsp), ]
      
      mod_dir = "./ModelBuild/ALL/"
      load(file = paste(mod_dir, 'plsglm_',j, sep="" ))
      
      mod.r2=rep(0,length(plsglm))
      mod.aic=rep(0,length(plsglm))
      #if want to use less than 100
      mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 5)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 5)])
      plsglm <- plsglm[mask[,1]]
      
      for(bb in 1: length(plsglm)){
        mod.aic[bb] <-plsglm[[bb]]$score$aic
      }
      delta.aic <- mod.aic - min(mod.aic)
      weights <- softmax(-0.5*delta.aic)
      
      output.daic = output.sd.daic = rep(0,dim(hsp)[1])
      w.daic <- rep(0, length(weights))
      w.daic <- weights
      rm(output, output.sd)
      for(bb in 1:10){#length(w.daic)){
        pred.val.data <- predict(plsglm[[bb]]$mod, newdata = hsp, ncomp=plsglm[[bb]]$ncomp, type='response',  se.fit = T)
        output.daic <- output.daic + pred.val.data$fit * w.daic[bb]
        output.sd.daic <-  output.sd.daic + as.vector(pred.val.data$se.fit) * weights[bb]
        #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
        if(!exists("output")){
          output <- as.vector(pred.val.data$fit)
        }else{
          output <- rbind(output,  as.vector(pred.val.data$fit))
        }
      }
      rm(plsglm)
      output.daic <- as.data.frame(as.matrix(output.daic))
      output.sd.daic <- as.data.frame(as.matrix(output.sd.daic))
      
      colnames(output.daic) <-  "yhat"
      colnames(output.sd.daic) <-  "yhat_unc"
      output.daic <- output.daic[!is.nan(output.daic$yhat), ]
      output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]
      out_dir <- args$path
      pdf(paste(paste(args$path, "products/", sep="/"), args$xmin, args$ymax, j, "hist.pdf", sep="_"))
      hist(output)
      dev.off()
      
      inc_ave <-apply(output,2, function(x) mean(x, na.rm=T))
      inc_ave_sd <- apply(output,2, function(x) sd(x, na.rm=T))
      #pix.ave <- pixelRaster(inc_ave,x.size, goodPix.pos, extnt, epsg, paste(j, "ave", sep="_"))
      warning(print("traitsToshp"))
      itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, args$epsg,itc_shp, paste(j, "ave", sep="_"), output.daic)
      itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, args$epsg,itc_shp, paste(j, "ster", sep="_"),output.sd.daic)
      itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, args$epsg,itc_shp,paste(j, "md_sd", sep="_"), inc_ave_sd)
    }
    #final <- getSpatialRegression(NeonSite = NeonSite, nm = c("LMA_g.m2", "P_pct"), hsp =hsp, epsg = epsg)
    
    print("getSpatialRegression regr ok")
    itc_shp <- itc_shp[!is.na(unlist(itc_shp@data[2])),]
    itc_shp <- itc_shp[!is.nan(unlist(itc_shp@data[2])),]
    itc_shp@data <- data.frame(sapply(itc_shp@data, function(x) as.numeric(as.character(x))))
    writeOGR(itc_shp, paste(args$path, "products/", sep=""), paste(c(unlist(strsplit(itc, "_"))[1:5],args$xmin,args$ymin, "shp_product"), collapse="_"), driver="ESRI Shapefile",overwrite_layer=T)
    
    warning(print(paste(itc, "ok!")))
    #delete tif
    system(paste("rm ", "./tmp/", hsp.str, sep=""))
    warning("did it!")
  }else{
    warning(paste("WARNING!:",args$path, args$f, args$xx, args$yy,"not found"))
  }
 }
 
 
 ## Define the function
gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                  pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}

crownIT <- function(in.dir = '.', method = 'silva', cores = 1, epsg=NULL){
  library(foreach)
  library(doParallel)
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  setwd("./OSBS")
  pt <-  "./inputs/Geofiles/RSData/Classified_point_cloud/"
  listFiles <- list.files(pt, pattern = "*.laz")[-(1:30)]
  i = listFiles[1]
  results <- foreach(i = listFiles) %dopar% {
  for(i in listFiles) { 
    library(raster)
    library(lidR)
    #source("../src/polygonize.R")
	    source("../src/sp_functions.R")
		
    pt <-  "./inputs/Geofiles/RSData/Classified_point_cloud/"
    las = readLAS(paste(pt, i, sep=""))

    # normalization
    lasnormalize(las, method = "knnidw", k = 10L)
    #lasnormalize(las, method = "kriging")
    
    # compute a canopy image
    chm = grid_canopy(las, res = 0.5, subcircle = 0.2, na.fill = "knnidw", k = 4)
    chm = as.raster(chm)
    kernel = matrix(1,3,3)
    chm = raster::focal(chm, w = kernel, fun = mean)
    chm = raster::focal(chm, w = kernel, fun = mean)
    
	if(method=='silva'){
    #silva 2016
    ttops = tree_detection(chm, 5, 2)
    crowns <-lastrees_silva(las, chm, ttops, max_cr_factor = 0.6, exclusion = 0.3, extra = T)
    x <- gdal_polygonizeR(crowns)
    writeOGR(x, dsn="../itcShp/", paste(i, "silva", sep="_"), overwrite_layer = T, check_exists = T, driver="ESRI Shapefile")
    } else if(method == "watershed"){
	# Watershed
    crowns <-lastrees_watershed(las, chm, th_tree = 5,extra = T)
    x <- gdal_polygonizeR(crowns)
    writeOGR(x, dsn="../itcShp/", paste(i, "watershed", sep="_"), overwrite_layer = T, check_exists = T, driver="ESRI Shapefile")
	}ελσε ιφ(μετηοδ = 'dalponte'){
		# #dalponte
    ttops = tree_detection(chm, 5, 2)
    crowns <- lastrees_dalponte(las, chm, ttops, max_cr = 7, extra = T)
    x <- gdal_polygonizeR(crowns)
    writeOGR(x, dsn="../itcShp/", paste(i, "dalponte", sep="_"), overwrite_layer = T, check_exists = T, driver="ESRI Shapefile")
	}
  }
  print(results)
  stopCluster(cl)
}
