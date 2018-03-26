
predict_tile <- function(bb, plsglm, hsp, w.daic){
  #output.daic = output.sd.daic = rep(0,dim(hsp)[1])
  out = list()
  #pred.val.data <- predict.withsd(plsglm[[bb]]$mod, newdata = hsp, ncomp=plsglm[[bb]]$ncomp, type='response',  se.fit = T)
  pred.val.data <- predict.withsd(plsglm[[bb]]$mod, newdata = hsp, ncomp=plsglm[[bb]]$ncomp, type='response',  wt = plsglm[[bb]]$mod$FinalModel$weights)
  out$output.daic <-  pred.val.data[,1] * w.daic[bb]
  out$pi.upper <- pred.val.data[,3]* w.daic[bb]
  out$pi.lower <- pred.val.data[,2]* w.daic[bb]

  #out$output.sd.daic <- as.vector(pred.val.data$se.fit) * w.daic[bb]
  out$output <- as.vector(pred.val.data[,1])
  return(out)
}

traitsToshp <- function(x.size, goodPix.pos, extnt, epsg,itc_shp, trName, inc_ave){
  #allocate new rasters
  #library(velox)
  tpt1 = tpt2 = tpt3 = tpt4 <- matrix(NA, nrow = x.size[1]*x.size[2])
  tpt1[goodPix.pos,] <- inc_ave[,1]
  tpt2[goodPix.pos] <- inc_ave[,2]
  tpt3[goodPix.pos] <- inc_ave[,3]
  tpt4[goodPix.pos] <- inc_ave[,4]


  dim(tpt1) <- c(x.size[1],x.size[2])
  dim(tpt2) <- c(x.size[1],x.size[2])
  dim(tpt3) <- c(x.size[1],x.size[2])
  dim(tpt4) <- c(x.size[1],x.size[2])


  r.temp1 <-raster(tpt1, xmn=extnt[1], xmx=extnt[2],
                   ymn=extnt[3], ymx=extnt[4], crs=CRS(paste("+init=epsg:", epsg, sep="")))

  r.temp2 <-raster(tpt2, xmn=extnt[1], xmx=extnt[2],
                   ymn=extnt[3], ymx=extnt[4], crs=CRS(paste("+init=epsg:", epsg, sep="")))

  r.temp3 <-raster(tpt3, xmn=extnt[1], xmx=extnt[2],
                   ymn=extnt[3], ymx=extnt[4], crs=CRS(paste("+init=epsg:", epsg, sep="")))
  r.temp4 <-raster(tpt4, xmn=extnt[1], xmx=extnt[2],
                   ymn=extnt[3], ymx=extnt[4], crs=CRS(paste("+init=epsg:", epsg, sep="")))

  r.temp <- raster::stack(r.temp1, r.temp2, r.temp3, r.temp4)
  rm (r.temp1, r.temp2, r.temp3, r.temp4, tpt1, tpt2, tpt3, tpt4)
  #r.temp <- mask(r.temp, itc_shp)
  # writeRaster(r.ave, paste("./outputs/average", j,i, sep="_"), overwrite=TRUE, format='GTiff')
  itc_shp <- crop(itc_shp,extent(r.temp))
  poly.id <-  sapply(slot(itc_shp, "polygons"), function(x) slot(x, "ID"))
  #tmp.val <- extract(r.temp, crop(itc_shp,extent(r.temp)))
  tmp.val <- raster::extract(r.temp, itc_shp)
  tmp.val <- Reduce(rbind, lapply(tmp.val,  function(x) apply(x,2,mean, na.rm=TRUE)))

  df.N = as.data.frame(cbind(tmp.val, poly.id))
  colnames(df.N) <- c(trName, "ID")
  #rownames(df.N) <- seq(0,dim(tmp.val)[1]-1)
  rownames(df.N) <- poly.id
  foo <- spCbind(itc_shp, df.N)
  # foo@data <- transform(foo@data, val = as.numeric(val),
  #                       ID = as.numeric(ID))
  # colnames(foo@data)[which(colnames(foo@data)=="ave")]= trName
  #final <- foo[!is.na(unlist(foo@data$ave)),]
  return(foo)
}
pixelRaster<- function(inc_ave,x.size, goodPix.pos, extnt, epsg, trName){
  #allocate new rasters
  template <- matrix(NA, nrow = x.size[1]*x.size[2], ncol = 1)
  template[goodPix.pos] <- inc_ave
  dim(template) <- c(x.size[1],x.size[2])
  r.temp <-raster(template, xmn=extnt[1], xmx=extnt[2],
                  ymn=extnt[3], ymx=extnt[4], crs=CRS(paste("+init=epsg:", epsg, sep="")))
  writeRaster(r.temp, filename=paste(trName, "tif", sep="."), format="GTiff",overwrite=TRUE )

}

normalizeImg<-function(hsp, scaled = F){
  if(scaled == F) {hsp <- hsp /10000}
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
  # Write vector normalized spectra to CSV
  return(normDF)
}

getSpatialRegression <- function(NeonSite = "OSBS", nm = c("LMA_g.m2","N_pct", "P_pct"), hsp = NULL, epsg = NULL){
  library(rgdal)
  library("tsensembler")
  library(plsRglm)
  library(maptools)
  library(data.table)
  x.size <- dim(hsp)
  extnt <- extent(hsp)
  hsp <- as.array(hsp)
  dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
  #hsp <- as.data.frame(hsp)
  hsp <- as.data.table(hsp)
  hsp <- normalizeImg(hsp)
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
  sp_path <- './spatialPhase'
  itc_shp <-readOGR(paste(sp_path, "segmentation", NeonSite, sep="/"), itc, stringsAsFactors = FALSE)
  proj4string(itc_shp) <-  CRS(paste("+init=epsg:", epsg, sep=""))
  #itc_shp <- spTransform(itc_shp, CRS(paste("+init=epsg:", epsg, sep="")))
  for(j in nm){
    if(j == "P_pct"){
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

    for(bb in 1: length(plsglm)){
      mod.aic[bb] <-plsglm[[bb]]$score$aic
    }
    delta.aic <- mod.aic - min(mod.aic)
    weights <- softmax(-0.5*delta.aic)

    output.daic = output.sd.daic = rep(0,dim(hsp)[1])
    w.daic <- rep(0, length(weights))
    w.daic <- weights
    rm(output, output.sd)
    for(bb in 1:length(w.daic)){
      pls.mod.train <- plsglm[[bb]]$mod
      optim.ncomps <- plsglm[[bb]]$ncomp
      pred.val.data <- predict(pls.mod.train, newdata = hsp, ncomp=optim.ncomps, type='response',  se.fit = T)

      output.daic <- output.daic + pred.val.data$fit * w.daic[bb]
      output.sd.daic <-  output.sd.daic + as.vector(pred.val.data$se.fit) * weights[bb]
      #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
      if(!exists("output")){
        output <- as.vector(pred.val.data$fit)
        #output.sd <-  as.vector(pred.val.data$se.fit)
      }else{
        output <- rbind(output,  as.vector(pred.val.data$fit))
        #output.sd <- rbind(output.sd,  as.vector(pred.val.data$se.fit))
      }
    }
    rm(plsglm)
    output.daic <- as.data.frame(as.matrix(output.daic))
    output.sd.daic <- as.data.frame(as.matrix(output.sd.daic))

    colnames(output.daic) <-  "yhat"
    colnames(output.sd.daic) <-  "yhat_unc"
    output.daic <- output.daic[!is.nan(output.daic$yhat), ]
    output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]

    pdf(paste("hist", j, ".pdf", sep=""))
    hist(output)
    dev.off()

    inc_ave <-apply(output,2, function(x) mean(x, na.rm=T))
    inc_ave_sd <- apply(output,2, function(x) sd(x, na.rm=T))
    #pix.ave <- pixelRaster(inc_ave,x.size, goodPix.pos, extnt, epsg, paste(j, "ave", sep="_"))

    itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, epsg,itc_shp, paste(j, "ave", sep="_"), output.daic)
    itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, epsg,itc_shp, paste(j, "ster", sep="_"),output.sd.daic)
    itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, epsg,itc_shp,paste(j, "md_sd", sep="_"), inc_ave_sd)
  }
  return(itc_shp)
}

gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath="./src/polygonize/", readpoly=TRUE, quiet=TRUE) {
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
    writeRaster(x, {f <- tempfile(fileext='.tif', tmpdir = "../tmp")})
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

crownIT <- function(las_id = NULL, epsg=NULL, pt = NULL, tileID = NULL, out_path = NULL,  NeonSite = "OSBS", method = 'silva',
                    max_cr_factor = 0.6, exclusion = 0.5, mv = 7, minh = 7){
  library(raster)
  library(lidR)
  source("./src/function_spatial.R")
  #pt <-  paste("./spatialPhase/", NeonSite, "/ClassifiedPointCloud/", sep = "")
  i <- list.files(pt, pattern = las_id)

  las = readLAS(paste(pt, i, sep=""))

  # normalization
  lasnormalize(las, method = "knnidw", k = 10L)

  # compute a canopy image
  chm = grid_canopy(las, res = 0.5, subcircle = 0.2, na.fill = "knnidw", k = 4, p = 1)
  chm = as.raster(chm)
  kernel = matrix(1,3,3)
  chm = raster::focal(chm, w = kernel, fun = mean)
  chm = raster::focal(chm, w = kernel, fun = mean)

  ttops = tree_detection(chm, mv, minh)
  crowns <-lastrees_silva(las, chm, ttops, max_cr_factor = max_cr_factor, exclusion = exclusion, extra = T)
  x <- gdal_polygonizeR(crowns, pypath="/ufrc/ewhite/s.marconi/Marconi2018/src/polygonize/")
  writeOGR(x, dsn=out_path, paste(tileID, "silva", sep="_"), overwrite_layer = T, check_exists = T, driver="ESRI Shapefile")
}
