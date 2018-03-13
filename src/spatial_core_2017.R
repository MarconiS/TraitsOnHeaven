#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)

predict_tile <- function(bb, plsglm, hsp, w.daic){
  #output.daic = output.sd.daic = rep(0,dim(hsp)[1])
  out = list()
  pred.val.data <- predict(plsglm[[bb]]$mod, newdata = hsp, ncomp=plsglm[[bb]]$ncomp, type='response',  se.fit = T)
  out$output.daic <-  pred.val.data$fit * w.daic[bb]
  out$output.sd.daic <- as.vector(pred.val.data$se.fit) * w.daic[bb]
  out$output <- as.vector(pred.val.data$fit)
  return(out)
}

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
  print(tile.args)
  
  #create layer of delineated crowns
  NeonSite <- unlist(strsplit(tile.args$f, "_"))[3]
  f_nm_1 = as.integer(tile.args$xmin)
  f_nm_2 = as.integer(tile.args$ymin)
  
  #python much faster: produce the ith hiperspectral tile
  #  tryCatch(
  #	  system2("python", args=(sprintf('"%1$s" "%2$s" "%3$s"',"./src/stripeToRaster.py",
  #                        tile.args$f, tile.args$path))), 
  #           error = function(e) {print("memory stuff") # or whatever error handling code you want
  #  })
  tryCatch({
    
    #produce the itc layer
    #    crownIT(las_id = paste(f_nm_1, f_nm_2, sep ="_"), epsg=tile.args$epsg, 
    #            NeonSite = NeonSite, method = 'silva', max_cr_factor = 0.6, exclusion = 0.5, mv = 7, minh = 7)
    #    print("ITC saved")
    hsp.str <- list.files("./spatialPhase/rasters/", pattern = paste(f_nm_1, f_nm_2, sep =""))
    hsp <- brick(paste("./spatialPhase/rasters/",hsp.str, sep=""))
    itc_path <-paste("./spatialPhase/", NeonSite, "ITCs", sep="/")
    itc <- list.files(itc_path, pattern = paste(f_nm_1, f_nm_2, sep ="_"))[1] # paste(f_nm_1, f_nm_2, sep ="_"))
    if(!is.na(itc) && !is.nan(itc) && length(itc)!=0){
      itc <- paste(c(head(unlist(strsplit(itc, "_")),-1), tile.args$itc_method), collapse="_")
      warning(print(itc))
      x.size <- dim(hsp)
      extnt <- extent(hsp)
      hsp <- as.array(hsp)
      dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
      hsp <- as.data.frame(hsp)
      #apply method on tif file
      if(mean(hsp$V90)==0 || mean(hsp$V90)==15000 || mean(hsp$V90)==-9999){
        warning('corrupted hps data')
        break
      }
      else{
        # Set bad bands to zero
        ndvi <- (hsp$V90- hsp$V58)/(hsp$V58 + hsp$V90) <0.7
        nir860 <- (hsp$V96 + hsp$V97)/20000 < 0.3
        naval = as.logical(ndvi * nir860)
        rm(ndvi, nir860)
        hsp[naval,] = NA
        #hsp[hsp <0 |hsp>1]=NA
        # Vector normalize spectra
        foo <- as.matrix(hsp) * as.matrix(hsp)
        foo <- foo /100000000
        
        #normMat<- )apply(as.matrix(hsp), function(x) x^2)
        #normMat=sqrt(apply(as.matrix(hsp)^2,FUN=sum,MAR=1, na.rm=TRUE))
        #normMat = rep(NA, dim(hsp)[1])
        normMat = apply(foo, FUN = sum, MAR = 1, na.rm=T)
        #		for(i in 1:dim(hsp)[1]){
        #			normMat[i] = (sum(foo[i,], na.rm=TRUE))
        #		}
        rm(foo)
        normMat <- sqrt(normMat)
        #normMat <- normMat + 10e-6
        hsp <- hsp /10000
        for(i in 1:dim(hsp)[2]){
          hsp[,i]=hsp[,i]/ normMat
        }
        #normMat=matrix(data=rep(normMat,ncol(hsp)),ncol=ncol(hsp))
        #normMat=hsp/normMat
        #hsp = as.data.table(hsp)
        rm(normMat)
        #if(tile.args$scaled == F) {hsp <- hsp /10000}
        
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
        itc_shp <-readOGR(itc_path, itc, stringsAsFactors = FALSE)
        proj4string(itc_shp) <-  CRS(paste("+init=epsg:", tile.args$epsg, sep=""))
        warning(print("for cicle"))
        #for(j in nm){
        if(nm == "P_pct"){
          hsp[hsp ==0] <- 0.000001
          hsp<- hsp[, -c(1:2)]
          hsp <-t(diff(t(log(hsp)),differences=1, lag=3))
          hsp<- cbind(site_bands, hsp)
        }
        hsp <- hsp[complete.cases(hsp), ]
        dim(hsp)
        if(dim(hsp)[1]!=0){
          mod_dir = "./ModelBuild/ALL/"
          load(file = paste(mod_dir, 'plsglm_',nm, sep="" ))
          
          mod.r2=rep(0,length(plsglm))
          mod.aic=rep(0,length(plsglm))
          for(bb in 1: length(plsglm)){
            #mod.aic[bb] <-plsglm[[bb]]$score$aic #plsglm[[bb]]$aic
            mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
              sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
          }
          #if want to use less than 100
          mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), ntop)), mod.r2[mod.r2 %in% tail(sort(mod.r2), ntop)])
          plsglm <- plsglm[mask[,1]]
          mod.r2=rep(0,length(plsglm))
          mod.aic=rep(0,length(plsglm))
          for(bb in 1: ntop){
            mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
              sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
            mod.aic[bb] <-plsglm[[bb]]$score$aic
          }
          delta.aic <- (max(mod.r2)-mod.r2)
          w.daic <- softmax(-0.5*delta.aic)
          
          #output.daic = output.sd.daic = rep(0,dim(hsp)[1])
          #rm(output, output.sd)
          
          preds_all <- lapply(1:length(w.daic), predict_tile, plsglm = plsglm, hsp = hsp, w.daic = w.daic)
          #output.daic <- lapply(function(x, outres){outres[[x]]$output.daic})
          output.daic = output.sd.daic =0
          output = NULL
          for(x in 1:ntop){
            output.daic <- output.daic + preds_all[[x]]$output.daic
            output.sd.daic <- output.sd.daic + preds_all[[x]]$output.sd.daic
            output <- rbind(output, preds_all[[x]]$output)
          }
          rm(preds_all)
          rm(plsglm)
          output.daic <- as.data.frame(as.matrix(output.daic))
          output.sd.daic <- as.data.frame(as.matrix(output.sd.daic))
          
          colnames(output.daic) <-  "yhat"
          colnames(output.sd.daic) <-  "yhat_unc"
          output.daic <- output.daic[!is.nan(output.daic$yhat), ]
          output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]
          out_dir <- tile.args$path
          #save(output, file = paste(paste(tile.args$path, "products/", sep="/"), tile.args$xmin, tile.args$ymax, nm, "hist.RData", sep="_"))
          
          inc_ave_sd <- apply(output,2, function(x) sd(x, na.rm=T))
          
          warning(print("traitsToshp"))
          id <- sub("_", "", nm)
          itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, epsg = tile.args$epsg, itc_shp= itc_shp, trName = c(id, "sd.fit", "md_sd"),
                                 inc_ave = cbind(output.daic, output.sd.daic, inc_ave_sd))
          print("getSpatialRegression regr ok")
          itc_shp <- itc_shp[!is.na(unlist(itc_shp@data[2])),]
          itc_shp <- itc_shp[!is.nan(unlist(itc_shp@data[2])),]
          itc_shp@data <- data.frame(sapply(itc_shp@data, function(x) as.numeric(as.character(x))))
          writeOGR(itc_shp, paste(tile.args$path, "products/", sep=""), paste(c(unlist(strsplit(itc, "_"))[1:5],tile.args$xmin,tile.args$ymin,nm, "shp_product"), collapse="~"), driver="ESRI Shapefile",overwrite_layer=T)
          
          warning("did it!")
        }
        warning(print(paste(itc, "hsp empty!")))
      }
    }else{
      warning(paste("WARNING!:",tile.args$path, tile.args$f, "not found"))
    }
  })
}


# pt = "/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/OSBS/Reflectance/"
# f = "NEON_D03_OSBS_DP3_394000_3280000_reflectance.h5"
# pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018"
# itc_method = "silva"
# scaled = F
getCrownpath <- function(pt =NULL, pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018", 
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
getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = "OSBS", nm = "LMA_g.m2")
getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = "OSBS", nm = "N_pct")
getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/", NeonSite = "OSBS", nm = "P_pct")
