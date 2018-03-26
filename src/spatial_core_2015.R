#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)

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

traitsFromITC <- function(tile.args, nm = NULL, old=F, ntop = 20){
  library(raster)
  source(paste(getwd(), '/src/function_spatial.R', sep=""))
  source(paste(getwd(), '/src/functions_build.R', sep=""))

  library(rgdal)
  library("tsensembler")
  library(plsRglm)
  library(maptools)
  library(data.table)
  rasterOptions(tmpdir="./tmp")
  tile.args <- data.frame((tile.args), stringsAsFactors = F)

  tile.args <- transform(tile.args,scaled = as.logical(scaled), epsg = as.integer(epsg))

  #create layer of delineated crowns
  hsp <- brick(paste(tile.args$path,tile.args$f,".tif", sep=""))

  itc <- list.files(tile.args$pt_itc, pattern = paste(tile.args$f, "_silva.shp", sep=""))
  if(!is.na(itc) && !is.nan(itc) && length(itc)!=0){
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
      naval[is.na(naval)] = FALSE
      rm(ndvi, nir860)
      hsp[naval,] = NA
      #hsp[hsp <0 |hsp>1]=NA
      # Vector normalize spectra
      foo <- as.matrix(hsp) * as.matrix(hsp)
      foo <- foo /100000000
      normMat = apply(foo, FUN = sum, MAR = 1, na.rm=T)
      rm(foo)
      normMat <- sqrt(normMat)
      hsp <- hsp /10000
      for(i in 1:dim(hsp)[2]){
        hsp[,i]=hsp[,i]/ normMat
      }
      rm(normMat)

      goodPix.pos <- which(!is.na(hsp[,1]))
      hsp=as.matrix(hsp)
      site_bands <- matrix(0, dim(hsp)[1],2)
      #change in the future
      if(tile.args$NeonSite=="OSBS"){
        site_bands[,1]<-1
      }else if(tile.args$NeonSite=="TALL"){
        site_bands[,2]<-1
      }
      hsp<- cbind(site_bands, hsp)

      #end of normalization of the spectra


      warning(print("for cicle"))
      #for(j in nm){
      hsp <- hsp[complete.cases(hsp), ]
      dim(hsp)
      if(dim(hsp)[1]!=0){
        for(j in c("LMA_g.m2", "N_pct", "P_pct")){
          if(j == "P_pct"){
            hsp[hsp ==0] <- 0.000001
            hsp<- hsp[, -c(1:2)]
            hsp <- log(hsp)
            hsp <- t(hsp)
            hsp <- diff(hsp,differences=1, lag=3)
            hsp <-t(hsp)
            hsp<- cbind(site_bands, hsp)
          }

          itc_shp <-readOGR(tile.args$pt_itc, paste(tile.args$f, "silva", sep="_"), stringsAsFactors = FALSE)
          proj4string(itc_shp) <-  CRS(paste("+init=epsg:", tile.args$epsg, sep=""))
          mod_dir = "./ModelBuild/ALL/"
          load(file = paste(mod_dir, 'plsglm_',j, sep="" ))

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
          output.daic = output.uppi = output.lopi = 0
          output = NULL
          for(x in 1:ntop){
            output.daic <- output.daic + preds_all[[x]]$output.daic
            output.uppi <- output.uppi + preds_all[[x]]$pi.upper
            output.lopi <- output.lopi + preds_all[[x]]$pi.lower
            #output.sd.daic <- output.sd.daic + preds_all[[x]]$output.sd.daic
            output <- rbind(output, preds_all[[x]]$output)
          }
          rm(preds_all)
          rm(plsglm)
          output.daic <- as.data.frame(as.matrix(output.daic))
          output.uppi <- as.data.frame(as.matrix(output.uppi))
          output.lopi <- as.data.frame(as.matrix(output.lopi))

          #output.sd.daic <- as.data.frame(as.matrix(output.sd.daic))

          colnames(output.daic) <-  "yhat"
          colnames(output.uppi) <-  "yPIup"
          colnames(output.lopi) <-  "yPIlo"

          #colnames(output.sd.daic) <-  "yhat_unc"
          output.daic <- output.daic[!is.nan(output.daic$yhat), ]
          output.uppi <- output.uppi[!is.nan(output.uppi$yPIup), ]
          output.lopi <- output.lopi[!is.nan(output.lopi$yPIlo), ]
          #output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]
          out_dir <- tile.args$path
          #save(output, file = paste(paste(tile.args$path, "products/", sep="/"), tile.args$xmin, tile.args$ymax, nm, "hist.RData", sep="_"))

          inc_ave_sd <- apply(output,2, function(x) sd(x, na.rm=T))

          warning(print("traitsToshp"))
          id <- sub("_", "", j)
          itc_shp <- traitsToshp(x.size, goodPix.pos, extnt, epsg = tile.args$epsg, itc_shp= itc_shp,
                                 trName = c(id, paste(id, "up", sep = "."), paste(id, "lo", sep = "."), paste(id, "sd", sep = ".")),
                                 inc_ave = cbind(output.daic, output.uppi, output.lopi, inc_ave_sd))
          print("getSpatialRegression regr ok")
          itc_shp <- itc_shp[!is.na(unlist(itc_shp@data[2])),]
          itc_shp <- itc_shp[!is.nan(unlist(itc_shp@data[2])),]
          itc_shp@data <- data.frame(sapply(itc_shp@data, function(x) as.numeric(as.character(x))))
          out_path = paste(gsub('.{7}$', '', tile.args$pt), "4/Predictions/", sep="")
          writeOGR(itc_shp, out_path, paste(tile.args$f,j, tile.args$NeonSite, sep="~"), driver="ESRI Shapefile",overwrite_layer=T)

          warning("did it!")
        }
      }
      warning(print(paste(itc, "hsp empty!")))
    }
  }else{
    warning(paste("WARNING!:",tile.args$path, tile.args$f, "not found"))
  }
}


# pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/"
# f = "IPLLS7044E.tif"
# pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018"
# itc_method = "silva"
# scaled = F
# NeonSite = "TALL"
getCrownpath <- function(pt =NULL, pathSpatial="/ufrc/ewhite/s.marconi/Marconi2018",
                         f = NULL, nm = NULL, NeonSite = "NULL", cores = 1, itc_method = "silva", scaled = F){
  setwd(pathSpatial)
  .libPaths(c("/usr/lib/R/site-library","/usr/local/lib/R/site-library","/usr/lib/R/library",
              "/home/s.marconi/R/x86_64-pc-linux-gnu-library/3.4",
              "/apps/gcc/6.3.0/R/3.4.0/lib64/R/library"))

  print(NeonSite)
  library(rgdal)
  if(NeonSite=="OSBS"){
    pt_ras = paste(pt, "D03/OSBS/L4/Rasters/", sep = "")
    pt_itc = paste(pt, "D03/OSBS/L4/ITCs/", sep = "")
    epsg <- "32617"
  }else if(NeonSite=="TALL"){
    pt_ras = paste(pt, "D08/TALL/L4/Rasters/", sep = "")
    pt_itc = paste(pt, "D08/TALL/L4/ITCs/", sep = "")
    epsg <- "32616"
  }

  f <- gsub('.{4}$', '', f)

  tile.args <- data.frame(epsg, NeonSite, pt_ras, pt_itc, f, scaled, itc_method, stringsAsFactors = F)
  colnames(tile.args) <- c("epsg", "NeonSite", "path", "pt_itc", "f", "scaled", "itc_method")
  traitsFromITC(tile.args, nm = nm, old = !newmeta)
}
getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/", NeonSite = "TALL", nm = "LMA_g.m2")
#getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/", NeonSite = "OSBS", nm = "N_pct")
#getCrownpath(f = args[7], pt ="/ufrc/ewhite/s.marconi/NeonData/2015_Campaign/", NeonSite = "OSBS", nm = "P_pct")
