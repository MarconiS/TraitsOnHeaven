#create an object carrying 
#mapping LMA as example
#get vector of best random models
rm(list=ls(all=TRUE))   # clear workspace

library(rgdal)
library(raster)
library(itcSegment)
require(MASS)
library(sm)
library(readr)
library(pls)
setwd('/Users/sergiomarconi/Documents/Projects/TraitsOnHeaven/')
source(paste(getwd(), '/src/itcSrc.R', sep=""))
source(paste(getwd(), '/src/functionsImage.R', sep=""))
source(paste(getwd(), '/src/src.R', sep=""))
source(paste(getwd(), '/src/modelSrc.R', sep=""))

scalar1 <- function(x) {x / sqrt(sum(x^2))}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


NeonSite <- "OSBS"
setwd('~/Documents/Projects/TraitsOnHeaven/')
out.dir = paste(getwd(), NeonSite,"outputs/", sep="/")
in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")
names <- c("name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
#debugging
laps =rounds=1
path = paste(getwd(), NeonSite, sep="/")
setwd(path)

proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
names <- c("name","LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
md.store1D = list()
md.store2D = list()
polys.df <- get.plot.extent(plots = read.centroids("OSBS_diversity_plot_centroids"), 80)

allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("./inputs/Spectra/CrownPix.csv")
all_Bands=as.character(allBand$BandName)
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
#for (i in as.character(polys.df$id)) {
plot.ave <- list()
plot.sd <- list()
listPlot <- (read.csv('inputs/FinalSet/plotsList.csv', header = F, stringsAsFactors = F))
names="name"
j = "name"
load(file = paste(out.dir, 'models_comps_',j, sep="" ))
load(file = paste(out.dir, 'models_out_',j, sep="" ))
load(file = paste(out.dir, 'models_stats_',j, sep="" ))

for (i in listPlot[['V1']]){#forToday[1]) {
  print(i)
  rasters <- loadIMG(in.dir, proj, img.lb = i)
  hsp <- as.array(rasters$hsp)
  x.size <- dim(hsp)
  dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
  hsp <- as.data.frame(hsp)
  colnames(hsp) <- names(allData)[3:428]
  hsp <- normalizeImg(hsp)
  hsp <- hsp[,colSums(is.na(hsp))<nrow(hsp)]
  goodPix.pos <- which(!is.na(hsp[,1]))
  # remove the first bands: 
  # hsp <- hsp[,-c(seq(1,8))]
  
  hsp=as.matrix(hsp)
  #hsp[,bad_Bands]=NA
  
  token = 0
  jj = 0
  weights <- list()
  for(j in names){
    jj = jj +1
    token = token +1
    
    # load(file = paste(out.dir, 'models_comps_',j, sep="" ))
    # load(file = paste(out.dir, 'models_out_',j, sep="" ))
    # load(file = paste(out.dir, 'models_stats_',j, sep="" ))
    mod.r2 <- rep(NA, length(mod.stats))
    for(bb in 1: length(mod.r2)){
      mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
    }
    mask <- which(mod.r2 %in% tail(sort(mod.r2), 20) & mod.r2 >0.0)
    
    if(j != "name"){weights[[j]] <- mod.r2[mask] /sum(mod.r2[mask])
    }else{
      #norm.R2 <- scalar1(mod.r2[mask])
      norm.R2 <- mod.r2[mask]/sum(mod.r2[mask])
      multiplier <- 10^3 #decimalplaces(min(norm.R2))
      weights[[j]] <- round(norm.R2 * multiplier)
    }
    
    tkn <- 0 
    for(k in mask){
      tkn <- tkn + 1
      #read the ith model
      pls.mod.train <- mod.out[[k]]
      optim.ncomps <- 6 #mod.comps[k]
      #standardize
      #pls.mod.train <- pls.calImg(train.data, 25,nm = names, jj, norm = F)
      
      #calculate number of components given min test PRESS or RMSEP
      #optim.ncomps <- opt.compsImg(pls.mod.train, Y)
      if(names != "name"){
        dat <- data.frame(X=I(hsp)) 
        md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = dat, ncomp=optim.ncomps, type='response')
        if(!exists("md.all")){
          md.all <- as.vector(md.plot)
        }else{
          md.all <- cbind(md.all, md.plot)
        }
      }else{
        dat <- data.frame(X=I(hsp)) 
        #md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = hsp, ncomp=optim.ncomps, type='class')
        prova <- as.vector(predict(eval(parse(text = paste('pls.mod.train$', j,sep=""))), 
                                   newdata = na.omit(hsp), ncomp=optim.ncomps, type="class"))
        for(donCare in 1: weights[[j]][tkn]){
          md.plot <- as.numeric(rep(NA, x.size[1]*x.size[2]))
          md.plot[goodPix.pos] <- as.numeric(prova)
          
          if(!exists("md.all")){
            md.all <- as.vector(md.plot)
          }else{
            md.all <- cbind(md.all, md.plot)
          }
        }
        #here get the one with highest
      }
      dim(md.plot) <- c(x.size[1],x.size[2])
      #image(md.plot)
      md.store2D[[names[j]]] <- md.plot
    }
    md.store1D[[j]] <- md.all
    rm(md.all)
  }
  
  #trees <- loadITC(in.dir, proj)
  xyz <- read.csv(paste(in.dir, "Geofiles/RSdata/pointCloud/ptcloud_", i, ".csv", sep=""))
  #colnames(xyz) <- c("X", "Y", "Z")
  lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg = 32617, resolution = 0.9, 
                     MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                     TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
  #trees$untg <- spTransform(trees$untg, CRS("+init=epsg:32617"))
  #trees$tg <- spTransform(trees$tg, CRS("+init=epsg:32617"))
  lasITC <- spTransform(lasITC, CRS("+init=epsg:32617"))
  #plot(lasITC)
  
  out.ave <- list()
  out.sd <-list()
  for(j in names){
    png(paste('./outputs/Maps_', i,'_', j, '.png',sep=""))
    par(pty="s")
    par(mfrow=c(1,2))
    #ave.out <-apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1, function(x) weighted.mean(x, weights[[j]], na.rm=T))
    #var.out <- apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,na.rm=TRUE, sd)
    if(j != "name"){
      ave.out <-apply(md.store1D[[j]],1, function(x) weighted.mean(x, weights[[j]], na.rm=T))
      var.out <- apply(md.store1D[[j]],1,na.rm=TRUE, sd)
      #ave.out <- apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,na.rm=TRUE, mean)
    }else{
      ave.out=rep(NA, dim(md.store1D[[j]])[1])
      var.out = rep(NA, dim(md.store1D[[j]])[1])
      for(ii in 1: dim(md.store1D[[j]])[1]) {
        temp.freq <- table(md.store1D[[j]][ii,], useNA = "ifany")
        ave.out[ii] <- as.numeric(names(temp.freq)[which(temp.freq == max(temp.freq))])
        if(!is.na(ave.out[ii])){var.out[ii] <- max(temp.freq)/sum(temp.freq)}
      }
      
    }
    dim(ave.out) <- c(x.size[1],x.size[2])
    dim(var.out) <- c(x.size[1],x.size[2])
    r.ave <-raster(ave.out, xmn=extent(rasters$hsp)[1], xmx=extent(rasters$hsp)[2],
                   ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS("+init=epsg:32617"))
    r.ave <- mask(r.ave, lasITC)
    writeRaster(r.ave, paste("./outputs/average", j,i, sep="_"), overwrite=TRUE, format='GTiff')
    if(j != "name"){
      plot(r.ave, col=terrain.colors(255))
    }else{
      plot(r.ave, col=bpy.colors(16))
    }
    plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
    
    r.sd <-raster(var.out, xmn=extent(rasters$hsp)[1], xmx=extent(rasters$hsp)[2],
                  ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS("+init=epsg:32617"))
    r.sd <- mask(r.sd, lasITC)
    writeRaster(r.sd, paste("./outputs/variance", j,i, sep="_"),overwrite=TRUE, format='GTiff')
    
    plot(r.sd, col=heat.colors(255))
    plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
    
    #plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 2)
    dev.off()
    
    foo <- extract(r.ave, lasITC)
    if(j != "name"){
      ave.itc <- lapply(foo, na.rm=TRUE,mean)
      sd.itc <- lapply (foo,na.rm=TRUE, sd)
    }else{
      ave.itc <- rep(NA, length(foo))
      sd.itc <- rep(NA, length(foo))
      for(ii in 1: length(foo)) {
        temp.freq <- table(foo[ii], useNA = "no")
        if(is.null(names(temp.freq))){
          ave.itc[ii] <- NA
        }else{
          ave.itc[ii] <- as.numeric(names(temp.freq)[which(temp.freq == max(temp.freq))])
        }
        if(!is.na(ave.out[ii])){sd.itc[ii] <- max(temp.freq)/sum(temp.freq)}
      }
    }
    out.ave[[j]] <- unlist(ave.itc)
    out.sd[[j]] <- unlist(sd.itc)
  }
  plot.ave[[i]] <- out.ave
  plot.sd[[i]] <- out.sd
  
  save(out.ave, file = paste(NeonSite,i, "species_aveOutputs", sep = "_"))
  save(out.sd, file = paste(NeonSite,i, "species_sdOutputs", sep = "_"))
}
