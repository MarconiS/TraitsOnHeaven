#create an object carrying 
#mapping LMA as example
#get vector of best random models
library(rgdal)
library(raster)
library(itcSegment)
require(MASS)
library(sm)
library(readr)
library(pls)

forToday <- c("OSBS_001", "OSBS_002", "OSBS_005", "OSBS_008", "OSBS_016", "OSBS_017", "OSBS_021", "OSBS_026", "OSBS_033", "OSBS_041")

source(paste(getwd(), '/src/itcSrc.R', sep=""))
loops <- 300
NeonSite <- "OSBS"
setwd('~/Documents/Projects/TraitsOnHeaven/')
out.dir = paste(getwd(), NeonSite,"outputs/", sep="/")
in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")
names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
#debugging
laps =rounds=1
path = paste(getwd(), NeonSite, sep="/")
setwd(path)

vect <- read_csv(paste(in.dir, 'FinalSet/listPix.csv', sep=""))
proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
md.store1D = list()
md.store2D = list()
cr.Traits <- read.csv(paste(in.dir, "/Spectra/trainCrownTraits.csv",sep=""))
polys.df <- get.plot.extent(plots = read.centroids("OSBS_diversity_plot_centroids"), 80)

allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("./inputs/Spectra/CrownPix.csv")
all_Bands=as.character(allBand$BandName)
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
#for (i in as.character(polys.df$id)) {
plot.ave <- list()
plot.sd <- list()
listPlot <- (read.csv('inputs/FinalSet/plotsList.csv', header = F, stringsAsFactors = F))
for (i in listPlot[['V1']]){#forToday[1]) {
  rasters <- loadIMG(in.dir, proj, img.lb = i)
  hsp <- as.array(rasters$hsp)
  x.size <- dim(hsp)
  dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
  hsp <- as.data.frame(hsp)
  colnames(hsp) <- names(allData)[3:428]
  hsp <- normalizeImg(hsp)
  hsp <- hsp[,colSums(is.na(hsp))<nrow(hsp)]
  hsp <- hsp[,-c(seq(1,8))]
  #hsp[,bad_Bands]=NA
  #hsp=as.matrix(hsp[, colSums(is.na(hsp)) == 0])
  hsp=as.matrix(hsp)
  #hsp[,bad_Bands]=NA
  
  token = 0
  jj = 0
  for(j in names){
    jj = jj +1
    token = token +1
    pixFile <- vect[,which(colnames(vect) %in% j)]
    pixFile <- pixFile[!is.na(pixFile),]
    pixFile <- eval(parse(text = paste('pixFile$',j,sep="")))
    for(k in pixFile){
      aug.spectra <- imp.spectra(k, paste(in.dir, 'FinalSet/chosenSets/', sep=""))
      aug.spectra$X.1 = NULL
      aug.spectra$X.2 = NULL
      
      aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
      X<- aug.spectra[grepl("band", names(aug.spectra))]
      Y <- aug.spectra[colnames(aug.spectra) %in% j]
      aug.X <- data.frame(aug.spectra$name, Y, X)
      # Subset data into cal/val by site
      eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
      train.data <- eval.set$train
      train.data$aug.spectra.name = NULL
      train.data <- train.data[ , apply(train.data, 2, function(x) !any(is.na(x)))]
      
      #standardize
      pls.mod.train <- pls.calImg(train.data, 25,nm = names, jj, norm = F)
      
      #calculate number of components given min test PRESS or RMSEP
      optim.ncomps <- opt.compsImg(pls.mod.train, Y)
      dat <- data.frame(X=I(hsp)) 
      md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = dat, ncomp=optim.ncomps, type='response')
      
      if(!exists("md.all")){
        md.all <- as.vector(md.plot)
      }else{
        md.all <- cbind(md.all, md.plot)
      }
      
      dim(md.plot) <- c(x.size[1],x.size[2])
      image(md.plot)
      md.store2D[[names[j]]] <- md.plot
    }
    md.store1D[[j]] <- md.all
    rm(md.all)
  }
  
  #trees <- loadITC(in.dir, proj)
  xyz <- read.csv(paste(in.dir, "Geofiles/LiDAR/plot_", i, ".txt", sep=""))
  colnames(xyz) <- c("X", "Y", "Z")
  lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg = 32617, resolution = 0.9, 
                     MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                     TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
  #trees$untg <- spTransform(trees$untg, CRS("+init=epsg:32617"))
  #trees$tg <- spTransform(trees$tg, CRS("+init=epsg:32617"))
  lasITC <- spTransform(lasITC, CRS("+init=epsg:32617"))
  plot(lasITC)
  
  out.ave <- list()
  out.sd <-list()
  for(j in names){
    png(paste('./outputs/Maps_', i,'_', j, '.png',sep=""))
    par(pty="s")
    par(mfrow=c(1,2))
    ave.out <-apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,na.rm=TRUE, mean)
    var.out <- apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,na.rm=TRUE, sd)
    # hist(ave.out, main = paste("Average", j, "pixel"), xlab = "gm^{-2}", cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7)
    # hist(var.out, main = paste("Average", j, "pixel"), xlab = "gm^{-2}", cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7)
    # plot(density(ave.out))
    # 
    # plot(density(ave.out), type="l", pch=0.01, main = j)
    # up.b <- ave.out+var.out
    # low.b <- ave.out-var.out
    # plot(density(up.b), type="l", pch=0.01,  col='red', add=T)
    #sm.density.compare(c(low.b, ave.out, up.b), c(rep(1, length(low.b)),rep(2, length(ave.out)),rep(3, length(up.b))))
    
    
    dim(ave.out) <- c(x.size[1],x.size[2])
    dim(var.out) <- c(x.size[1],x.size[2])
    r <-raster(ave.out, xmn=extent(rasters$hsp)[1], xmx=extent(rasters$hsp)[2],
               ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS("+init=epsg:32617"))
    r <- mask(r, lasITC)
    plot(r, col=terrain.colors(255))
    plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
    
    r <-raster(var.out, xmn=extent(rasters$hsp)[1], xmx=extent(rasters$hsp)[2],
               ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS("+init=epsg:32617"))
    r <- mask(r, lasITC)
    
    plot(r, col=heat.colors(255))
    plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
    
    #plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 2)
    dev.off()
    
    foo <- extract(r, lasITC)
    ave.itc <- lapply(foo, na.rm=TRUE,mean)
    sd.itc <- lapply (foo,na.rm=TRUE, sd)
    out.ave[[j]] <- unlist(ave.itc)
    out.sd[[j]] <- unlist(sd.itc)
    
    #hist(unlist(ave.itc), xlab = 'gm-2', main = 'Average ITC P',cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3)
    #hist(unlist(sd.itc), xlab = 'gm-2', main = 'St. deviation ITC P',cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3)
  }
  plot.ave[[i]] <- out.ave
  plot.sd[[i]] <- out.sd
}
