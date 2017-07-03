#mapping LMA as example
#get vector of best random models
library(rgdal)
library(raster)
library(itcSegment)

vect <- read.csv(paste(in.dir, 'Top10/first100.csv', sep=""))
proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")

pl <-c(F,T,F,T,T,F)
md.store1D = list()
md.store2D = list()
in.dir = paste(getwd(), "/inputs/", sep="")
cr.Traits <- read.csv(paste(in.dir, "spectra/CrownTraits.csv",sep=""))
polys.df <- get.plot.extent(plots = read.centroids("OSBS_diversity_plot_centroids"), 80)


allBand=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/MappingTraits/inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/MappingTraits/inputs/Spectra/CrownPix.csv")
all_Bands=as.character(allBand$BandName)
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
i <- as.character(polys.df$id[27])
for (i in as.character(polys.df$id[27])) {
  #fill in later
  rasters <- loadIMG(in.dir, proj, img.lb = i)
  require(MASS)
  hsp <- as.array(rasters$hsp)
  x.size <- dim(hsp)
  dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
  hsp <- as.data.frame(hsp)
  colnames(hsp) <- names(allData)[3:428]
  hsp[,bad_Bands]=NA
  hsp=as.matrix(hsp[, colSums(is.na(hsp)) == 0])
}
token = 0
rasters <- loadIMG(in.dir, proj, img.lb = i)
for(j in names){
  token = token +1
  for(k in vect[,token]){
    aug.spectra <- imp.spectra(paste('Bootstrap_2','/onePix1Crown_',j, k, '.csv', sep = ''), in.dir)
    aug.spectra$X.1 = NULL
    aug.spectra$X.2 = NULL
    
    aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
    X <- aug.spectra[,20:length(aug.spectra[1,])]
    filter <- apply(X, 2, max)
    if(any(filter > 1)){
      X <- X[-which(filter > 1)]
    }
    X=X[, colSums(is.na(X)) == 0]
    Y <- aug.spectra[,10:15]
    X_corr = cor(as.matrix(X), Y)
    aug.X <- data.frame(aug.spectra$name, Y, X)
    # Subset data into cal/val by site
    eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
    train.data <- eval.set$train
    if(!pl[token]){
      pls.mod.train <- pls.cal(train.data, 40, token, norm = F)
      #calculate number of components given min test PRESS or RMSEP
      optim.ncomps <- opt.comps(pls.mod.train, token)
      #plot(pls.mod.train$LMA_g.m2, "loadings", comps = 1:2, legendpos = "topleft", xlab = "nm")
      abline(h = 0)
      
      dat <- data.frame(X=I(hsp)) 
      md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = dat, ncomp=optim.ncomps)
      
    }else{
      grid =10^ seq (10,-2, length =100)
      traits <- train.data[,2:7]
      leaf.spec <- train.data[,8:length(train.data[1,])]
      leaf.trait <- traits[,j]
      lasso.mod =cv.glmnet(as.matrix(leaf.spec),as.vector(leaf.trait),alpha =1, lambda =grid)
      bestlam =lasso.mod$lambda.min
      md.plot=predict(lasso.mod ,s=bestlam ,newx=hsp)
    }
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

trees <- loadITC(in.dir, proj)
xyz <- read.csv("/Users/sergiomarconi/Documents/PhD/Projects/MappingTraits/inputs/LiDAR/xyz.txt")
colnames(xyz) <- c("X", "Y", "Z")
lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg = 32617, resolution = 0.9, 
                   MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                   TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
trees$untg <- spTransform(trees$untg, CRS("+init=epsg:32617"))
trees$tg <- spTransform(trees$tg, CRS("+init=epsg:32617"))
lasITC <- spTransform(lasITC, CRS("+init=epsg:32617"))
eval.extent <-  union(trees$tg, trees$untg)
eval.extent <- eval.extent[lasITC,]
plot(eval.extent)
plot(lasITC)

library(sm)


for(j in names){
  png(paste('./outputs/Maps_', j,'.png',sep=""))
  par(pty="s")
  par(mfrow=c(1,2))
  ave.out <-apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,mean)
  var.out <- apply(eval(parse(text = paste('md.store1D$',j,sep=""))),1,sd)
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
  r <-raster(ave.out, xmn=extent(rasters$cmh)[1], xmx=extent(rasters$cmh)[2],
             ymn=extent(rasters$cmh)[3], ymx=extent(rasters$cmh)[4], crs=CRS("+init=epsg:32617"))
  r <- mask(r, lasITC)
  plot(r, col=terrain.colors(255))
  plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
  
  r <-raster(var.out, xmn=extent(rasters$cmh)[1], xmx=extent(rasters$cmh)[2],
             ymn=extent(rasters$cmh)[3], ymx=extent(rasters$cmh)[4], crs=CRS("+init=epsg:32617"))
  r <- mask(r, lasITC)
  
  plot(r, col=heat.colors(255))
  plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
  
  #plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 2)
  dev.off()
  
  foo <- extract(r, lasITC)
  ave.itc <- lapply(foo, mean)
  sd.itc <- lapply (foo, sd)
  #hist(unlist(ave.itc), xlab = 'gm-2', main = 'Average ITC P',cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3)
  #hist(unlist(sd.itc), xlab = 'gm-2', main = 'St. deviation ITC P',cex.axis = 1.3, cex.lab = 1.3,cex.main = 1.3)
  
}

