#-----getSpatialRegression------------------------------------------------------------------------------------------

getSpatialRegression <- function(NeonSite = "OSBS",
                                 names = c("name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                                 tile = 80,
                                 out.dir = paste(getwd(), NeonSite,"outputs/", sep="/"),
                                 in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/"),
                                 proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"){
  md.store1D = list()
  md.store2D = list()
  polys.df <- get.plot.extent(plots = read.centroids(paste(NeonSite,"diversity_plot_centroids", sep="_"), tile))
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  plot.ave <- list()
  plot.sd <- list()
  listPlot <- (read.csv('inputs/plotsList.csv', header = F, stringsAsFactors = F))
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
    hsp=as.matrix(hsp)

    token = 0
    jj = 0
    weights <- list()
    for(j in names){
      jj = jj +1
      token = token +1
      
      load(file = paste(out.dir, 'models_comps_',j, sep="" ))
      load(file = paste(out.dir, 'models_out_',j, sep="" ))
      load(file = paste(out.dir, 'models_stats_',j, sep="" ))
      mod.r2 <- rep(NA, length(mod.stats))
      for(bb in 1: length(mod.r2)){
        mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
      }
      mask <- which(mod.r2 %in% tail(sort(mod.r2), 100) & mod.r2 >0.0)
      if(j != "name"){weights[[j]] <- mod.r2[mask] /sum(mod.r2[mask])}
      
      for(k in mask){
        #read the ith model
        pls.mod.train <- mod.out[[k]]
        optim.ncomps <- mod.comps[k]

        if(names != "name"){
          dat <- data.frame(X=I(hsp)) 
          md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = dat, ncomp=optim.ncomps, type='response')
        }else{
          dat <- data.frame(X=I(hsp)) 
          prova <- as.vector(predict(eval(parse(text = paste('pls.mod.train$', j,sep=""))), 
                                     newdata = na.omit(hsp), ncomp=optim.ncomps, type="class"))
          md.plot <- as.numeric(rep(NA, x.size[1]*x.size[2]))
          md.plot[goodPix.pos] <- as.numeric(prova)
        }
        
        if(!exists("md.all")){
          md.all <- as.vector(md.plot)
        }else{
          md.all <- cbind(md.all, md.plot)
        }
        
        dim(md.plot) <- c(x.size[1],x.size[2])
        md.store2D[[names[j]]] <- md.plot
      }
      md.store1D[[j]] <- md.all
      rm(md.all)
    }
    
    xyz <- read.csv(paste(in.dir, "Geofiles/LiDAR/ptcloud_", i, ".csv", sep=""))
    lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg, resolution = 0.9, 
                       MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                       TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
    lasITC <- spTransform(lasITC, CRS(proj))
    
    out.ave <- list()
    out.sd <-list()
    
    for(j in names){
      png(paste('./outputs/Maps_', i,'_', j, '.png',sep=""))
      par(pty="s")
      par(mfrow=c(1,2))
      if(j != "name"){
        ave.out <-apply(md.store1D[[j]],1, function(x) weighted.mean(x, weights[[j]], na.rm=T))
        var.out <- apply(md.store1D[[j]],1,na.rm=TRUE, sd)
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
                     ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS(proj))
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
      
      foo <- raster::extract(r.ave, lasITC)
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
}
#-----getTreeRegression------------------------------------------------------------------------------------------

getTreeRegression<- function(NeonSite = "OSBS", 
                             names =  c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                             out.dir = paste(getwd(), NeonSite,"outputs/", sep="/"),
                             in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")){
  
  listPlot <- (read.csv('OSBS/inputs/Plot_class.csv', header = T, stringsAsFactors = F))
  for (i in listPlot[['Plot_ID']]){#forToday[1]) {
    f <- paste("~/Documents/Projects/TraitsOnHeaven/", "/traits_regression_csv/", NeonSite, "_", i, "_", "aveOutputs" ,sep="" )
    f.sp <-  paste("~/Documents/Projects/TraitsOnHeaven/", "/species_classification_csv/", NeonSite, "_", i, "_", "species_aveOutputs" ,sep="" )
    load(f)
    foo <- unlist(out.ave)
    crowns <- length(out.ave[[1]])
    categories <- rep(names[1], crowns)
    
    for(tkn in 2: length(names)){
      categories <- c(categories, rep(names[tkn], crowns))
    }
    
    plot.id <- rep(i, length(names) * crowns)
    foo <- data.frame(foo, factor(categories), as.character(plot.id))
    colnames(foo) <- c("value", "trait", "Plot_ID")
    foo <- inner_join(foo, listPlot, by = "Plot_ID")
    ####
    load(f.sp)
    foo.sp <- rep(unlist(out.ave), length(categories) / crowns)
    species.id <- read_csv(paste(in.dir, "Spectra/speciesName.csv", sep=""))
    foo <- data.frame(foo, (foo.sp))
    colnames(foo) <- c("value", "trait", "Plot_ID", "class","species")
    foo <- inner_join(foo, species.id, by = "species")
    
    ####
    if(!exists("allPlotsOut")){
      allPlotsOut <- foo
    }else{
      allPlotsOut <- rbind(allPlotsOut, foo)
    }
  }
  write_csv(allPlotsOut, paste(out.dir, "individualTreesOut.csv", sep=""))
  
  for(trName in names){
    p <- ggplot(subset(allPlotsOut, trait == trName), aes(factor(Plot_ID), value))
    p + theme_bw(base_size = 16) + geom_violin(aes(fill = factor(class)), draw_quantiles = c(0.25, 0.5, 0.75)) + 
      labs(title = paste(trName, "distribution in OSBS\n"), alpha = "",x = "Plot ID", y = "value", color = "Species\n", fill = "Ecosystem type\n") +
      geom_jitter(aes(colour = factor(name), alpha = 0.5), height = 0, width = 0.1) + scale_colour_manual(values = c("yellow", "black", "red", "gray", "orange", "brown")) +
      theme(axis.text.x=element_text(angle=90,hjust=1)) 
    ggsave(paste(out.dir, trName, "_dist.png", sep=""), width=30, height = 20, units="in")
  }
}
