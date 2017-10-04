#-----getSpatialRegression------------------------------------------------------------------------------------------


getSpatialRegression <- function(NeonSite = "OSBS", names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"), in.dir = NULL, out.dir = NULL,
                                 tile = 80, epsg = NULL, proje = NULL, normalz = T){
  md.store1D = list()
  md.store2D = list()
  allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
  allData=read.csv("./inputs/Spectra/CrownPix.csv")
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  plot.ave <- list()
  plot.sd <- list()
  listPlot <- (read.csv('inputs/plotsList.csv', header = F, stringsAsFactors = F))
  for (i in listPlot[['V1']]){#forToday[1]) {
    if(exists("lasITC")){rm(lasITC)}
    print(i)
    rasters <- loadIMG(in.dir, proje, img.lb = i)
    hsp <- as.array(rasters$hsp)
    x.size <- dim(hsp)
    dim(hsp) <- c(x.size[1]*x.size[2], x.size[3])
    hsp <- as.data.frame(hsp)
    colnames(hsp) <- names(allData)[3:428]
    hsp <- normalizeImg(hsp)
    goodPix.pos <- which(!is.na(hsp[,1]))
    hsp=as.matrix(hsp)
    if(normalz){hsp <-t(diff(t(log(hsp)),differences=1, lag=3))}
    hsp <- hsp[,colSums(is.na(hsp))<nrow(hsp)]

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
        if(j == "name"){
          mod.r2[bb] <- mean(mod.stats[[bb]]$name$R2)
        }else{
          mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
        }
      }
      mask <- which(mod.r2 %in% tail(sort(mod.r2), 100))
      if(j != "name"){weights[[j]] <- mod.r2[mask] /sum(mod.r2[mask])}
      
      for(k in mask){
        #read the ith model
        pls.mod.train <- mod.out[[k]]
        optim.ncomps <- mod.comps[k]
        
        if(j != "name"){
          dat <- data.frame(X=I(hsp)) 
          md.plot <- predict(eval(parse(text = paste('pls.mod.train$',j,sep=""))), newdata = dat, ncomp=optim.ncomps, type='response')
        }else{
          dat <- data.frame(X=I(hsp)) 
          prova <- as.vector(predict(eval(parse(text = paste('pls.mod.train$', j,sep=""))), 
                                     newdata = na.omit(hsp), ncomp=optim.ncomps, type="response"))
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
    if(!exists("lasITC")){
      xyz <- read.csv(paste(in.dir, "Geofiles/RSData/pointCloud/ptcloud_", i, ".csv", sep=""))
      lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg, resolution = 0.5, 
                         MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                         TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
      lasITC <- spTransform(lasITC, CRS(proje))
      writeOGR(lasITC, paste(out.dir, "crownDelim", sep=""), paste(i,"_itc", by=""), overwrite_layer = T, 
               check_exists = T, driver="ESRI Shapefile")
    }
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
                     ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS(proje))
      r.ave <- mask(r.ave, lasITC)
      writeRaster(r.ave, paste("./outputs/average", j,i, sep="_"), overwrite=TRUE, format='GTiff')
      if(j != "name"){
        plot(r.ave, col=terrain.colors(255))
      }else{
        plot(r.ave, col=bpy.colors(16))
      }
      plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
      
      r.sd <-raster(var.out, xmn=extent(rasters$hsp)[1], xmx=extent(rasters$hsp)[2],
                    ymn=extent(rasters$hsp)[3], ymx=extent(rasters$hsp)[4], crs=CRS(proje))
      r.sd <- mask(r.sd, lasITC)
      writeRaster(r.sd, paste("./outputs/rasters/variance", j,i, sep="_"),overwrite=TRUE, format='GTiff')
      
      plot(r.sd, col=heat.colors(255))
      plot(lasITC,axes=T, border="blue", add=TRUE, lwd = 1.5)
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
    if(j != "name"){
      setwd(paste(out.dir, "traits_regression_csv/", sep=""))
      save(out.ave, file = paste(out.dir, "traits_regression_csv/", NeonSite,"_",i, "_aveOutputs", sep = ""))
      save(out.sd, file = paste(out.dir, "species_classification_csv/", NeonSite,"_",i, "_sdOutputs", sep = ""))
    } else{
      setwd(paste(out.dir, "species_classification_csv/", sep=""))
      save(out.ave, file = paste(NeonSite,i, "species_aveOutputs", sep = "_"))
      save(out.sd, file = paste(NeonSite,i, "species_sdOutputs", sep = "_"))
    }
    rm(lasITC)
    setwd("../..")
  }
  
}
#-----getTreeRegression------------------------------------------------------------------------------------------

getTreeRegression<- function(NeonSite = "OSBS", 
                             names =  c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
                             in.dir = NULL, out.dir = NULL){
  rm(allPlotsOut)
  listPlot <- (read.csv(paste(getwd(), '/inputs/Plot_class.csv', sep=""), header = T, stringsAsFactors = F))
  for (i in listPlot[['Plot_ID']]){#forToday[1]) {
    f <- paste(out.dir, "/traits_regression_csv/", NeonSite, "_", i, "_", "aveOutputs" ,sep="" )
    f.sp <-  paste(out.dir, "/species_classification_csv/", NeonSite, "_", i, "_", "species_aveOutputs" ,sep="" )
    load(f)
    foo <- unlist(out.ave)
    crowns <- length(out.ave[[1]])
    categories <- rep(names[1], crowns)
    
    for(tkn in 2: length(names)){
      categories <- c(categories, rep(names[tkn], crowns))
    }
    
    plot.id <- rep(i, length(names) * crowns)
    foo <- data.frame(foo, factor(categories), as.character(plot.id), stringsAsFactors = F)
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
  write_csv(allPlotsOut, paste(out.dir, "tree_out/individualTreesOut.csv", sep=""))
  
  allPlotsOut <- cbind(allPlotsOut, rep("predicted", dim(allPlotsOut)[1]))
  colnames(allPlotsOut)[7] <- "origin"
  cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
  cr.Traits <- cr.Traits[,colnames(cr.Traits) %in% c("LMA_g.m2", "d13C", "d15N","C_pct", "N_pct", "P_pct","name")]
  cr.Traits <- melt(cr.Traits)
  cr.Traits <- inner_join(cr.Traits, species.id, by = "name")
  cr.Traits <- cbind(cr.Traits, rep("observed", dim(cr.Traits)[1]))
  colnames(cr.Traits) <- c("name", "trait", "value", "species","origin")
  final.df <- rbind(cr.Traits, allPlotsOut[,c("name", "trait", "value", "species","origin")])
  
  for(trName in names){
    #p <- ggplot(subset(allPlotsOut, trait == trName), aes(factor(Plot_ID), value))
    p <- ggplot(subset(final.df, trait == trName), aes(factor(name), value))
    p + theme_bw(base_size = 16) + geom_violin(aes(colour = factor(origin), alpha = 0.5), draw_quantiles = c(0.25, 0.5, 0.75)) + #+ 
      labs(title = paste(trName, "distribution in", NeonSite,"\n"), alpha = "",x = "species", y = trName) + 
      theme(axis.text.x=element_text(angle=90,hjust=1)) 
    
    
    # p + theme_bw(base_size = 16) + geom_violin(aes(fill = factor(class)), draw_quantiles = c(0.25, 0.5, 0.75)) + 
    #   labs(title = paste(trName, "distribution in OSBS\n"), alpha = "",x = "Plot ID", y = "value", color = "Species\n", fill = "Ecosystem type\n") +
    #   geom_jitter(aes(colour = factor(name), alpha = 0.5), height = 0, width = 0.1) + #scale_colour_manual(values = c("yellow", "black", "red", "gray", "orange", "brown")) +
    #   theme(axis.text.x=element_text(angle=90,hjust=1)) 
    ggsave(paste(out.dir, "tree_out/", trName, "_dist.png", sep=""), width=30, height = 20, units="in")
  }
}


getPointCloud <- function(NeonString = "2014_OSBS_", in.dir = NULL, out.dir = NULL){
  path.json <- "inputs/Geofiles/RSdata/pointCloud/extract_ptCloud.json"
  myjsonList = fromJSON(file = path.json)
  listPlot <- (read.csv(paste(getwd(), '/inputs/Plot_class.csv', sep=""), header = T, stringsAsFactors = F))
  for (i in listPlot[['Plot_ID']]){#forToday[1]) {
    #load the tif file
    r <- raster(paste(in.dir, 'Geofiles/RSdata/hs/', i, '_hyper.tif', sep=""))
    f1 <- paste(NeonString, "1_", floor(extent(r)[1]/1000)*1000, "_", floor(extent(r)[3]/1000)*1000, "_", "colorized.laz", sep ="")
    f2 <- paste(NeonString, "1_", floor(extent(r)[1]/1000)*1000, "_", floor(extent(r)[3]/1000)*1000, ".laz", sep="")
    if(file.exists(paste(in.dir,"Geofiles/RSdata/Classified_point_cloud/", f1, sep = ""))){
      myjsonList$pipeline[[1]] <- f1
    }else if(file.exists(paste(in.dir,"Geofiles/RSdata/Classified_point_cloud/", f2, sep = ""))){
      myjsonList$pipeline[[1]] <- f2
    }else{
      warning("no laz file found for this plot")
    }
    poly.fill <- paste(extent(r)[1], extent(r)[3], ",", 
                       extent(r)[2], extent(r)[3], ",",
                       extent(r)[2], extent(r)[4], ",",
                       extent(r)[1], extent(r)[4], ",",
                       extent(r)[1], extent(r)[3],sep = " ")
    myjsonList$pipeline[[2]]$polygon <- paste("POLYGON((",poly.fill, "))", sep = "" )
    myjsonList$pipeline[[3]]$filename <- paste(NeonString, i, '.las', sep="")
    newjson= toJSON(myjsonList)
    write(newjson,paste(getwd(), '/inputs/Geofiles/RSdata/Classified_point_cloud/extract_ptCloud.json', sep=""))
    setwd(paste(in.dir, "/Geofiles/RSdata/Classified_point_cloud",sep="" ))
    system("pdal pipeline extract_ptCloud.json")
    system(paste("pdal translate ",myjsonList$pipeline[[3]]$filename, " ", 
                 paste('../pointCloud/ptcloud_', i, '.csv', sep=""), sep = ""))
    setwd("../../../..")
  }
}
