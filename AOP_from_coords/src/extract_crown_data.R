
extract_crown_data <- function(centroids, hps_f, f_path, chm_f, epsg, wd,NeonSites =NULL,  buffer = 20, cores = 2){
  library(foreach)
  library(doParallel)
  library(rhdf5)
  
  source(paste(wd, "src/polygonize.R", sep=""))
  source(paste(wd, "src/extract_data.R", sep=""))
  source(paste(wd, "src/get_itcs_in_tile.R", sep=""))
  
  cr_per_path <- rep(NA, length(hps_f))
  for(ii in 1:length(hps_f)){
    cr_per_path[ii] <- dim(get_itcs_in_tile(hps_f[ii], centroids=centroids, NeonSites = NeonSites, f_path=f_path))[1]
  }
  #cr_per_path <- lapply(hps_f, get_number_itcs_in_tile, centroids=centroids, f_path=f_path)
  #cr_per_path<-unlist(cr_per_path)
  #remove paths with no data
  clean_hps <- hps_f[cr_per_path>0]
  cr_per_path <- cr_per_path[cr_per_path>0]  
  
  #now, you want to loop only trough the ith in the clean_hps/cr_per_path
  
  for(z in 1:length(clean_hps)){ 
    
    registerDoSEQ()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())

    #results <- foreach(z = clean_hps, .combine = 'cbind', .verbose = T) %:%
    results <- foreach(mm = 1:cr_per_path[z], .verbose = T) %dopar% {
    #for(mm in 1:cr_per_path[z]){
      source(paste(wd, "src/polygonize.R", sep=""))
      source(paste(wd, "src/extract_data.R", sep=""))
      source(paste(wd, "src/get_itcs_in_tile.R", sep=""))
      
      # print(paste(mm, z, clean_hps[z]))
      # print(search())
      #"/Users/sergiomarconi/anaconda3/bin/"
      token = unlist(strsplit(clean_hps[z], split = "_"))[6] #3 for old tiles
      itcextract <- get_itcs_in_tile(clean_hps[z], centroids, NeonSites= NeonSites, f_path)
      tryCatch({chm_pt <- list.files(chm_f, pattern = paste(as.integer(itcextract$easting[mm]/1000)*1000, "_",
                                                            as.integer(itcextract$northing[mm]/1000)*1000, sep=""))
      chm <- (paste(chm_f, chm_pt, sep="/"))
      itcExtract(x = itcextract[mm,], f= paste(f_path,clean_hps[z],sep = "/"), chm = chm, buffer = buffer,
                 epsg = epsg, token = token, wd = wd, pybin = "/home/s.marconi/.conda/envs/quetzal3/bin")
      print(paste(mm, z, clean_hps[z]))
      },error=function(e){})
    }
    stopCluster(cl)
  }
  #return(results)
} 
