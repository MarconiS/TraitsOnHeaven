# once extracted csv data, make it a unique model
fetch_aop_csv <- function(in_dir = "./AOP_Retriever/outputs/itc_aop_spectra/", 
                          out_dir = "./Model_builder/inputs/raw" ){
  library("readr")
  library(dplyr)
  source("./Model_builder/src/utility.R")
  AOPsites = list.files(in_dir)
  CrownPix <- NULL
  for(NeonSite in AOPsites){
    data_files <- list.files(paste(in_dir, NeonSite, sep="/"), pattern = ".csv")
    for(ii in data_files){
      id <- unlist(strsplit(ii, "_"))
      tree_id <- paste(id[3], unlist(strsplit(id[4], ".c"))[1], sep="_")
      X = read_csv(paste(in_dir,NeonSite, ii, sep="/"))
      X = cbind(id[2],id[3], unlist(strsplit(id[4], ".c"))[1], X)
      col_bands = rep(NA, 373)
      col_bands[1:3] = c("site", "species", "crown_ID")
      for(nid in 4: 373){
        col_bands[nid] = paste("band_", unlist(strsplit(colnames(X[nid]), paste(tree_id, ".", sep="")))[2], sep="")
      }
      colnames(X) <- col_bands
      CrownPix <- rbind(CrownPix, X)
    }
  }
  #clean data
  CrownPix <- CrownPix %>%
    na.rm %>%
    optical.filter %>%
    na.rm %>%
    na.omit %>%
    no.val.rm 
  
  #normalize spectra
  CrownPix[grepl("band", names(CrownPix))] <- hiper_features(CrownPix)
  
  CrownPix <- cbind(CrownPix$crown_ID, category_to_dummy(CrownPix$site), CrownPix[grepl("band", names(CrownPix))])
  colnames(CrownPix)[1]<- "pixel_crownID"
  write_csv(CrownPix, paste(out_dir,"CrownPix.csv", sep="/"))
  return(unique(CrownPix$pixel_crownID))
}
