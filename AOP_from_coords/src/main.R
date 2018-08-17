
# pt <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP1.30003.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L1/DiscreteLidar/Classified_point_cloud/", sep="")
# f_path <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP1.30006.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L1/Spectrometer/H5/", sep="")
# chm_f <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP3.30003.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L3/DiscreteLidar/CHM/", sep="")
# 

rm(list=ls())

get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32617, 32616, 32616, 32612, 32613, 32617) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", "JORN", "OSBS")
  return(dictionary[colnames(dictionary)==utm])
}

wd = "/ufrc/ewhite/s.marconi/Chapter1/AOP_from_coords/"
inputs <- "/ufrc/ewhite/s.marconi/Chapter1/TOS_Retriever/out/utm_dataset.csv"
options(scipen=999)

library(readr)
library(dplyr)
source(paste(wd, "src/itcSegmentFast.R", sep=""))
source(paste(wd, "src/extract_data.R", sep=""))
source(paste(wd, "src/extract_crown_data.R", sep=""))

dataset <- read_csv(inputs) %>%
  dplyr::select(individualID, taxonID.x, collectDate, siteID.x, plotID.x, domainID.x, utmZone, UTM_E, UTM_N) %>%
  unique
colnames(dataset) <-  c("treeID","taxaID","collectDate","siteID","plotID","domainID", "utmZone","easting", "northing")
dataset[which(dataset$siteID=="STEI"),] <- convert_stei(dataset[which(dataset$siteID=="STEI"),])

for(NeonSites in unique(dataset$siteID)){
  tryCatch({
    centroids <- dataset[dataset$siteID %in% NeonSites,] 
    year <- unlist(strsplit(as.character(unique(dataset$collectDate)), split = "-"))[1]
    
    if(NeonSites %in% c( "GRSM", "ORNL")){
      year <- 2016
    }
    tileID <- unique(cbind(as.character(as.integer(centroids$easting/1000)*1000), 
                           as.character(as.integer(centroids$northing/1000)*1000)))
    
    epsg <- get_epsg_from_utm(unique(centroids$siteID))

    pt <- paste("//orange/ewhite/NeonData/", NeonSites, "/DP1.30003.001/", year, "/FullSite/", unique(centroids$domainID), "/", 
                "/", year,"_", NeonSites, "/", "L1/DiscreteLidar/ClassifiedPointCloud/", sep="")
    f_path <- paste("//orange/ewhite/NeonData/", NeonSites, "/DP1.30006.001/", year, "/FullSite/", unique(centroids$domainID), "/", 
                    "/", year,"_", NeonSites,"/", "L1/Spectrometer/H5/", sep="")
    
    chm_f <- paste("//orange/ewhite/NeonData/", NeonSites, "/DP1.30003.001/", year, "/FullSite/", unique(centroids$domainID), "/",
                   "/", year,"_", NeonSites, "/", "L3/CHM/", sep="")
    
    crownITC(pt, wd = wd, pattern = paste(tileID[,1], "_", tileID[,2], sep=""), epsg = epsg, cores = 32, chm_f = chm_f)
    hps_f = list.files(f_path)
    
    #hps_f = NULL, f_path = NULL, chm_f = NULL, epsg=NULL, buffer = 20, cores = 2
    extract_crown_data(centroids = centroids, hps_f = hps_f, f_path = f_path, chm_f = chm_f, epsg=epsg, wd = wd,NeonSites=NeonSites, cores = 32)
  },error=function(e){})
  
}
