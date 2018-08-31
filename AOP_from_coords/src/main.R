#!/usr/bin/env Rscript
argument_values = commandArgs(trailingOnly = FALSE)
print(argument_values)
# pt <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP1.30003.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L1/DiscreteLidar/Classified_point_cloud/", sep="")
# f_path <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP1.30006.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L1/Spectrometer/H5/", sep="")
# chm_f <- paste("//orange/ewhite/NeonData/", unique(centroids$siteID),"/DP3.30003.001/",year, "/FullSite/", unique(centroids$domainID), "/", 
#             year, "_", unique(centroids$siteID), "/L3/DiscreteLidar/CHM/", sep="")
# 

get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32617, 32616, 32616, 32612, 32613, 32617, 32617) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", "JORN", "OSBS", "MLBS")
  return(dictionary[colnames(dictionary)==utm])
}

wd = "./AOP_from_coords/"
inputs <- "./TOS_retriever/out/utm_dataset.csv"

#inputs <- "./AOP_from_coords/inputs/Dimensions_centroids.csv"
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
#dataset <- dataset[!dataset$siteID=="MLBS",]

NeonSites = argument_values[6]

# for(NeonSites in unique(dataset$siteID)){
#   #tryCatch({
centroids <- dataset[dataset$siteID %in% NeonSites,] %>%
  unique
year <- unlist(strsplit(as.character(unique(dataset$collectDate)), split = "-"))[1]

if(NeonSites %in% c( "GRSM", "ORNL")){
  year <- 2016
}else{
  year <- 2017
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

crownITC(pt, wd = wd, pttrn = paste(tileID[,1], "_", tileID[,2], sep=""),
         epsg = epsg, cores = 64, chm_f = chm_f,
         pybin = "/home/s.marconi/.conda/envs/quetzal3/bin")
#pybin = "/Users/sergiomarconi/anaconda3/bin/")
hps_f = list.files(f_path)
extract_crown_data(centroids = centroids, hps_f = hps_f, f_path = f_path, 
                   chm_f = chm_f, epsg=epsg, wd = wd,NeonSites=NeonSites, cores = 32)

