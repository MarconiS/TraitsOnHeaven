#extract data from 2015 GRSM
get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32617, 32616, 32616, 32612, 32613, 32617, 32617) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", "JORN", "OSBS", "MLBS")
  return(dictionary[colnames(dictionary)==utm])
}

wd = "./AOP_from_coords/"
inputs <- "./TOS_retriever/out/utm_dataset.csv"
#inputs <- "./TOS_retriever/out/missing_utm.csv"
#inputs <- "./AOP_from_coords/inputs/Dimensions_centroids.csv"
options(scipen=999)

library(readr)
library(dplyr)
source(paste(wd, "src/itcSegmentFast.R", sep=""))
source(paste(wd, "src/extract_data.R", sep=""))
source(paste(wd, "src/extract_crown_data.R", sep=""))

GRSM_2015_extract <- function(year = 2015){
  pt_h5 <- "/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D7/GRSM/2015/GRSM_L1/GRSM_Spectrometer/Reflectance/2015080313/"
  files = list.files(pt_h5, pattern = ".h5")
  inputs <- "./TOS_retriever/out/utm_dataset.csv"
  data_pool <- read_csv(inputs) %>%
    dplyr::select(individualID, eventID, taxonID, siteID,  domainID,utmZone, stemDiameter, height, maxCrownDiameter, UTM_E, UTM_N) %>%
    unique
  
  smplY <- str_sub(data_pool$eventID, start= -4)
  centroids <- data_pool[smplY==year,] %>%
    filter(siteID == "GRSM") %>% unique
  colnames(centroids)[colnames(centroids) %in% c("UTM_E", "UTM_N")] <-  c("easting", "northing")
  tileID <- unique(cbind(as.character(as.integer(centroids$easting/1000)*1000),
                         as.character(as.integer(centroids$northing/1000)*1000)))
  epsg <- get_epsg_from_utm(unique(centroids$siteID))
  
  pt <- paste("/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D7/GRSM/2015/GRSM_L1/GRSM_Lidar/Classified_point_cloud/")
  f_path <- paste("/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D7/GRSM/2015/GRSM_L1/GRSM_Spectrometer/Reflectance/d1/", sep="")
  
  chm_f <- paste("/Volumes/groups/lab-white-ernest/NEON_AOP/AOP_1.3a_w_WF_v1.1a/1.3a/D7/GRSM/2015/GRSM_L3/GRSM_Lidar/CHM/")
  hps_f = list.files(f_path, pattern = ".h5")
  extract_crown_data_old(centroids = centroids, hps_f = hps_f, f_path = f_path, 
                     chm_f = chm_f, epsg=epsg, wd = wd,NeonSites="GRSM", cores = 2)
  
  }