library(rgdal)
library(sp)
fp = "/Users/sergiomarconi/Documents/Data/OSBS_crown_polygons/"
up <- readOGR(fp, "OSBS_sample_polygons_edits_Feb2018")
centroids <- data.frame(getSpPPolygonsLabptSlots(up))
colnames(centroids) <- c("UTM_N", "UTM_E")
centroids$site <- "OSBS"
centroids$individualID <- up@data$ID


fp = "/Users/sergiomarconi/Documents/Data/TALL_crown_polygons/"
up <- readOGR(fp, "TALL_sample_crowns_edits_May2018")
centroids2 <- data.frame(getSpPPolygonsLabptSlots(up))
colnames(centroids2) <- c("UTM_N", "UTM_E")
centroids2$site <- "TALL"
centroids2$individualID <- up@data$tree

fp = "/Users/sergiomarconi/Documents/Data/MLBS_crown_polygons/"
up <- readOGR(fp, "MLBS_sample_polygons")
centroids3 <- data.frame(getSpPPolygonsLabptSlots(up))
colnames(centroids3) <- c("UTM_N", "UTM_E")
centroids3$site <- "MLBS"
centroids3$individualID <- up@data$tag

final <- rbind(centroids, centroids2, centroids3)
readr::write_csv(final, "~/Documents/Data/Dimensions_centroids.csv")
