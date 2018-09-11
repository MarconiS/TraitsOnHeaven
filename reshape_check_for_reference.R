library(rgdal)
library(raster)
example <- stack("/Users/sergiomarconi/Dropbox (UFL)/NEON_tiles/outputs/itcTiff/140311_OSBS_DIOSP_NEON.PLA.D03.OSBS.01018.tif")
plotRGB(example, 17, 55, 113, stretch="lin")
shpITC <- readOGR("/Users/sergiomarconi/Documents/tmp/", "140311_OSBS_DIOSP_NEON.PLA.D03.OSBS.01018")
spectra <- readr::read_csv("/Users/sergiomarconi/Dropbox (UFL)/NEON_tiles/outputs/Spectra/140311_OSBS_DIOSP_NEON.PLA.D03.OSBS.01018.csv") 
plot(shpITC, add = T, border = "red")
sp2 <- raster::extract(example, shpITC)
sp2

test <- as.matrix(example)
dim(test)

a <- apply(test, 1, sum)
b <- apply(spectra, 1, sum)
c <-match(b, a)
remake <- rep(0, 1600)
remake[c] <- 1
dim(remake) <- c(40,40)
remake <- raster(t(remake))
plot(remake)
plot(shpITC)

#test if is a reshape problem
clip_raster<- raster::crop(example, extent(shpITC))
plotRGB(clip_raster, 15,55,113, stretch="lin")
plotRGB(example, 15,55,113, stretch="lin")
plot(shpITC, add = T, border = "red")
plotRGB(clip_raster, 15,55,113, stretch="lin")
plot(shpITC, add = T, border = "red")

result <- raster::extract(clip_raster, shpITC)
dim(result[[1]])
