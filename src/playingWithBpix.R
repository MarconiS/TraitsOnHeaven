allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("./inputs/Spectra/CrownPix.csv")
all.pixPerCrown <- allData %>%
  group_by(pixel_crownID) %>%
  summarise(n = n())
# Find bad bands
all_Bands=as.character(allBand$BandName)
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
#allData[,-c(1,2)] <- allData[,-c(1,2)] /10000
# Set bad bands to zero
ndvi <- (allData$band_90- allData$band_58)/(allData$band_58 + allData$band_90)
nir860 <- (allData$band_96 + allData$band_97)/2 
bad.pixs <- allData[which(ndvi < 0.7 | nir860 < .3),]
allData[which(ndvi < 0.7 | nir860 < .3),]=NA
#allData[which(nir860 < .3),] = NA
allData <- allData[complete.cases(allData), ]
allData[,bad_Bands]=NA
foo<- allData[complete.cases(allData), ]
bad.pixPerCrown <- allData %>%
  group_by(pixel_crownID) %>%
  summarise(n = n())
#extract crowns with no pixels
noCrData <- anti_join(all.pixPerCrown, bad.pixPerCrown, by = "pixel_crownID")

#remove any reflectance bigger than 1
pixel_crownID <- allData[,2]
allData <- allData[,-c(1,2)]
allData[allData>1]=NA
allData <- cbind(pixel_crownID,allData)

allData$pixel_crownID <- factor(allData$pixel_crownID)
foo <- melt(allData)
foo <- foo[, colnames(foo) %in% c("pixel_crownID", "variable", "value")]
ggplot(foo, aes(as.integer(variable), value)) + geom_smooth()

bad.pixs <- bad.pixs[complete.cases(bad.pixs), ]
bad.pixs[,bad_Bands]=NA
bad.pixs$pixel_crownID <- factor(bad.pixs$pixel_crownID)
foo.bad <- melt(bad.pixs)
foo.bad <- foo.bad[, colnames(foo.bad) %in% c("pixel_crownID", "variable", "value")]
ggplot(foo.bad, aes(x = as.integer(variable), y = value)) + geom_smooth() 

write_csv(bad.pixs, "badPix.csv")
write_csv(allData, "goodPix.csv")
