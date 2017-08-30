rm(list=ls(all=TRUE))   # clear workspace
setwd('/Users/sergiomarconi/Documents/Projects/TraitsOnHeaven/')
#---------------- Load required libraries ---------------------------------------------------------#
library(pls)
#library(plotrix)
#library(prospectr)
library (leaps)
library (glmnet)
library (boot)
library(ggplot2)
library(reshape2)
library(readr)
#library(neonAOP)
library(parallel)
library(doMC)
library(doParallel)

source(paste(getwd(), '/src/modelSrc.R', sep=""))
source(paste(getwd(), '/src/src.R', sep=""))

# arguments
#parameters (to put in the arguments form to automize)
loops <- 300
NeonSite <- "OSBS"
setwd('~/Documents/Projects/TraitsOnHeaven/')
out.dir = paste(getwd(), NeonSite,"outputs/", sep="/")
in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")
names <- c('abs_2yr',	'abs_5yr',	'rel_2yr',	'rel_5yr')

#names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
#debugging
laps =rounds=1
path = paste(getwd(), NeonSite, sep="/")
setwd(path)
main <- function(loops=100, out.dir = paste(getwd(), "/outputs/", sep=""),in.dir = paste(getwd(), 
          "/inputs/", sep=""),  names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")){
  
  CrownIDS <- read_csv(paste(path, "inputs", "Spectra", 'CrownIDS.csv', sep="/"))
  cid <- eval(parse(text = paste("CrownIDS", NeonSite, sep="$")))
  cid <- cid[!is.na(cid)]
  set.seed(14)
  cid.train <- sample(cid, length(cid)*1)
  cid.test <- cid[!(cid %in% cid.train)]
  write.csv(cid.train, paste(path,'inputs/trainID.csv', sep="/"))
  write.csv(cid.test, paste(path,'inputs/testID.csv', sep="/"))
  #getTrainCrownsx
  alltratis <- read.csv(paste(path, "/inputs/Spectra/CrownTraits.csv",sep=""))
  alltratis.tr <- alltratis[alltratis$pixel_crownID %in% cid.train,]  
  alltratis.ts <- alltratis[alltratis$pixel_crownID %in% cid.test,]  
  write.csv(alltratis.tr, paste(path,'inputs/Spectra/trainCrownTraits.csv', sep="/"))
  write.csv(alltratis.ts, paste(path,'inputs/Spectra/testCrownTraits.csv', sep="/"))
  
  nCrowns = length(cid.train) #85
  
  setwd(path)
  normalize()
  for(rounds in 1:3){
    setwd(path)
    rPix(rounds, loops, unqCrown = cid.train, names, path=path)
    PLSandLasso(rounds,loops, names)
    setwd(path) #was in.dir
    DiscardAndRerun(names,rounds, loops, nCrowns)
  }
  #return best 100 of last filtering
  storeFinalSet(in.dir, out.dir, rounds = 10, names, nentries = 40)
}