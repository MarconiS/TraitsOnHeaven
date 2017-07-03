rm(list=ls(all=TRUE))   # clear workspace
#---------------- Load required libraries ---------------------------------------------------------#
library(pls)
library(plotrix)
library(prospectr)
library (leaps)
library (glmnet)
library (boot)
library(ggplot2)
library(reshape2)
library(readr)
library(neonAOP)


setwd("/Users/sergiomarconi/OneDrive/Projects/OSBS")
source(paste(getwd(), '/src/modelSrc.R', sep=""))
#source(paste(getwd(), '/src/washSrc.R', sep=""))
source(paste(getwd(), '/src/src.R', sep=""))
#source(paste(getwd(), '/src/plotSrc.R', sep=""))

# arguments
#parameters (to put in the arguments form to automize)

loops <- 300
#Set ouput directory
out.dir = paste(getwd(), "/outputs/", sep="")
in.dir = paste(getwd(), "/inputs/", sep="")
names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
laps =rounds=1
nCrowns = 76 #85

main <- function(loops=100, out.dir = paste(getwd(), "/outputs/", sep=""),in.dir = paste(getwd(), 
          "/inputs/", sep=""),  names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")){
  
  normalize()
  for(rounds in 1:10){
    rPix(rounds, loops, names)
    PLSandLasso(rounds,loops, names)
    setwd(in.dir)
    DiscardAndRerun(names,rounds, loops, nCrowns)
  }
  
  #return best 100 of last filtering
  
}