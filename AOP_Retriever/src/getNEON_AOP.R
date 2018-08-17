#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = FALSE)


library(devtools)

install_github("NEONScience/NEON-utilities/neonUtilities", dependencies=TRUE)
library(neonUtilities)

print(args)

prd = args[7]
ste = args[8]
yr = args[9]
byFileAOP(prd, site = ste, year = yr, check.size=F, savepath = paste("/orange/ewhite/NeonData/",  ste, sep=""))

