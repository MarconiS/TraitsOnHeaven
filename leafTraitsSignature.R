# script for calibrating the leaf lv spectr PLSR models and ncertainty estimates.
#---------------- Load required libraries ---------------------------------------------------------#

rm(list=ls(all=TRUE))   # clear workspace
library(pls)
library(plotrix)
library(prospectr)
source(paste(getwd(), '/src.R', sep=""))
closeAll()

# Script options
pls.options(plsralg = "oscorespls")
Create.fig <- FALSE

# Set ouput directory
out.dir = paste(getwd(), "/Outputs/LMA/", sep="") 
in.dir = paste(getwd(), "/Inputs/", sep="") 
aug.spectra <- imp.spectra('ASD_Chem.csv', in.dir)
X <- aug.spectra[,11:length(aug.spectra[1,])]
Y <- aug.spectra[,6:10]
#diag.spectra()
X_corr = cor(as.matrix(X), Y)
sp.corr(X, Y, out.dir)

# Build full PLSR dataset: from wawelength 500 to 2048
aug.X <- data.frame(aug.spectra$Site,Y,X[,51:1050])

#--------------------------------------------------------------------------------------------------#
# Subset data into cal/val by site
eval.set <- cut.set(aug.X,out.dir)
train.data <- eval.set$train
test.data <- eval.set$test
### Output cal data for VIP models
write.csv(train.data,file=paste(out.dir,"Splitted_train_Dataset.csv",sep=""),row.names=FALSE)
write.csv(test.data,file=paste(out.dir,"Splitted_test_Dataset.csv",sep=""),row.names=FALSE)
rm(eval.set)
#--------------------------------------------------------------------------------------------------#
# Run calibration PLSR analysis to select optimal number of components
pls.mod.train <- pls.cal(train.data, rep(15,5))
#--------------------------------------------------------------------------------------------------#
test.comp <- opt.comps(pls.mod.train, test.data,use.press = FALSE, h0 = TRUE)
#calculate number of components given min test PRESS or RMSEP
optim.ncomps <- opt.comps(pls.mod.train, train.datak)
#--------------------------------------------------------------------------------------------------#
pred.val.data <- predict.pls(pls.mod.train, test.data, optim.ncomps)
out.data <- res.out(pred.val.data, test.data)

