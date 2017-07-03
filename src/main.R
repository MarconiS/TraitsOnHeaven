# # script for calibrating the leaf lv spectr PLSR models and ncertainty estimates.
# #---------------- Load required libraries ---------------------------------------------------------#
# 
# rm(list=ls(all=TRUE))   # clear workspace
# setwd("/Users/sergiomarconi/Documents/PhD/Projects/MappingTraits")
# library(pls)
# library(plotrix)
# library(prospectr)
# library (leaps)
# library (glmnet)
# library (boot)
# library(ggplot2)
# library(reshape2)
# library(readr)
# # source(paste(getwd(), '/src/DiscardandRerun.R', sep=""))
# # source(paste(getwd(), '/src/extractRandom.R', sep=""))
# # source(paste(getwd(), '/src/normAditya.R', sep=""))
# source(paste(getwd(), '/src/src.R', sep=""))
# source(paste(getwd(), '/src/itcSrc.R', sep=""))
# 
# 
# varPlot <- function(dat){
#   dat <- as.data.frame(dat)
#   colnames(dat) <- c("LMA_g.m2", ".LMA_g.m2","d13C",".d13C","d15N", ".d15N","C_pct", ".C_pct","N_pct",".N_pct", "P_pct", ".P_pct")
#   colors = rep(c("red","blue"),6)
#   boxplot(dat, col=colors, ylim=c(-1,1))
#   abline(h = T, b = 0,a=0)
# }
# closeAll()
# 
# #parameters (to put in the arguments form to automize)
# Create.fig <- T
# iiround <- T
# loops <- 1000
# site <- "OSBS"
# #Set ouput directory
# out.dir = paste(getwd(), "/outputs/", sep="")
# in.dir = paste(getwd(), "/inputs/", sep="")
# names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
# laps =1
# 
# # PLSandLasso <- function(out.dir = paste(getwd(), "/outputs/", sep=""),in.dir = paste(getwd(), "/inputs/", sep=""),
# #                Create.fig = T,iiround = F, site = "OSBS",loops = 1000){
# PLSandLasso <- function(rounds, loops, out.dir = paste(getwd(), "/outputs/", sep=""),
#                         names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),
#                         site = "TALL", in.dir = paste(getwd(), "/inputs/", sep=""),
#                         Create.fig = F){
#   setwd(in.dir)
#   cr.Traits <- read.csv(paste(in.dir, "spectra/CrownTraits.csv",sep=""))
# 
#   nCrowns <- dim(cr.Traits)[1]
#   
#   for (j in 1:6){
#     R = matrix(NA, ncol=4, nrow=loops)
#     
#     for(laps in 1:loops){
#       aug.spectra <- imp.spectra(paste('Bootstrap_',rounds,'/onePix1Crown_',names[j], laps, '.csv', sep = ''), in.dir)
#       
#       aug.spectra$X.1 = NULL
#       aug.spectra$X.2 = NULL
#       aug.spectra$X.3 = NULL
#       aug.spectra$X.4 = NULL
#       aug.spectra$X.5 = NULL
#       aug.spectra$X.6 = NULL
#       aug.spectra$X.7 = NULL
#       aug.spectra$X.8 = NULL
#       aug.spectra$X.9 = NULL
#       aug.spectra$X.10 = NULL
#       
#       aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
#       X <- aug.spectra[,20:length(aug.spectra[1,])]
#       filter <- apply(X, 2, max)
#       if(any(filter > 1)){
#         X <- X[-which(filter > 1)]
#       }
#       X=X[, colSums(is.na(X)) == 0]
#       Y <- aug.spectra[,10:15]
#       
#       X_corr = cor(as.matrix(X), Y)
#       if(Create.fig){
#         
#         wl <- read.csv(paste(in.dir, "spectra/bands.csv", sep=''))
#         X_corr.wl = cbind(wl, X_corr)
#         X_corr.wl = melt(X_corr.wl, "wl")
#         ggplot(X_corr.wl, aes(wl, value)) +  geom_line() + facet_wrap(~variable, scales = "free")
#         
#         sp.trait <- cbind(aug.spectra$name,  Y)
#         tr.melt = melt(sp.trait, "aug.spectra$name")
#         colnames(tr.melt) = c("species", "trait", "value")
#         ggplot(tr.melt, aes(species, value)) +  geom_boxplot() + 
#           facet_wrap(~trait, scales = "free") + theme(plot.subtitle = element_text(vjust = 1), 
#                                                       plot.caption = element_text(vjust = 1), axis.text.x = element_text(angle = 90))
#         ggplot(tr.melt, aes(x=value)) +  geom_density() + 
#           facet_wrap(~trait, scales = "free") + theme(plot.subtitle = element_text(vjust = 1), 
#                                                       plot.caption = element_text(vjust = 1), axis.text.x = element_text(angle = 90))
#         
#         X.temp = X;
#         colnames(X.temp) <- t(wl[,1])
#         sp.spectra <- (cbind(aug.spectra$name, X.temp))
#         tr.melt = melt(sp.spectra, "aug.spectra$name")
#         colnames(tr.melt) = c("species", "band", "value")
#         
#       }
#       aug.X <- data.frame(aug.spectra$name, Y, X)
#       # Subset data into cal/val by site
#       eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
#       train.data <- eval.set$train
#       test.data <- eval.set$test
#       
#       print(paste(laps, names(Y[j])))
#       
#       # Run calibration PLSR analysis to select optimal number of components
#       pls.mod.train <- pls.cal(train.data, 40, j, norm = F)
#       #calculate number of components given min test PRESS or RMSEP
#       optim.ncomps <- opt.comps(pls.mod.train, j)
# 
#       pred.val.data <- predict.pls(pls.mod.train, test.data, optim.ncomps,j, norm = F)
#       out.data <- res.out(pred.val.data, test.data, j)
#       
#       #-- run lasso -------------------------------------------------------------------------------------#
#       lasso.reg <- lasso.asd(test.data,train.data, j, norm = F)
#       print.coeffs(pls.mod.train,optim.ncomps, lasso.reg, j, out.dir) 
#       
#       setwd(out.dir)
#       
#       #I apology for the lack of organization in the plotting scheme
#       R[laps,1:2] <- c(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$r2",  sep = ""))),
#                        mean(eval(parse(text = paste("out.data$", names(Y[j]), "$R2", sep = "")))))
#       R[laps,3:4] <- c(sqrt(eval(parse(text = paste("lasso.reg$", names(Y[j]),"$mse",  sep = "")))),
#                        sqrt(mean(eval(parse(text = paste("out.data$", names(Y[j]), "$MSE", sep = ""))))))
#     } 
#     assign(paste("R.", names(Y[j]),  sep = ""), R)
#   }
#   
#   R2est <- as.data.frame(cbind(R.LMA_g.m2[,1:2],R.d13C[,1:2], R.d15N[,1:2], R.C_pct[,1:2], R.N_pct[,1:2],R.P_pct[,1:2]))
#   colnames(R2est) <- c("LMA Lasso", "LMA PLS","δ13C Lasso","δ13C PLS","δ15N Lasso", "δ15N PLS","C Lasso", "C PLS","N Lasso","N PLS", "P Lasso", ".P_pct")
#   
#   MSEest <- as.data.frame(cbind(R.LMA_g.m2[,3:4], R.d13C[,3:4], R.d15N[,3:4],  R.C_pct[,3:4], R.N_pct[,3:4],R.P_pct[,3:4]))
#   colnames(MSEest) <- c("LMA_g.m2", ".LMA_g.m2","d13C",".d13C","d15N", ".d15N","C_pct", ".C_pct","N_pct",".N_pct", "P_pct", ".P_pct")
#   varPlot(R2est)
#   write.csv(R2est, file = paste(names[j],"DataBoost2pixR.csv",sep="_"))
#   write.csv(MSEest, file = paste(names[j],"Databoost2pixMSE.csv",sep="_"))
# }
# 
# 
# 
# 
# 
