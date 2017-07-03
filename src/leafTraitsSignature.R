# script for calibrating the leaf lv spectr PLSR models and ncertainty estimates.
#---------------- Load required libraries ---------------------------------------------------------#

rm(list=ls(all=TRUE))   # clear workspace
library(pls)
library(plotrix)
library(prospectr)
library (leaps)
library (glmnet)
library (boot)
setwd("/Users/sergiomarconi/Documents/Classes/Dropbox/FinalProject_Marconi/FinalProject_Marconi")
source(paste(getwd(), '/src.R', sep=""))
closeAll()

set.seed(1)
# Script options
pls.options(plsralg = "oscorespls")
Create.fig <- T

# Set ouput directory
out.dir = paste(getwd(), "/Outputs/", sep="") 
in.dir = paste(getwd(), "/Inputs/", sep="") 
aug.spectra <- imp.spectra('ASD_Chem.csv', in.dir)
norm_wholeC <- imp.spectra('CrownPix_norm.csv', in.dir)



X <- aug.spectra[,19:length(aug.spectra[1,])]
filter <- apply(X, 2, max)
which(filter > 1)
X <- X[-which(filter > 1)]
Y <- aug.spectra[,7:18]
#diag.spectra()
X_corr = cor(as.matrix(X), Y)
# sp.corr(X, Y, out.dir)

# # Build full PLSR dataset: from wawelength 500 to 2048
# bb <- read_csv(paste(in.dir, "bad_bands.csv", sep=""), col_names = F)
# aug.X <- data.frame(aug.spectra$Site,Y,X[,51:1050])
aug.X <- data.frame(aug.spectra$site, Y, X)
#--------------------------------------------------------------------------------------------------#
# Subset data into cal/val by site
eval.set <- cut.set(aug.X,out.dir)
train.data <- eval.set$train
test.data <- eval.set$test

#--------------------------------------------------------------------------------------------------#
# Run calibration PLSR analysis to select optimal number of components
pls.mod.train <- pls.cal(train.data, rep(15,12))
#--------------------------------------------------------------------------------------------------#
test.comp <- opt.comps(pls.mod.train, test.data,use.press = TRUE)
#calculate number of components given min test PRESS or RMSEP
optim.ncomps <- opt.comps(pls.mod.train, train.data)
optim.ncomps = rep(5,12)
#--------------------------------------------------------------------------------------------------#
pred.val.data <- predict.pls(pls.mod.train, test.data, optim.ncomps)
out.data <- res.out(pred.val.data, test.data)
train.val.data <- predict.pls(pls.mod.train, train.data, optim.ncomps)
train.out <- res.out(train.val.data, train.data)

#-- run a generalized linear model with all the features: clearly not a realistic option, ---------#
#-- since there are more features than observations, we are in front of the " p > n issue ---------#
#--   Hence, this model is not an option, and will not be used for the actual analyses    ---------#

full.reg<- full.asd(test.data,train.data)
#-- run lasso -------------------------------------------------------------------------------------#
lasso.reg <- lasso.asd(test.data,train.data)
#-- run step forward selection -------------------------------------------------------------------------------------#
fwd.reg <- step.fwd(test.data,train.data)

#-------------------------------------------------------------------------------------------------------------------#
# ------------------##############              PLOTTING AND TABLEING               ##############------------------#
#-------------------------------------------------------------------------------------------------------------------#
setwd(out.dir)
#I apology for the lack of organization in the plotting scheme
R.dry_mass_g <- c(lasso.reg$dry_mass_g$r2, mean(out.data$dry_mass_g$R2))
R.leaf_area_cm2 <- c(lasso.reg$leaf_area_cm2$r2,mean(out.data$leaf_area_cm2$R2))
R.SLA_cm2.g <- c(lasso.reg$SLA_cm2.g$r2, mean(out.data$SLA_cm2.g$R2))
R.LMA_g.m2 <- c( lasso.reg$LMA_g.m2$r2, mean(out.data$LMA_g.m2$R2))
R.d13C <- c(lasso.reg$d13C$r2, mean(out.data$d13C$R2))
R.d15N <- c(lasso.reg$d15N$r2, mean(out.data$d15N$R2))
R.C_pct <- c(lasso.reg$C_pct$r2,mean(out.data$C_pct$R2))
R.N_pct <- c(lasso.reg$N_pct$r2, mean(out.data$N_pct$R2))
R.P_pct <- c(lasso.reg$P_pct$r2, mean(out.data$P_pct$R2))
R.C_N <- c(lasso.reg$C_N$r2, mean(out.data$C_N$R2))
R.C_P <- c(lasso.reg$C_P$r2, mean(out.data$C_P$R2))
R.N_P <- c(lasso.reg$N_P$r2, mean(out.data$N_P$R2))

mse.N <- c(sqrt(fwd.reg$totalN$mse), sqrt(lasso.reg$totalN$mse), sqrt(mean(out.data$totalN$MSE)))
mse.C <- c(sqrt(fwd.reg$totalC$mse), sqrt(lasso.reg$totalC$mse), sqrt(mean(out.data$totalC$MSE)))
mse.Cell <- c(sqrt(fwd.reg$Cellulose$mse), sqrt(lasso.reg$Cellulose$mse), sqrt(mean(out.data$Cellulose$MSE)))
mse.Lig <- c(sqrt(fwd.reg$Lignin$mse), sqrt(lasso.reg$Lignin$mse), sqrt(mean(out.data$Lignin$MSE)))
mse.DW <- c(sqrt(fwd.reg$dry_wgt$mse), sqrt(lasso.reg$dry_wgt$mse), sqrt(mean(out.data$dry_wgt$MSE)))

#train stats
tr.R.N <- c(fwd.reg$totalN$tr.R2, lasso.reg$totalN$tr.R2, mean(train.out$totalN$R2))
tr.R.C <- c(fwd.reg$totalC$tr.R2, lasso.reg$totalC$tr.R2,mean(train.out$totalC$R2))
tr.R.Cell <- c(fwd.reg$Cellulose$tr.R2, lasso.reg$Cellulose$tr.R2, mean(train.out$Cellulose$R2))
tr.R.Lig <- c(fwd.reg$Lignin$tr.R2, lasso.reg$Lignin$tr.R2, mean(train.out$Lignin$R2))
tr.R.DW <- c(fwd.reg$dry_wgt$tr.R2, lasso.reg$dry_wgt$tr.R2, mean(train.out$dry_wgt$R2))

tr.mse.N <- c(sqrt(fwd.reg$totalN$tr.MSE), sqrt(lasso.reg$totalN$tr.MSE), sqrt(mean(train.out$totalN$MSE)))
tr.mse.C <- c(sqrt(fwd.reg$totalC$tr.MSE), sqrt(lasso.reg$totalC$tr.MSE), sqrt(mean(train.out$totalC$MSE)))
tr.mse.Cell <- c(sqrt(fwd.reg$Cellulose$tr.MSE), sqrt(lasso.reg$Cellulose$tr.MSE), sqrt(mean(train.out$Cellulose$MSE)))
tr.mse.Lig <- c(sqrt(fwd.reg$Lignin$tr.MSE), sqrt(lasso.reg$Lignin$tr.MSE), sqrt(mean(train.out$Lignin$MSE)))
tr.mse.DW <- c(sqrt(fwd.reg$dry_wgt$tr.MSE), sqrt(lasso.reg$dry_wgt$tr.MSE), sqrt(mean(train.out$dry_wgt$MSE)))


#euclidian length of each coefficient vector, as proposed by Zhao et al., 2012
beta_ecDist.N <- c(sqrt(sum(fwd.reg$totalN$coeff^2)), sqrt(sum(lasso.reg$totalN$coeff^2)), 
                   sqrt(sum(coef(pls.mod.train$totalN, ncomp=optim.ncomps[j])^2)))
beta_ecDist.C <- c(sqrt(sum(fwd.reg$totalC$coeff^2)), sqrt(sum(lasso.reg$totalC$coeff^2)), 
                   sqrt(sum(coef(pls.mod.train$totalC, ncomp=optim.ncomps[j])^2)))
beta_ecDist.Cell <- c(sqrt(sum(fwd.reg$Cellulose$coeff^2)), sqrt(sum(lasso.reg$Cellulose$coeff^2)), 
                      sqrt(sum(coef(pls.mod.train$Cellulose, ncomp=optim.ncomps[j])^2)))
beta_ecDist.Lig <- c(sqrt(sum(fwd.reg$Lignin$coeff^2)), sqrt(sum(lasso.reg$Lignin$coeff^2)), 
                     sqrt(sum(coef(pls.mod.train$Lignin, ncomp=optim.ncomps[j])^2)))
beta_ecDist.DW <- c(sqrt(sum(fwd.reg$dry_wgt$coeff^2)), sqrt(sum(lasso.reg$dry_wgt$coeff^2)), 
                    sqrt(sum(coef(pls.mod.train$dry_wgt, ncomp=optim.ncomps[j])^2)))

#ncomp/pars

vars.N <- c(fwd.reg$totalN$npar, length(summary(lasso.reg$totalN$coeff)[,1]), optim.ncomps[1])
vars.C <- c(fwd.reg$totalC$npar, length(summary(lasso.reg$totalC$coeff)[,1]), optim.ncomps[2])
vars.Cell <- c(fwd.reg$Cellulose$npar, length(summary(lasso.reg$Cellulose$coeff)[,1]), optim.ncomps[3])
vars.Lig <- c(fwd.reg$Lignin$npar, length(summary(lasso.reg$Lignin$coeff)[,1]), optim.ncomps[4])
vars.DW <- c(fwd.reg$dry_wgt$npar, length(summary(lasso.reg$dry_wgt$coeff)[,1]), optim.ncomps[5])

#which coeffs
tt.SR_Las <-  list()
tt.Las_PLS<-  list()
for (j in 1:12) {
  traits <- train.data[,2:13]
  leaf.spec <- train.data[,14:length(train.data[1,])]
  leaf.spec <- standardNormalVariate(X = t(diff(t(log(leaf.spec)),differences=1, lag=3)))
  data = data.frame(traits[,j], leaf.spec)
  colnames(data)[1] <- "Y"
  test.traits <- test.data[,2:13]
  test.spec <- test.data[,14:length(test.data[1,])]
  test.spec <- standardNormalVariate(X = t(diff(t(log(test.spec)),differences=1, lag=3)))
  data.test <- data.frame(test.traits[j], test.spec)
  colnames(data.test)[1] <- "Y"
  
  # pdf(paste(names(Y)[j],'_coeffs.pdf',sep=""),height=12,width=8)
  # par(mfrow =c(3,1))
  # sfs.coef.plot <- as.numeric(gsub("^.*?X","",names(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$coeff[-1]",sep=""))))))
  # sfs.coef.plot<- data.frame(cbind(as.numeric(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$coeff[-1]",sep="")))), sfs.coef.plot))
  # colnames(sfs.coef.plot)=c("coeff", "band")
  # plot(sfs.coef.plot$band ,sfs.coef.plot$coeff,ylim = c(min(sfs.coef.plot$coeff),max(sfs.coef.plot$coeff)), xlim = c(500,2500),main=names(Y)[j],ylab="SF coef",xlab="WAVELENGTH (nm)", cex=2)
  # abline(h = 0, col='red')
  # 
  # Lasso <- data.frame(cbind(summary(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$coeff",sep=""))))$x, 
  #                           summary(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$coeff",sep=""))))$i*2 + 498))
  # colnames(Lasso)=c("coeff", "band")
  # plot(Lasso$band,Lasso$coeff, ylim = c(min(Lasso$coeff),max(Lasso$coeff)), xlim = c(500,2500),ylab="Lasso coef",xlab="WAVELENGTH (nm)", cex=2)
  # abline(h = 0, col='red')
  # 
  # waves <-seq(506,2498,2)
  # coefs.pls <- coef(eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep=""))), ncomp=optim.ncomps[j], intercept=FALSE)
  # plot(waves,coefs.pls,cex=2,xlab="WAVELENGTH (nm)",ylab="Plsr coef",  ylim = c(min(coefs.pls),max(coefs.pls)))
  # lines(waves,coefs.pls,lwd=0.5)
  # abline(h=0,lty=2,col="red")
  # dev.off()
  
  # y_i ~ y_hat -------------------------------------------------------------
  pdf(paste(names(Y)[j],'_corr.pdf',sep=""),height=12,width=8)
  par(mfrow =c(3,2))
  # #SFS
  # plot(train.data[,j+1], predict.regsubsets(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$mod",sep=""))), data, 
  #                                           id=as.integer(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$npar",sep=""))))), xlab = "Measured", ylab="Predicted", main = "SFS train", cex=2)
  # abline(0,1, col="red", lty=2)
  # 
  # plot(test.data[,j+1], predict.regsubsets(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$mod",sep=""))), data.test, 
  #                                          id=as.integer(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$npar",sep=""))))), xlab = "Measured", ylab="Predicted", main = "SFS test", cex=2)
  # abline(0,1, col="red", lty=2)
  #LASSO
  plot(train.data[,j+1], predict(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$mod",sep=""))), 
                                 s=eval(parse(text = paste('lasso.reg$',names(Y)[j],"$lambda",sep=""))) ,newx=as.matrix(leaf.spec)), xlab = "Measured", ylab="Predicted", main = "LASSO train", cex=2)
  abline(0,1, col="red", lty=2)
  plot(test.data[,j+1], predict(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$mod",sep=""))), 
                                s=eval(parse(text = paste('lasso.reg$',names(Y)[j],"$lambda",sep=""))) ,newx=as.matrix(test.spec)), xlab = "Measured", ylab="Predicted", main = "LASSO test", cex=2)
  abline(0,1, col="red", lty=2)
  # PLSR 
  plot(eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep=""))),ncomp=optim.ncomps[j], asp=1, xlab = "Measured", ylab="Predicted",main="PLS train", cex=2)
  abline(0,1, col="red", lty=2)
  plot(test.data[,j+1],eval(parse(text = paste('pred.val.data$',names(Y)[j],sep=""))), xlab = "Measured", ylab="Predicted", main = "PLS test", cex=2)
  abline(0,1, col="red", lty=2)
  dev.off()
  
  # #paired ttest
  # tt.SR_Las[[names(train.data[,2:6])[j]]]  <- t.test(test.data[,j+1]-predict.regsubsets(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$mod",sep=""))), data.test, 
  #                                        id=as.integer(eval(parse(text = paste('fwd.reg$',names(Y)[j],"$npar",sep=""))))), 
  #                                        test.data[,j+1]-predict(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$mod",sep=""))), 
  #                             s=eval(parse(text = paste('lasso.reg$',names(Y)[j],"$lambda",sep=""))),
  #                             newx=as.matrix(test.spec)), paired=TRUE)
  # print(tt.SR_Las[[names(train.data[,2:6])[j]]]$p.value)
  # tt.Las_PLS[[names(train.data[,2:6])[j]]]  <- t.test(test.data[,j+1]-predict(eval(parse(text = paste('lasso.reg$',names(Y)[j],"$mod",sep=""))), 
  #                               s=eval(parse(text = paste('lasso.reg$',names(Y)[j],"$lambda",sep=""))),
  #                               newx=as.matrix(test.spec)), test.data[,j+1]-eval(parse(text = paste('pred.val.data$',names(Y)[j],sep=""))), paired=TRUE)
  # print(tt.Las_PLS[[names(train.data[,2:6])[j]]]$p.value)
}




