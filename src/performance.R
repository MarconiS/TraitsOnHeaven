rm(list=ls(all=TRUE))   # clear workspace

library("tsensembler")
library(readr)
library(pls)
#library(plsRglm)

library(tidyverse)
NeonSite = "ALL"
wd = "/ufrc/ewhite/s.marconi/Marconi2018/ModelBuild/"
setwd(wd)
source('../src/functions_build.R')
source('../src/src_build.R')
source( '../src/functions_infer.R')
source('../src/src_infer.R')
#source(paste(getwd(), '/src/pls_par.R', sep=""))


path = paste(wd, NeonSite, sep="/")
setwd(path)
out.dir = paste(getwd(), "/outputs/", sep="")
in.dir = paste(getwd(), "/inputs/", sep="")

nm = "N"
j="N_pct"
get100=T 
method = 'press'
normz =F
LOG=F
if(get100){
  loops = 1000
  if(normz){
    load(file = paste('./pls_',j, sep="" ))
  }else{
    load(file = paste('./pls_',j, sep="" ))
  }
  mod.aic <- rep(NA, loops)
  mod.r2 <- rep(NA, loops)
  weights = rep(0,loops)
  for(bb in 1: length(mod.aic)){
    #mod.aic[bb] <-foo[[bb]]$score$aic
    mod.r2[bb] <- 1 - sum((foo[[bb]]$pred$fit - (foo[[bb]]$score$data$Y.test))^2) /
      sum((foo[[bb]]$score$data$Y.test - mean(foo[[bb]]$score$data$Y.test))^2)
  }
  mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 100)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 100)])
  #mask <- cbind(which(mod.aic %in% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])
  #mask <- mask[order(mask[,2]),]
  plsglm <- foo[mask[,1]]
  save(plsglm, file = paste("pls_", nm,  sep = ""))
}else{
  load(file=paste("./pls_", nm,  sep = ""))
  
}
mod.r2=rep(0,length(plsglm))
mod.aic=rep(0,length(plsglm))

for(bb in 1: length(plsglm)){
  #mod.aic[bb] <-plsglm[[bb]]$score$aic #plsglm[[bb]]$aic
  mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
    sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
}
#mod.aic
mod.r2

if(method=='aic'){
  delta.aic <- mod.aic - min(mod.aic)
  weights <- softmax(-0.5*delta.aic)
}else if(method =='press'){
  delta.aic <- abs(mod.r2 - max(mod.r2))
  #weights <- delta.aic/sum(delta.aic)
  weights <- softmax(-0.5*delta.aic)
  
}
#
#weights <- softmax(-0.5*delta.aic)
#weights[mask[,1]] <- softmax(0.5*delta.press)
weights
#weights <- softmax(0.5*delta.press)

#out of bag
test.data.y <- read.csv(paste(in.dir, "Spectra/CrownTraits_outBag.csv", sep=""))
test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
if(LOG){test.data.y[names(test.data.y)==j] <- log(test.data.y[names(test.data.y)==j])}
test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag.csv", sep=""))
crownID <- as.data.frame(test.data.x$pixel_crownID)
test.data.x <- test.data.x[grepl("band", names(test.data.x))]
spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
if(normz){  
  spectra_log_dif_snv[spectra_log_dif_snv==0] <- 0.0000001
  spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv))),differences=1, lag=3))
}    
#test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
test.PLS = (as.matrix(spectra_log_dif_snv))

rm(output)
rm(output.sd)
output.daic =  output.sd.daic <- crownID
output.daic$yhat = output.sd.daic$yhat = 0
w.daic <- rep(0, length(weights))
w.daic <- weights
out <- list()
pred.val.data <- list()
for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  optim.ncomps <- plsglm[[bb]]$ncomp
  pred.val.data$fit<- predict(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, type='response')

  output.daic$yhat <- output.daic$yhat + pred.val.data$fit * w.daic[bb]
  #output.sd.daic$yhat <-  output.sd.daic$yhat + as.vector(pred.val.data$se.fit) * weights[bb]
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
    #output.sd <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))
  }else{
    output <- rbind(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))))
    #output.sd <- rbind(output.sd, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))))
  }
}
colnames(output) <- c("pixel_crownID", "modelID", "yhat")
#colnames(output.sd) <- c("pixel_crownID", "modelID", "yhat_unc")
output.daic <- as.data.frame(as.matrix(output.daic))
colnames(output.daic) <- c("pixel_crownID", "yhat")
#colnames(output.sd.daic) <- c("pixel_crownID", "yhat_unc")


output <- output[!is.infinite(output$yhat), ]
output.sd <- output.sd[!is.infinite(output.sd$yhat), ]
output <- output[!is.nan(output$yhat), ]
output.sd <- output.sd[!is.nan(output.sd$yhat), ]
output.daic <- output.daic[!is.nan(output.daic$yhat), ]
output.sd.daic <- output.sd.daic[!is.nan(output.sd.daic$yhat), ]


#output$yhat <- output$yhat * weights
pixel.mat <- inner_join(output, test.data.y, by = "pixel_crownID")
pixel.daic.mat <- inner_join(output.daic, test.data.y, by = "pixel_crownID")
colnames(pixel.daic.mat)[3] <- "y"
colnames(pixel.mat)[4] <- "y"
pixel.based.r2 <-  1 - sum((pixel.mat$yhat - (pixel.mat$y))^2) / sum((pixel.mat$y - mean(pixel.mat$y))^2)
pixel.daic.r2 <-  1 - sum((pixel.daic.mat$yhat - (pixel.daic.mat$y))^2) / sum((pixel.daic.mat$y - mean(pixel.daic.mat$y))^2)

#pixel.based <- cor(pixel.mat$yhat, pixel.mat$y)^2
crown.based <- pixel.mat %>%
  group_by(pixel_crownID) %>%
  summarise(yhat = mean(yhat))
crown.based <- inner_join(crown.based, test.data.y, by = "pixel_crownID")
colnames(crown.based) <- c("pixel_crownID", "yhat", "y")
crown.based.r2 <-  1 - sum((crown.based$yhat - (crown.based$y))^2) / sum((crown.based$y - mean(crown.based$y))^2)
#crown.based.r2 <- cor(crown.based$yhat, crown.based$y)^2

crown.based.daic <- pixel.daic.mat %>%
  group_by(pixel_crownID) %>%
  summarise(yhat = mean(yhat))
crown.based.daic <- inner_join(crown.based.daic, test.data.y, by = "pixel_crownID")
colnames(crown.based.daic) <- c("pixel_crownID", "yhat", "y")
crown.daic.r2 <-  1 - sum((crown.based.daic$yhat - (crown.based.daic$y))^2) / sum((crown.based.daic$y - mean(crown.based.daic$y))^2)


cor(crown.based.daic$yhat,crown.based.daic$y)^2

crown.daic.r2
crown.based.r2
pixel.based.r2
pixel.daic.r2

print("sd")
mean(output.sd.daic$yhat_unc, na.rm = T)


out.sd <- output.sd.daic %>%
  group_by(pixel_crownID) %>%
  summarise(yhat = mean(yhat_unc, na.rm = T))

mean(out.sd$yhat)

print("rmse")
rmse(pixel.daic.mat$yhat, pixel.daic.mat$y)
rmse(crown.based.daic$yhat, crown.based.daic$y)



mlim = c(min(min(pixel.daic.mat$yhat), min(pixel.daic.mat$y)))
Mlim = c(max(min(pixel.daic.mat$yhat), max(pixel.daic.mat$y)))
plot(pixel.daic.mat$yhat, pixel.daic.mat$y, xlim = c(mlim,Mlim),ylim = c(mlim,Mlim))
par(new=TRUE)
plot(crown.based.daic$yhat, crown.based.daic$y,xlim = c(mlim,Mlim),ylim = c(mlim,Mlim), col="red")
abline(0, 1) 

write.csv(crown.based.daic, paste(out.dir, nm, "_crown_res.csv", sep=""))
write.csv(pixel.daic.mat, paste(out.dir, nm, "_pix_res.csv", sep=""))
