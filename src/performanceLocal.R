rm(list=ls(all=TRUE))   # clear workspace

library("tsensembler")
library(readr)
library(plsRglm)
library(tidyverse)
NeonSite = "ALL"
wd = "/Users/sergiomarconi/Documents/Projects/TraitsOnHeaven/"
mddir <- "/Users/sergiomarconi/Dropbox (UFL)/Marconi2018/readableModels/ALL/"
setwd(wd)
source(paste(getwd(), '/src/functions_build.R', sep=""))
source(paste(getwd(), '/src/src_build.R', sep=""))
source(paste(getwd(), '/src/functions_spatial.R', sep=""))
source(paste(getwd(), '/src/src_spatial.R', sep=""))
#source(paste(getwd(), '/src/pls_par.R', sep=""))


path = paste(wd, NeonSite, sep="/")
setwd(path)
out.dir = paste(getwd(), "/outputs/", sep="")
in.dir = paste(getwd(), "/inputs/", sep="")

nm = "N"
j="N_pct"
get100=F
method = 'aic'
normz =F
LOG = F
cutoff = 10
loops = 100

if(get100){
  
  if(normz){
    load(file = paste(mddir, 'pls_glm',nm, sep="" ))
  }else{
    load(file = paste(mddir,'pls_glm_',j, sep="" ))
  }
  mod.aic <- rep(NA, loops)
  mod.r2 <- rep(NA, loops)
  weights = rep(0,loops)
  for(bb in 1: length(mod.aic)){
    mod.aic[bb] <-foo[[bb]]$score$aic
    mod.r2[bb] <- 1 - sum((foo[[bb]]$pred$fit - (foo[[bb]]$score$data$Y.test))^2) /
      sum((foo[[bb]]$score$data$Y.test - mean(foo[[bb]]$score$data$Y.test))^2)
  }
  mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), cutoff)), mod.r2[mod.r2 %in% tail(sort(mod.r2), cutoff)])[1:cutoff, ]
  #mask <- cbind(which(mod.aic %in% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])
  #mask <- mask[order(mask[,2]),]
  plsglm <- foo[mask[,1]]
}else{
  load(file=paste(mddir, "./plsglm_", j,  sep = ""))
  mod.aic <- rep(NA, loops)
  mod.r2 <- rep(NA, loops)
  weights = rep(0,loops)
  for(bb in 1: length(mod.aic)){
    mod.aic[bb] <-plsglm[[bb]]$score$aic
    mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
      sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
  }
  mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), cutoff)), mod.r2[mod.r2 %in% tail(sort(mod.r2), cutoff)])[1:cutoff, ]
  plsglm <- plsglm[mask[,1]]
}
# mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), cutoff)), mod.r2[mod.r2 %in% tail(sort(mod.r2), cutoff)])[1:cutoff, ]
# plsglm <- plsglm[mask[,1]]

mod.r2=rep(0,length(plsglm))
mod.aic=rep(0,length(plsglm))

for(bb in 1: length(plsglm)){
  mod.aic[bb] <-plsglm[[bb]]$score$aic #plsglm[[bb]]$aic
  mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
    sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
}
mod.aic
mod.r2

if(method=='aic'){
  delta.aic <- mod.aic - min(mod.aic)
  weights <- softmax(-0.5*delta.aic)
}else if(method =='press'){
  delta.aic <- abs(mod.r2 - min(mod.r2))
  weights <- delta.aic/sum(delta.aic)
  
}
#
#weights <- softmax(-0.5*delta.aic)
#weights[mask[,1]] <- softmax(0.5*delta.press)
weights
#weights <- softmax(0.5*delta.press)

#out of bag
test.data.y <- read.csv(paste(in.dir, "Spectra/CrownTraits_outBag.csv", sep=""))
test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
#test.data.y <- test.data.y[test.data.y$pixel_crownID >2000,]
if(LOG){test.data.y[names(test.data.y)==j] <- log(test.data.y[names(test.data.y)==j])}
test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag.csv", sep=""))
crownID <- as.data.frame(test.data.x$pixel_crownID)
siteBands <- test.data.x[,2:3]
test.data.x <- test.data.x[grepl("band", names(test.data.x))]
spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
if(normz){  
  spectra_log_dif_snv[spectra_log_dif_snv==0] <- 0.0000001
  spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv[,-c(1:2)]))),differences=1, lag=3))
  spectra_log_dif_snv <- cbind(siteBands, spectra_log_dif_snv)
}    
#test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
test.PLS = (as.matrix(spectra_log_dif_snv))

rm(output)
rm(output.sd)
output.daic =  output.sd.daic <- crownID
out.validation = length(plsglm[[1]]$pred$fit)
output.daic$yhat = output.sd.daic$yhat = out.validation$yhat = 0

w.daic <- rep(0, length(weights))
w.daic <- weights
for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  optim.ncomps <- plsglm[[bb]]$ncomp
  pred.val.data <- predict(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, type='response',  se.fit = T)
  output.daic$yhat <- output.daic$yhat + pred.val.data$fit * w.daic[bb]
  output.sd.daic$yhat <-  output.sd.daic$yhat + as.vector(pred.val.data$se.fit) * weights[bb]
  out.validation$yhat <- out.validation$yhat + plsglm[[bb]]$pred$fit * weights[bb]
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
    output.sd <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))
  }else{
    output <- rbind(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))))
    output.sd <- rbind(output.sd, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))))
  }
}
valid.ensamble <- 1 - sum((out.validation$yhat - (plsglm[[1]]$score$data$Y.test))^2) /
  sum((plsglm[[1]]$score$data$Y.test - mean(plsglm[[1]]$score$data$Y.test))^2)
colnames(output) <- c("pixel_crownID", "modelID", "yhat")
colnames(output.sd) <- c("pixel_crownID", "modelID", "yhat_unc")
output.daic <- as.data.frame(as.matrix(output.daic))
colnames(output.daic) <- c("pixel_crownID", "yhat")
colnames(output.sd.daic) <- c("pixel_crownID", "yhat_unc")


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
pix.r2.vect <- pixel.mat %>%
  group_by(modelID) %>%
  summarise(r2 = 1 - sum((yhat - (y))^2) / sum((y - mean(y))^2))
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

pix.r2.vect$bpix <- pixel.daic.r2
pix.r2.vect$bcrown <-crown.daic.r2

cor(crown.based.daic$yhat,crown.based.daic$y)^2

crown.daic.r2
#crown.based.r2
#pixel.based.r2
pixel.daic.r2
valid.ensamble
summary(mod.r2)
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
write.csv(pix.r2.vect, paste(out.dir, nm, "_r2_vec.csv", sep=""))



# ------------------
#
#   need more organize
#-------------------
# get uncertainty for single bands (for each model)
out <- list()
rm(unc100, avePar, signific100)
par.boot.out <- list()
for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  out["ncomp"] <- plsglm[[bb]]$ncomp
  set.seed(123)
  parUnc=bootplsglm(pls.mod.train,R=1000)
  #boxplots.bootpls(parUnc,las=2,mar=c(5,2,1,1)+0.1)
  par.ci=confints.bootpls(parUnc)
  par.boot.out[[bb]] <- parUnc
  significance <- (par.ci[,7]<0&par.ci[,8]<0)|(par.ci[,7]>0&par.ci[,8]>0)
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("avePar")){
    unc100 <- cbind(as.vector(rownames(par.ci)), par.ci[,7:8])
    signific100 <- (significance)
    avePar <- cbind(1:371, apply(parUnc$t, 2, median) )
  }else{
    unc100 <- cbind(unc100, par.ci[,7:8])
    signific100 <- cbind(signific100, significance)
    avePar <- cbind(avePar, apply(parUnc$t, 2, median))
    
  }
}
matind=(signific100)
colnames(matind) <- paste("md", 1:100, sep="_")
as.matrix(matind)

library(ggplot2)
library(reshape2)

melted <- melt(t(matind))
ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  scale_fill_manual(values = c("white", "black"))

pr.band <- as.data.frame(apply(matind,1, table)) /100
pr.band <- pr.band[grepl("Freq", colnames(pr.band))]
pr.band
plot(1:368, pr.band[2,], type="l")
dd <- as.data.frame(avePar)
d <- melt(dd, "V1")
p <-  ggplot(d, aes(x = (V1), y = value)) 
p   + theme_bw() + geom_point()

# ------------------
#
#   need more organize
#-------------------


# #test set reflectance
# refl.oob <- cbind(crownID, test.data.x[,-c(1,2)])
# colnames(refl.oob)[1]<-"crID"
# d <- refl.oob %>% tidyr::gather("Object","value", names(refl.oob)[-1])
# d$value[d$value<0.001]=NA
# p <-  ggplot(d, aes(x = Object, y = value)) 
# p  + facet_grid(crID~.)  + theme_bw() + geom_line()+
#   stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1)
# 
# d2 <- refl.oob %>%
#   group_by(crID) %>%
#   summarise_all(.funs = c(sd, mean))
# names(d2)[2:427] <- paste("b",1:426,sep="_")
# names(d2)[428:853]<- paste("sd",1:426,sep="_")
# d <- d2 %>% tidyr::gather("Object","value", names(d2)[-1])
# 
# d <-  melt(refl.oob, id.vars =names(refl.oob)[-1])
# d <- refl.oob %>% tidyr::gather("Object","value", names(refl.oob)[-1])

pix.sd <- inner_join(output.sd, test.data.y, by = "pixel_crownID")
colnames(pix.sd)[3] <- "y_sd"

cr.sd <- pix.sd %>%
  group_by(pixel_crownID, modelID) %>%
  summarise(yhat = mean(y_sd))
cr.sd <- inner_join(cr.sd, test.data.y, by = "pixel_crownID")
colnames(cr.sd) <- c("pixel_crownID", "modelID","y_sd", "y")

library(reshape2)

sd.compare <- inner_join(cr.sd, pix.sd, by = "pixel_crownID")
sd.compare <- sd.compare[names(sd.compare) %in% 
                           c("pixel_crownID", "modelID.x","y_sd.x", "y_sd.y", j)]
names(sd.compare) <- c("pixel_crownID", "modelID","Crown", "Pixel", "y")
d <- sd.compare %>% tidyr::gather("Object","value", c("Crown", "Pixel"))
d <- unique(d)

d$perc <- d$value#/d$y 
p <-  ggplot(d, aes(x = (modelID), y = perc)) 
p  + facet_grid(. ~ pixel_crownID)  + theme_bw()+
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+ 
  geom_boxplot(aes(fill=Object),notch = T, alpha = 0.4, width = 0.04) + 
  ylab(bquote(C*'( '*se['md']*')')) + theme(text = element_text(size=18)) + xlab("Crown ID")

ggsave(paste(out.dir, nm, "_sd_compare.png", sep=""), width = 11.8, height = 4.93)


p <-  ggplot(cr.sd, aes(x = (y), y = y_sd)) 
names(cr.sd)[4] <- "y"
cr.sd$perc <- cr.sd$y_sd #/ cr.sd$y 
p <-  ggplot(cr.sd, aes(x = (y), y = perc)) 
p +  theme_bw()+geom_smooth(alpha = 0.2, method = "lm") + geom_boxplot(aes(x=y, y=perc,group=y), width=0.008,notch=T, fill="grey", alpha=0.2) + 
  ylab(bquote(P*'( '*se['md']*')'))+ xlab(bquote(P['ob']*'(%)'))+ theme(text = element_text(size=20))

#xlab(bquote(C*'(g '*m^2*')'))

ggsave(paste(out.dir, nm, "_sd_per_unit.png", sep=""), width = 4.93, height = 4.93,  units="in")

prova <- output %>%
  group_by(modelID, pixel_crownID) %>%
  summarise(y_md = mean(yhat))
prova <- inner_join(prova, test.data.y, by = "pixel_crownID")
names(prova)[4] <- "y"
prova$dY <-  prova$y_md - prova$y
#prova$dY <- prova$y - prova$y_md
prova$dY_perc <- prova$dY # / prova$y

p <-  ggplot(prova, aes(x = (y), y = dY_perc)) 
p +  geom_hline(yintercept = 0)  + geom_smooth(alpha = 0.2, method = "lm", formula = y ~ splines::bs(x, 2))+  theme_bw() + 
  geom_boxplot(aes(x=y, y=dY_perc,group=y), notch=T, fill="grey", width=0.008,  alpha=0.2) +
  ylab(bquote(P*'( '*res['(ob-md)']*')'))+ xlab(bquote(P*'(%)')) + theme(text = element_text(size=20))

ggsave(paste(out.dir, nm, "_residuals_per_unit.png", sep=""), units="in", width = 4.93, height = 4.93)

#residuals compare crown against pixel
res.compare <- inner_join(prova, pixel.mat, by = "pixel_crownID")
res.compare <- res.compare[names(res.compare) %in%  
                             c("pixel_crownID",  "modelID.x","y_md","y.y", "dY_perc", "yhat")]
names(res.compare) <- c("modelID","pixel_crownID", "Crown","dCrown", "Pixel", "y")
res.compare$dPixel <- (res.compare$Pixel - res.compare$y)#/res.compare$y
res.compare <- unique(res.compare)
d <- res.compare %>% 
  tidyr::gather("Object","value", c("dCrown", "dPixel"))
d <- unique(d)
p <-  ggplot(d, aes(x = (modelID), y = value)) 
p  + facet_grid(. ~ pixel_crownID)  + theme_bw()+
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+ 
  geom_boxplot(aes(fill=Object),notch = T, alpha = 0.4, width = 0.08) + 
  ylab(bquote(LMA*'( '*res['(ob-md)']*')')) + theme(text = element_text(size=18)) + xlab("Crown ID") + geom_hline(yintercept = 0)

ggsave(paste(out.dir, nm, "_residuals_compare.png", sep=""), units="in", width = 11.8, height = 4.93)




# #load accuracy
# d <- accuracy %>% tidyr::gather("Object","value", c("perm", "bag_pix", "bag_crown"))
# d <- unique(d)
# p <-  ggplot(d, aes(x = (modelID), y = value)) 
# p  +  geom_boxplot(data=subset(d,Object=="perm"), notch = T, alpha = 0.4) +  
#   geom_jitter(data=subset(d,Object=="perm"), alpha = 0.2, width = 0.2)+
#   geom_point(data=subset(d,Object=="bag_pix"), col= "blue",shape =20, size = 2) + 
#   geom_point(data=subset(d,Object=="bag_crown"), col= "red", shape =24, size = 2 ) +
#   theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+ theme_bw()+
#   ylab(bquote(R['OOB']^2)) + theme(text = element_text(size=18)) + xlab("")
# 
