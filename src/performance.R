rm(list=ls(all=TRUE))   # clear workspace


NeonSite = "ALL"
wd = "/ufrc/ewhite/s.marconi/Marconi2018/ModelBuild/"
setwd(wd)
source('../src/functions_build.R')
source('../src/src_build.R')
source( '../src/functions_infer.R')
source('../src/src_infer.R')
#source(paste(getwd(), '/src/pls_par.R', sep=""))

library(tidyverse)
library(ggplot2)
library(fifer)
library(reshape2)
library(stargazer)
library(tsensembler)



path = paste(wd, NeonSite, sep="/")
setwd(path)
out.dir = paste(getwd(), "/outputs/", sep="")
in.dir = paste(getwd(), "/inputs/", sep="")

nm = "P"
j="P_pct"
get100=F
method = 'press'
normz = T
LOG=F
glm = T
if(get100){
  loops = 1000
  if(glm){
    load(file = paste('./pls_glm_',j, sep="" ))
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
  #mask <- cbind(which(mod.aic %iwen% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])
  #mask <- mask[order(mask[,2]),]
  plsglm <- foo[mask[,1]]
  save(plsglm, file = paste("plsglm_", j,  sep = ""))
}else{
  load(file=paste("./plsglm_", j,  sep = ""))

}
mod.r2=rep(0,length(plsglm))
mod.aic=rep(0,length(plsglm))

for(bb in 1: length(plsglm)){
  #mod.aic[bb] <-plsglm[[bb]]$score$aic #plsglm[[bb]]$aic
  mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred$fit - (plsglm[[bb]]$score$data$Y.test))^2) /
    sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
}
mod.aic
summary(mod.r2)

if(method=='aic'){
  delta.aic <- mod.aic - min(mod.aic)
  weights <- softmax(-0.5*delta.aic)
}else if(method =='press'){
  delta.aic <- (max(mod.r2)-mod.r2)
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
test.data.y <- test.data.y[-5,]

if(LOG){test.data.y[names(test.data.y)==j] <- log(test.data.y[names(test.data.y)==j])}
test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag_new.csv", sep=""))
x.id <- test.data.x[,1]
ndvi <- (test.data.x$band_90- test.data.x$band_58)/(test.data.x$band_58 + test.data.x$band_90)
nir860 <- (test.data.x$band_96 + test.data.x$band_97)/2
test.data.x[which(ndvi < 0.7 | nir860 < .3),]=NA
#allData[which(nir860 < .3),] = NA
test.data.x <- test.data.x[complete.cases(test.data.x), ]
pName<- test.data.x$pixel_crownID

allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
test.data.x[,bad_Bands]=NA

test.data.x <- test.data.x[grepl("band", names(test.data.x))]

test.data.x[test.data.x>1]=NA
#test.data.x$pixel_crownID <- pName
specMat=as.matrix(test.data.x[,-c(1:2)])
# Vector normalize spectra
# equalizer=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
# norm_par <- read.csv("./outputs/Equal_norm_vectors.csv")
# mu_train = norm_par[,1]
# sd_train =  norm_par[,2]
# normMat=matrix(data=rep(equalizer,ncol(specMat)),ncol=ncol(specMat))
# for(i in 1:length(sd_train)){
#   normMat[,i]= (specMat[,i] - mu_train[i])/(sd_train[i])
# }
normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
normMat=specMat/normMat
test.data.x[,-(1:2)] <- normMat

crownID <- as.data.frame(pName)
spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
if(normz){
  spectra_log_dif_snv[spectra_log_dif_snv==0] <- 0.0000001
  if(NeonSite =="ALL"){
    spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv[,-(1:2)]))),differences=1, lag=3))
    spectra_log_dif_snv<- cbind(test.data.x[,1:2], spectra_log_dif_snv)
  }else{
    spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv))),differences=1, lag=3))
  }
}
#test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
test.PLS = (as.matrix(spectra_log_dif_snv))

rm(output)
rm(output.sd)
output.daic =  output.up.daic = output.lw.daic <- crownID
output.daic$yhat = output.up.daic$yhat = output.lw.daic$yhat = 0
w.daic <- rep(0, length(weights))
w.daic <- weights
out <- list()
pred.val.data <- list()
r2.out <- rep(NA, length(w.daic))

for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  optim.ncomps <- plsglm[[bb]]$ncomp
  if(glm){
    foo_pi <- predict.withsd(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, wt = pls.mod.train$FinalModel$weights,  type='response')
    pred.val.data$fit <- foo_pi[,1]
    pred.val.data$upper <- foo_pi[,3]
    pred.val.data$lower <- foo_pi[,2]

    out$output.daic <-  pred.val.data$fit
    out$upper <-  pred.val.data$upper
    out$lower <-  pred.val.data$lower

    #output.sd.daic$yhat <-  output.sd.daic$yhat + as.vector(pred.val.data$se.fit) * weights[bb]
    ppercr <- as.data.frame(cbind(pName, pred.val.data$fit))
    colnames(ppercr) <- c("pixel_crownID", "yhat")
    appu <- stratified(ppercr, c("pixel_crownID"), 1)
    #appu <- inner_join(ppercr, test.data.y, by="pixel_crownID")
    r2.out[bb] <- 1 - sum((appu$yhat - test.data.y[,2])^2) /
       sum((test.data.y[,2] - mean(test.data.y[,2]))^2)
   #r2.out[bb] <- 1 - sum((appu$yhat - appu[,3])^2) /
    #  sum((appu[,3] - mean(appu[,3]))^2)
  }else{
    pred.val.data$fit <- predict(pls.mod.train, newdata = test.PLS, ncomp=optim.ncomps, type='response', se.fit = T)
  }
  output.daic$yhat <- output.daic$yhat + pred.val.data$fit * w.daic[bb]
  output.up.daic$yhat <- output.up.daic$yhat + pred.val.data$upper * w.daic[bb]
  output.lw.daic$yhat <- output.lw.daic$yhat + pred.val.data$lower * w.daic[bb]

  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
    if(glm){
      #output.sd <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))
    }
  }else{
    output <- rbind(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))))
    if(glm){
      #output.sd <- rbind(output.sd, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$se.fit))))
    }
  }
}
colnames(output) <- c("pixel_crownID", "modelID", "yhat")
output.daic <- as.data.frame(as.matrix(output.daic))
colnames(output.daic) <- c("pixel_crownID", "yhat")
if(glm){
  colnames(output.up.daic) <- c("pixel_crownID", "yhat_up")
  colnames(output.lw.daic) <- c("pixel_crownID", "yhat_lw")

#  output.sd <- output.sd[!is.infinite(output.sd$yhat), ]
#  output.sd <- output.sd[!is.nan(output.sd$yhat), ]
  output.up.daic <- output.up.daic[!is.nan(output.daic$yhat), ]
  output.lw.daic <- output.lw.daic[!is.nan(output.daic$yhat), ]

}

output <- output[!is.infinite(output$yhat), ]
output <- output[!is.nan(output$yhat), ]
output.daic <- output.daic[!is.nan(output.daic$yhat), ]


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
write.csv(crown.based.daic, paste(out.dir, nm, "_crown_res.csv", sep=""))
write.csv(pixel.daic.mat, paste(out.dir, nm, "_pix_res.csv", sep=""))

pix.up <- inner_join(output.up.daic, test.data.y, by = "pixel_crownID")
colnames(pix.up)[2] <- "y_up"
cr.up <- pix.up %>%
  group_by(pixel_crownID) %>%
  summarise(y95pi = mean(y_up))
cr.up <- inner_join(cr.up, test.data.y, by = "pixel_crownID")
colnames(cr.up) <- c("pixel_crownID", "y_up", "y")

pix.lw <- inner_join(output.lw.daic, test.data.y, by = "pixel_crownID")
colnames(pix.lw)[2] <- "y_lw"
cr.lw <- pix.lw %>%
  group_by(pixel_crownID) %>%
  summarise(y5pi = mean(y_lw))
cr.lw <- inner_join(cr.lw, test.data.y, by = "pixel_crownID")
colnames(cr.lw) <- c("pixel_crownID", "y_lw", "y")



sd_compare_crown <- cbind(inner_join( cr.up, cr.lw, by = "pixel_crownID"), "Crown")[,-3]
sd_compare_crown <- inner_join(crown.based.daic, sd_compare_crown, by = "pixel_crownID")[,-6]

sd_compare_pix <- cbind(inner_join(pix.up, pix.lw, by = "pixel_crownID"), "Pix")[,-3]
sd_compare_pix <- inner_join(pixel.daic.mat,sd_compare_pix, by = "pixel_crownID")[,-6]
colnames(sd_compare_crown) = c("pixel_crownID","yhat", "y", "PI95", "PI5",  "sim")
colnames(sd_compare_pix) = colnames(sd_compare_crown)

#sd.compare <- inner_join(sd_compare_crown, sd_compare_pix, by = "pixel_crownID")
sd.compare <- rbind(sd_compare_crown, sd_compare_pix)
#names(sd.compare) <- c("pixel_crownID","yhat",PI95c", "PI95p","PI5c", "y","PI5p")
#d <- sd.compare %>% tidyr::gather("Object","value", c("PI95c", "PI95p","PI5c", "PI5p", "y"))
df <- unique(sd.compare)

#d$perc <- d$value/d$y
# p  +  theme_bw()+
#   theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#   geom_line(aes(fill=Object),notch = T, alpha = 0.4, width = 0.04) +
#   ylab(bquote('( '*sigma['(ob-md)']/ mu['ob']*')')) + theme(text = element_text(size=18)) + xlab("Crown ID")

p <- ggplot(df, aes(factor(pixel_crownID), y, colour = sim))
p + theme_bw()+geom_ribbon(aes(ymin = PI5, ymax = PI95), position = position_dodge(0.5))+
  geom_point(aes(y = yhat), position = position_dodge(0.5)) + geom_point(aes(y = y), colour = "black")+
theme(legend.title=element_blank())

ggsave(paste(out.dir, nm, "_sd_compare.png", sep=""), width = 11.8, height = 4.93)


# p <-  ggplot(cr.sd, aes(x = (y), y = y_sd))
# names(cr.sd)[4] <- "y"
# cr.sd$perc <- cr.sd$y_sd / cr.sd$y
# p <-  ggplot(cr.sd, aes(x = (y), y = perc))
# p +  theme_bw()+geom_smooth(alpha = 0.2) + geom_boxplot(aes(x=y, y=perc,group=y), width=0.008,notch=T, fill="grey", alpha=0.2) +
#   ylab(bquote('( '*sigma['(ob-md)']/ mu['ob']*')'))+ xlab(nm) + theme(text = element_text(size=20))

#xlab(bquote(C*'(g '*m^2*')'))

# ggsave(paste(out.dir, nm, "_sd_per_unit.png", sep=""), width = 4.93, height = 4.93,  units="in")

prova <- output %>%
  group_by(modelID, pixel_crownID) %>%
  summarise(y_md = mean(yhat))
prova <- inner_join(prova, test.data.y, by = "pixel_crownID")
names(prova)[4] <- "y"
prova$dY <-  prova$y_md - prova$y
#prova$dY <- prova$y - prova$y_md
prova$dY_perc <- prova$dY / prova$y

p <-  ggplot(prova, aes(x = (y), y = dY_perc))
p +  geom_hline(yintercept = 0)  +   theme_bw() + geom_boxplot(aes(x=y, y=dY_perc,group=y), notch=T, fill="grey", width=0.008,  alpha=0.2) +
  ylab(bquote('( '*delta['(ob-md)']/ mu['ob']*')'))+ xlab(nm) + theme(text = element_text(size=20))

ggsave(paste(out.dir, nm, "_residuals_per_unit.png", sep=""), units="in", width = 4.93, height = 4.93)

#residuals compare crown against pixel
res.compare <- inner_join(prova, pixel.mat, by = "pixel_crownID")
res.compare <- res.compare[names(res.compare) %in%
                             c("pixel_crownID",  "modelID.x","y_md","y.y", "dY_perc", "yhat")]
names(res.compare) <- c("modelID","pixel_crownID", "Crown","dCrown", "Pixel", "y")
res.compare$dPixel <- (res.compare$Pixel - res.compare$y)/res.compare$y
res.compare <- unique(res.compare)
d <- res.compare %>%
  tidyr::gather("Object","value", c("dCrown", "dPixel"))
d <- unique(d)
p <-  ggplot(d, aes(x = (modelID), y = value))
p  + facet_grid(. ~ pixel_crownID)  + theme_bw()+
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  geom_boxplot(aes(fill=Object),notch = T, alpha = 0.4, width = 0.04) +
  ylab(bquote('( '*delta['(ob-md)']/ mu['ob']*')')) + theme(text = element_text(size=18)) + xlab("Crown ID") + geom_hline(yintercept = 0)

ggsave(paste(out.dir, nm, "_residuals_compare.png", sep=""), units="in", width = 11.8, height = 4.93)


# crown.daic.r2
# crown.based.r2
# pixel.based.r2
# pixel.daic.r2

#print("sd")
#mean(output.sd.daic$yhat_unc, na.rm = T)
out.sd <- output.sd.daic %>%
  group_by(pixel_crownID) %>%
  summarise(yhat = mean(yhat_unc, na.rm = T))

#mean(out.sd$yhat)

#print("rmse")
# rmse(pixel.daic.mat$yhat, pixel.daic.mat$y)
# rmse(crown.based.daic$yhat, crown.based.daic$y)

tab_stats <- data.frame(crown.daic.r2,pixel.daic.r2,
                        mean(output.sd.daic$yhat_unc, na.rm = T),mean(out.sd$yhat),
                        rmse(pixel.daic.mat$yhat, pixel.daic.mat$y), rmse(crown.based.daic$yhat, crown.based.daic$y))
colnames(tab_stats) <- c("R2 Crown", "R2 pixel", "sd Crown", "sd pixel","rmse Crown", "rmse pixel")
#your tables
#titles
pdf(paste(out.dir, nm, "_modeltype_compare.pdf", sep=""))
plot(density(r2.out), xlim = c(-2,1))
polygon(density(r2.out), col="green", border="black")
#hist(r2.out, xlim = c(-10,1), breaks = 20, col='grey')
abline(v = crown.daic.r2, col = 'blue', lwd=3,cex=1.5, lty=2)
abline(v = pixel.daic.r2, col = 'magenta',lwd=3,cex=1.5, lty=2)
abline(v = 0, col = 'black',lwd=3,cex=1.5)
dev.off()

out <- list()
rm(avePar)
par.boot.out <- list()
for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  out["ncomp"] <- plsglm[[bb]]$ncomp
  set.seed(123)
  parUnc=bootplsglm(pls.mod.train,R=1000)
   par.boot.out[[bb]] <- parUnc
  if(!exists("avePar")){
     avePar <- cbind(1:371, apply(parUnc$t, 2, median) )
  }else{
    avePar <- cbind(avePar, apply(parUnc$t, 2, median))

  }
}
dd <- as.data.frame(avePar)
write.csv(summary(t(dd)), paste(nm,"paramUnc.csv", sep="_"))


stargazer(tab_stats[1,], summary=FALSE, rownames=FALSE)
