rm(list=ls(all=TRUE))   # clear workspace


library(tidyverse)
library(ggplot2)
library(plsRglm)
library(reshape2)
library(stargazer)
library(tsensembler)
library(fifer)
#/ufrc/ewhite/s.marconi/Marconi2018/ModelBuild/ALL/
wd = "/ufrc/ewhite/s.marconi/Chapter1/Model_builder/"
source(paste(wd, './src/src_build.R', sep=""))
out_dir = paste(wd, "/outputs/", sep="")
in_dir = paste(wd, "/inputs/", sep="")

nm = j="P_pct"
get100=F
method = 'aic'
normz = T
exp_itc = ""
pls_dir = "//orange/ewhite/s.marconi/traitsOnHeaven/original_pls_glm/"
if(get100){
  loops = 1000
  load(file = paste(out_dir, '/pls_glm_',j, sep="" ))

  mod.aic <- rep(NA, loops)
  mod.r2 <- rep(NA, loops)
  weights = rep(0,loops)
  for(bb in 1: length(mod.aic)){
    mod.aic[bb] <- foo[[bb]]$score$aic
    mod.r2[bb] <- 1 - sum((foo[[bb]]$pred - (foo[[bb]]$score$data$Y.test))^2) /
      sum((foo[[bb]]$score$data$Y.test - mean(foo[[bb]]$score$data$Y.test))^2)
  }
  # mask <- cbind(which(mod.r2 >0), mod.r2[mod.r2 >0])
  # foo <- foo[mask[,1]]
  mask <- cbind(which(mod.aic %in% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])
  plsglm <- foo[mask[,1]]
  save(plsglm, file = paste(pls_dir, "plsglm_daic", j,  exp_itc, sep = ""))
}else{
  #load(file=paste("//orange/ewhite/s.marconi/traitsOnHeaven/original_pls_glm/plsglm_daic", j,  sep = ""))
  load(file=paste(pls_dir, "/plsglm_", j, exp_itc,  sep = ""))
}
mod.r2=rep(0,length(plsglm))
mod.aic=rep(0,length(plsglm))

for(bb in 1: length(plsglm)){
  mod.aic[bb] <- plsglm[[bb]]$mod$FinalModel$aic
  #mod.r2[bb] <- 1 - sum((plsglm[[bb]]$pred - unlist(plsglm[[bb]]$score$data$Y.test))^2) /
  #  sum((plsglm[[bb]]$score$data$Y.test - mean(plsglm[[bb]]$score$data$Y.test))^2)
}
mod.aic <- scale(mod.aic)
summary(mod.aic)
summary(mod.r2)

if(method=='aic'){
  delta.aic <- mod.aic[,1] - min(mod.aic)
  weights <- softmax(-0.5*delta.aic)
}else if(method =='press'){
  delta.aic <- (max(mod.r2)-mod.r2)
  weights <- softmax(-0.5*delta.aic)

}
#to calculate residuals on training/validation set
train.data.y <- read.csv(paste(in_dir, "Spectra/CrownTraits.csv", sep=""))
train.data.y <-train.data.y[colnames(train.data.y) %in% c("pixel_crownID",j)]
#out of bag
test.data.y <- read.csv(paste(in_dir, "Spectra/CrownTraits_outBag.csv", sep=""))
nsites <- length(unique(test.data.y$site))

test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
test.data.y <- test.data.y[-5,]

#test.data.x <- read.csv("///ufrc/ewhite/s.marconi/Chapter1/Model_builder/inputs/CrownPix_outBag.csv")
test.data.x <- read.csv(paste(in_dir, "Spectra/CrownPix_outBag.csv", sep=""))
head(test.data.x)
test.data.x <- test.data.x[-(which(colnames(test.data.x) %in% c("band_1", "band_370")))]
#test.data.x[test.data.x==0] = 0.0000001

#x.id <- test.data.x[,1]
# ndvi <- (test.data.x$band_90- test.data.x$band_58)/(test.data.x$band_58 + test.data.x$band_90)
# nir860 <- (test.data.x$band_96 + test.data.x$band_97)/2
# test.data.x[which(ndvi < 0.7 | nir860 < .3),]=NA

test.data.x <- test.data.x[complete.cases(test.data.x), ]
pName<- test.data.x$pixel_crownID
test.data.x <- test.data.x[grepl("band", names(test.data.x))]

#allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
#bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
#test.data.x[,bad_Bands]=NA


# test.data.x[test.data.x>1]=NA
# specMat=as.matrix(test.data.x[,-c(1:2)])

# normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
# normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
# normMat=specMat/normMat
# test.data.x[,-(1:2)] <- normMat

crownID <- as.data.frame(pName)
spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
if(normz==T){
  foot <-t(diff(t(log(spectra_log_dif_snv[,-c(1:nsites)])),differences=1, lag=3))
  spectra_log_dif_snv <- cbind(spectra_log_dif_snv[,c(1:nsites)], foot)
}

test.PLS = (as.matrix(spectra_log_dif_snv))

rm(output)
rm(output.sd)
output.daic =  output.up.daic = output.lw.daic =  crownID
output.residuals <- train.data.y
output.daic$yhat = output.up.daic$yhat =  output.residuals$yhat = output.lw.daic$yhat = 0
w.daic <- rep(0, length(weights))
w.daic <- weights
out <- list()
pred.val.data <- list()
r2.out = r2.appu <- rep(NA, length(w.daic))

for(bb in 1:length(w.daic)){
  pls.mod.train <- plsglm[[bb]]$mod
  optim.ncomps <- plsglm[[bb]]$ncomp
  #if(glm){
  foo_pi <- predict.withsd(pls.mod.train, newdata = test.PLS,
                           ncomp=optim.ncomps, wt = pls.mod.train$FinalModel$weights,  type='response')
  pred.val.data$fit <- foo_pi[,1]
  pred.val.data$upper <- foo_pi[,3]
  pred.val.data$lower <- foo_pi[,2]
  out$output.daic <-  pred.val.data$fit
  out$upper <-  pred.val.data$upper
  out$lower <-  pred.val.data$lower

  ppercr <- as.data.frame(cbind(pName, pred.val.data$fit))
  colnames(ppercr) <- c("pixel_crownID", "yhat")
  #appu <- stratified(ppercr, c("pixel_crownID"), 1)
  #ppercrown <- inner_join(ppercr, test.data.y, by ="pixel_crownID")
  #r2.appu[bb] <- 1 - sum((appu$yhat - test.data.y[,2])^2) /
  #  sum((test.data.y[,2] - mean(test.data.y[,2]))^2)
  #r2.out[bb] <- 1 - sum((ppercrown$yhat - ppercrown[,3])^2) /
  #  sum((ppercrown[,3] - mean(ppercrown[,3]))^2)

  output.daic$yhat <- output.daic$yhat + pred.val.data$fit * w.daic[bb]
  output.up.daic$yhat <- output.up.daic$yhat + pred.val.data$upper * w.daic[bb]
  output.lw.daic$yhat <- output.lw.daic$yhat + pred.val.data$lower * w.daic[bb]
  output.residuals$yhat <- output.residuals$yhat + pls.mod.train$Yresidus * w.daic[bb]
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
  }else{
    output <- rbind(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))))
  }
}
colnames(output) <- c("pixel_crownID", "modelID", "yhat")
output.daic <- as.data.frame(as.matrix(output.daic))
colnames(output.daic) <- c("pixel_crownID", "yhat")

colnames(output.residuals)[1] <- "pixel_crownID"
#output.residuals <- inner_join(output.residuals, test.data.y, by ="pixel_crownID")
colnames(output.residuals) <- c("pixel_crownID", "y","yhat")

colnames(output.up.daic) <- c("pixel_crownID", "yhat_up")
colnames(output.lw.daic) <- c("pixel_crownID", "yhat_lw")

output.up.daic <- output.up.daic[!is.nan(output.daic$yhat), ]
output.lw.daic <- output.lw.daic[!is.nan(output.daic$yhat), ]

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
  summarise(yhat = median(yhat))
crown.based <- inner_join(crown.based, test.data.y, by = "pixel_crownID")
colnames(crown.based) <- c("pixel_crownID", "yhat", "y")
crown.based.r2 <-  1 - sum((crown.based$yhat - (crown.based$y))^2) / sum((crown.based$y - mean(crown.based$y))^2)
#crown.based.r2 <- cor(crown.based$yhat, crown.based$y)^2

crown.based.daic <- pixel.daic.mat %>%
  group_by(pixel_crownID) %>%
  summarise(yhat = median(yhat))
crown.based.daic <- inner_join(crown.based.daic, test.data.y, by = "pixel_crownID")
colnames(crown.based.daic) <- c("pixel_crownID", "yhat", "y")
crown.daic.r2 <-  1 - sum((crown.based.daic$yhat - (crown.based.daic$y))^2) / sum((crown.based.daic$y - mean(crown.based.daic$y))^2)
write.csv(crown.based.daic, paste(out_dir, nm, "_crown_res.csv", sep=""))
write.csv(pixel.daic.mat, paste(out_dir, nm, "_pix_res.csv", sep=""))

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
sd.compare <- rbind(sd_compare_crown, sd_compare_pix)

tags <- read_csv(paste(in_dir, "Spectra/CrownTraits_tags.csv", sep=""))
sd.compare <- inner_join(sd.compare, tags, by = "pixel_crownID")
df <- unique(sd.compare)
p <- ggplot(df, aes(factor(tag), y, colour = sim))
p + theme_bw()+geom_ribbon(aes(ymin = PI5, ymax = PI95), position = position_dodge(0.5))+ coord_flip()+
  geom_point(aes(y = yhat), position = position_dodge(0.5)) + geom_point(aes(y = y), colour = "black")+
  theme(legend.title=element_blank()) + guides(colour = guide_legend(override.aes = list(size=10)))

ggsave(paste(out_dir, nm, "_sd_compare.png", sep=""), width = 4.93, height = 11.5)

p <-  ggplot(output.residuals, aes(yhat)) + theme_bw()
p +  geom_histogram(binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)))

ggsave(paste(out_dir, nm, "_res_distribution.png", sep=""), width = 11.8, height = 4.93)


# p <-  ggplot(output.residuals, aes(x = factor(pixel_crownID), y = yhat))
# p +  geom_hline(yintercept = 0)  + theme_bw() + geom_boxplot()
# ggsave(paste(out_dir, nm, "_res_distribution_perCrown.png", sep=""), width = 11.8, height = 4.93)

p <-  ggplot(output.residuals, aes(x = (y), y = yhat))
p +  geom_hline(yintercept = 0)  + theme_bw() + geom_point()

ggsave(paste(out_dir, nm, "_res_distribution_perSize.png", sep=""), width = 11.8, height = 4.93)


# prova <- output %>%
#   group_by(modelID, pixel_crownID) %>%
#   summarise(y_md = mean(yhat))
# prova <- inner_join(prova, test.data.y, by = "pixel_crownID")
# names(prova)[4] <- "y"
# prova$dY <-  prova$y_md - prova$y
# prova$dY_perc <- prova$dY / prova$y
#
# p <-  ggplot(prova, aes(x = (y), y = dY_perc))
# p +  geom_hline(yintercept = 0)  +   theme_bw() + geom_boxplot(aes(x=y, y=dY_perc,group=y), notch=T, fill="grey", width=0.008,  alpha=0.2) +
#   ylab(bquote('( '*delta['(ob-md)']/ mu['ob']*')'))+ xlab(nm) + theme(text = element_text(size=20))
#
# ggsave(paste(out_dir, nm, "_residuals_per_unit.png", sep=""), units="in", width = 4.93, height = 4.93)

#residuals compare crown against pixel
# res.compare <- inner_join(prova, pixel.mat, by = "pixel_crownID")
# res.compare <- res.compare[names(res.compare) %in%
#                              c("pixel_crownID",  "modelID.x","y_md","y.y", "dY_perc", "yhat")]
# names(res.compare) <- c("modelID","pixel_crownID", "Crown","dCrown", "Pixel", "y")
# res.compare$dPixel <- (res.compare$Pixel - res.compare$y)/res.compare$y
# res.compare <- unique(res.compare)
# d <- res.compare %>%
#   tidyr::gather("Object","value", c("dCrown", "dPixel"))
# d <- unique(d)
# p <-  ggplot(d, aes(x = (modelID), y = value))
# p  + facet_grid(. ~ pixel_crownID)  + theme_bw()+
#   theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#   geom_boxplot(aes(fill=Object),notch = T, alpha = 0.4, width = 0.04) +
#   ylab(bquote('( '*delta['(ob-md)']/ mu['ob']*')')) + theme(text = element_text(size=18)) + xlab("Crown ID") + geom_hline(yintercept = 0)
#
# ggsave(paste(out_dir, nm, "_residuals_compare.png", sep=""), units="in", width = 11.8, height = 4.93)

# out.sd <- output.sd.daic %>%
#   group_by(pixel_crownID) %>%
#   summarise(yhat = mean(yhat_unc, na.rm = T))

#mean(out.sd$yhat)

#print("rmse")
# rmse(pixel.daic.mat$yhat, pixel.daic.mat$y)
# rmse(crown.based.daic$yhat, crown.based.daic$y)

tab_stats <- data.frame(crown.daic.r2,pixel.daic.r2,
                        mean(cr.up$y_up), mean(cr.lw$y_lw),mean(cr.up$y_up-cr.lw$y_lw),
                        mean(pix.up$y_up-pix.lw$y_lw),  mean(pix.up$y_up), mean(pix.lw$y_lw),
                        rmse(pixel.daic.mat$yhat, pixel.daic.mat$y), rmse(crown.based.daic$yhat, crown.based.daic$y))
colnames(tab_stats) <- c("R2 Crown", "R2 pixel", "PI95 Crown","PI5 Crown","PI95int Crown", "PI95int Pix", "PI95 Pix","PI5 Pix","rmse Crown", "rmse pixel")
#your tables
#titles
pdf(paste(out_dir, nm, "_modeltype_compare.pdf", sep=""))
plot(density(r2.out, bw = 0.015), xlim = c(-0.5,1))
polygon(density(r2.out, bw = 0.015), col="green", border="black")
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
rmse_crown <- rmse(crown.based.daic$yhat, crown.based.daic$y)
rmse_pix_ensamble <- rmse(pixel.daic.mat$yhat, pixel.daic.mat$y)
rmse_pix <- rmse(pixel.mat$yhat, pixel.mat$y)
write.csv(summary(t(dd)), paste(nm,"paramUnc.csv", sep="_"))


stargazer(tab_stats[1,], summary=FALSE, rownames=FALSE)
