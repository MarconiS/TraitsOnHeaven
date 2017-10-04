# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

plot_upper_cor <- function(cormat){
  
  #this is when uncertainty associated to species is 0
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
}

library(ggplot2)
library(reshape2)

#is it better predict? or better know species?
library(tidyverse)
wd = '~/Documents/Projects/TraitsOnHeaven/'
NeonSite = "OSBS"
names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
setwd(wd)
path = paste(wd, NeonSite, sep="/")
setwd(path)
out.dir = paste(getwd(), "/outputs/", sep="")
in.dir = paste(getwd(), "/inputs/", sep="")
cr.Traits <- read.csv("train_for_baseline.csv", stringsAsFactors = F)
cr.Traits <- cr.Traits[,colnames(cr.Traits) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]
lab.sp <- read.csv("spLabels.csv", stringsAsFactors = F)
lab.for.weights <- read.csv("labTrain.csv", stringsAsFactors = F)
cr.Traits <- inner_join(lab.for.weights, cr.Traits, by = "name")
sp.pred <- as.matrix(read_csv(paste(out.dir, "predicted_name.csv",sep="")))
sp.freq <- apply(sp.pred, 1, function(x){table(x)/length(x)})
obs <- read.csv("test_for_baseline.csv", stringsAsFactors = F)
obs <- obs[,colnames(obs) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]
traits.pred <- read.csv(paste(out.dir, "bagged_preds.csv",sep=""), stringsAsFactors = F)

ave.Sp <- cr.Traits %>%
  group_by(name) %>%
  summarize(link = mean(id), LMA_g.m2 = mean(LMA_g.m2),d13C = mean(d13C),d15N = mean(d15N),C_pct = mean(C_pct),N_pct = mean(N_pct),P_pct = mean(P_pct))

sd.Sp <- cr.Traits %>%
  group_by(name) %>%
  summarize( link = mean(id), LMA_g.m2 = sd(LMA_g.m2), d13C = sd(d13C),d15N = sd(d15N),C_pct = sd(C_pct),N_pct = sd(N_pct),P_pct = sd(P_pct))

# assign average value given the relative probability (frequency) of each  individual species to an ITC
trait.from.class <- as.data.frame(matrix(NA, nrow = length(obs[,1]), ncol = length(names)+1))
colnames(trait.from.class) <- c("name", names)
for(ii in 1: length(obs[,1])){
  pix.id <- obs$name[ii]
  load.sp <- as.numeric(names(sp.freq[[ii]]))
  sp.trait.test <- ave.Sp[ave.Sp$link %in% load.sp, ]
  sp.trait.test <- sp.trait.test[order(sp.trait.test$link),]
  sp.trait.test[colnames(sp.trait.test) %in% names] <- sp.trait.test[colnames(sp.trait.test) %in% names] * as.numeric(sp.freq[[ii]])
  trait.from.class[ii, ] <- c(pix.id,  (round(apply(sp.trait.test[colnames(sp.trait.test) %in% names], 2, sum), 2)))
}

#see why only separately!
for(ii in 2:7){
  trait.from.class[,ii]<-as.numeric(trait.from.class[,ii])
}

#calculate averages from global dataset
globNet <- read_csv("~/Documents/Projects/TraitsOnHeaven/globNet.csv")
average.traits <- globNet %>%
  group_by(SpID) %>%
  summarize(LMA_g.m2 = 10^(mean(`log LMA`,na.rm=TRUE)), N_pct= 10^(mean(`log Nmass`, na.rm=TRUE))) #, P_pct = mean(`log Pmass`, na.rm=TRUE))
colnames(average.traits)[1]<- "name"

colnames(obs)[1]<-"name"
obs.sorted <- as.data.frame((obs$name), stringsAsFactors = F)
colnames(obs.sorted)<-"name"
species_outs <- inner_join(obs.sorted, ave.Sp, by="name")
global.species_outs <- as.data.frame(obs[obs$name %in% average.traits$name, colnames(obs) %in% c("name", "LMA_g.m2", "N_pct")], stringsAsFactors = F)
colnames(species_outs)[1] <- "name"
global.species_outs <- merge(global.species_outs,average.traits, by="name")
global.species_outs <- unique(global.species_outs)
colnames(global.species_outs)[4:5]<- c("pred LMA_g.m2", "pred N_pct")
sp.globnet.cor <- round(cor(global.species_outs[,2:3], global.species_outs[,4:5]), 2)
traits_out <- cbind(traits.pred, obs)

#ensamble
ensamble <- (trait.from.class[,-1] + traits_out[,1:6])/2
ensamble.pure <- (species_outs[,-c(1,2)] + traits_out[,1:6])/2


sp.cor <- round(cor(species_outs[,-c(1,2)], obs[,2:7]), 2)
tr.sp.cor <- round(cor(trait.from.class[,-1], obs[,2:7]), 2)
tr.cor <- round(cor(traits_out[,-7])[1:6, 7:12], 2)
en.cor <- round(cor(ensamble, obs[,2:7]), 2)
en.pure.cor <- round(cor(ensamble.pure, obs[,2:7]), 2)


sp.cor
tr.cor
tr.sp.cor
en.cor
en.pure.cor
plot_upper_cor(sp.cor)
plot_upper_cor(tr.cor)
plot_upper_cor(en.pure.cor)
plot_upper_cor(round(en.pure.cor-sp.cor, 2))
plot_upper_cor(round(en.pure.cor-tr.cor, 2))
plot_upper_cor(round(sp.cor-tr.cor, 2))
plot_upper_cor(round(en.cor - tr.cor, 2))

plot_upper_cor(tr.sp.cor)

plot_upper_cor(tr.cor)
plot_upper_cor(round(sp.cor - tr.cor, 2))
tr.glob <- round(cor(traits_out[,c(1,5)], traits_out[,c(8,12)]), 2)
plot_upper_cor(round(sp.globnet.cor - tr.glob, 2))



