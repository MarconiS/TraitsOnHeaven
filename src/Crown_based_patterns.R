rm(list=ls(all=TRUE))   # clear workspace
NeonSite = "OSBS"
library(tidyverse)
setwd('/Users/sergiomarconi/Documents/Projects/TraitsOnHeaven/')
names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
listPlot <- (read.csv('OSBS/inputs/Plot_class.csv', header = T, stringsAsFactors = F))
out.dir = paste(getwd(), NeonSite,"outputs/", sep="/")
in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")

for (i in listPlot[['Plot_ID']]){#forToday[1]) {
  f <- paste("~/Documents/Projects/TraitsOnHeaven/", NeonSite, "/traits_regression_csv/", NeonSite, "_", i, "_", "aveOutputs" ,sep="" )
  f.sp <-  paste("~/Documents/Projects/TraitsOnHeaven/", NeonSite, "/species_classification_csv/", NeonSite, "_", i, "_", "species_aveOutputs" ,sep="" )
  load(f)
  foo <- unlist(out.ave)
  crowns <- length(out.ave[[1]])
  categories <- c( rep("LMA_g.m2", crowns),rep("d13C", crowns),rep("d15N", crowns),
                  rep("C_pct", crowns),rep("N_pct", crowns), rep("P_pct", crowns))
  
  plot.id <- rep(i, length(names) * crowns)
  foo <- data.frame(foo, factor(categories), as.character(plot.id))
  colnames(foo) <- c("value", "trait", "Plot_ID")
  foo <- inner_join(foo, listPlot, by = "Plot_ID")
  ####
  load(f.sp)
  foo.sp <- rep(unlist(out.ave), length(categories) / crowns)
  species.id <- read_csv(paste(in.dir, "Spectra/speciesName.csv", sep=""))
  foo <- data.frame(foo, (foo.sp))
  colnames(foo) <- c("value", "trait", "Plot_ID", "class","species")
  foo <- inner_join(foo, species.id, by = "species")
  
  ####
  if(!exists("allPlotsOut")){
    allPlotsOut <- foo
  }else{
    allPlotsOut <- rbind(allPlotsOut, foo)
  }
}

for(trName in names){
  p <- ggplot(subset(allPlotsOut, trait == trName), aes(factor(Plot_ID), value))
  p + theme_bw(base_size = 16) + geom_violin(aes(fill = factor(class)), draw_quantiles = c(0.25, 0.5, 0.75)) + 
    labs(title = paste(trName, "distribution in OSBS\n"), alpha = "",x = "Plot ID", y = "value", color = "Species\n", fill = "Ecosystem type\n") +
    geom_jitter(aes(colour = factor(name), alpha = 0.5), height = 0, width = 0.1) + scale_colour_manual(values = c("yellow", "black", "red", "gray", "orange", "brown")) +
    theme(axis.text.x=element_text(angle=90,hjust=1)) 
  ggsave(paste(out.dir, trName, "_dist.png", sep=""), width=30, height = 20, units="in")
}

#facet_grid(. ~ plot, scales = "free")


