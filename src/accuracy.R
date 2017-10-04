load(file = paste(out.dir, 'models_comps_name', sep="" ))
load(file = paste(out.dir, 'models_out_name', sep="" ))
load(file = paste(out.dir, 'models_stats_name', sep="" ))
mod.r2 <- rep(NA, length(mod.stats))
test.data <- read_csv("test_classification.csv")
test.data <- test.data[colnames(test.data) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]

lab.sp <- read.csv("spLabels.csv", stringsAsFactors = F)
for(bb in 1: length(mod.r2)){
  mod.r2[bb] <- mean(mod.stats[[bb]]$name$R2)
}
mask <- which(mod.r2 %in% tail(sort(mod.r2), 100))
allBand=read.csv("./inputs/Spectra/neon_aop_bands.csv")
allData=read.csv("./inputs/Spectra/CrownPix_norm.csv")
crownID = test.data$pixel_crownID
allData = allData[allData$pixel_crownID %in% crownID,]
X.all <- allData[grepl("band", names(allData))]
normalz = T
if(normalz){X.all <-t(diff(t(log(as.matrix(X.all))),differences=1, lag=3))}

X.all=X.all[, colSums(is.na(X.all)) == 0]
test.PLS = data.frame(X=I(X.all))

rm(output)
for(jj in mask){
  pls.mod.train <- mod.out[[jj]]
  optim.ncomps <- mod.comps[jj]
  pred.val.data <- round(predict(pls.mod.train$name,newdata = test.PLS, ncomp=optim.ncomps, type="response"))
  
  if(!exists("output")){
    output <- cbind(crownID, jj, pred.val.data)
  }else{
    output <- rbind(output, cbind(crownID, jj, pred.val.data))
  }
}
write.csv(output, "species_classification_object.csv", row.names=FALSE)


foo <- as.data.frame(output)
rm(test.out.crown)
for(index in mask){
  tmp <- foo[foo$jj %in% index, ]
  for(Cid in unique(crownID)){
    maj.vote <- tmp[tmp$crownID %in% Cid, ]
    maj.vote <- names(which(table(maj.vote$pred.val.data)== (max(table(maj.vote$pred.val.data)))))
    token = 0
    while(token < length(maj.vote)){
      token = token +1
      if(!exists("test.out.crown")){
        test.out.crown <- as.integer(cbind(Cid, index, maj.vote[token]))
      }else{
        test.out.crown <- rbind(test.out.crown, as.integer(cbind(Cid, index, maj.vote[token])))
      }
    }
  }
}
rm(maj.vote)
test.out.crown <- as.data.frame(test.out.crown)
colnames(test.out.crown) <- c("Cid", "index", "maj.vote")
#test.out.crown <- as.data.frame(cbind(as.integer(test.out.crown[,1]), as.integer(test.out.crown[,2]), as.integer(test.out.crown[,3])))
for(Cid in unique(crownID)){
  junkie <- test.out.crown[test.out.crown$Cid %in% Cid, ]
  token = 0
  foo.vote <- names(which(table(junkie$maj.vote)== (max(table(junkie$maj.vote)))))
  while(token < length(foo.vote)){
    token = token + 1
    if(!exists("maj.vote")){
      maj.vote <-  as.integer(cbind(Cid, foo.vote[token]))
    }else{
      maj.vote <- rbind(maj.vote, as.integer(cbind(Cid, foo.vote[token])))
    }
  }
}
maj.vote <- as.data.frame(maj.vote, row.names=FALSE)
colnames(maj.vote) <- c("pixel_crownID", "class")
test.data <- read_csv("test_classification.csv")
test.data <- test.data[colnames(test.data) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]

#maj.vote$pixel_crownID = as.numeric(labels(maj.vote$pixel_crownID))
test.data <- inner_join(test.data, lab.sp[,c(1,3)], by = "name")
test.data <- unique(test.data)
get_test <- inner_join(test.data, maj.vote, by = "pixel_crownID")
get_test$realLabel <- as.numeric(factor(get_test$name))
accuracy = sum(get_test$spID_train == get_test$class) / length(get_test$pixel_crownID)
accuracy
baseline <- read.csv("baseline.csv")
get_test$baseline <- baseline$X0
accuracy = sum(get_test$baseline == get_test$name) / length(get_test$pixel_crownID)
accuracy
output <- read_csv("species_classification_object.csv")
colnames(output) <- c("pixel_crownID", "modelID", "class")
pixel.lv <- inner_join(test.data, output, by = "pixel_crownID")
accuracy = sum(pixel.lv$spID_train == pixel.lv$class) / length(pixel.lv$pixel_crownID)
accuracy




