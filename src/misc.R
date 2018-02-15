#-----perform_summary------------------------------------------------------------------------------------------
perform_summary <- function(names=NULL,out.name = NULL, in.dir = NULL, out.dir = NULL, out.of.bag = T, normlz = F, weighted = T){
  for(j in names){
    load(file = paste(out.dir, 'models_comps_',j, sep="" ))
    load(file = paste(out.dir, 'models_out_',j, sep="" ))
    load(file = paste(out.dir, 'models_stats_',j, sep="" ))
    mod.r2 <- rep(NA, length(mod.stats))
    if(j == "name"){
      for(bb in 1: length(mod.r2)){
        mod.r2[bb] <- mean(mod.stats[[bb]]$name$R2)
      }
      mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 100)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 100)])
      mask <- mask[order(mask[,2]),]
      pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
      pred.weighted <- rep(0, length(mod.stats[[1]]$name$pred))
      norm.R2 <- scalar1(mask[,2])
      multiplier <- 10^2 #decimalplaces(min(norm.R2))
      weights <- floor(norm.R2 * multiplier)
      
      test.data <- read_csv("test_classification.csv")
      test.data <- test.data[colnames(test.data) %in% c("pixel_crownID", "name", "LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")]
      
      lab.sp <- read.csv("spLabels.csv", stringsAsFactors = F)
      
    }else{
      for(bb in 1: length(mod.r2)){
        mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
      }
      mask <- cbind(which(mod.r2 %in% tail(sort(mod.r2), 100)), mod.r2[mod.r2 %in% tail(sort(mod.r2), 100)])
      mask <- mask[order(mask[,2]),]
      pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
      if(out.of.bag){
        test.data.y <- read_csv(paste(in.dir, "Spectra/CrownTraits_outBag.csv", sep=""))
        test.data.y <- test.data.y[colnames(test.data.y) %in% c("pixel_crownID",j)]
        test.data.x <- read.csv(paste(in.dir, "Spectra/CrownPix_outBag.csv", sep=""))
        crownID <- as.data.frame(test.data.x$pixel_crownID)
        test.data.x <- test.data.x[grepl("band", names(test.data.x))]
        spectra_log_dif_snv=test.data.x[, colSums(is.na(test.data.x)) == 0]
        if(normlz){  spectra_log_dif_snv <-t(diff(t(log(as.matrix(spectra_log_dif_snv))),differences=1, lag=3))}    
        
      }else{
        predictions <- sapply(mod.stats, "[[", 2)
      }
    }
    
    if(j == "name"){
      rm(pred.weighted, pred.weight, pred, predictions)
      tkn <- 0
      for(bb in (mask[,1])){
        tkn <- tkn + 1
        if(exists("pred.weighted")){
          predictions <- cbind(predictions, matrix(mod.stats[[bb]]$name$pred, nrow(test.data)))
          pred.weighted <- cbind(pred.weighted, matrix(mod.stats[[bb]]$name$pred, nrow(test.data), as.integer(mask[tkn,2])))
        }else{
          predictions <- matrix(mod.stats[[bb]]$name$pred, nrow(test.data))
          pred.weighted <-  matrix(mod.stats[[bb]]$name$pred, nrow(test.data), as.integer(mask[tkn,2]))
        }
      }
      pred=rep(NA, length(mod.stats[[1]]$name$pred))
      pred.weight=rep(NA, length(mod.stats[[1]]$name$pred))
      for(ii in 1: length(mod.stats[[bb]]$name$pred)) {
        temp.freq <- table(predictions[ii,])
        temp.freq.weight <- table(pred.weighted[ii,])
        pred[ii] <- names(temp.freq)[which(temp.freq == max(temp.freq))]
        pred.weight[ii] <- names(temp.freq.weight)[which(temp.freq.weight == max(temp.freq.weight))]
      }
      pred <- as.data.frame(cbind(test.data$pixel_crownID, pred))
      pred.weight <- as.data.frame(cbind(test.data$pixel_crownID, pred.weight))
      colnames(pred) = c("pixel_crownID","spID")
      colnames(pred.weight) = c("pixel_crownID","spID_w")
    }else{
      if(out.of.bag){
        test.PLS = data.frame( X=I(as.matrix(spectra_log_dif_snv)))
        rm(output)
        for(jj in mask[,1]){
          pls.mod.train <- mod.out[[jj]]
          optim.ncomps <- mod.comps[jj]
          pred.val.data <- predict(pls.mod.train[[j]], newdata = test.PLS, ncomp=optim.ncomps, type='response')
          
          #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
          if(!exists("output")){
            output <- cbind(crownID, rep(jj, dim(test.data.x)[1]), as.vector(pred.val.data))
          }else{
            output <- rbind(output, cbind(crownID, rep(jj, dim(test.data.x)[1]), as.vector(pred.val.data)))
          }
        }
      }else{
        #????????????????
        pred=rep(NA, length(mod.stats[[1]]$pred))
        pred.weight=rep(NA, length(mod.stats[[1]]$pred))
        if(!exists("pred.weighted")){
          pred.weighted <- mod.stats[[bb]]$pred * mod.r2[bb]
          pred <- mod.stats[[bb]]$pred
        }else{
          pred.weighted <- pred.weighted + mod.stats[[bb]]$pred * mod.r2[bb]
          pred <- pred + mod.stats[[bb]]$pred
        }
      }
    }
    
    if(j == "name"){
      rm(evalFinal)
      evalFinal <- merge(lab.sp, test.data, by = "name")
      evalFinal <-unique(evalFinal)
      evalFinal <- merge(evalFinal, pred.weight, by = "pixel_crownID")
      evalFinal <- merge(evalFinal, pred, by = "pixel_crownID")
      correct <- evalFinal$spID_train == evalFinal$spID
      accuracy <- sum(as.numeric(correct))/length(correct)
      correct.weigth <- evalFinal$spID_train == evalFinal$spID_w
      accuracy.weight <- sum(as.numeric(correct.weigth))/length(correct.weigth)
      print(accuracy)
      print(accuracy.weight)
      
      if(!exists("out")){
        out <- c(j, accuracy, accuracy.weight)
      }else{
        out <- rbind(out, c(j, accuracy, accuracy.weight))
      }
    }else{
      if(out.of.bag){
        colnames(output) <- c("pixel_crownID", "modelID", "yhat")
        if(weighted){
          foo <- mask[,2]/sum(mask[,2])
          weights <- NULL
          for(temp in 1:length(foo)){
            weights<- c(weights, rep(foo[temp], dim(test.data.x)[1]))
          }
        }else{
          weights <- rep(1, length(output$yhat))
        }
        #output$yhat <- output$yhat * weights
        pixel.mat <- inner_join(output, test.data.y, by = "pixel_crownID")
        pixel.based <- cor(pixel.mat$yhat, eval(parse(text = paste('pixel.mat$',j,sep=""))))^2
        crown.based <- pixel.mat %>%
          group_by(pixel_crownID) %>%
          summarise(yhat = mean(yhat))
        crown.based <- inner_join(crown.based, test.data.y, by = "pixel_crownID")
        colnames(crown.based) <- c("pixel_crownID", "yhat", "y")
        #crown.based <- cor(crown.based$yhat, crown.based$y)^2
        #plot(lm(crown.based$yhat ~ crown.based$y))
        par(mfrow=c(2,1))
        plot(crown.based$yhat ~ crown.based$y, main = "1:1 crown object leaf N (%)", cex = 2, xlab = "OOB observed N (%)",ylab = "OOB predicted N (%)", xlim = c(0.5, 2.5), ylim = c(0.5, 2.5), pch=21, bg="blue")
        abline(0,1, col="red")
        plot(pixel.mat$yhat ~ eval(parse(text = paste('pixel.mat$',j,sep=""))), main = "1:1 pixel leaf N (%)", cex = 2, xlab = "OOB observed N (%)",ylab = "OOB predicted N (%)", xlim = c(0.5, 2.5), ylim = c(0.5, 2.5), pch=21, bg="blue")
        abline(0,1, col="red")
        print(paste(pixel.based, "|", crown.based))
        if(!exists("out")){
          out <- c(j, pixel.based, crown.based)
        }else{
          out <- rbind(out, c(j, pixel.based, crown.based))
        }
      }else{
        predictions <- predictions / length(mod.stats[[1]]$pred)
        pred.weighted <- pred.weighted /sum(mod.r2[mask])
        if(!exists("out")){
          out <- c(j, cor(pred, mod.stats[[1]]$obs)^2, cor(pred.weighted, mod.stats[[1]]$obs)^2)
        }else{
          out <- rbind(out, c(j, cor(pred, mod.stats[[1]]$obs)^2, cor(pred.weighted, mod.stats[[1]]$obs)^2))
        }
      }
    }
    write_csv(data.frame(pred.weighted), paste(out.dir, "predicted_",j, ".csv", sep=""))
  }     
  write_csv(data.frame(out), paste(out.dir, out.name, sep=""))
  return(out)
}

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



PLS_DA <- function(loops = 1000, names = c("name"), norm = T, j = 1, proportions = 0.7,
                   out.dir = NULL, in.dir = NULL){
  cr.Traits <- read.csv(paste(in.dir, "Spectra/crownTraits.csv",sep=""), stringsAsFactors = F)
  nCrowns <- dim(cr.Traits)[1]
  names = "name"
  #aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_1.csv', sep = ''), in.dir)
  
  #matrix to store performances in
  mod.out = vector("list", loops)
  mod.stats = vector("list", loops)
  mod.comps = rep(NA, loops)
  for(laps in 1:loops) {
    aug.spectra <- imp.spectra(paste('Permutations/onePix1Crown_', laps, '.csv', sep = ''), in.dir)
    
    aug.spectra$X= NULL
    aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
    X<- aug.spectra[grepl("band", names(aug.spectra))]
    X=X[, colSums(is.na(X)) == 0]
    Y <- aug.spectra[,names(aug.spectra) %in% names]
    names(Y) <- "name"
    aug.X <- data.frame(aug.spectra$name, Y, X)
    # Subset data into cal/val by site
    prova <- data.frame(unique(cbind(aug.spectra$pixel_crownID, aug.spectra$name)))
    colnames(prova) <- c("pixel_crownID", "species")
    set.seed(14)
    out <- prova %>% 
      group_by(species) %>%
      filter(pixel_crownID %in% sample(pixel_crownID, ceiling(proportions*length(pixel_crownID))))
    
    pixID <- data.frame(as.integer(levels(out$pixel_crownID)[out$pixel_crownID]))
    names(pixID) <- "pixel_crownID"
    train.data <- inner_join(pixID, aug.spectra, by = "pixel_crownID")
    test.data <- anti_join(aug.spectra, by = "pixel_crownID", pixID)
    train.data <- train.data[order(train.data$pixel_crownID),]
    test.data <- test.data[order(test.data$pixel_crownID),]
    
    write.csv(test.data, "test_classification.csv")
    write.csv(train.data, "train_classification.csv")
    
    print(paste(laps, names(Y[j])))
    # Run calibration PLSR analysis to select optimal number of components
    labelsSp <- as.data.frame(cbind(unique(Y), as.numeric(factor(unique(Y)))), stringsAsFactors = F)
    colnames(labelsSp) <- c("name", "id")
    lab.test <- merge(test.data[,colnames(test.data)%in%c("pixel_crownID", "name")], labelsSp[labelsSp$name %in% test.data$name, ], by = "name")
    lab.test <- lab.test[order(lab.test$pixel_crownID),]
    write.csv(lab.test, "labtest.csv")
    
    lab.train <- merge(train.data[,colnames(train.data)%in%c("pixel_crownID", "name")], labelsSp[labelsSp$name %in% train.data$name, ], by = "name")
    lab.train <- lab.train[order(lab.train$pixel_crownID),]
    
    pls.mod.train <- pls.cal(train.data, numcomps=35, j=1,nm = names, normalz = norm, lab.train)
    #calculate number of components given min test PRESS or RMSEP
    #directly take ncomps
    # foo <- unlist(lapply(pls.mod.train$name$RMSEP, min))
    # nOcomps <- which(foo == min(foo, na.rm=T))
    if(!is.nan(pls.mod.train$name$validation$PRESS[1])){
      optim.ncomps <- which(pls.mod.train$name$validation$PRESS==min(pls.mod.train$name$validation$PRESS))
      #no.comps <- which(tmp.pls$RMSEP[[ncomps]]== min(tmp.pls$RMSEP[[ncomps]]))
      best.pred <- round(pls.mod.train$name$fitted.values)[,,optim.ncomps]
      sum(best.pred==lab.train$id)/length(lab.train$id)
      pred.val.data <- predict.pls(pls.mod.train, test.data, nm = names,norm = norm, optim.ncomps=optim.ncomps,j, lab.test)
      mod.out[[laps]] <- pls.mod.train 
      mod.stats[[laps]] <- pred.val.data
      mod.comps[laps] <- optim.ncomps
    }else{
      warning("pls didn't converge!")
      mod.out[[laps]] <- pls.mod.train 
      out.data = data.frame(lab.test$name,rep(NA, length(lab.test$name)),rep(NA, length(lab.test$name)), -9999, -9999)
      names(out.data) = c("obs","pred","res", "R2", "MSE")
      mod.stats[[laps]] <- out.data
      mod.comps[laps] <- 0
    }
    setwd(out.dir)
  }
  save(mod.out, file = paste("models_out_", names[j],  sep = ""))
  save(mod.stats,  file = paste("models_stats_", names[j],  sep = ""))
  save(mod.comps, file = paste("models_comps_", names[j],  sep = ""))
  
}

#-----scale to 0:1------------------------------------------------------------------------------------------
scalar1 <- function(x) {x / sqrt(sum(x^2))}
#-----decimalplaces------------------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}



#----
#relies on performaGator
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
}