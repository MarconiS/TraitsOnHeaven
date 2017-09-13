perform_summary <- function(names=NULL,out.name = NULL){
  
  for(j in names){
    load(file = paste(out.dir, 'models_comps_',j, sep="" ))
    load(file = paste(out.dir, 'models_out_',j, sep="" ))
    load(file = paste(out.dir, 'models_stats_',j, sep="" ))
    mod.r2 <- rep(NA, length(mod.stats))
    
    for(bb in 1: length(mod.r2)){
      mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
    }
    mask <- which(mod.r2 %in% tail(sort(mod.r2), 100) & mod.r2 > 0.0)
    pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
    if(j == "name"){
      norm.R2 <- scalar1(mod.r2[mask])
      multiplier <- 10^4 #decimalplaces(min(norm.R2))
      weights <- floor(norm.R2 * multiplier)
    }  
    predictions <- sapply(mod.stats, "[[", 2)
    rm(pred.weighted, pred.weighted, pred)
    if(j == "name"){
      tkn <- 0
      for(bb in (mask)){
        tkn <- tkn + 1
        if(exists("pred.weighted")){
          pred.weighted <- cbind(pred.weighted, matrix(predictions[, bb], nrow(predictions), as.integer(weights[tkn])))
        }else{
          pred.weighted <- matrix(predictions[, bb], nrow(predictions), as.integer(weights[tkn]))
        }
      }
    }
    for(bb in mask){
      if(j == "name"){
        pred=rep(NA, length(mod.stats[[1]]$pred))
        pred.weight=rep(NA, length(mod.stats[[1]]$pred))
        
        for(ii in 1: length(mod.stats[[1]]$pred)) {
          temp.freq <- table(predictions[ii,])
          temp.freq.weight <- table(pred.weighted[ii,])
          pred[ii] <- names(temp.freq)[which(temp.freq == max(temp.freq))]
          pred.weight[ii] <- names(temp.freq.weight)[which(temp.freq.weight == max(temp.freq.weight))]
        }
      }else{
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
      correct <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred)
      accuracy <- sum(as.numeric(correct))/length(correct)
      correct.weigth <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred.weight)
      accuracy.weight <- sum(as.numeric(correct.weigth))/length(correct.weigth)
      out[[j]] <- c(j, accuracy, accuracy.weight)
    }else{
      predictions <- predictions / length(mod.stats[[1]]$pred)
      pred.weighted <- pred.weighted /sum(mod.r2[mask])
      out[[j]] <- c(j, cor(pred, mod.stats[[1]]$obs)^2, cor(pred.weighted, mod.stats[[1]]$obs)^2)
    }
    write_csv(out, paste(out.dir, out.name, sep="/"), append = T)
  }
}