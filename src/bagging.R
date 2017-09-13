# load("~/Documents/Projects/TraitsOnHeaven/OSBS/outputs/models_stats_name")
# load("~/Documents/Projects/TraitsOnHeaven/ALL/outputs/models_out_name")
# load("~/Documents/Projects/TraitsOnHeaven/ALL/outputs/models_comps_name")
load(file = paste(out.dir, 'models_comps_',j, sep="" ))
load(file = paste(out.dir, 'models_out_',j, sep="" ))
load(file = paste(out.dir, 'models_stats_',j, sep="" ))
mod.r2 <- rep(NA, length(mod.stats))
deltaPRESS <- rep(NA, length(mod.stats))
names <- "name"
scalar1 <- function(x) {x / sqrt(sum(x^2))}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

for(bb in 1: length(mod.r2)){
  if(names == "name"){
    mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
  }else{
    mod.r2[bb] <- mean(mod.stats[[bb]]$R2)
    deltaPRESS[bb] <- eval(parse(text = paste("mod.out[[bb]]$",names, sep="")))$validation$PRESS[mod.comps[bb]]
  }
}

deltaPRESS <- max(deltaPRESS) - deltaPRESS
boxplot(mod.r2)
mask <- which(mod.r2 > 0.15)
mask <- which(mod.r2 %in% tail(sort(mod.r2), 100))
predictions <- rep(0, length(mod.stats[[1]]$pred))
pred.weighted <- rep(0, length(mod.stats[[1]]$pred))
pred.deltaPRESS <- rep(0, length(mod.stats[[1]]$pred))
sum.deltaPRESS <- 0

if(names == "name"){
  norm.R2 <- scalar1(mod.r2[mask])
  multiplier <- 10^4 #decimalplaces(min(norm.R2))
  weights.class <- floor(norm.R2 * multiplier)
}  
predictions <- sapply(mod.stats, "[[", 2)
rm(predictions.weighted)
tkn <- 0
for(bb in (mask)){
  tkn <- tkn + 1
  if(exists("predictions.weighted")){
    predictions.weighted <- cbind(predictions.weighted , 
                                  matrix(predictions[, bb], nrow(predictions), as.integer(weights.class[tkn])))
  }else{
    predictions.weighted <- matrix(predictions[, bb], nrow(predictions), as.integer(weights.class[tkn]))
  }
}

for(bb in mask){
  if(names == "name"){
    #predictions <- sapply(mod.stats, "[[", 2)
    pred=rep(NA, length(mod.stats[[1]]$pred))
    pred.weight=rep(NA, length(mod.stats[[1]]$pred))
    
    for(ii in 1: length(mod.stats[[1]]$pred)) {
      temp.freq <- table(predictions[ii,])
      temp.freq.weight <- table(predictions.weighted[ii,])
      pred[ii] <- names(temp.freq)[which(temp.freq == max(temp.freq))]
      pred.weight[ii] <- names(temp.freq.weight)[which(temp.freq.weight == max(temp.freq.weight))]
    }
  }else{
    predictions <- predictions + mod.stats[[bb]]$pred
    pred.weighted <- pred.weighted + mod.stats[[bb]]$pred * mod.r2[bb]
    pred.deltaPRESS <- pred.deltaPRESS + mod.stats[[bb]]$pred * exp(-0.5 *deltaPRESS[bb])
    sum.deltaPRESS <- sum.deltaPRESS + exp(-0.5 *deltaPRESS[bb])
  }
}

pred
pred.weight
if(names == "name"){
  correct <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred)
  accuracy <- sum(as.numeric(correct))/length(correct)
  
  correct.weigth <- as.numeric(factor(mod.stats[[1]]$obs)) == as.numeric(pred.weight)
  accuracy.weight <- sum(as.numeric(correct.weigth))/length(correct.weigth)
}
accuracy
accuracy.weight

predictions <- predictions / length(mod.stats[[1]]$pred)
pred.weighted <- pred.weighted /sum(mod.r2[mask])
pred.deltaPRESS <- pred.deltaPRESS /sum.deltaPRESS

cor(predictions, mod.stats[[1]]$obs)^2
cor(pred.weighted, mod.stats[[1]]$obs)^2
cor(pred.deltaPRESS, mod.stats[[1]]$obs)^2

if(names = "name"){
  table(mod.stats)
}

