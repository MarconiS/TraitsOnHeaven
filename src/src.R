#--------------------------------------------------------------------------------------------------#
closeAll <- function()
{
  # Close all devices and delete all variables.
  
  graphics.off()          # close any open graphics
  closeAllConnections()   # close any open connections to files
}
#--------------------------------------------------------------------------------------------------#

imp.spectra <- function(f, pwd)
{
  # Import dry spectra dataset
  aug.data <- read.table(paste(pwd,f,sep=""), header=TRUE,sep=",")
  return(aug.data)
}
#--------------------------------------------------------------------------------------------------#
sp.corr <- function(X,Y, pwd)
{
  spec_corr <- data.frame(cor(X, Y))
  waves <- data.frame(seq(1,1050,1),seq(400,2498,2))
  mean.spec <- colMeans(X)
  spec.quant <- apply(X,2,quantile,probs=c(0.05,0.95))
  # Output correlation data
  pdf(paste(pwd, '/','FFT_Spectra_Correlations.pdf',sep=""),height=12,width=8)
  par(mfrow=c(6,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
  matplot(mean.spec, type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n",ylim=c(0,1.4))
  ind <- pretty(seq(from = 400, to = 2498, by = 2)) # Using pretty to standardize the axis
  ind <- ind[ind >= 400 & ind <= 2498]
  ind <- (ind - 399) / 1
  axis(1, ind, colnames(X)[ind]) # add column names to wavelengths
  # CIs
  lines(spec.quant[1,],lty=1,col="dark grey")
  lines(spec.quant[2,],lty=1,col="dark grey")
  legend("topleft",legend=c("Mean","95% CI"),col=c("black","dark grey"),lwd=3)
  
  for (j in 1:length(spec_corr)) {
    plot(waves[,2],spec_corr[,j],xlab="WAVELENGTH (nm)",ylab="CORRELATION", main=names(spec_corr)[j], cex=0.01)
    abline(h=0,lty=2,lwd=1.5,col="grey80")
  }
  dev.off()
}

predict.regsubsets =function (object ,newdata ,id){
  form=as.formula(Y~.)
  mat=model.matrix (form ,newdata )
  coefi =coef(object ,id=id)
  xvars =names (coefi )
  return(mat[,xvars ]%*% coefi)
}

cut.set<-function(aug.X,out.dir, c.id, prop = 0.7){
  species <- unique(aug.X$aug.spectra.name)
  
  # Sample proportion for cal data
  train.data <- 0
  test.data <- 0
  j <- 1
  cr.id <- NULL
  for (i in as.character(species)){
    set.seed(1)
    #subset by species
    temp.data <- aug.X[which(aug.X$aug.spectra.name==i),]
    rows <- sample(1:nrow(temp.data),floor(prop*nrow(temp.data)))
    foo <- c.id[which(aug.X$aug.spectra.name==i)]
    # set.seed(1)
    # foo <- sample(1:length(foo),floor(prop*length(foo)))
    
    cr.id <- c(cr.id, foo[rows])
    
    cal.data = droplevels(temp.data[rows,])
    val.data = droplevels(temp.data[-rows,])
    
    if(j==1){
      train.data <- cal.data
      test.data <- val.data
    } else {
      train.data <- rbind(train.data,cal.data)
      test.data <- rbind(test.data,val.data)
    }
    
    j <- j+1
  }
  return(list(train=train.data, test=test.data, cr.id=cr.id))
}


# Calibration PLS model ---------------------------------------------------
pls.cal <- function(train.data, comps, j, nm, normalz = F){
  spectra <- as.matrix(train.data[grepl("band", names(train.data))])
  traits <- as.matrix(train.data[,names(train.data) %in% nm])
  pls.summ <- list()
  spectra_log_dif_snv <- spectra
  #spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  if(normalz){spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))}
  
  leaf.trait <- traits[,j]
  if(is.character(leaf.trait)){
    #clean
    spectra_log_dif_snv=spectra_log_dif_snv[, colSums(is.na(spectra_log_dif_snv)) == 0]
    tmp.y<-matrix(as.numeric(factor(leaf.trait)),ncol=1) # make numeric matrix
    #train.PLS = data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    tmp.pls<-caret:::plsda(y = factor(tmp.y), x = spectra_log_dif_snv,  ncomp = comps,  probMethod = "softmax", type = "class")

  }else{
    train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = plsr(Y ~ X,scale=F, ncomp=comps,validation="LOO", trace=TRUE, method = "oscorespls", data = train.PLS) 
    
  }
  pls.summ[[nm[j]]] <- tmp.pls 
  
  return(pls.summ)
}

pls.calImg <- function(train.data, comps, j, nm, normalz = F){
  spectra <- as.matrix(train.data[grepl("band", names(train.data))])
  spectra<- spectra[,-c(seq(1,8))]
  #train.data <- train.data[ , apply(train.data, 2, function(x) !any(is.na(x)))]
  traits <- as.matrix(train.data[,1])
  pls.summ <- list()
  spectra_log_dif_snv <- spectra
  #spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  if(normalz){  spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))
  }
  
  #leaf.trait <- traits[,j]
  leaf.trait <- traits
  #get only 80% training, 20% to figure optimum components out?    
  train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
  tmp.pls = plsr(Y ~ X,scale=F, ncomp=comps,validation="LOO", trace=TRUE, method = "oscorespls", data = train.PLS) 
  pls.summ[[nm[j]]] <- tmp.pls 
  
  return(pls.summ)
}


#--------------------------------------------------------------------------------------------------#
opt.comps <- function(pls.mod.train, Y, j){
  ncomps <- NA
  tmp.pls = eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep="")))
  if(j != 1){
    ncomps <- which(tmp.pls$validation$PRESS==min(tmp.pls$validation$PRESS[3:length(tmp.pls$validation$PRESS)]))
  }else{
    ncomps <- tmp.pls$ncomp
  }
  return(ncomps)
}
#--------------------------------------------------------------------------------------------------#


Q2 <- function(test.PLS, tmp.pls)
{
  dims <- dim(test.PLS)
  PRESS = tmp.pls$validation$PRESS
  SS = sum((tmp.pls$Y)^2)
  TSS = sum((test.PLS$Y-mean(test.PLS$Y))^2)
  Q2 = 1-(PRESS/TSS)  
  
  # Calculate RMSECV
  RMSECV = PRESS/dims[1]
  return(RMSECV)
}

predict.pls <- function(pls.mod.train, test.data, optim.ncomps,j, nm, norm = F){
  spectra <- as.matrix(test.data[grepl("band", names(test.data))])
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  pred <- list()
  spectra_log_dif_snv <- spectra
  #spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  if(norm){  spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))
  }    
  leaf.trait <- traits[,j]
  if(is.character(leaf.trait)){
    spectra_log_dif_snv=spectra_log_dif_snv[, colSums(is.na(spectra_log_dif_snv)) == 0]
    
    tmp.y<-matrix(as.numeric(factor(leaf.trait)),ncol=1) # make numeric matrix
    test.PLS = data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    pred[[nm[j]]] <- as.vector(predict(eval(parse(text = paste('pls.mod.train$', nm[j],sep=""))), 
                                       newdata = spectra_log_dif_snv, ncomp=optim.ncomps, type="class"))
    
  }else{
    test.PLS <- data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    pred[[nm[j]]] <- as.vector(predict(eval(parse(text = paste('pls.mod.train$',
                                                               nm[j],sep=""))), newdata = test.PLS, ncomp=optim.ncomps, type="response"))
  }
  return(pred)
}

res.out <- function(pred.val.data, train.data, test.data, j, nm)
{
  out <- list()
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  # Build output dataset
  # PLSR Summary statistics
  pred.data <- as.vector(eval(parse(text = paste('pred.val.data$',nm[j],sep=""))))
  if(is.character(traits[,j])){
    correct <- as.numeric(factor(traits[,j])) == as.numeric(pred.data)
    accuracy <- sum(as.numeric(correct))/length(correct)
    temp = data.frame(as.numeric(factor(traits[,j])),as.numeric(pred.data),correct, accuracy)
    names(temp) = c("obs","pred","res", "R2")
  }else{
    res <- traits[,j]- pred.data
    MSE.test <- mean(res^2)
    RMSE.test <- sqrt(MSE.test)
    R2 <- 1- length(pred.data) * MSE.test / sum((traits[,j]- mean(traits[,j]))^2)
    ### Output val dataset
    temp = data.frame(traits[,j],pred.data,res, R2, MSE.test)
    names(temp) = c("obs","pred","res", "R2", "MSE")
  }

  out[[names(train.data[,2:7])[j]]] <- temp 
  
  return(temp)
}


# LASSO SELECTION ----------------------------------------------

lasso.asd <- function(test.data,train.data, j, nm, norm = F)
{
  grid =10^ seq (10,-2, length =100)
  traits <- train.data[,names(train.data) %in% nm]
  leaf.spec <- train.data[grepl("band", names(train.data))]
  #leaf.spec <-t(diff(t(log(leaf.spec)),differences=1, lag=3))
  if(norm){    
    leaf.spec <- standardNormalVariate(X = t(diff(t(log(leaf.spec)),differences=1, lag=3)))
  }
  test.traits <- test.data[,names(test.data) %in% nm]
  test.spec <- test.data[grepl("band", names(test.data))]
  #test.spec <-t(diff(t(log(test.spec)),differences=1, lag=3))
  if(norm){    
    test.spec <- standardNormalVariate(X = t(diff(t(log(test.spec)),differences=1, lag=3)))
  }
  out <- list()
  leaf.trait <- traits[,j]
  lasso.mod =cv.glmnet(as.matrix(leaf.spec),as.vector(leaf.trait),alpha =1, lambda =grid)
  
  bestlam =lasso.mod$lambda.min
  tr.MSE = sum((traits[,j] - predict(lasso.mod ,s=bestlam ,newx=as.matrix(leaf.spec)))^2) / length(traits[,j])
  tr.R2 = 1 - length(traits[,j]) / sum((traits[,j]- mean(traits[,j]))^2) * tr.MSE
  
  lasso.pred=predict(lasso.mod ,s=bestlam ,newx=as.matrix(test.spec))
  lasso.mse <- mean((lasso.pred - test.traits[,j])^2)
  lasso.r2 <- 1-length(test.traits[,1]) * lasso.mse / sum((test.traits[,j]-mean(test.traits[,j]))^2)
  lasso.coef=predict(lasso.mod ,type ="coefficients" ,s=bestlam ) #[1:41,]
  bands = names(leaf.spec)[lasso.coef@i]
  out[[nm[j]]] <- list(mod= lasso.mod, bands = bands, coeff = lasso.coef, mse = lasso.mse, r2=lasso.r2, lambda = bestlam, tr.R2 = tr.R2, tr.MSE = tr.MSE)
  
  return(out)
}

print.coeffs<-function(pls.mod.opt, optim.ncomps, rounds, lasso.reg, Y,  j, out.dir, token=NULL)
{
  tmp.pls = eval(parse(text = paste('pls.mod.opt$',names(Y)[j],sep="")))
  FF <- as.matrix(t(c(tmp.pls$validation$PRESS[optim.ncomps], tmp.pls$coefficients[,,1])))
  write.table(FF, file = paste(out.dir,"pls_coeffs_",rounds, names(Y)[j],token,".csv", sep=""), sep = ",", 
              col.names = FALSE, append=TRUE, row.names = FALSE)   
  
  tmp.las = eval(parse(text = paste('lasso.reg$',names(Y)[j],sep="")))
  FF <- as.matrix(t(c(tmp.las$r2, tmp.las$coeff@i)))
  write.table(FF, file = paste(out.dir, "las_coeffs_",rounds,names(Y)[j],token,".csv", sep=""), sep = ",", 
              col.names = FALSE, append=TRUE, row.names = FALSE)  
}




#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#' 
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}

#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
#'              the PRESS statistic.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}


extractPerm <- function(sz=1000){
  # Read data
  allBand=imp.spectra( "Spectra/neon_aop_bands.csv",in.dir)
  allData=imp.spectra("Spectra/CrownPix.csv",in.dir)   
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  
  # Set bad bands to zero
  allData[,bad_Bands]=NA
  allData=allData[, colSums(is.na(allData)) == 0]
  
  # Find unique crowns (only for plotting purposes)
  unqCrown=unique(allData$pixel_crownID)
  bootDat <- allData[1:length(unqCrown), ]
  bootDat[] <- NA
  for (laps in 1:sz){
    tk = 1
    for(i in unqCrown){
      bootDat[tk,] <- allData[sample(which(allData$pixel_crownID==i), 1),]
      tk = tk +1
    }
    write.csv(bootDat, paste('inputs/Bootstrap_2/onePix1Crown_', laps, '.csv', sep = ''))
  }
}



normAditya <- function(){
  # Read data
  allBand=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/neon_aop_bands.csv")
  allData=read.csv("/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/CrownPix.csv")
  # Find bad bands
  all_Bands=as.character(allBand$BandName)
  bad_Bands=as.character(allBand[which(allBand$noise==1),"BandName"])
  # Set bad bands to zero
  allData[,bad_Bands]=NA
  # Find unique crowns (only for plotting purposes)
  unqCrown=unique(allData$pixel_crownID)
  # Extract spectra into matrix
  specMat=as.matrix(allData[,all_Bands])
  # Vector normalize spectra
  normMat=sqrt(apply(specMat^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat=matrix(data=rep(normMat,ncol(specMat)),ncol=ncol(specMat))
  normMat=specMat/normMat
  # Write vector normalized spectra back into dataframe
  normDF=allData
  normDF[,all_Bands]=normMat
  # Write vector normalized spectra to CSV
  write.csv(normDF, "/Users/sergiomarconi/Documents/PhD/Projects/JTDM/TraitsOnHeaven/Inputs/CrownPix_norm.csv",row.names=FALSE)
}

