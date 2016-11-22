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
  waves <- data.frame(seq(1,2151,1),seq(350,2500,1))
  mean.spec <- colMeans(X)
  spec.quant <- apply(X,2,quantile,probs=c(0.05,0.95))
  # Output correlation data
  pdf(paste(pwd, '/','FFT_Spectra_Correlations.pdf',sep=""),height=12,width=8)
  par(mfrow=c(6,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
  matplot(mean.spec, type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n",ylim=c(0,0.9))
  ind <- pretty(seq(from = 350, to = 2500, by = 1)) # Using pretty to standardize the axis
  ind <- ind[ind >= 350 & ind <= 2500]
  ind <- (ind - 349) / 1
  axis(1, ind, colnames(X)[ind]) # add column names to wavelengths
  # CIs
  lines(spec.quant[1,],lty=1,col="dark grey")
  lines(spec.quant[2,],lty=1,col="dark grey")
  legend("topleft",legend=c("Mean","95% CI"),col=c("black","dark grey"),lwd=3)
  
  for (j in 1:length(spec_corr)) {
    plot(waves[,2],spec_corr[,j],xlab="WAVELENGTH (nm)",ylab="CORRELATION", main=names(spec_corr)[j], cex=0.01)
    #lines(waves[,2],spec_corr[,j],lwd=4)
    abline(h=0,lty=2,lwd=1.5,col="grey80")
    #box(lwd=2)
  }
  dev.off()
}


cut.set<-function(aug.X,out.dir){
  sites <- unique(aug.X$aug.spectra.Site)
  # Sample proportion for cal data
  prop <- 0.7
  
  # Random seed
  create.seed <- FALSE  #TRUE/FALSE
  if (create.seed){
    set.seed(as.vector(round(runif(5,min=5,max=9))))
    ### Write out seed
    .Random.seed[1:6]
    seed.save <- .Random.seed
    write.table(seed.save,paste(out.dir,"random.seed",sep=""));
  }
  ### Read in previous random seed
  seed <- read.table(paste(out.dir,"random.seed",sep=""))[,1];
  .Random.seed <- seed
  
  train.data <- 0
  test.data <- 0
  j <- 1
  for (i in sites){
    temp.data <- aug.X[which(aug.X$aug.spectra.Site==i),]
    rows <- sample(1:nrow(temp.data),floor(prop*nrow(temp.data)))
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
  return(list(train=train.data, test=test.data))
}


# Calibration PLS model ---------------------------------------------------
pls.cal <- function(train.data, comps, scaling = FALSE, pl = FALSE){
  spectra <- as.matrix(train.data[,7:length(train.data[1,])])
  traits <- as.matrix(train.data[,2:6])
  pls.summ <- list()
  spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  
  for (j in 1:5) {
    leaf.trait <- traits[,j]
    #standard normal variate transform [log(first derivative)]
    train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = plsr(Y ~ X ,scale=scaling, ncomp=comps[j],validation="LOO", trace=TRUE, method = "oscorespls", data = train.PLS) 
    if (pl) {
      predplot(tmp.pls, ncomp = 8:14, asp = 1, line = TRUE,which = c("train","validation"),
               xlim=c(5,300),ylim=c(5,300))
    }
    pls.summ[[names(train.data[,2:6])[j]]] <- tmp.pls 
  }
  return(pls.summ)
}

#--------------------------------------------------------------------------------------------------#
opt.comps <- function(pls.mod.train, test.data, plots.ok = FALSE, use.press = TRUE)
{
  ncomps = rep(0,5)
  spectra <- as.matrix(test.data[,7:length(test.data[1,])])
  traits <- as.matrix(test.data[,2:6])
  spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))

  for (j in 1:5) {
    leaf.trait <- traits[,j]
    test.PLS <- data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep="")))
    #use RMSE as discrimant to determine components optimization
    if(!use.press){
      rms <- RMSEP(eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep=""))), newdata = test.PLS)$val
      if(plots.ok){    
        plot(RMSEP(eval(parse(text = paste('pls$',names(Y)[j],sep=""))),estimate=c("test"),newdata = test.PLS), main="MODEL RMSEP",
             xlab="NUM OF COMPONENTS",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
        
        r2 <- R2(eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep=""))), newdata = test.PLS)
        plot(R2(eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep=""))),estimate=c("test"),newdata = test.PLS), main="MODEL R2",
             xlab="NUM OF COMPONENTS",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
      }
    }else{
      #--------------------------------------------------------------------------------------------------#
      # Calculate Q2 statistic
      q2<- Q2(test.PLS, eval(parse(text = paste('pls.mod.train$',names(Y)[j],sep="")))) 
    }
    if(!use.press){
      rms= rms[seq(2,length(rms),2)]
      ncomps[j] <- which(rms==min(rms))
    }else{
      #q2 = q2[seq(2,length(q2),2)]
      ncomps[j] <- which(q2==min(q2))
    }
  }
  ncomps <- pmin(0.8 * length(test.data[,1]), ncomps)
  return(ncomps)
}

Q2 <- function(test.PLS, tmp.pls)
{
  dims <- dim(test.PLS)
  PRESS = tmp.pls$validation$PRESS
  SS = sum((tmp.pls$Y)^2)
  TSS = sum((test.PLS$Y-mean(test.PLS$Y))^2)
  Q2 = 1-(PRESS/TSS)  # using SS not TSS
  
  # Calculate RMSECV
  RMSECV = PRESS/dims[1]
  return(RMSECV)
}