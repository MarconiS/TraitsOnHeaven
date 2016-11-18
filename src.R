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
  prop <- 0.8
  
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
pls.cal <- function(train.data, comps, scaling){
  spectra <- as.matrix(train.data[,7:length(train.data[1,])])
  traits <- as.matrix(train.data[,2:6])
  pls.summ <- list()
  spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
  
  for (j in 1:5) {
    leaf.trait <- traits[,j]
    #standard normal variate transform [log(first derivative)]
    train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = plsr(Y ~ X ,scale=scaling, ncomp=comps[j],validation="LOO", trace=TRUE, method = "oscorespls", data = train.PLS) 
    pls.summ[[names(train.data[,2:6])[j]]] <- tmp.pls  
  }
  return(pls.summ)
}

# Jackknife test here.  Determine optimal number of components
opt.comps <- function(train.data, k, comps, i)
{  
  dims = dim(train.data)
  i <- 30
  k <- 20
  comps <- 15
  jk.out <- array(data=NA,dim=c(i,comps+1)) # num of components plus intercept
  for (j in 1:i) {
    LeafLMA.pls = plsr(Y~X,scale=FALSE,ncomp=comps,validation="LOO",
                       trace=TRUE, method = "oscorespls", data=train.data)
    vec <- as.vector(RMSEP(LeafLMA.pls)$val)
    vec.sub.adj <- vec[seq(2,length(vec),2)]
    jk.out[i,] <- vec.sub.adj
  }
  
  # Boxplot of results
  boxplot(jk.out, xaxt="n",xlab="NUM OF COMPONENTS",ylab="RMSEP") 
  numcomps <- comps+1
  axis(1,at=1:numcomps,0:15)
  box(lwd=2.2)
  
  ttest <- t.test(jk.out[,12],jk.out[,13],alternative = "two.sided",paired=F,var.equal = FALSE)
  ttest
  
  rm(LeafLMA.pls,comps,i,iterations,vec,vec.sub.adj,ttest,segs,numcomps)
  dims <- dim(cal.plsr.data)
  k <- round(dims[1]/10)
  segs = cvsegments(dims[1],k = k, type="random")
  LeafLMA.pls = plsr(LMA_g_DW_m2~Spectra,scale=FALSE,ncomp=15,validation="CV",segments=segs,
                     trace=TRUE, method = "oscorespls", data=cal.plsr.data)
  
  # Examine raw PLSR output
  summary(LeafLMA.pls)
  plot(RMSEP(LeafLMA.pls), legendpos = "topright")
  names(LeafLMA.pls)
  
  predplot(LeafLMA.pls, ncomp = 9:11, asp = 1, line = TRUE,which = c("train","validation"),
           xlim=c(5,300),ylim=c(5,300))
}
#--------------------------------------------------------------------------------------------------#



