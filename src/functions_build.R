#-------imp.spectra-------------------------------------------------------------------------------------------

imp.spectra <- function(f, pwd)
{
  # Import dry spectra dataset
  aug.data <- read.table(paste(pwd,f,sep=""), header=TRUE,sep=",")
  return(aug.data)
}
#-------cut.set-------------------------------------------------------------------------------------------

cut.set<-function(aug.X,out.dir, c.id, prop = 0.7){
  species <- unique(aug.X$aug.spectra.name)
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

# pls.cal ---------------------------------------------------
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

#-------opt.comps-------------------------------------------------------------------------------------------
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
#-----predict.pls------------------------------------------------------------------------------------------
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
#-----res.out------------------------------------------------------------------------------------------
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