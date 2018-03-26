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
    rows <- sample(1:nrow(temp.data),ceiling(prop*nrow(temp.data)))
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
pls.cal <- function(train.data, numcomps = 15, nm,j,  normalz = F, lab.train = NULL){

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
    tmp.y<-matrix(as.numeric(lab.train$id,ncol=1)) # make numeric matrix
    train.PLS = data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    #tmp.pls<-plsda(y = factor(tmp.y), x = spectra_log_dif_snv,  ncomp = comps,  probMethod = "softmax", type = "class")
    tmp.pls<- plsr(Y ~ X, ncomp=numcomps,validation="LOO", probMethod = "softmax", trace=TRUE, method = "oscorespls", data = train.PLS, type = "class")
  }else{
    train.PLS = data.frame(Y = I(leaf.trait), X=I(spectra_log_dif_snv))
    tmp.pls = plsr(Y ~ X,scale=F, ncomp=numcomps,validation="LOO", trace=TRUE, method = "oscorespls", probMethod = "softmax", data = train.PLS)
  }
  pls.summ[[nm[j]]] <- tmp.pls

  return(pls.summ)
}

#-------opt.comps-------------------------------------------------------------------------------------------
opt.comps <- function(pls.mod.train, Y, j, Class = F){
  ncomps <- NA
  tmp.pls = pls.mod.train[[1]]
  if(!Class){
    ncomps <- which(tmp.pls$validation$PRESS==min(tmp.pls$validation$PRESS[3:length(tmp.pls$validation$PRESS)]))
  }
  return(c(ncomps))
}
#-----predict.pls------------------------------------------------------------------------------------------
predict.pls <- function(pls.mod.train, test.data, optim.ncomps=NULL,j, nm, norm = F, lab.test = NULL){
  spectra <- as.matrix(test.data[grepl("band", names(test.data))])
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  out <- list()
  spectra_log_dif_snv <- spectra
  if(norm){  spectra_log_dif_snv <-t(diff(t(log(spectra)),differences=1, lag=3))}
  leaf.trait <- traits[,j]
  if(is.character(leaf.trait)){
    spectra_log_dif_snv=spectra_log_dif_snv[, colSums(is.na(spectra_log_dif_snv)) == 0]

    tmp.y<-matrix(as.numeric(lab.test$id, ncol=1)) # make numeric matrix
    test.PLS <- data.frame(Y = I(tmp.y), X=I(spectra_log_dif_snv))
    pred<-predict(pls.mod.train$name,newdata=test.PLS, comps = optim.ncomps)
    accu = (sum(round(pred)==lab.test$id)/length(lab.test$id))
    print(accu)
    #out.data = data.frame(lab.test$name,pred.val.data$,rep(NA, length(lab.test$name)), -9999, -9999)
    temp = data.frame(traits[,j],round(pred), rep(accu,length(traits[,j])))
    names(temp) = c("obs","pred", "R2")

    out[[nm[j]]] <- temp

  }else{
    test.PLS <- data.frame(X=I(spectra_log_dif_snv))
    out[[nm[j]]] <- as.vector(predict(pls.mod.train[[1]], newdata = test.PLS, ncomp=optim.ncomps, type="response"))
  }
  return(out)
}
#-----res.out------------------------------------------------------------------------------------------
res.out <- function(pred.val.data, train.data, test.data, j, nm, labelsSp = NULL)
{
  out <- list()
  traits <- as.matrix(test.data[,names(test.data) %in% nm])
  # Build output dataset
  # PLSR Summary statistics
  if(is.character(traits[,j])){
    pred.data <- round(as.vector(pred.val.data))
    correct <- labelsSp$id == as.numeric(pred.data)
    accuracy <- sum(as.numeric(correct))/length(correct)
    temp = data.frame(as.numeric(factor(traits[,j])),as.numeric(pred.data),correct, accuracy)
    names(temp) = c("obs","pred","res", "R2")
  }else{
    pred.data <- as.vector(eval(parse(text = paste('pred.val.data$',nm[j],sep=""))))
    res <- traits[,j]- pred.data
    MSE.test <- mean(res^2)
    RMSE.test <- sqrt(MSE.test)
    R2 <- 1- length(pred.data) * MSE.test / sum((traits[,j]- mean(traits[,j]))^2)
    ### Output val dataset
    temp = data.frame(traits[,j],pred.data,res, R2, MSE.test)
    names(temp) = c("obs","pred","res", "R2", "MSE")
  }

  out[[names(train.data[,names(train.data) %in% names])[j]]] <- temp

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

#-----pls_glm--------------------------------------------------------------------------------------$

pls_glm <- function(j = "P_pct", loops = NULL,cores = 2, nrmlz = F, out.dir = NULL, in.dir = NULL, nsites = 2){
  library(foreach)
  library(doParallel)

  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())

  ll = as.character(1:loops)
  mod.out <- foreach(laps = ll, .verbose = T, .multicombine = TRUE) %dopar% {

    library('plsRglm')
    library(readr)
    library(gnm)

    source(paste('../../src/functions_build.R', sep=""))
    source(paste('../../src/src_build.R', sep=""))
    out.dir = paste(getwd(), "/outputs/", sep="")
    in.dir = paste(getwd(), "/inputs/", sep="")

    cr.Traits <- read.csv(paste(in.dir, "Spectra/trainCrownTraits.csv",sep=""), stringsAsFactors = F)
    nCrowns <- dim(cr.Traits)[1]
    aug.spectra <- read.csv(paste(in.dir, 'Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    aug.spectra$X= NULL
    aug.spectra <- merge(cr.Traits, aug.spectra, by.x = "pixel_crownID", by.y = "pixel_crownID")
    X<- aug.spectra[grepl("band", names(aug.spectra))]
    X=X[, colSums(is.na(X)) == 0]
    Y <- as.data.frame(aug.spectra[,names(aug.spectra) %in% j])
    Y <- (round(Y,3))
    if(dim(Y)[2]==1){names(Y) = j}
    aug.X <- data.frame(aug.spectra$name, Y, X)

    # Subset data into cal/val by site
    eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
    train.data <- eval.set$train
    test.data <- eval.set$test
    colnames(train.data)[1] <- "name"
    colnames(test.data)[1] <- "name"
    # write_csv(train.data, paste(in.dir, "datSets/train_for_baseline_", laps,".csv", sep=""))
    # write_csv(test.data, paste(in.dir, "datSets/test_for_baseline_", laps,".csv", sep=""))

    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    Y <- as.vector(train.data[,names(train.data) %in% j])
    #K=nrow(Y)
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }

    train.PLS<- cv.plsRglm(dataY=log(Y),dataX=X, nt=15,NK=1, K=10,modele="pls-glm-family",family=gaussian())
    out <- list()
    #out["summary"] <- summary(train.PLS)
    press <- unlist(kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))

    #chosen to be built
    mod <- plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]),modele="pls-glm-gaussian")
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }

    Y.test <- as.vector(test.data[,names(test.data) %in% j])
    out[["pred"]] <- predict(mod,newdata=X.tst,type="response",comps=as.integer(out["ncomp"]), se.fit = T)
    out[["mod"]] <- mod
    out[["score"]] <- gnm(out$pred$fit ~ (Y.test), family = gaussian())
    out["aic"]<- AIC(out$score)
    #r2 <- summary(lm(pred$fit~ Y.test))
    out
  }
  #print(paste(j, "done"))
  #save(mod.out, file = paste("pls_glm_", j,  sep = ""))
  stopCluster(cl)
  return(mod.out)
}

predict.withsd <- function(object,newdata,comps=object$computed_nt,type=c("link", "response", "terms", "scores", "class", "probs"),se.fit=FALSE, wt = NULL, dispersion = NULL,methodNA="adaptative",verbose=TRUE,...)
{
  if (!inherits(object, "plsRglmmodel"))
    stop("Primary argument much be a plsRglmmodel object")
  if(missing(type)){type="link"}
  if (!(type %in% c("link",  "response", "terms", "scores", "class", "probs")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="link"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="response"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="terms"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="class"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}
    }
    if (type=="probs"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}
    }
    if (type=="scores"){return(object$tt[,1:comps])}
  } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$call$formula)){
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("subset", "weights"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$data <- newdata
      mf$formula <- object$call$formula[-2]
      mf$drop.unused.levels <- TRUE
      mf$na.action <- na.pass
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      mt <- attr(mf, "terms")
      #    attr(mt,"intercept")<-0L
      newdata.frame <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)[,-1]
      else matrix(, nrnd, 0L)
      weights <- as.vector(model.weights(mf))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
    newdataNA <- !is.na(newdata)
    newdata.scaledwotNA <- as.matrix(newdata.scaled)
    newdata.scaledwotNA[!newdataNA] <- 0
    ttpredY <- NULL
    if (methodNA=="adaptative") {
      for(ii in 1:nrnd){
        if (all(newdataNA[ii,])){
          ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps)))
        }
        else {
          if(verbose){cat("Missing value in row ",ii,".\n")}
          ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
        }}}
    if(methodNA=="missingdata") {
      if(verbose){cat("Prediction as if missing values in every row.\n")}
      for (ii in 1:nrnd) {
        ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
      }
    }
    colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY)
    if("glm" %in% class(object$FinalModel)){
      if (type=="link"){return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit,dispersion = dispersion,...))
      }
      #if (type=="response"){return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit,dispersion = dispersion,...))
      #}
      if (type=="terms"){return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit,dispersion = dispersion,...))
      }
      if (type=="response"){return(predict.lm(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit,interval = 'prediction',weights = wt, dispersion = dispersion,...))
      }
    }
    if("polr" %in% class(object$FinalModel)){
      if (type=="class"){return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}
      if (type=="probs"){return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}
    }
    if (type=="scores"){colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="");
    return(ttpredY[,1:comps,drop=FALSE])}
  }
}
