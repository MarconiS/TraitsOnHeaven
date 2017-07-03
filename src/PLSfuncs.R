
################### GET JACOBIAN

getJ=function(X,y,h)
	{
	X=as.matrix(X)
	#print(dim(X))
	y=as.matrix(y,ncol=1)
	#print(dim(y))
	p=ncol(X)
	S=t(X)%*%(X)
	#print(dim(S))
	s=t(X)%*%y
	#print(dim(s))
	ai=s
	atildei=as.matrix(0,nrow=p)
	Hi=matrix(0,nrow=p,ncol=p)
	
	for (i in seq(h))
		{
		if (i==1)
			{
			dahdy=t(X)
			atildei=s
			datildehdy=dahdy
			qi=(t(atildei)%*%S%*%atildei)
			Hi=(atildei%*%t(atildei))/qi[1]
			dbhdy=Hi%*%t(X)+(kronecker(t(s),atildei)%*%datildehdy+(t(s)%*%atildei)[1]*datildehdy)/qi[1]-(2*(t(s)%*%atildei)[1]*Hi%*%S%*%datildehdy)/qi[1]
			}
		else
			{
		      dahdy=dahdy-S%*%Hi%*%t(X)-S%*%((kronecker(t(s),atildei)%*%datildehdy+(t(s)%*%atildei)[1]*diag(p)%*%datildehdy)-2*(t(s)%*%atildei)[1]*Hi%*%S%*%datildehdy)/qi[1]
		      ai=ai-S%*%Hi%*%s
		      datildehdy=(diag(p)-Hi%*%S)%*%dahdy-((kronecker(t(ai)%*%S,atildei)%*%datildehdy+(t(ai)%*%S%*%atildei)[1]*diag(p)%*%datildehdy)-2*(t(ai)%*%S%*%atildei)[1]*Hi%*%S%*%datildehdy)/qi[1]
		      atildei=(diag(p)-Hi%*%S)%*%ai
		      qi=(t(atildei)%*%S%*%atildei)
		      Hi=atildei%*%t(atildei)/qi[1]
		      dbhdy=dbhdy+Hi%*%t(X)+((kronecker(t(s),atildei)%*%datildehdy+(t(s)%*%atildei)[1]*datildehdy)-2*(t(s)%*%atildei)[1]*Hi%*%S%*%datildehdy)/qi[1]
			}
		}

	return(dbhdy)
	}

################### GET ANALYTICAL CIs

getACI=function(allData,inVar,iLst,predY)
	{
	N=nrow(allData)
	Y=as.matrix(allData[,inVar],ncol=1)
	X=as.matrix(allData[,iLst])
	X=scale(X,center=T)
	J=getJ(X,Y,nComp)
	r=Y-predY

	rPrimeZero=(diag(N)-X%*%J) 				# Denham (1997) Eqn 18, Serneels (2004) Eqn. 10
	df=sum(diag(t(rPrimeZero)%*%(rPrimeZero))) 	# Serneels (2004) Eqn. 10
	#rZero=scale(r,center=T)
	rZero=r-rPrimeZero%*%(Y-scale(Y,center=T))
	rTilde=r-rZero
	sigma=sqrt((t(rTilde)%*%rTilde)/df)[1] 		# Denham (1997) Eqn. 18
	cInt=c()
	for (x in seq(nrow(allData)))
		{
		newVec=(X[x,])
		int95=dt(0.95,df)*sigma*sqrt(((N+1)/N)+t(newVec)%*%J%*%t(J)%*%newVec)
		cInt=c(cInt,int95)
		}
	return(list(cInt,N,df,sigma,J))
	}

################### GET VIP

VIPjh <- function(object, j, h) 
      { 
      if (object$method != "oscorespls") stop("Only implemented for orthogonal scores algorithm. Refit with 'method = \"oscorespls\"'")

      if (nrow(object$Yloadings) > 1) stop("Only implemented for single-response models")

      b <- c(object$Yloadings)[1:h]
      T <- object$scores[,1:h, drop = FALSE]
      SS <- b^2 * colSums(T^2)
      W <- object$loading.weights[,1:h, drop = FALSE]
      Wnorm2 <- colSums(W^2) 
	sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS)) 
      }


