#' Confidence Region for Mean Rotation
#'
#' Find the radius of a 100(1-a)% confidence region for the projected mean
#'
#' @param Qs A n-by-4 matrix where each row corresponds to a random rotation in matrix form
#' @param method Character string specifying which type of interval is required
#' @param alpha The alhpa level desired
#' @return radius of the confidence region centered at the projected mean
#' @cite rancourt2000
#' @S3method region Q4
#' @method region Q4
#' @examples
#' Qs<-ruars(20,rcayley,space='Q4',kappa=100)
#' region(Qs,method='rancourt',alpha=0.9)

region.Q4<-function(Qs,method,alpha){
	
	Qs<-formatQ4(Qs)
	
	if(method%in%c('Rancourt','rancourt')){
		
		r<-rancourtCR.Q4(Qs=Qs,a=alpha)
		return(r)
		
	}else{
		stop("Only the Rancourt option is available now.")
	}
	
}

rancourtCR.SO3<-function(Rs,a){
	Qs<-Q4(Rs)
	r<-rancourtCR.Q4(Qs,a)
	return(r)
}

rancourtCR.Q4<-function(Qs,a){
	#This takes a sample qs and returns the radius of the confidence region
	#centered at the projected mean
	n<-nrow(Qs)
	Shat<-mean(Qs)
	Phat<-pMat(Shat)
	
	Rhat<-Qs%*%Phat
	resids<-matrix(0,n,3)
	VarShat<-matrix(0,3,3)
	
	resids<-2*Rhat[,1]*matrix(Rhat[,2:4],n,3)
	
	VarShat<-t(resids)%*%resids/(n-1)
	
	RtR<-t(Rhat)%*%Rhat
	Ahat<-(diag(RtR[1,1],3,3)-RtR[-1,-1])/n
	
	avec<-matrix(axis2(Shat),1,3)
	
	Tm<-n*avec%*%Ahat%*%solve(VarShat)%*%Ahat%*%t(avec)
	r<-sqrt(qchisq(alpha,3)/Tm)
	return(r)
}



cdfuns<-function(Rs,Shat){
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	c<-0
	d<-0
	
	StO<-centeringSO3(Rs,Shat)
	
	for(i in 1:n){
		StOi<-matrix(StO[i,],3,3)
		c<-c+(3-sum(diag(StOi%*%StOi)))
		d<-d+sum(diag(StOi))
	}
	
	c<-c/(6*n)
	d<-d/(3*n)
	return(list(c=c,d=d))
}

ZhangCI<-function(Rs,m,alpha){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#m is the number of resamples to find q_1-alpha
	#alpha is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	Shat<-mean(Rs)
	
	tstar<-rep(0,m)
	tstarPivot<-rep(0,m)
	
	for(i in 1:m){
		Ostar<-as.SO3(Rs[sample(1:n,replace=T),])
		ShatStar<-mean(Ostar)
		tstar[i]<-3-sum(diag(t(Shat)%*%ShatStar))
		
		cd<-cdfuns(Ostar,ShatStar)
		
		tstarPivot[i]<-tstar[i]*2*n*cd$d^2/cd$c
		
	}
	
	return(list(ts=as.numeric(quantile(tstar,alpha)),tsPivot=as.numeric(quantile(tstarPivot,alpha))))
}


fisherAxis<-function(Qs,Shat){
	
	n<-nrow(Qs)
	svdQs<-svd(t(Qs)%*%Qs/n)
	mhat<-svdQs$v[,1]
	Mhat<-t(svdQs$v[,-1])
	etad<-svdQs$d[1]
	etas<-svdQs$d[-1]
	G<-matrix(0,3,3)
	
	for(j in 1:3){
		for(k in j:3){
			denom<-1/(n*(etad-etas[j])*(etad-etas[k]))
			
			for(i in 1:n){
				G[j,k]<-G[k,j]<-G[j,k]+(Mhat[j,]%*%Qs[i,])*(Mhat[k,]%*%Qs[i,])*(mhat%*%Qs[i,])^2*denom
			}
		}
	}
	
	Ginv<-try(solve(G),silent=T)
	
	if(class(Ginv)!='matrix'){
		Ginv<-diag(1/diag(G))
		print('Diagonal used')
	}
	
	Tm<-n*Shat%*%t(Mhat)%*%Ginv%*%Mhat%*%t(Shat)
	return(Tm)
}