#' Confidence Region for Mean Rotation
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the projected mean
#'
#' @param Qs A \eqn{n\times 4}{n-by-4} matrix where each row corresponds to a random rotation in matrix form
#' @param method Character string specifying which type of interval is required
#' @param alpha The alpha level desired
#' @param ... Additional arguments
#' @return radius of the confidence region centered at the projected mean
#' @cite rancourt2000
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Qs,method='rancourt',alpha=0.9)

region<-function(Qs,method,alpha,...){
	UseMethod("region")
}


#' @rdname region
#' @method region Q4
#' @S3method region Q4

region.Q4<-function(Qs,method,alpha,...){
	
	Qs<-formatQ4(Qs)
	
	if(method%in%c('Rancourt','rancourt')){
		
		r<-rancourtCR.Q4(Qs=Qs,a=alpha)
		
		return(r)
		
	}else	if(method%in%c('Zhang','zhang')){
		
		r<-zhangCR.Q4(Qs=Qs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Fisher','fisher')){
		
		r<-fisherCR.Q4(Qs=Qs,a=alpha)
		
		return(r)
		
	}else{
		
		stop("Only the Rancourt, Zhang and Fisher options are currently available")
		
	}
	
}


#' @rdname region
#' @method region SO3
#' @S3method region SO3

region.SO3<-function(Rs,method,alpha,...){
	
	Rs<-formatSO3(Rs)
	
	if(method%in%c('Rancourt','rancourt')){
		
		r<-rancourtCR.SO3(Rs=Rs,a=alpha)
		return(r)
		
	}else	if(method%in%c('Zhang','zhang')){
		
		r<-zhangCR.SO3(Rs=Rs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Fisher','fisher')){
		
		r<-fisherCR.SO3(Rs=Rs,a=alpha)
		
		return(r)
		
	}else{
		
		stop("Only the Rancourt, Zhang and Fisher options are currently available")
		
	}
	
}

#' Rancourt CR Method
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the projected mean \cite{rancourt2000}
#'
#' This works in the same way as done in \cite{bingham09} which assumes rotational 
#' symmetry and is therefore conservative.
#'
#' @param Qs A \eqn{n\times 4}{n-by-4} matrix where each row corresponds to a random rotation in matrix form
#' @param a The alpha level desired
#' @return radius of the confidence region centered at the projected mean
#' @cite rancourt2000
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='rancourt',alpha=0.9)

rancourtCR<-function(Qs,a){
	UseMethod("rancourtCR")
}


#' @rdname rancourtCR
#' @method rancourtCR Q4
#' @S3method rancourtCR Q4

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
	
	Tm<-min(diag(n*Ahat%*%solve(VarShat)%*%Ahat))
	
	r<-sqrt(qchisq(a,3)/Tm)
	return(r)
}


#' @rdname rancourtCR
#' @method rancourtCR SO3
#' @S3method rancourtCR SO3


rancourtCR.SO3<-function(Rs,a){
	Qs<-Q4(Rs)
	r<-rancourtCR.Q4(Qs,a)
	return(r)
}

#' Zhang CR Method
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the projected mean
#'
#' @param Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix or quaternion form
#' @param a The alpha level desired
#' @param m Number of replicates to use to estiamte cut point
#' @param pivot should the pivotal (T) or non-pivotal (F) method be used
#' @param estimator Mean or median
#' @return radius of the confidence region centered at the projected mean
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='zhang',alpha=0.9)

zhangCR<-function(Qs,a,m,pivot,estimator){
	UseMethod("zhangCR")
}


#' @rdname zhangCR
#' @method zhangCR SO3
#' @S3method zhangCR SO3

zhangCR.SO3<-function(Rs,a,m=300,pivot=T,estimator='mean'){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#m is the number of resamples to find q_1-a
	#a is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
  
  if(estimator=='mean'){
	  Shat<-mean(Rs)
	}else{
    Shat<-median(Rs,type='intrinsic')
  }

	tstar<-rep(0,m)
	
	if(pivot){
		
		tstarPivot<-rep(0,m)
	
		for(i in 1:m){
			
			Ostar<-as.SO3(Rs[sample(n,replace=T),])
      
			if(estimator=='mean'){
			  ShatStar<-mean(Ostar)
			}else{
			  ShatStar<-median(Ostar,type='intrinsic')
			}
      
			tstar[i]<-3-sum(diag(t(Shat)%*%ShatStar))
		
			cd<-cdfuns(Ostar,ShatStar)
		
			tstarPivot[i]<-tstar[i]*2*n*cd$d^2/cd$c
		
		}
		
		cdhat<-cdfuns(Rs,Shat)
		
		qhat<-as.numeric(quantile(tstarPivot,a))*cdhat$c/(2*n*cdhat$d)
		
		return(acos(1-qhat/2))
		
	}else{
		
		for(i in 1:m){
			
			Ostar<-as.SO3(Rs[sample(n,replace=T),])
			
			if(estimator=='mean'){
			  ShatStar<-mean(Ostar)
			}else{
			  ShatStar<-median(Ostar,type='intrinsic')
			}
			
			tstar[i]<-3-sum(diag(t(Shat)%*%ShatStar))
			
		}
		
		qhat<-as.numeric(quantile(tstar,a))
		
		return(acos(1-qhat/2))
	}
	
}

#' @rdname zhangCR
#' @method zhangCR Q4
#' @S3method zhangCR Q4

zhangCR.Q4<-function(Qs,alpha,m=300,pivot=T,estimator='mean'){
	
	Rs<-SO3(Qs)
	
	r<-zhangCR.SO3(Rs,alpha,m,pivot,estimator)
	
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




#' Fisher Mean Polax Axis CR Method
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the projected mean \cite{fisher1996}
#'
#' This works in the same way as done in \cite{bingham09} which assumes rotational 
#' symmetry and is therefore conservative.
#'
#' @param Qs A \eqn{n\times 4}{n-by-4} matrix where each row corresponds to a random rotation in matrix form
#' @param a The alpha level desired
#' @param boot Should the bootstrap or normal theory critical value be used
#' @param m number of bootstrap replicates to use to estimate critical value
#' @return radius of the confidence region centered at the projected mean
#' @cite fisher1996
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='fisher',alpha=0.9)

fisherCR<-function(Qs,a,boot,m){
	UseMethod("fisherCR")
}


#' @rdname fisherCR
#' @method fisherCR Q4
#' @S3method fisherCR Q4

fisherCR.Q4<-function(Qs,a,boot=T,m=300){
	
	Qs<-formatQ4(Qs)
	n<-nrow(Qs)
	mhat<-mean(Qs)
	
	if(boot){
		Tstats<-rep(0,m)
	
		for(i in 1:m){
			Qsi<-as.Q4(Qs[sample(n,replace=T),])
			Tstats[i]<-fisherAxis(Qsi,mhat)
		}
	
		qhat<-as.numeric(quantile(Tstats,a))
		
	}else{
		
		qhat<-qchisq(a,3)
		
	}
	
	rsym<-optim(.05,optimAxis,Qs=Qs,cut=qhat,method='Brent',lower=0,upper=pi)$par
	
	return(rsym)
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
	
	Tm<-n*Shat%*%t(Mhat)%*%solve(G)%*%Mhat%*%t(Shat)
	return(Tm)
}

optimAxis<-function(r,Qs,cut){
	
	Shat<-Q4(axis2(mean(Qs)),r)
	
	Tm<-fisherAxis(Qs,Shat)
	
	return((Tm-cut)^2)
}


#' @rdname fisherCR
#' @method fisherCR SO3
#' @S3method fisherCR SO3

fisherCR.SO3<-function(Rs,a,boot=T,m=300){
	
	Qs<-Q4(Rs)
	r<-fisherCR.Q4(Qs,a,boot,m)
	
	return(r)
}