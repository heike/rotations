

RivestCI<-function(qs,S=id.Q4){
	#This takes as input the dataset and true central direction S
	n<-nrow(qs)
	Shat<-mean(qs)
	Phat<-pMat(Shat)
	
	Rhat<-qs%*%Phat
	resids<-matrix(0,n,3)
	VarShat<-matrix(0,3,3)
	
	resids<-2*Rhat[,1]*matrix(Rhat[,2:4],n,3)
	
	VarShat<-t(resids)%*%resids/(n-1)
	
	RtR<-t(Rhat)%*%Rhat
	Ahat<-(diag(RtR[1,1],3,3)-RtR[-1,-1])/n
	
	St<-as.Q4(matrix(c(S[1],-S[2:4]),1,4))
	StShat<-qMult(St,Shat)
	avec<-matrix(axis2(StShat)*angle(StShat),1,3)
	
	Tm<-n*avec%*%Ahat%*%solve(VarShat)%*%Ahat%*%t(avec)
	
	return(Tm)
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