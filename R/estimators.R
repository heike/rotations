#' Mean Rotation
#'
#' Compute the intrinsic or projected mean of a sample of rotations
#'
#' This function takes a sample of \eqn{3\times 3}{3-by-3} rotations (in the form of a \eqn{n\times 9}{n-by-9} matrix where \eqn{n>1} is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})}{argmin d^2(bar(R),S)} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum Ri/n} and the distance metric \eqn{d_D}{d}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param Rs A \eqn{n\times 9}{n-by-9} matrix where each row corresponds to a random rotation in matrix form
#' @param type String indicating 'projected' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @param ... additional arguments passed to mean
#' @return Estimate of the projected or intrinsic mean of the sample
#' @seealso \code{\link{median.SO3}}
#' @cite moakher02, manton04
#' @S3method mean SO3
#' @method mean SO3
#' @examples
#' Rs<-ruars(20,rvmises,kappa=0.01)
#' mean(Rs)

mean.SO3 <- function(Rs, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(Rs)	
	
	if(nrow(Rs)==1)
		return(Rs)
	
  if (!(type %in% c("projected", "intrinsic")))
    stop("type needs to be one of 'projected' or 'intrinsic'.")
  
  R <- project.SO3(matrix(colMeans(Rs), 3, 3))
  
  if (type == "intrinsic") {
    n <- nrow(Rs)
    d <- 1
    iter <- 0
    s <- matrix(0, 3, 3)
    
    while (d >= epsilon) {
      
      R <- R %*% exp.skew(s)
      
      s <- matrix(colMeans(t(apply(Rs, 1, tLogMat, S = R))), 3, 3)
      
      d <- norm(s, type = "F")
      
      iter <- iter + 1
      
      if (iter >= maxIter) {
        warning(paste("No convergence in ", iter, " iterations."))
        return(as.SO3(R))
      }
    }
    
  }
  class(R)<-"SO3"
  return(R)
}

#' Mean Rotation
#' 
#' Compute the projected or intrinsic mean of a sample of rotations
#'
#' This function takes a sample of \eqn{n} unit quaternions and approximates the mean rotation.  If the projected mean
#' is called for then the according to \cite{tyler1981} an estimate of the mean is the eigenvector corresponding to the largest
#' eigen value of \eqn{\frac{1}{n}\sum_{i=1}^nq_i^\top q_i}{Q`Q/n}.  If the intrinsic
#' mean is called then the quaternions are transformed into \eqn{3\times 3}{3-by-3} matrices and the \code{mean.SO3}
#' function is called.
#'
#' 
#' @param Qs A \eqn{n\times 4}{n-by-4} matrix where each row corresponds to a random rotation in unit quaternion
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{mean.SO3}}
#' @cite moakher02, manton04
#' @S3method mean Q4
#' @method mean Q4
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,space="Q4")
#' mean(Qs,type='intrinsic')

mean.Q4 <- function(Qs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
	
	Qs<-formatQ4(Qs)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	if(type=='projected'){
		
		QtQ<-t(Qs)%*%Qs
		R<-formatQ4(svd(QtQ)$u[,1])
		if(sign(R[1])==-1)
			R<--R
		
	}else{
		
		Rs<-SO3(Qs)
  	R<-mean(Rs,type,epsilon,maxIter)
		R<-Q4.SO3(R)
	}
  return(R)
  
}


#' Median Rotation
#' 
#' Compute the projected or intrinsic median of a sample of rotations
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}{argmin\sum d(Ri,S)}.  If the choice of distance metrid, \eqn{d_D}{d}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#'
#' @param x A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion form (\eqn{p=4})
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @param ... additional arguments
#' @return an estimate of the projected or intrinsic mean
#' @seealso \code{\link{mean.SO3}}
#' @cite hartley11
#' @export

median<-function(x,...){
  UseMethod("median")
}

#' @rdname median
#' @method median SO3
#' @S3method median SO3

median.SO3 <- function(Rs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	
	if(nrow(Rs)==1)
		return(Rs)
	
  stopifnot(type %in% c("projected", "intrinsic"))
  
  S <- mean.SO3(Rs)
  d <- 1
  iter <- 1
  delta <- matrix(0, 3, 3)
  
  while (d >= epsilon) {
    
    if (type == "projected") {
    	
    	cRs<-Rs-matrix(as.vector(S),n,9,byrow=T)
    	vn<-sqrt(rowSums(cRs^2))
      
      delta <- matrix(colSums(Rs/vn)/sum(1/vn), 3, 3)
      
      Snew <- project.SO3(delta)
      
      d <- norm(Snew - S, type = "F")
      S <- Snew
      
    } else if (type == "intrinsic") {
      
      S <- exp.skew(delta) %*% S 
      
      v <- apply(Rs, 1, tLogMat2, S = S)
      vn <- sqrt(colSums(v^2))
      vn <- pmax(epsilon, vn) # make sure we don't dividde by zero
 #     if (iter ==25) browser()
      delta <- matrix(colSums(v/vn)/sum(1/vn), 3, 3)
      d <- norm(delta, type = "F")
#      print(d)
#      print(delta)
      
    }
    
    iter <- iter + 1
    
    if (iter >= maxIter) {
      warning(paste("Unique solution wasn't found after ", iter, " iterations."))
      
      return(as.SO3(S))
    }
  }
  class(S)<-"SO3"
  return(S)
}


#' @rdname median
#' @method median Q4
#' @S3method median Q4

median.Q4 <- function(Qs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
	
	Qs<-formatQ4(Qs)
	
	if(length(Qs)==4)
		return(Qs)

  Rs<-SO3(Qs)
  
  R<-median.SO3(Rs,type,epsilon,maxIter)
  
  return(Q4.SO3(R))
}


#' Weighted Mean Rotation
#'
#' Compute the weighted intrinsic or projected mean of a sample of rotations
#'
#' This function takes a sample of \eqn{3\times 3}{3-by-3} rotations (in the form of a \eqn{n\times 9}{n-by-9} matrix where \eqn{n>1} is the sample size) and returns the weighted projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})}{argmin d(bar(R),S)} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum R_i/n} and the distance metric \eqn{d_D}{d}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param Rs A \eqn{n\times 9}{n-by-9} matrix where each row corresponds to a random rotation in matrix form
#' @param w a numerical vector of weights the same length as the number of rows in Rs giving the weights to use for elements of Rs
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @param ... only used for consistency with mean.default
#' @return weighted projected mean of the sample
#' @seealso \code{\link{median.SO3}} \code{\link{mean.SO3}}
#' @cite moakher02
#' @S3method weighted.mean SO3
#' @method weighted.mean SO3
#' @examples
#' Rs<-ruars(20,rvmises,kappa=0.01)
#' wt<-abs(1/angle(Rs))
#' weighted.mean(Rs,wt)

weighted.mean.SO3 <- function(Rs, w, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(Rs)
	
	if(nrow(Rs)==1)
		return(Rs)
	
	if(length(w)!=nrow(Rs))
		stop("'Rs' and 'w' must have same length")
	
	if (!(type %in% c("projected", "intrinsic")))
		stop("type needs to be one of 'projected' or 'intrinsic'.")
	
	if(any(w<0))
		warning("Negative weights were given.  Their absolute value is used.")
	
	w<-abs(w/sum(w))
	
	wRs<-w*Rs
	
	R <- as.SO3(project.SO3(matrix(colSums(wRs), 3, 3)))
	
	if (type == "intrinsic") {
		n <- nrow(Rs)
		d <- 1
		iter <- 0
		s <- matrix(0, 3, 3)
		
		while (d >= epsilon) {
			
			R <- R %*% exp.skew(s)
			
			s <- matrix(colSums(w*t(apply(Rs, 1, tLogMat, S = R))), 3, 3)
			
			d <- norm(s, type = "F")
			
			iter <- iter + 1
			
			if (iter >= maxIter) {
				warning(paste("No convergence in ", iter, " iterations."))
				return(as.SO3(R))
			}
		}
		R<-as.SO3(R)	
	}
	
	return(R)
}

#' Weighted Rotation Median
#' 
#' Compute the weighted projected or intrinsic mean of a sample of rotations
#'
#' This function takes a sample of n unit quaternions and approximates the mean rotation.  If the projected mean
#' is called for then the quaternions are turned reparameterized to matrices and mean.SO3 is called.  If the intrinsic
#' mean is called then according to \cite{gramkow01} a better approximation is achieved by taking average quaternion
#' and normalizing.  Our simulations don't match this claim.
#'
#' 
#' @param Qs A \eqn{n\times 4} matrix where each row corresponds to a random rotation in unit quaternion
#' @param w a numerical vector of weights the same length as Rs giving the weights to use for elements of Rs
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return weighted projected or intrinsic mean of the sample
#' @seealso \code{\link{mean.SO3}}
#' @cite moakher02, manton04
#' @S3method weighted.mean Q4
#' @method weighted.mean Q4
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' wt<-abs(1/r)
#' Qs<-genR(r,space="Q4")
#' weighted.mean(Qs,wt)

weighted.mean.Q4 <- function(Qs, w, type = "projected", epsilon = 1e-05, maxIter = 2000) {
	
	Qs<-formatQ4(Qs)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	Rs<-as.SO3(t(apply(Qs,1,SO3.Q4)))
	
	R<-weighted.mean(Rs,w,type,epsilon,maxIter)
	
	return(Q4.SO3(R))
	
}
