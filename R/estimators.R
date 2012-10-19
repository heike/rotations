#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of \eqn{3\times 3} rotations (in the form of a \eqn{n\times 9} matrix where n is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and the distance metric \eqn{d_D}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param x A \eqn{n\times 9} matrix where each row corresponds to a random rotation in matrix form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @param ... only used for consistency with mean.default
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{median.SO3}}
#' @cite moakher02, manton04
#' @S3method mean SO3
#' @method mean SO3
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' mean(Rs)

mean.SO3 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
  if (!all(apply(x, 1, is.SO3))) 
    warning("At least one of the observations is not in SO(3).  Use result with caution.")
  
  if (!(type %in% c("projected", "intrinsic")))
    stop("type needs to be one of 'projected' or 'intrinsic'.")
  
  R <- project.SO3(matrix(colMeans(x), 3, 3))
  
  
  if (type == "intrinsic") {
    n <- nrow(x)
    d <- 1
    iter <- 0
    s <- matrix(0, 3, 3)
    
    while (d >= epsilon) {
      
      R <- R %*% exp.skew(s)
      
      s <- matrix(colMeans(t(apply(x, 1, tLogMat, S = R))), 3, 3)
      
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

#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of n unit quaternions and approximates the mean rotation.  If the projected mean
#' is called for then the quaternions are turned reparameterized to matrices and mean.SO3 is called.  If the intrinsic
#' mean is called then according to \cite{gramkow01} a better approximation is achieved by taking average quaternion
#' and normalizing.  Our simulations don't match this claim.
#'
#' 
#' @param x A \eqn{n\times 4} matrix where each row corresponds to a random rotation in unit quaternion
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
#' r<-rvmises(20,0.01)
#' Qs<-genR(r,space="Q4")
#' mean(Qs,type='intrinsic')

mean.Q4 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  Qs <- x
  Rs<-t(apply(Qs,1,SO3.Q4))
  
  R<-mean.SO3(Rs,type,epsilon,maxIter)
  
  return(Q4.SO3(R))
  
}

#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of \eqn{3\times 3} rotations (in the form of a \eqn{n\times 9} matrix where n is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and the distance metric \eqn{d_D}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param x A \eqn{n\times 3} matrix where each row corresponds to a random rotation in Euler angle form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{mean.SO3}}
#' @cite moakher02, manton04
#' @S3method mean EA
#' @method mean EA
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' EAs<-genR(r,space="EA")
#' mean(EAs)

mean.EA <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  EAs <- x
  Rs<-as.SO3(t(apply(EAs,1,SO3.EA)))
  
  R<-mean.SO3(Rs,type,epsilon,maxIter)
  
  return(EA.SO3(R))
}


#' Compute the projected or intrinsic median estimate of the central direction
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice of distance metrid, \eqn{d_D}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#'
#' @param x A \eqn{n\times 9} matrix where each row corresponds to a random rotation in matrix form
#' @param ... additional arguments
#' @return the median
#' @seealso \code{\link{mean.SO3}}
#' @cite hartley11
#' @export

median<-function(x,...){
  UseMethod("median")
}

#' @return \code{NULL}
#' 
#' @rdname median
#' @method median SO3
#' @S3method median SO3

median.SO3 <- function(Rs, type = "projected", epsilon = 1e-05, maxIter = 2000, na.rm=FALSE) {
  
  if (!all(apply(Rs, 1, is.SO3))) 
    warning("At least one of the given observations is not in SO(3).  Use result with caution.")
  
  stopifnot(type %in% c("projected", "intrinsic"))
  
  S <- mean.SO3(Rs)
  d <- 1
  iter <- 1
  delta <- matrix(0, 3, 3)
  
  while (d >= epsilon) {
    
    if (type == "projected") {
      vn <- apply(Rs, 1, vecNorm, type = "F", S = S)
      
      delta <- matrix(colSums(Rs/vn)/sum(1/vn), 3, 3)
      
      Snew <- project.SO3(delta)
      
      d <- norm(Snew - S, type = "F")
      S <- Snew
      
    } else if (type == "intrinsic") {
      
      #      S <- exp.skew(delta) %*% S ## needs to be multiplied the other way round
      S <- S %*% exp.skew(delta)
      
      v <- t(apply(Rs, 1, tLogMat2, S = S))
      vn <- apply(v, 1, vecNorm, S = diag(0, 3, 3), type = "F")
      
      delta <- matrix(colSums(v/vn)/sum(1/vn), 3, 3)
      d <- norm(delta, type = "F")
      
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

#' @return \code{NULL}
#' 
#' @rdname median
#' @method median Q4
#' @S3method median Q4

median.Q4 <- function(Qs, type = "projected", epsilon = 1e-05, maxIter = 2000, na.rm=FALSE) {

  Rs<-t(apply(Qs,1,SO3.Q4))
  
  R<-median.SO3(Rs,type,epsilon,maxIter)
  
  return(Q4.SO3(R))
}


#' @return \code{NULL}
#' 
#' @rdname median
#' @method median EA
#' @S3method median EA

median.EA <- function(EAs, type = "projected", epsilon = 1e-05, maxIter = 2000, na.rm=FALSE) {

  Rs<-t(apply(EAs,1,SO3.EA))
  
  R<-median.SO3(Rs,type,epsilon,maxIter)
  
  return(EA.SO3(R))
}
