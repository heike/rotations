#' Find confidence interval radius for central orientation estimate
#' 
#' Given a sample or rotations and estimator for the central orientation, a confidence interval radius will
#' be estimated based on a bootstrap procedure.  The level of confidence, bootstrap iterations and bootstrap
#' sample sizes are all options for the user to set.  The procedure is as follows:
#'           0) Estimate Shat from (R_1,...,R_n) in Rs
#'           1) Sample R_1,...,R_m from Rs with replacement
#'           2) Estimate ShatStar from bootstrap sample
#'           3) Compute That=d_r(Shat,ShatStar)
#'           4) Repeat steps 1-3 B times
#'           5) Report q% percentile of That to be the radius of the confidence 'cone'
#'  @param Rs sample of size n rotations in matrix format
#'  @param main The method of central orientation estimation; e.g. mean or median, other methods can be used, but will result in a warning
#'  @param B bootstrap iterations
#'  @param m bootstrap sample size
#'  @param alpha level of confidence
#'  @param ... additional arguments passed to 'fun' such as 'type' and 'epsilon'
#'  @return the radius of the confidence cone
#'  @references bingham10
#'  @export
#'  @examples
#'  Rs<-genR(rvmises(20))
#'  CIradius(Rs,mean,type='projected')  ## calls CIradius.SO3
#'  Qs<-genR(rvmises(20),space='Q4')
#'  CIradius(Qs,mean,type='projected')  ## calls CIradius.Q4


CIradius <- function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...)
{
  UseMethod( "CIradius" )
}


#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius SO3
#' @S3method CIradius SO3

CIradius.SO3 <- function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...){
  
  # This is more conservative than Bingham's method since d_r > max abs angle always
  args <- as.list(match.call())[-1]
  fname <- as.character(args$main)
  
  #  if (!(fname %in% c("mean", "median")))
  #    stop(sprintf("'%s' is not valid function for estimating main direction. Use 'mean' or 'median'.", fname))
  n<-nrow(Rs)
  Shat<-main(Rs,...)
  
  That <- rep(NA, B)
  for (i in 1:B) {  
    samp<-sample(1:n,m,replace=T)
    ShatStar<-main(as.SO3(Rs[samp,]),...)
    That[i] <- dist.SO3(Shat,ShatStar,method='intrinsic',p=1)
  }
  
  return(quantile(That,alpha))  
}

#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius Q4
#' @S3method CIradius Q4

CIradius.Q4<-function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...){  
  Rs<-as.SO3(t(apply(Rs,1,SO3.Q4)))
  #  args <- as.list(match.call())[-1]
  #  fname <- as.character(args$main)
  
  return(CIradius(Rs,main,B,m,alpha,...))
}

#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius EA
#' @S3method CIradius EA

CIradius.EA<-function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...){  
  Rs<-as.SO3(t(apply(Rs,1,SO3.EA)))
  #  args <- as.list(match.call())[-1]
  #  fname <- as.character(args$main)
  
  return(CIradius(Rs,main,B,m,alpha,...))
}
