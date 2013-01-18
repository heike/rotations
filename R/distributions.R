#' Angular distributions in the rotations packages
#' 
#' Density and random variate generation for symmetric probability distributions in the rotations package
#' 
#' The functions for the density function and random variate generation are named in the usual form dxxxx and rxxxx 
#' respectively.        
#' \itemize{
#' 	\item See \code{\link{dcayley}} for the Cayley distribution.
#' 	\item See \code{\link{dfisher}}for the matrix Fisher distribution.
#' 	\item See \code{\link{dhaar}} for the uniform distribution on the circle.
#' 	\item See \code{\link{dvmises}} for the von Mises-Fisher distribution.
#' }
#' 
#' @name Angular-distributions

NULL

#' Sample of size n from target density f
#'
#' @author Heike Hofmann
#' @param n number of sample wanted
#' @param f target density
#' @param M maximum number in uniform proposal density
#' @param ... additional arguments sent to arsample
#' @return a vector of size n of observations from target density

rar <- function(n, f, M, ...) {
  res <- vector("numeric", length = n)
  for (i in 1:n) res[i] <- arsample.unif(f, M, ...)
  return(res)
}


#' The Symmetric Cayley Distribution
#'
#' Density and random generation for the Cayley distribution with concentration kappa
#'
#' The symmetric Cayley distribution with concentration kappa (or circular variance nu) had density 
#' \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}
#'
#' @name Cayley
#' @aliases Cayley rcayley dcayley
#' @usage dcayley(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE)
#' @usage rcayley(n, kappa = 1, nu = NULL)
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa Concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return \code{dcayley} gives the density,\code{pcayley} gives the distribution function, \code{rcayley} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package
#' @cite Schaeben97 leon06

NULL


#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

dcayley <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
 	den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}

#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

pcayley<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dcayley,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail) 
    return(cdf) else return((1-cdf))
}

#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
  return(theta)
}

#' The Matrix Fisher Distribution
#'
#' Density and random generation for the matrix Fisher distribution with concentration kappa
#'
#' The matrix Fisher distribution with concentration kappa (or circular variance nu) has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#'
#' @name Fisher
#' @aliases Fisher dfisher rfisher
#' @usage dfisher(r, kappa = 1, nu = NULL, Haar = TRUE)
#' @usage pfisher(r, kappa = 1, nu = NULL, lower.tail=TRUE)
#' @usage rfisher(n, kappa = 1, nu = NULL)
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return \code{dfisher} gives the density, \code{pfisher} gives the distribution function, \code{rfisher} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

dfisher <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  n<-length(r)
  den<-rep(0,n)
  
 	den <- exp(2 * kappa * cos(r)) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
  
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

pfisher<-function(q,kappa=1, nu= NULL, lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dfisher,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  M <- max(dfisher(seq(-pi, pi, length = 1000), kappa,Haar=F))
  return(rar(n, dfisher, M, kappa = kappa, Haar=F))
}

#' Haar Measure
#'
#' Uniform density on the circle
#' 
#' The uniform density on the circle  (also referred to as Haar measure)
#' has the density \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}
#'
#' @name Haar
#' @aliases Haar dhaar phaar rhaar
#' @usage dhaar(r)
#' @usage phaar(q, lower.tail=T)
#' @usage rhaar(n)
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @return \code{dhaar} gives the density, \code{phaar} gives the distribution function, \code{rhaar} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

dhaar <- function(r){
	
	den <-(1 - cos(r))/(2 * pi)
	
	return(den)
} 

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

phaar<-function(q,lower.tail=TRUE){
  
  cdf<-(q-sin(q)+pi)/(2*pi)
  
  ind<-which(cdf>1)
  cdf[ind]<-1
  
  ind2<-which(cdf<0)
  cdf[ind2]<-0
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export


rhaar<-function(n){
	
	lenn<-length(n)
	if(lenn>1)
		n<-lenn
	
  return(rar(n, dhaar, 1/pi))
}

#' The circular-von Mises distribution
#'
#' Density for the the circular von Mises-based distribution with concentration kappa
#' 
#' The circular von Mises-based distribution has the density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa cos(r)}.}
#'
#' @name Mises
#' @aliases Mises dvmises rvmises
#' @usage dvmises(r, kappa = 1, nu = NULL, Haar = TRUE)
#' @usage pvmises(q, kappa = 1, nu = NULL, lower.tail=TRUE)
#' @usage rvmises(n, kappa = 1, nu = NULL)
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @param lower.tail logica; if TRUE probabilites are \eqn{P(X\le x)} otherwise, \eqn{P(X>x)}
#' @return \code{dvmises} gives the density, \code{pvmises} gives the distribution function, \code{rvmises} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

dvmises <- function(r, kappa = 1, nu = NULL, Haar = T) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

pvmises<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dvmises,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail) 
    return(cdf) else return((1-cdf))
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  u <- runif(3, 0, 1)
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a))/(2 * kappa)
  r <- (1 + b^2)/(2 * b)
  theta <- rep(10, n)
  
  for (i in 1:n) {
    
    while (theta[i] == 10) {
      # Step 1
      u <- runif(3, 0, 1)
      z <- cos(pi * u[1])
      f <- (1 + r * z)/(r + z)
      c <- kappa * (r - f)
      
      # Step 2
      u <- runif(3, 0, 1)
      if ((c * (2 - c) - u[2]) > 0) {
        
        theta[i] = sign(u[3] - 0.5) * acos(f)
        
      } else {
        
        if (log(c/u[2]) + 1 - c < 0) {
          u <- runif(3, 0, 1)
        } else {
          u <- runif(3, 0, 1)
          theta[i] = sign(u[3] - 0.5) * acos(f)
        }
      }
    }
  }
  return(theta)
}



#' UARS density function
#' 
#' Evaluate the UARS density with a given angular distribution.
#' 
#' @param os Value at which to evaluate the UARS density
#' @param S principal direction of the distribution
#' @param kappa concentration of the distribution
#' @param dangle The function to evaulate the angles from: e.g. dcayley, dvmises, dfisher, dhaar
#' @param ... additional arguments passed to the angular distribution
#' @return density value at o
#' @export

duars<-function(os,S=diag(3),kappa=1,dangle,...){
	
	os<-formatSO3(os)
	rs<-angle(os)
	cr<-dangle(rs,kappa,...)	
	trStO<-2*cos(rs)+1
	
	den<-4*pi*cr/(3-trStO)
	
	return(den)
}

#' UARS distribution function
#' 
#' Evaluate the UARS distributions fuction with a given angular distribution
#' 
#' @param os Value at which to evaluate the UARS density
#' @param S principal direction of the distribution
#' @param kappa concentration of the distribution
#' @param dangle The function to evaulate the angles from: e.g. dcayley, dvmises, dfisher, dhaar.  If left at NULL, the empirical CDF is used
#' @param ... additional arguments passed to the angular distribution
#' @return cdf evaulated at each os value

puars<-function(os,S=diag(3),kappa=1,pangle=NULL,...){
	
	#This is not a true CDF, but it will work for now
	os<-formatSO3(os)
	rs<-angle(os)
	
	if(is.null(pangle)){
		
		n<-length(rs)
		cr<-rep(0,n)
		
		for(i in 1:length(rs))
			cr[i]<-length(which(rs<=rs[i]))/n
		
	}else{		
		cr<-pangle(rs,kappa,...)
	}
	
	#trStO<-2*cos(rs)+1
	
	#den<-4*pi*cr/(3-trStO)
	
	return(cr)
	
}
#' UARS random deviates
#' 
#' Produce random deviates from a chosen UARS distribution.
#' 
#' @param n number of observations. If \code{length(n)>1}, the length is taken to be n
#' @param rangle The function from which to simulate angles: e.g. rcayley, rvmises, rhaar, rfisher
#' @param S principal direction of the distribution
#' @param kappa concentration of the distribution
#' @param space Indicates the desired representation: matrix (SO3), quaternion (Q4) or Euler angles (EA)
#' @param ... additional arguments passed to the angular function
#' @return random deviates from the specified UARS distribution
#' @export

ruars<-function(n,rangle,S=NULL,kappa=1,space="SO3",...){
  
  r<-rangle(n,kappa,...)
  Rs<-genR(r,S,space)
  
  return(Rs)
}